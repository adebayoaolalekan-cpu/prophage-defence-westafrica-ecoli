#!/usr/bin/env python3
"""
Stratified selection of the genomic pilot set.

Prophage Burden, Prophage-Associated Virulence Determinants and Anti-Phage
Defence Systems in West African Escherichia coli.

The pilot has two jobs, and they need different sampling frames, so the
selection is built as two blocks.

Block A, the inference block. Drawn only from the two originating studies that
satisfy the eligibility rule for the human versus poultry contrast
(PRJNA293225-A and PRJNA1019343), restricted to the high-contiguity collection,
balanced across study and ecological source. This block supplies the observed
prophage count distribution and the negative binomial dispersion estimate that
Section 10.1 needs. Dispersion estimated anywhere else would not describe the
primary inference set.

Block B, the catalogue block. Drawn from the remainder of the retrieval,
stratified across ecological source, country and quality band. This block
supplies runtime, memory and disk figures, workflow failure rates, output
completeness, boundary truncation rates and the catalogue-wide burden
distribution. It deliberately includes genomes that fail the broad-quality
threshold, because the truncation rate among excluded genomes is what tells you
whether the quality filter is earning its place.

Assembly statistics exist for the two qualifying studies but not for the whole
retrieval, so block B is drawn in two stages. Stage 1 writes an oversampled
candidate frame. Assembly statistics are then fetched for that frame with
fetch_candidate_stats.sh. Stage 2 assigns quality bands and draws the final
block.

Usage
    python3 pilot_select.py --stage 1
    bash fetch_candidate_stats.sh          # or the documented API route
    python3 pilot_select.py --stage 2

Every draw uses a fixed seed. Re-running reproduces the same selection.
"""

import argparse
import csv
import json
import os
import random
import sys
from collections import Counter, defaultdict

SEED = 20261125

# Quality band thresholds, Table 9 of the protocol. Completeness and
# contamination come from CheckM2 inside the workflow and are therefore not
# available at selection time. Selection uses the assembly-statistic half of
# each definition. Genomes that later fail CheckM2 are reported as such rather
# than replaced, because replacing them after the fact would make the pilot
# a quality-selected sample of a quality-selected sample.
BROAD_MIN_N50 = 50_000
BROAD_MAX_CONTIGS = 300
HIGHC_MIN_N50 = 100_000
HIGHC_MAX_CONTIGS = 200

# Block A: cells are (originating_study, ecological_category). Twelve genomes
# per cell gives 24 human and 24 poultry, 24 per study.
BLOCK_A_PER_CELL = 12
BLOCK_A_CELLS = [
    ("PRJNA293225-A", "Human"),
    ("PRJNA293225-A", "Poultry"),
    ("PRJNA1019343", "Human"),
    ("PRJNA1019343", "Poultry"),
]

# Block B: target per ecological source. Human is the largest part of the
# retrieval and gets the largest allocation. Environmental is allocated above
# its share because it carries the catalogue description on its own now that
# the environmental contrasts have been dropped from inference.
BLOCK_B_TARGETS = {
    "Human": 15,
    "Environmental": 12,
    "Poultry": 9,
    "Other animal": 6,
    "Food": 6,
}

# Within each source, the split across quality bands.
BAND_SPLIT = {
    "high_contiguity": 0.50,
    "broad_only": 0.30,
    "below_broad": 0.20,
}

OVERSAMPLE = 3  # candidate frame size per final slot
MIN_DISTINCT_STUDIES = 8  # hard floor on block B study spread

QUALIFYING_STUDIES = {"PRJNA293225", "PRJNA1019343"}


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def assign_band(n50, contigs):
    """Return the quality band from assembly statistics alone."""
    if n50 is None or contigs is None:
        return "unknown"
    if n50 > HIGHC_MIN_N50 and contigs <= HIGHC_MAX_CONTIGS:
        return "high_contiguity"
    if n50 >= BROAD_MIN_N50 and contigs <= BROAD_MAX_CONTIGS:
        return "broad_only"
    return "below_broad"


def spread_draw(pool, n, key_funcs, rng):
    """
    Draw n records from pool, spreading them across the levels of each key in
    key_funcs in turn. Round-robin over the levels of the first key, then the
    second within each level, so that no single country or study swallows a
    cell. Falls back to a plain random draw once the spread keys are exhausted.
    """
    if n <= 0 or not pool:
        return []
    if not key_funcs:
        return rng.sample(pool, min(n, len(pool)))

    head, rest = key_funcs[0], key_funcs[1:]
    groups = defaultdict(list)
    for r in pool:
        groups[head(r)].append(r)

    # Visit the levels in order of ascending availability so scarce levels are
    # served before the common ones absorb the quota.
    levels = sorted(groups, key=lambda k: (len(groups[k]), str(k)))
    taken, want = [], n
    while want > 0:
        progress = False
        for lv in levels:
            if want == 0:
                break
            g = groups[lv]
            if not g:
                continue
            picked = spread_draw(g, 1, rest, rng)
            if not picked:
                continue
            rec = picked[0]
            g.remove(rec)
            taken.append(rec)
            want -= 1
            progress = True
        if not progress:
            break
    return taken


def to_int(v):
    try:
        return int(float(v))
    except (TypeError, ValueError):
        return None


# ----------------------------------------------------------------------------
# Stage 1
# ----------------------------------------------------------------------------

def stage1(args):
    rng = random.Random(SEED)

    qc = read_tsv(args.qc_table)
    cat = read_tsv(args.catalogue)

    # --- Block A ---------------------------------------------------------
    block_a, a_log = [], []
    for study, source in BLOCK_A_CELLS:
        pool = [r for r in qc
                if r["originating_study"] == study
                and r["ecological_category"] == source
                and r["passes_high_contiguity"].strip().upper() == "YES"]
        picked = spread_draw(pool, BLOCK_A_PER_CELL,
                             [lambda r: r.get("isolation_source", "")], rng)
        for r in picked:
            block_a.append({
                "assembly_accession": r["assembly_accession"],
                "block": "A_inference",
                "originating_study": study,
                "ecological_category": source,
                "country": r["country"],
                "collection_date": r["collection_date"],
                "isolation_source": r["isolation_source"],
                "assembly_level": r["assembly_level"],
                "contig_n50": r["contig_n50"],
                "number_of_contigs": r["number_of_contigs"],
                "total_length": r["total_length"],
                "quality_band": "high_contiguity",
                "stratum": f"{study}|{source}|high_contiguity",
            })
        a_log.append({
            "cell": f"{study}|{source}",
            "available": len(pool),
            "requested": BLOCK_A_PER_CELL,
            "drawn": len(picked),
            "shortfall": max(0, BLOCK_A_PER_CELL - len(picked)),
        })

    # --- Block B candidate frame -----------------------------------------
    # Exclude every accession already in block A, and every accession belonging
    # to the two qualifying studies. Block B must not borrow power from the
    # inference set, and a genome cannot sit in both blocks.
    in_a = {r["assembly_accession"] for r in block_a}
    qualifying_accs = {r["assembly_accession"] for r in qc}

    def acc_stem(a):
        return a.split(".")[0]

    excluded = in_a | {acc_stem(a) for a in qualifying_accs}

    eligible = [r for r in cat
                if acc_stem(r["assembly_accession"]) not in excluded
                and r["ecological_category"] in BLOCK_B_TARGETS]

    cand, b_log = [], []
    for source, target in BLOCK_B_TARGETS.items():
        pool = [r for r in eligible if r["ecological_category"] == source]
        want = target * OVERSAMPLE
        picked = spread_draw(
            pool, want,
            [lambda r: r["country_harmonised"], lambda r: r["bioproject"]],
            rng)
        for r in picked:
            cand.append({
                "assembly_accession": r["assembly_accession"],
                "block": "B_catalogue_candidate",
                "originating_study": r["bioproject"],
                "ecological_category": source,
                "country": r["country_harmonised"],
                "collection_date": r["collection_date"],
                "isolation_source": r["isolation_source_raw"],
                "umbrella_project_flag": r.get("umbrella_project_flag", ""),
            })
        b_log.append({
            "cell": source,
            "available": len(pool),
            "requested": want,
            "drawn": len(picked),
            "shortfall": max(0, want - len(picked)),
            "final_target": target,
        })

    write_tsv(args.out_block_a, block_a, list(block_a[0].keys()))
    write_tsv(args.out_candidates, cand, list(cand[0].keys()))

    with open(args.out_log, "w") as fh:
        json.dump({
            "seed": SEED,
            "stage": 1,
            "block_a": {
                "n": len(block_a),
                "cells": a_log,
                "source_counts": dict(Counter(r["ecological_category"] for r in block_a)),
                "study_counts": dict(Counter(r["originating_study"] for r in block_a)),
                "country_counts": dict(Counter(r["country"] for r in block_a)),
            },
            "block_b_candidate_frame": {
                "n": len(cand),
                "cells": b_log,
                "source_counts": dict(Counter(r["ecological_category"] for r in cand)),
                "country_counts": dict(Counter(r["country"] for r in cand)),
                "distinct_bioprojects": len({r["originating_study"] for r in cand}),
            },
            "thresholds": {
                "broad": {"min_n50": BROAD_MIN_N50, "max_contigs": BROAD_MAX_CONTIGS},
                "high_contiguity": {"min_n50": HIGHC_MIN_N50, "max_contigs": HIGHC_MAX_CONTIGS},
            },
        }, fh, indent=2)

    print(f"block A: {len(block_a)} genomes", file=sys.stderr)
    print(f"block B candidate frame: {len(cand)} accessions "
          f"({len({r['originating_study'] for r in cand})} bioprojects)",
          file=sys.stderr)
    for e in a_log + b_log:
        if e["shortfall"]:
            print(f"  shortfall in {e['cell']}: {e['shortfall']}", file=sys.stderr)


# ----------------------------------------------------------------------------
# Stage 2
# ----------------------------------------------------------------------------

def stage2(args):
    rng = random.Random(SEED + 1)

    cand = read_tsv(args.out_candidates)
    stats = {}
    with open(args.candidate_stats, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            stats[r["assembly_accession"].split(".")[0]] = r

    unresolved = []
    for r in cand:
        s = stats.get(r["assembly_accession"].split(".")[0])
        if s is None:
            unresolved.append(r["assembly_accession"])
            r["contig_n50"] = r["number_of_contigs"] = r["total_length"] = ""
            r["assembly_level"] = ""
            r["quality_band"] = "unknown"
            continue
        n50 = to_int(s.get("contig_n50"))
        ctg = to_int(s.get("number_of_contigs"))
        r["contig_n50"] = "" if n50 is None else n50
        r["number_of_contigs"] = "" if ctg is None else ctg
        r["total_length"] = s.get("total_length", "")
        r["assembly_level"] = s.get("assembly_level", "")
        r["quality_band"] = assign_band(n50, ctg)

    usable = [r for r in cand if r["quality_band"] != "unknown"]

    block_b, b_log = [], []
    for source, target in BLOCK_B_TARGETS.items():
        src_pool = [r for r in usable if r["ecological_category"] == source]
        # Integer allocation across bands, remainder to the largest share.
        alloc = {b: int(target * f) for b, f in BAND_SPLIT.items()}
        rem = target - sum(alloc.values())
        for b in sorted(BAND_SPLIT, key=lambda k: -BAND_SPLIT[k]):
            if rem == 0:
                break
            alloc[b] += 1
            rem -= 1

        drawn_source = []
        deficit = 0
        for band, k in alloc.items():
            pool = [r for r in src_pool if r["quality_band"] == band]
            picked = spread_draw(
                pool, k,
                [lambda r: r["country"], lambda r: r["originating_study"]],
                rng)
            drawn_source.extend(picked)
            b_log.append({
                "cell": f"{source}|{band}",
                "available": len(pool),
                "requested": k,
                "drawn": len(picked),
                "shortfall": max(0, k - len(picked)),
            })
            deficit += max(0, k - len(picked))

        # Backfill a shortfall from the other bands of the same source, so the
        # source total is held and the loss falls on the band balance rather
        # than on the ecological balance. Band balance is a convenience.
        # Source balance is what the catalogue description rests on.
        if deficit:
            taken = {r["assembly_accession"] for r in drawn_source}
            spare = [r for r in src_pool if r["assembly_accession"] not in taken]
            extra = spread_draw(
                spare, deficit,
                [lambda r: r["quality_band"], lambda r: r["country"]], rng)
            for r in extra:
                r["backfilled"] = "YES"
            drawn_source.extend(extra)
            b_log.append({
                "cell": f"{source}|backfill",
                "available": len(spare),
                "requested": deficit,
                "drawn": len(extra),
                "shortfall": max(0, deficit - len(extra)),
            })

        for r in drawn_source:
            r["block"] = "B_catalogue"
            r["stratum"] = (f"{r['ecological_category']}|{r['country']}|"
                            f"{r['quality_band']}")
            r.setdefault("backfilled", "NO")
        block_b.extend(drawn_source)

    n_studies = len({r["originating_study"] for r in block_b})

    block_a = read_tsv(args.out_block_a)
    for r in block_a:
        r.setdefault("umbrella_project_flag", "")
        r.setdefault("backfilled", "NO")

    fields = ["assembly_accession", "block", "stratum", "originating_study",
              "ecological_category", "country", "collection_date",
              "isolation_source", "assembly_level", "contig_n50",
              "number_of_contigs", "total_length", "quality_band",
              "umbrella_project_flag", "backfilled"]
    combined = block_a + block_b
    write_tsv(args.out_manifest, combined, fields)

    # The sample sheet the workflow consumes. Kept separate and minimal so a
    # metadata correction does not force the workflow to re-run.
    write_tsv(args.out_samplesheet,
              [{"sample": r["assembly_accession"],
                "accession": r["assembly_accession"],
                "block": r["block"],
                "quality_band": r["quality_band"],
                "ecological_category": r["ecological_category"],
                "originating_study": r["originating_study"]}
               for r in combined],
              ["sample", "accession", "block", "quality_band",
               "ecological_category", "originating_study"])

    summary = {
        "seed": SEED + 1,
        "stage": 2,
        "n_total": len(combined),
        "n_block_a": len(block_a),
        "n_block_b": len(block_b),
        "block_b_cells": b_log,
        "block_b_distinct_studies": n_studies,
        "block_b_study_floor": MIN_DISTINCT_STUDIES,
        "block_b_study_floor_met": n_studies >= MIN_DISTINCT_STUDIES,
        "unresolved_candidate_accessions": len(unresolved),
        "composition": {
            "by_source": dict(Counter(r["ecological_category"] for r in combined)),
            "by_band": dict(Counter(r["quality_band"] for r in combined)),
            "by_country": dict(Counter(r["country"] for r in combined)),
            "by_assembly_level": dict(Counter(r["assembly_level"] for r in combined)),
            "by_source_and_band": {
                f"{a}|{b}": c for (a, b), c in
                Counter((r["ecological_category"], r["quality_band"])
                        for r in combined).items()},
        },
        "backfilled": sum(1 for r in combined if r.get("backfilled") == "YES"),
    }
    with open(args.out_summary, "w") as fh:
        json.dump(summary, fh, indent=2)

    print(f"pilot set: {len(combined)} genomes "
          f"(A={len(block_a)}, B={len(block_b)})", file=sys.stderr)
    print(f"block B distinct originating studies: {n_studies} "
          f"(floor {MIN_DISTINCT_STUDIES}: "
          f"{'met' if n_studies >= MIN_DISTINCT_STUDIES else 'NOT MET'})",
          file=sys.stderr)
    if unresolved:
        print(f"WARNING {len(unresolved)} candidate accessions had no assembly "
              f"statistics and were dropped from the frame", file=sys.stderr)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stage", type=int, choices=(1, 2), required=True)
    p.add_argument("--catalogue",
                   default="resources/ena_westafrica_ecoli_metadata.tsv")
    p.add_argument("--qc-table",
                   default="resources/qualifying_studies_assembly_qc.tsv")
    p.add_argument("--candidate-stats",
                   default="resources/candidate_assembly_stats.tsv")
    p.add_argument("--out-block-a", default="config/pilot_block_a.tsv")
    p.add_argument("--out-candidates", default="config/pilot_block_b_candidates.tsv")
    p.add_argument("--out-manifest", default="config/pilot_manifest.tsv")
    p.add_argument("--out-samplesheet", default="config/samples.tsv")
    p.add_argument("--out-log", default="config/pilot_selection_stage1.json")
    p.add_argument("--out-summary", default="config/pilot_selection_summary.json")
    a = p.parse_args()
    os.makedirs(os.path.dirname(a.out_manifest) or ".", exist_ok=True)
    (stage1 if a.stage == 1 else stage2)(a)


if __name__ == "__main__":
    main()
