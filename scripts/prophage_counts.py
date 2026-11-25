"""
The burden table: one row per genome, with the counts the analysis uses.

Protocol Step 9. The primary burden measure is the count of internally
resolved geNomad predictions per genome. The count restricted to predictions
supported by both geNomad and the independent predictor is the sensitivity
analysis. Boundary-truncated candidates are carried in their own columns and
never folded into the primary count.

Every covariate the Section 10.1 model needs is attached here, so the
re-simulation reads one file: source, study, country, N50, contig count,
completeness, quality band and the PopPUNK lineage.
"""

import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv, write_tsv, num

FIELDS = ["sample", "block", "ecological_category", "originating_study",
          "country", "collection_date", "lineage", "quality_band",
          "passes_broad", "passes_high_contiguity",
          "contig_n50", "number_of_contigs", "total_length",
          "completeness", "contamination",
          "n_prophage_primary", "n_prophage_all_predictions",
          "n_truncated_left", "n_truncated_right", "n_truncated_both",
          "n_prophage_intersection", "n_prophage_high_concordance",
          "n_refined", "n_no_breakpoint", "n_unstable_boundary",
          "n_shared_region", "n_clipped_flank",
          "total_prophage_bp", "median_prophage_length",
          "prophage_fraction_of_genome"]


def median(vals):
    v = sorted(vals)
    if not v:
        return ""
    n = len(v)
    return v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2


def main(samples, outdir, bands_path, lineage_path, manifest_path, out_path):
    bands = {r["sample"]: r for r in read_tsv(bands_path)}
    lineage = {r["sample"]: r["lineage"] for r in read_tsv(lineage_path)}
    man = {r["assembly_accession"]: r for r in read_tsv(manifest_path)}
    pairs = read_tsv(f"{outdir}/reports/concordance_pairs.tsv")

    matched = {}
    for p in pairs:
        matched.setdefault(p["sample"], {})[p["genomad_id"]] = \
            p["concordance_class"]

    rows = []
    for s in samples:
        ref = read_tsv(f"{outdir}/prophage/refined/{s}.tsv")
        b = bands.get(s, {})
        m = man.get(s, {})
        mt = matched.get(s, {})

        primary = [r for r in ref
                   if r["boundary_class_refined"] == "internally_resolved"]
        bc = Counter(r["boundary_class_refined"] for r in ref)
        lengths = [int(r["refined_length"]) for r in primary
                   if num(r["refined_length"]) is not None]
        total_len = num(b.get("total_length"))

        rows.append({
            "sample": s,
            "block": b.get("block") or m.get("block", ""),
            "ecological_category": b.get("ecological_category")
                                   or m.get("ecological_category", ""),
            "originating_study": b.get("originating_study")
                                  or m.get("originating_study", ""),
            "country": b.get("country") or m.get("country", ""),
            "collection_date": m.get("collection_date", ""),
            "lineage": lineage.get(s, ""),
            "quality_band": b.get("quality_band", ""),
            "passes_broad": b.get("passes_broad", ""),
            "passes_high_contiguity": b.get("passes_high_contiguity", ""),
            "contig_n50": b.get("contig_n50", ""),
            "number_of_contigs": b.get("number_of_contigs", ""),
            "total_length": b.get("total_length", ""),
            "completeness": b.get("completeness", ""),
            "contamination": b.get("contamination", ""),
            "n_prophage_primary": len(primary),
            "n_prophage_all_predictions": len(ref),
            "n_truncated_left": bc.get("truncated_left", 0),
            "n_truncated_right": bc.get("truncated_right", 0),
            "n_truncated_both": bc.get("truncated_both", 0),
            "n_prophage_intersection": sum(
                1 for r in primary if r["prediction_id"] in mt),
            "n_prophage_high_concordance": sum(
                1 for r in primary if mt.get(r["prediction_id"]) == "high"),
            "n_refined": sum(1 for r in ref
                             if r["refinement_status"].startswith("refined")),
            "n_no_breakpoint": sum(
                1 for r in ref
                if r["refinement_status"] == "no_breakpoint_geNomad_retained"),
            "n_unstable_boundary": sum(1 for r in ref
                                       if r["boundary_stability"] == "unstable"),
            "n_shared_region": sum(1 for r in ref
                                   if r["shared_region"] == "YES"),
            "n_clipped_flank": sum(1 for r in ref
                                   if r["flank_left_clipped"] == "YES"
                                   or r["flank_right_clipped"] == "YES"),
            "total_prophage_bp": sum(lengths),
            "median_prophage_length": median(lengths),
            "prophage_fraction_of_genome": (
                round(sum(lengths) / total_len, 5)
                if total_len else ""),
        })

    write_tsv(out_path, rows, FIELDS)

    counts = [r["n_prophage_primary"] for r in rows]
    sys.stderr.write(
        f"burden: {len(rows)} genomes, primary count mean "
        f"{sum(counts)/len(counts):.3f}, range {min(counts)}-{max(counts)}\n"
        if counts else "burden: no genomes\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(list(sm.params.samples), sm.params.outdir, sm.input.bands,
             sm.input.lineage, sm.input.manifest, sm.output.tsv)
    else:
        main(sys.argv[1].split(","), sys.argv[2], sys.argv[3], sys.argv[4],
             sys.argv[5], sys.argv[6])
