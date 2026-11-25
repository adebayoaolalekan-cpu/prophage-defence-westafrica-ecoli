#!/usr/bin/env python3
"""
Integration test for the pilot workflow's coordinate handling and statistics.

The tool calls themselves cannot be tested without the tools, but every piece
of arithmetic between them can be, and that is where the errors hide. This
builds a synthetic genome with proviruses at known coordinates, drives the
chain, and asserts the coordinates that come out the far end.

Run: python3 -I tests/test_chain.py
"""
import csv
import json
import os
import random
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
SCRIPTS = os.path.join(ROOT, "scripts")
sys.path.insert(0, SCRIPTS)

from _common import read_fasta, read_tsv  # noqa: E402

FAILS = []


def check(cond, msg):
    if cond:
        print(f"  pass  {msg}")
    else:
        print(f"  FAIL  {msg}")
        FAILS.append(msg)


def rand_dna(n, rng):
    return "".join(rng.choice("ACGT") for _ in range(n))


def run(script, *args):
    r = subprocess.run([sys.executable, "-I", os.path.join(SCRIPTS, script)]
                       + [str(a) for a in args],
                       capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout)
        print(r.stderr)
        raise SystemExit(f"{script} exited {r.returncode}")
    return r.stderr


def write_tsv(path, rows, fields):
    os.makedirs(os.path.dirname(path) or ".", exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)


def main():
    rng = random.Random(7)
    tmp = tempfile.mkdtemp(prefix="pilot_test_")
    S = "GCA_TEST001"

    # ---------------------------------------------------------------- genome
    # contig1 200 kb: an internal provirus at 100000-140000, and a second at
    # 155000-170000. Two windows merge when the gap between the predictions is
    # no more than twice the flank, so a 15 kb gap merges at the 10 kb flank
    # and separates at 5 kb. That is the case worth testing, because it is the
    # one where the primary and sensitivity extractions differ in structure
    # rather than only in length.
    # contig2  60 kb: a provirus starting at position 1, truncated left.
    # contig3  30 kb: no prediction.
    contigs = {
        f"{S}|contig1": rand_dna(200000, rng),
        f"{S}|contig2": rand_dna(60000, rng),
        f"{S}|contig3": rand_dna(30000, rng),
    }
    fna = os.path.join(tmp, "genome.fna")
    with open(fna, "w") as fh:
        for n, s in contigs.items():
            fh.write(f">{n}\n{s}\n")

    # ------------------------------------------------- synthetic geNomad out
    summary = os.path.join(tmp, "virus_summary.tsv")
    write_tsv(summary, [
        {"seq_name": f"{S}|contig1|provirus_100000_140000", "length": 40001,
         "topology": "Provirus", "coordinates": "100000-140000", "n_genes": 45,
         "genetic_code": 11, "virus_score": 0.98, "fdr": "NA",
         "n_hallmarks": 6, "marker_enrichment": 12.1,
         "taxonomy": "Viruses;Duplodnaviria;Caudoviricetes"},
        {"seq_name": f"{S}|contig1|provirus_155000_170000", "length": 15001,
         "topology": "Provirus", "coordinates": "155000-170000", "n_genes": 20,
         "genetic_code": 11, "virus_score": 0.91, "fdr": "NA",
         "n_hallmarks": 3, "marker_enrichment": 7.4, "taxonomy": "Viruses"},
        {"seq_name": f"{S}|contig2|provirus_1_25000", "length": 25000,
         "topology": "Provirus", "coordinates": "1-25000", "n_genes": 30,
         "genetic_code": 11, "virus_score": 0.95, "fdr": "NA",
         "n_hallmarks": 4, "marker_enrichment": 9.0, "taxonomy": "Viruses"},
        # A free virus and a plasmid-like call, which must not be counted.
        {"seq_name": f"{S}|contig3", "length": 30000, "topology": "DTR",
         "coordinates": "NA", "n_genes": 35, "genetic_code": 11,
         "virus_score": 0.99, "fdr": "NA", "n_hallmarks": 8,
         "marker_enrichment": 15.0, "taxonomy": "Viruses"},
    ], ["seq_name", "length", "topology", "coordinates", "n_genes",
        "genetic_code", "virus_score", "fdr", "n_hallmarks",
        "marker_enrichment", "taxonomy"])

    # ------------------------------------------------------- provirus_table
    print("\nprovirus_table.py")
    prov = os.path.join(tmp, "prov.tsv")
    run("provirus_table.py", summary, fna, prov, S)
    P = read_tsv(prov)
    check(len(P) == 3, "three integrated proviruses kept, free virus dropped")
    excl = read_tsv(prov.replace(".tsv", ".excluded.tsv"))
    check(len(excl) == 1, "the non-provirus prediction is recorded, not lost")
    byid = {r["prediction_id"]: r for r in P}
    p1, p2, p3 = P[0], P[1], P[2]
    check(p1["contig"] == f"{S}|contig1" and p1["genomad_start"] == "100000",
          "provirus suffix stripped, contig name with its sample prefix intact")
    check(p1["contig_length"] == "200000", "contig length read from the FASTA")
    check(p1["boundary_class"] == "internally_resolved",
          "internal prediction classified as internally resolved")
    check(p3["boundary_class"] == "truncated_left",
          "prediction starting at position 1 classified as truncated left")

    # ------------------------------------------------------- checkv_extract
    print("\ncheckv_extract.py, 10 kb flank")
    f10 = os.path.join(tmp, "ex.primary.fna")
    m10 = os.path.join(tmp, "ex.primary.map.tsv")
    run("checkv_extract.py", fna, prov, f10, m10, 10000, "true", S)
    seqs10 = dict(read_fasta(f10))
    M10 = read_tsv(m10)
    check(len(seqs10) == 2,
          "two regions at 10 kb: the overlapping pair merged, contig2 separate")
    shared = [m for m in M10 if m["shared_region"] == "YES"]
    check(len(shared) == 2, "both contig1 predictions share one region")
    r1 = [m for m in M10 if m["prediction_id"] == p1["prediction_id"]][0]
    # window1 = 90000-150000, window2 = 145000-180000, merged = 90000-180000
    check(r1["region_start"] == "90000" and r1["region_end"] == "180000",
          "merged region spans 90000 to 180000")
    check(len(seqs10[r1["region_id"]]) == 180000 - 90000 + 1,
          "extracted sequence length equals the inclusive region length")
    check(seqs10[r1["region_id"]] == contigs[f"{S}|contig1"][89999:180000],
          "extracted sequence is the right slice of the contig")
    check(r1["offset_start_in_region"] == "10001",
          "prediction start sits at offset 10001 inside the region")
    r3 = [m for m in M10 if m["prediction_id"] == p3["prediction_id"]][0]
    check(r3["region_start"] == "1" and r3["flank_left_clipped"] == "YES",
          "left flank clipped at the contig start and flagged")
    check(r3["region_end"] == "35000" and r3["flank_right_clipped"] == "NO",
          "right flank applied in full on contig2")

    print("\ncheckv_extract.py, 5 kb flank")
    f5 = os.path.join(tmp, "ex.sens.fna")
    m5 = os.path.join(tmp, "ex.sens.map.tsv")
    run("checkv_extract.py", fna, prov, f5, m5, 5000, "true", S)
    M5 = read_tsv(m5)
    check(len(dict(read_fasta(f5))) == 3,
          "three regions at 5 kb: a 15 kb gap separates the contig1 windows")
    check(all(m["shared_region"] == "NO" for m in M5),
          "no shared regions at the 5 kb flank")

    # The merge rule's consequence, asserted rather than left implicit: two
    # predictions closer together than twice the smaller flank share a region
    # at both flank lengths, so their boundaries cannot be refined
    # independently at either. The shared_region column is what marks those
    # cases for the analyst.
    near = os.path.join(tmp, "near.tsv")
    write_tsv(near, [
        dict(P[0], prediction_id="NEAR_A"),
        dict(P[0], prediction_id="NEAR_B", genomad_start="143000",
             genomad_end="150000"),
    ], list(P[0].keys()))
    fn5 = os.path.join(tmp, "near.5.fna")
    mn5 = os.path.join(tmp, "near.5.map.tsv")
    run("checkv_extract.py", fna, near, fn5, mn5, 5000, "true", S)
    check(len(dict(read_fasta(fn5))) == 1
          and all(m["shared_region"] == "YES" for m in read_tsv(mn5)),
          "predictions 3 kb apart share one region even at the 5 kb flank")

    # --------------------------------------------- synthetic CheckV results
    # On the merged 10 kb region (90000-170000, length 80001) CheckV is told
    # the viral parts are at region offsets 10001-50001 and 55001-70001, which
    # are exactly the two predictions. Refined contig coordinates must come
    # back as 100000-140000 and 145000-160000.
    def checkv_files(outdir, maps, viral_by_region, provirus_flag="Yes"):
        rows_c, rows_k = [], []
        for rid in sorted({m["region_id"] for m in maps}):
            rl = int([m for m in maps if m["region_id"] == rid][0]["region_length"])
            segs = viral_by_region.get(rid)
            if segs is None:
                rows_c.append({"contig_id": rid, "contig_length": rl,
                               "viral_length": rl, "host_length": 0,
                               "provirus": "No", "region_types": "viral",
                               "region_lengths": rl,
                               "region_coords_bp": f"1-{rl}"})
            else:
                types, lens, coords, hl, vl = [], [], [], 0, 0
                cur = 1
                for (a, b) in segs:
                    if a > cur:
                        types.append("host"); lens.append(a - cur)
                        coords.append(f"{cur}-{a-1}"); hl += a - cur
                    types.append("viral"); lens.append(b - a + 1)
                    coords.append(f"{a}-{b}"); vl += b - a + 1
                    cur = b + 1
                if cur <= rl:
                    types.append("host"); lens.append(rl - cur + 1)
                    coords.append(f"{cur}-{rl}"); hl += rl - cur + 1
                rows_c.append({
                    "contig_id": rid, "contig_length": rl,
                    "viral_length": vl, "host_length": hl,
                    "provirus": provirus_flag,
                    "region_types": ",".join(types),
                    "region_lengths": ",".join(str(x) for x in lens),
                    "region_coords_bp": ",".join(coords)})
            rows_k.append({"contig_id": rid, "contig_length": rl,
                           "provirus": provirus_flag, "proviral_length": rl,
                           "gene_count": 50, "viral_genes": 30,
                           "host_genes": 12, "checkv_quality": "High-quality",
                           "miuvig_quality": "High-quality",
                           "completeness": 94.2,
                           "completeness_method": "AAI-based",
                           "contamination": 0.0, "kmer_freq": 1.0,
                           "warnings": ""})
        os.makedirs(outdir, exist_ok=True)
        c = os.path.join(outdir, "contamination.tsv")
        k = os.path.join(outdir, "completeness.tsv")
        write_tsv(c, rows_c, ["contig_id", "contig_length", "viral_length",
                              "host_length", "provirus", "region_types",
                              "region_lengths", "region_coords_bp"])
        write_tsv(k, rows_k, list(rows_k[0].keys()))
        return c, k

    rid1 = r1["region_id"]
    rid3 = r3["region_id"]
    pc, pk = checkv_files(os.path.join(tmp, "cvp"), M10, {
        rid1: [(10001, 50001), (55001, 70001)],
        rid3: [(1, 25000)],
    })
    m5_by_pred = {m["prediction_id"]: m for m in M5}
    # At 5 kb the first prediction's region is 95000-145000, and CheckV is
    # given a viral segment shifted 300 bp, which is inside the tolerance.
    sens_viral = {}
    for pid, m in m5_by_pred.items():
        rs = int(m["region_start"])
        os_, oe_ = int(m["offset_start_in_region"]), int(m["offset_end_in_region"])
        shift = 300 if pid == p1["prediction_id"] else 0
        sens_viral[m["region_id"]] = [(max(1, os_ + shift), oe_ + shift)]
    sc, _ = checkv_files(os.path.join(tmp, "cvs"), M5, sens_viral)

    print("\nrefine_boundaries.py")
    refined = os.path.join(tmp, "refined.tsv")
    run("refine_boundaries.py", prov, pc, pk, m10, sc, m5, refined, S)
    R = {r["prediction_id"]: r for r in read_tsv(refined)}
    a = R[p1["prediction_id"]]
    b = R[p2["prediction_id"]]
    c = R[p3["prediction_id"]]
    check(a["refined_start"] == "100000" and a["refined_end"] == "140000",
          "first prediction refines back to its original contig coordinates")
    check(b["refined_start"] == "145000" and b["refined_end"] == "160000",
          "second prediction in the shared region gets the right viral segment")
    check(a["refinement_status"] == "refined_shared_region",
          "shared-region refinement is labelled as such")
    check(a["boundary_class_refined"] == "internally_resolved",
          "refined boundary class recomputed")
    check(c["boundary_class_refined"] == "truncated_left",
          "contig-start truncation survives refinement")
    check(a["start_shift"] == "300" and a["boundary_stability"] == "stable",
          "boundary shift between the two flanks measured, 300 bp is stable")

    # A region CheckV declines to call a provirus must retain geNomad.
    pc2, pk2 = checkv_files(os.path.join(tmp, "cvp2"), M10, {}, "No")
    refined2 = os.path.join(tmp, "refined2.tsv")
    run("refine_boundaries.py", prov, pc2, pk2, m10, sc, m5, refined2, S)
    R2 = {r["prediction_id"]: r for r in read_tsv(refined2)}
    check(R2[p1["prediction_id"]]["refinement_status"]
          == "no_breakpoint_geNomad_retained"
          and R2[p1["prediction_id"]]["refined_start"] == "100000",
          "no CheckV breakpoint keeps the geNomad coordinates and says so")

    # ------------------------------------------------------- concordance
    print("\nconcordance.py")
    out = os.path.join(tmp, "res")
    for sub in ("prophage/genomad", "prophage/independent", "reports"):
        os.makedirs(os.path.join(out, sub), exist_ok=True)
    import shutil
    shutil.copy(prov, os.path.join(out, "prophage/genomad", f"{S}.tsv"))
    write_tsv(os.path.join(out, "prophage/independent", f"{S}.tsv"), [
        # Near-identical to prediction 1: high concordance.
        {"sample": S, "prediction_id": f"{S}_ind0001", "tool": "pide",
         "contig": f"{S}|contig1", "contig_length": 200000,
         "start": 100500, "end": 139600, "length": 39101, "score": 0.9},
        # Engulfs prediction 2: overlaps, but reciprocal overlap fails.
        {"sample": S, "prediction_id": f"{S}_ind0002", "tool": "pide",
         "contig": f"{S}|contig1", "contig_length": 200000,
         "start": 141000, "end": 199000, "length": 58001, "score": 0.7},
        # Specific to the independent tool.
        {"sample": S, "prediction_id": f"{S}_ind0003", "tool": "pide",
         "contig": f"{S}|contig3", "contig_length": 30000,
         "start": 5000, "end": 20000, "length": 15001, "score": 0.8},
    ], ["sample", "prediction_id", "tool", "contig", "contig_length",
        "start", "end", "length", "score"])
    cs = os.path.join(out, "reports/concordance.tsv")
    cp = os.path.join(out, "reports/concordance_pairs.tsv")
    run("concordance.py", S, out, cs, cp, 0.80, 2000, "pide")
    CS = {r["sample"]: r for r in read_tsv(cs)}["ALL"]
    CP = {r["genomad_id"]: r for r in read_tsv(cp)}
    check(CS["n_matched"] == "2", "two overlapping pairs found")
    check(CS["n_high_concordance"] == "1",
          "only the near-identical pair meets the prespecified definition")
    check(CP[p1["prediction_id"]]["concordance_class"] == "high",
          "prediction 1 matched at high concordance")
    check(CP[p2["prediction_id"]]["concordance_class"] == "overlapping_only",
          "an engulfing prediction fails reciprocal overlap, as intended")
    check(CS["n_genomad_only"] == "1" and CS["n_independent_only"] == "1",
          "tool-specific predictions counted on both sides")

    # ------------------------------------------------------ dispersion
    print("\ndispersion.py, theta recovery from a known distribution")
    nrng = __import__("numpy").random.default_rng(11)
    TRUE_THETA = 6.0
    rows = []
    for i in range(96):
        blk = "A_inference" if i < 48 else "B_catalogue"
        src = ("Human" if i % 2 == 0 else "Poultry") if i < 48 else \
              ["Human", "Poultry", "Environmental", "Food"][i % 4]
        mu = 1.4 if src == "Poultry" else 0.8
        n = int(nrng.poisson(nrng.gamma(shape=TRUE_THETA,
                                        scale=mu / TRUE_THETA)))
        rows.append({
            "sample": f"G{i:03d}", "block": blk, "ecological_category": src,
            "originating_study": "PRJNA293225-A" if i % 3 else "PRJNA1019343",
            "country": "Nigeria" if i % 2 else "Ghana",
            "collection_date": "2019" if i % 2 else "2022",
            "lineage": str(i % 12), "quality_band": "high_contiguity",
            "passes_broad": "YES", "passes_high_contiguity": "YES",
            "contig_n50": 120000 + i * 500, "number_of_contigs": 150 + i,
            "total_length": 5000000, "completeness": 98.0 + (i % 5) * 0.2,
            "contamination": 0.5,
            "n_prophage_primary": n, "n_prophage_all_predictions": n + 1,
            "n_truncated_left": 0, "n_truncated_right": 1,
            "n_truncated_both": 0,
            "n_prophage_intersection": max(0, n - (i % 2)),
            "n_prophage_high_concordance": max(0, n - 1),
            "n_refined": n, "n_no_breakpoint": 0, "n_unstable_boundary": 0,
            "n_shared_region": 0, "n_clipped_flank": 0,
            "total_prophage_bp": n * 40000,
            "median_prophage_length": 40000,
            "prophage_fraction_of_genome": round(n * 40000 / 5000000, 5),
        })
    cnt = os.path.join(tmp, "counts.tsv")
    write_tsv(cnt, rows, list(rows[0].keys()))
    dj = os.path.join(tmp, "disp.json")
    dt = os.path.join(tmp, "disp.tsv")
    run("dispersion.py", cnt, dj, dt)
    D = json.load(open(dj))
    inf = next(e for e in D["estimates"]
               if e["set"] == "inference_block_high_contiguity"
               and e.get("measure") == "n_prophage_primary")
    print(f"        true theta {TRUE_THETA}, mom "
          f"{inf['theta_method_of_moments']}, mle "
          f"{inf['theta_mle_intercept_only']}, adjusted "
          f"{inf['theta_mle_covariate_adjusted']}, "
          f"ci [{inf['theta_bootstrap_ci_low']}, "
          f"{inf['theta_bootstrap_ci_high']}]")
    check(inf["n"] == 48, "inference block isolated correctly")
    check(inf["theta_for_resimulation"] is not None,
          "a theta is reported for the re-simulation")
    check(inf["theta_bootstrap_ci_low"] is not None,
          "a lineage-clustered bootstrap interval is reported")
    check(len(D["estimates"]) >= 8,
          "every set reported for both the primary and sensitivity measure")
    check(D["verdict"]["status"] in
          ("rr_1_5_supported", "rr_1_5_borderline", "rr_2_0_only",
           "poisson_adequate"),
          f"a verdict is issued: {D['verdict']['status']}")

    # ---------------------------------------------- selection reproducibility
    print("\npilot_select.py reproducibility")
    r1_ = subprocess.run([sys.executable, "-I",
                          os.path.join(SCRIPTS, "pilot_select.py"),
                          "--stage", "1",
                          "--catalogue", os.path.join(ROOT, "resources/ena_westafrica_ecoli_metadata.tsv"),
                          "--qc-table", os.path.join(ROOT, "resources/qualifying_studies_assembly_qc.tsv"),
                          "--out-block-a", os.path.join(tmp, "a1.tsv"),
                          "--out-candidates", os.path.join(tmp, "c1.tsv"),
                          "--out-manifest", os.path.join(tmp, "m1.tsv"),
                          "--out-samplesheet", os.path.join(tmp, "s1.tsv"),
                          "--out-log", os.path.join(tmp, "l1.json"),
                          "--out-summary", os.path.join(tmp, "y1.json")],
                         capture_output=True, text=True)
    check(r1_.returncode == 0, "stage 1 runs against the real catalogue")
    a_rows = read_tsv(os.path.join(tmp, "a1.tsv"))
    check(len(a_rows) == 48, "block A holds 48 genomes")
    from collections import Counter as C
    check(dict(C(r["ecological_category"] for r in a_rows))
          == {"Human": 24, "Poultry": 24},
          "block A is balanced 24 human against 24 poultry")
    check(all(r["quality_band"] == "high_contiguity" for r in a_rows),
          "block A is entirely high-contiguity")
    c_rows = read_tsv(os.path.join(tmp, "c1.tsv"))
    overlap = {r["assembly_accession"].split(".")[0] for r in a_rows} & \
              {r["assembly_accession"].split(".")[0] for r in c_rows}
    check(not overlap, "no accession appears in both blocks")
    qual = {r["assembly_accession"].split(".")[0]
            for r in read_tsv(os.path.join(ROOT, "resources/qualifying_studies_assembly_qc.tsv"))}
    leak = {r["assembly_accession"].split(".")[0] for r in c_rows} & qual
    check(not leak,
          "the catalogue block borrows nothing from the two qualifying studies")

    r2_ = subprocess.run([sys.executable, "-I",
                          os.path.join(SCRIPTS, "pilot_select.py"),
                          "--stage", "1",
                          "--catalogue", os.path.join(ROOT, "resources/ena_westafrica_ecoli_metadata.tsv"),
                          "--qc-table", os.path.join(ROOT, "resources/qualifying_studies_assembly_qc.tsv"),
                          "--out-block-a", os.path.join(tmp, "a2.tsv"),
                          "--out-candidates", os.path.join(tmp, "c2.tsv"),
                          "--out-manifest", os.path.join(tmp, "m2.tsv"),
                          "--out-samplesheet", os.path.join(tmp, "s2.tsv"),
                          "--out-log", os.path.join(tmp, "l2.json"),
                          "--out-summary", os.path.join(tmp, "y2.json")],
                         capture_output=True, text=True)
    same = (open(os.path.join(tmp, "a1.tsv")).read()
            == open(os.path.join(tmp, "a2.tsv")).read()
            and open(os.path.join(tmp, "c1.tsv")).read()
            == open(os.path.join(tmp, "c2.tsv")).read())
    check(r2_.returncode == 0 and same,
          "re-running the selection reproduces it exactly")

    print()
    if FAILS:
        print(f"{len(FAILS)} FAILED:")
        for f in FAILS:
            print(f"  - {f}")
        raise SystemExit(1)
    print("all checks passed")
    print(f"(artefacts under {tmp})")


if __name__ == "__main__":
    main()
