"""
Write the pilot report: the four things the pilot exists to establish, with the
numbers taken from the generated tables rather than restated by hand.
"""

import json
import os
import statistics
import sys
from collections import Counter

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv


def fmt(x, nd=3):
    if x is None or x == "":
        return "not available"
    try:
        f = float(x)
    except (TypeError, ValueError):
        return str(x)
    return f"{f:,.0f}" if f.is_integer() else f"{f:,.{nd}f}"


def tab(rows, cols, headers=None):
    headers = headers or cols
    out = ["| " + " | ".join(headers) + " |",
           "| " + " | ".join("---" for _ in headers) + " |"]
    for r in rows:
        out.append("| " + " | ".join(str(r.get(c, "")) for c in cols) + " |")
    return "\n".join(out)


def main(outdir, n_requested, out_md):
    counts = pd.read_csv(f"{outdir}/reports/prophage_counts.tsv", sep="\t",
                         dtype=str)
    counts["n"] = pd.to_numeric(counts["n_prophage_primary"], errors="coerce")
    disp = json.load(open(f"{outdir}/reports/dispersion_estimates.json"))
    per_tool = read_tsv(f"{outdir}/reports/resource_per_tool.tsv")
    conc = read_tsv(f"{outdir}/reports/concordance.tsv")
    defence = read_tsv(f"{outdir}/reports/defence_systems.tsv")
    crispr = read_tsv(f"{outdir}/reports/crispr_arrays.tsv")
    bands = read_tsv(f"{outdir}/qc/quality_bands.tsv")
    changed = read_tsv(f"{outdir}/qc/band_changes.tsv")
    prov = read_tsv(f"{outdir}/reports/provenance.tsv")

    all_row = next((r for r in per_tool if r["tool"] == "ALL"), {})
    conc_all = next((r for r in conc if r["sample"] == "ALL"), {})
    inf = next((r for r in disp["estimates"]
                if r["set"] == "inference_block_high_contiguity"
                and r.get("measure") == "n_prophage_primary"), {})
    verdict = disp.get("verdict", {})

    n_done = len(counts)
    y = counts["n"].dropna()
    band_counts = Counter(r["quality_band"] for r in bands)

    L = []
    A = L.append

    A("# Genomic pilot report")
    A("")
    A("Prophage Burden, Prophage-Associated Virulence Determinants and "
      "Anti-Phage Defence Systems in West African *Escherichia coli*")
    A("")
    A(f"Genomes requested: {n_requested}. Genomes completing the workflow: "
      f"{n_done}. Workflow failure rate across all jobs: "
      f"{fmt(all_row.get('failure_rate'), 4)}.")
    A("")

    # --- 1. dispersion -----------------------------------------------------
    A("## The dispersion estimate Section 10.1 was waiting on")
    A("")
    if not inf or inf.get("n", 0) == 0:
        A("The inference block returned no high-contiguity genomes, so the "
          "dispersion parameter could not be estimated and Section 10.1 "
          "cannot be re-run. This is a workflow failure rather than a "
          "finding, and the cause is in the quality band table.")
    else:
        A(f"The primary inference set holds {inf['n']} genomes. The observed "
          f"prophage count has a mean of {fmt(inf.get('mean_count'))} and a "
          f"variance of {fmt(inf.get('variance'))}, giving a variance to mean "
          f"ratio of {fmt(inf.get('variance_to_mean_ratio'))}. "
          f"{inf.get('zero_count_genomes', 0)} genomes carried no internally "
          f"resolved prophage and the maximum count was "
          f"{fmt(inf.get('max_count'))}.")
        A("")
        A(f"Dispersion was estimated three ways. Method of moments returned "
          f"theta of {fmt(inf.get('theta_method_of_moments'))}. The "
          f"intercept-only maximum likelihood fit returned "
          f"{fmt(inf.get('theta_mle_intercept_only'))}. The covariate-adjusted "
          f"fit, which is the one the re-simulation should use, returned "
          f"{fmt(inf.get('theta_mle_covariate_adjusted'))}, with a lineage-"
          f"clustered bootstrap interval of "
          f"{fmt(inf.get('theta_bootstrap_ci_low'))} to "
          f"{fmt(inf.get('theta_bootstrap_ci_high'))}.")
        A("")
        A(f"**Verdict: {verdict.get('status', 'unknown')}.** "
          f"{verdict.get('note', '')}")
        if inf.get("covariates_unavailable"):
            A("")
            A("Covariates the pilot could not supply to the adjusted fit: "
              + ", ".join(inf["covariates_unavailable"])
              + ". The adjusted theta is conditional on fewer covariates than "
                "the analysis model, which makes it conservative rather than "
                "optimistic.")
    A("")
    A("Burden distribution across the whole pilot:")
    A("")
    if len(y):
        A(tab([{"stat": k, "value": v} for k, v in [
            ("n", len(y)), ("mean", fmt(y.mean())),
            ("median", fmt(y.median())), ("variance", fmt(y.var(ddof=1))),
            ("minimum", fmt(y.min())), ("maximum", fmt(y.max())),
            ("genomes with zero", int((y == 0).sum())),
        ]], ["stat", "value"], ["Statistic", "Value"]))
    A("")

    # --- 2. lineage --------------------------------------------------------
    A("## Lineage partition")
    A("")
    lin = Counter(counts["lineage"].fillna("").tolist())
    lin.pop("", None)
    if lin:
        sizes = sorted(lin.values(), reverse=True)
        A(f"PopPUNK assigned {n_done} genomes to {len(lin)} lineages. The "
          f"largest holds {sizes[0]} genomes, the median lineage holds "
          f"{fmt(statistics.median(sizes))}, and {sum(1 for s in sizes if s == 1)} "
          f"lineages hold a single genome. The clustered standard errors in "
          f"the re-simulation should use this partition rather than the "
          f"generic one the feasibility screen assumed.")
    else:
        A("No lineage assignment was produced. The re-simulation would fall "
          "back to the generic partition, which is the weakness the pilot was "
          "meant to remove.")
    A("")

    # --- 3. resources ------------------------------------------------------
    A("## Computational cost")
    A("")
    A(f"Total wall time across the pilot: "
      f"{fmt(float(all_row.get('wall_seconds_total') or 0) / 3600, 2)} hours. "
      f"Total CPU time: "
      f"{fmt(float(all_row.get('cpu_seconds_total') or 0) / 3600, 2)} "
      f"core-hours. Output footprint: "
      f"{fmt(float(all_row.get('output_bytes_total') or 0) / 1e9, 2)} GB, "
      f"{fmt(all_row.get('output_mb_per_genome'))} MB per genome.")
    A("")
    A(f"Projected to the full retrieval of 3,243 genomes: "
      f"{fmt(all_row.get('projected_hours_3243_genomes'), 1)} hours of wall "
      f"time at the pilot's thread allocation and "
      f"{fmt(all_row.get('projected_gb_3243_genomes'), 1)} GB of output. The "
      f"projection is linear in genome count, which is reasonable for the "
      f"per-genome tools and optimistic for CheckM2 and PopPUNK, both of "
      f"which scale worse than linearly.")
    A("")
    A(tab([r for r in per_tool if r["tool"] != "ALL"],
          ["tool", "scope", "n_completed", "n_failed", "wall_seconds_median",
           "max_rss_mb_max", "projected_hours_3243_genomes",
           "projected_gb_3243_genomes"],
          ["Tool", "Scope", "Done", "Failed", "Median wall (s)",
           "Peak RSS (MB)", "Projected hours", "Projected GB"]))
    A("")

    # --- 4. concordance ----------------------------------------------------
    tool_name = next((r["value"] for r in prov
                      if r["item"] == "independent_predictor_used"), "unknown")
    A("## Boundary concordance and the independent predictor")
    A("")
    if conc_all:
        n_g = float(conc_all.get("n_genomad") or 0)
        n_h = float(conc_all.get("n_high_concordance") or 0)
        A(f"The independent predictor was {tool_name}. geNomad made "
          f"{fmt(conc_all.get('n_genomad'))} integrated-provirus predictions "
          f"and {tool_name} made {fmt(conc_all.get('n_independent'))}. "
          f"{fmt(conc_all.get('n_matched'))} overlapped, of which "
          f"{fmt(conc_all.get('n_high_concordance'))} met the prespecified "
          f"definition of high boundary concordance. That is "
          f"{fmt(100 * n_h / n_g if n_g else 0, 1)}% of geNomad's predictions.")
        A("")
        A(f"Median reciprocal overlap among matched pairs: "
          f"{fmt(conc_all.get('median_reciprocal_overlap'))}. Median start "
          f"difference: {fmt(conc_all.get('median_start_difference'))} bp. "
          f"Median end difference: "
          f"{fmt(conc_all.get('median_end_difference'))} bp. "
          f"{fmt(conc_all.get('n_genomad_only'))} geNomad predictions had no "
          f"counterpart and {fmt(conc_all.get('n_independent_only'))} "
          f"predictions were specific to {tool_name}.")
    A("")
    pide_row = next((r for r in per_tool if r["tool"] == "pide"), {})
    if pide_row:
        A(f"Practicality at collection scale. The independent predictor took "
          f"a median of {fmt(pide_row.get('wall_seconds_median'))} seconds "
          f"per genome with a peak resident set of "
          f"{fmt(pide_row.get('max_rss_mb_max'))} MB, projecting to "
          f"{fmt(pide_row.get('projected_hours_3243_genomes'), 1)} hours over "
          f"the full retrieval. Section 11.0 substitutes PhiSpy if that is "
          f"impractical, and the decision belongs with these numbers rather "
          f"than with an impression of how long it felt.")
    A("")

    # --- quality and truncation -------------------------------------------
    A("## Assembly quality and boundary truncation")
    A("")
    A("Quality bands as finally assigned, after CheckM2:")
    A("")
    A(tab([{"band": k, "n": v} for k, v in sorted(band_counts.items())],
          ["band", "n"], ["Band", "Genomes"]))
    A("")
    if changed:
        A(f"{len(changed)} genomes changed band between selection and the "
          f"definitive assignment. They are retained and reported rather than "
          f"replaced, because replacing them would make the pilot a "
          f"quality-selected sample of a quality-selected sample and would "
          f"bias the dispersion estimate downward.")
    else:
        A("No genome changed band between selection and the definitive "
          "assignment.")
    A("")
    tl = pd.to_numeric(counts["n_truncated_left"], errors="coerce").sum()
    tr = pd.to_numeric(counts["n_truncated_right"], errors="coerce").sum()
    tb = pd.to_numeric(counts["n_truncated_both"], errors="coerce").sum()
    allp = pd.to_numeric(counts["n_prophage_all_predictions"],
                         errors="coerce").sum()
    nb = pd.to_numeric(counts["n_no_breakpoint"], errors="coerce").sum()
    ns = pd.to_numeric(counts["n_unstable_boundary"], errors="coerce").sum()
    A(f"Of {fmt(allp)} predictions, {fmt(tl + tr + tb)} terminated at a "
      f"contig boundary and are excluded from the primary burden measure: "
      f"{fmt(tl)} truncated on the left, {fmt(tr)} on the right and "
      f"{fmt(tb)} at both ends. CheckV found no host junction for "
      f"{fmt(nb)} predictions despite the 10 kb flank, and those retain their "
      f"geNomad coordinates. {fmt(ns)} predictions moved by more than the "
      f"stability tolerance between the 10 kb and 5 kb extractions.")
    A("")

    # --- defence and crispr ------------------------------------------------
    A("## Defence systems and CRISPR arrays")
    A("")
    if defence:
        tot = sum(int(r["n_systems_total"] or 0) for r in defence)
        edge = sum(int(r["n_systems_near_contig_end"] or 0) for r in defence)
        types = sum(int(r["n_distinct_types"] or 0) for r in defence)
        A(f"DefenseFinder recovered {tot} systems across {len(defence)} "
          f"genomes, a mean of {fmt(tot / len(defence))} per genome, with a "
          f"mean of {fmt(types / len(defence))} distinct types. {edge} systems "
          f"sat within 5 kb of a contig end and are flagged, since a missing "
          f"component cannot be told apart from a component that was never "
          f"sequenced.")
    if crispr:
        arrays = len(crispr)
        sp = sum(int(r["n_spacers"] or 0) for r in crispr)
        orph = sum(1 for r in crispr
                   if r["association_class"] == "orphan_no_cas_on_contig")
        A("")
        A(f"CRISPRCasFinder returned {arrays} arrays carrying {sp} spacers. "
          f"{orph} arrays sat on a contig with no Cas locus anywhere on it. "
          f"Arrays with a Cas locus on the contig but beyond the 30 kb "
          f"association window are reported separately rather than called "
          f"orphans, because in *E. coli* CRISPR2 lies about 20 kb from "
          f"CRISPR1 and carries no adjacent *cas* genes.")
    A("")

    # --- what happens next -------------------------------------------------
    A("## What this changes in the protocol")
    A("")
    A("1. Section 10.1 is re-run with the theta above and the PopPUNK "
      "partition in `poppunk/lineage_assignments.tsv`, replacing the fixed "
      "dispersion and generic partition of the feasibility screen.")
    A("2. The minimum detectable effect for the human and poultry contrast is "
      "restated from that re-run. Condition 4 of Section 8.0 is then either "
      "satisfied at 1.5 or left at 2.0.")
    A("3. The independent-predictor decision is recorded, with the runtime "
      "figures as its justification.")
    A("4. Any parameter the pilot shows to need adjustment goes in the "
      "amendment log before the production run, not after it.")
    A("")
    A(f"Provenance for this run, including every parameter and database "
      f"release, is in `reports/provenance.tsv` ({len(prov)} entries).")
    A("")

    os.makedirs(os.path.dirname(out_md) or ".", exist_ok=True)
    with open(out_md, "w") as fh:
        fh.write("\n".join(L) + "\n")
    sys.stderr.write(f"wrote {out_md}\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.params.outdir, sm.params.n_requested, sm.output.md)
    else:
        main(sys.argv[1], int(sys.argv[2]), sys.argv[3])
