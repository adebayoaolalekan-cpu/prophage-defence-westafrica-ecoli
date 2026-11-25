"""
Collate Snakemake's benchmark files into the resource figures Section 11.0 asks
for: CPU, memory and disk requirements, runtime per genome and workflow
failure rates.

Snakemake writes one benchmark per job with wall seconds, CPU time, peak RSS,
peak VMS, peak USS, peak PSS, io_in, io_out and mean load. Those are read
rather than estimated, and a missing benchmark is counted as a failed job
rather than skipped, because the failure rate is one of the pilot's outputs.
"""

import glob
import json
import os
import sys
from collections import defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv, write_tsv, num

PER_JOB_FIELDS = ["tool", "sample", "wall_seconds", "cpu_seconds",
                  "max_rss_mb", "max_vms_mb", "io_in_mb", "io_out_mb",
                  "mean_load"]

PER_TOOL_FIELDS = ["tool", "scope", "n_expected", "n_completed", "n_failed",
                   "failure_rate", "wall_seconds_median",
                   "wall_seconds_mean", "wall_seconds_max",
                   "wall_seconds_total", "cpu_seconds_total",
                   "max_rss_mb_median", "max_rss_mb_max",
                   "output_bytes_total", "output_mb_per_genome",
                   "projected_hours_3243_genomes",
                   "projected_gb_3243_genomes"]

# Full retrieval after the NCBI completeness diff. Used only to project the
# production run from the pilot's observed per-genome cost.
PRODUCTION_N = 3243

# Rules that run once over the whole set rather than once per genome.
COLLECTION_SCOPE = {"checkm2", "poppunk"}

OUTPUT_DIRS = {
    "download": "genomes",
    "seqkit": "qc/seqkit",
    "checkm2": "qc/checkm2",
    "bakta": "annotation",
    "genomad": "genomad",
    "pide": "prophage/pide_raw",
    "checkv_primary": "checkv/primary",
    "checkv_sensitivity": "checkv/sensitivity",
    "defensefinder": "defence",
    "crisprcasfinder": "crispr",
    "poppunk": "poppunk",
}


def dir_bytes(path):
    total = 0
    for root, _, files in os.walk(path):
        for f in files:
            try:
                total += os.path.getsize(os.path.join(root, f))
            except OSError:
                pass
    return total


def median(v):
    v = sorted(x for x in v if x is not None)
    if not v:
        return ""
    n = len(v)
    return round(v[n // 2] if n % 2 else (v[n // 2 - 1] + v[n // 2]) / 2, 3)


def main(bench_dir, outdir, samples, out_jobs, out_tools):
    jobs = defaultdict(list)
    for path in sorted(glob.glob(os.path.join(bench_dir, "*", "*.tsv"))):
        tool = os.path.basename(os.path.dirname(path))
        sample = os.path.basename(path)[:-4]
        rows = read_tsv(path)
        if not rows:
            continue
        r = rows[0]
        jobs[tool].append({
            "tool": tool, "sample": sample,
            "wall_seconds": num(r.get("s")),
            "cpu_seconds": num(r.get("cpu_time")),
            "max_rss_mb": num(r.get("max_rss")),
            "max_vms_mb": num(r.get("max_vms")),
            "io_in_mb": num(r.get("io_in")),
            "io_out_mb": num(r.get("io_out")),
            "mean_load": num(r.get("mean_load")),
        })

    per_job = [j for t in sorted(jobs) for j in jobs[t]]
    write_tsv(out_jobs, per_job, PER_JOB_FIELDS)

    n_samples = len(samples)
    tool_rows = []
    for tool in sorted(set(list(jobs) + list(OUTPUT_DIRS))):
        js = jobs.get(tool, [])
        scope = "collection" if tool in COLLECTION_SCOPE else "per_genome"
        expected = 1 if scope == "collection" else n_samples
        completed = len(js)
        failed = max(0, expected - completed)

        walls = [j["wall_seconds"] for j in js]
        cpus = [j["cpu_seconds"] for j in js]
        rss = [j["max_rss_mb"] for j in js]
        wall_total = sum(w for w in walls if w is not None)

        od = OUTPUT_DIRS.get(tool)
        ob = dir_bytes(os.path.join(outdir, od)) if od else 0

        if scope == "per_genome" and completed:
            per_genome_wall = wall_total / completed
            per_genome_bytes = ob / completed
        elif completed:
            per_genome_wall = wall_total / max(1, n_samples)
            per_genome_bytes = ob / max(1, n_samples)
        else:
            per_genome_wall = per_genome_bytes = 0

        tool_rows.append({
            "tool": tool, "scope": scope,
            "n_expected": expected, "n_completed": completed,
            "n_failed": failed,
            "failure_rate": round(failed / expected, 4) if expected else "",
            "wall_seconds_median": median(walls),
            "wall_seconds_mean": round(wall_total / completed, 3) if completed else "",
            "wall_seconds_max": round(max([w for w in walls if w is not None], default=0), 3),
            "wall_seconds_total": round(wall_total, 3),
            "cpu_seconds_total": round(sum(c for c in cpus if c is not None), 3),
            "max_rss_mb_median": median(rss),
            "max_rss_mb_max": round(max([r for r in rss if r is not None], default=0), 2),
            "output_bytes_total": ob,
            "output_mb_per_genome": round(per_genome_bytes / 1e6, 3),
            "projected_hours_3243_genomes": round(
                per_genome_wall * PRODUCTION_N / 3600.0, 2),
            "projected_gb_3243_genomes": round(
                per_genome_bytes * PRODUCTION_N / 1e9, 2),
        })

    # Whole-workflow row.
    tool_rows.append({
        "tool": "ALL", "scope": "workflow",
        "n_expected": sum(r["n_expected"] for r in tool_rows),
        "n_completed": sum(r["n_completed"] for r in tool_rows),
        "n_failed": sum(r["n_failed"] for r in tool_rows),
        "failure_rate": round(
            sum(r["n_failed"] for r in tool_rows)
            / max(1, sum(r["n_expected"] for r in tool_rows)), 4),
        "wall_seconds_total": round(
            sum(r["wall_seconds_total"] or 0 for r in tool_rows), 3),
        "cpu_seconds_total": round(
            sum(r["cpu_seconds_total"] or 0 for r in tool_rows), 3),
        "output_bytes_total": sum(r["output_bytes_total"] for r in tool_rows),
        "output_mb_per_genome": round(
            sum(r["output_mb_per_genome"] or 0 for r in tool_rows), 3),
        "projected_hours_3243_genomes": round(
            sum(r["projected_hours_3243_genomes"] or 0 for r in tool_rows), 2),
        "projected_gb_3243_genomes": round(
            sum(r["projected_gb_3243_genomes"] or 0 for r in tool_rows), 2),
        "wall_seconds_median": "", "wall_seconds_mean": "",
        "wall_seconds_max": "", "max_rss_mb_median": "", "max_rss_mb_max": "",
    })

    write_tsv(out_tools, tool_rows, PER_TOOL_FIELDS)
    sys.stderr.write(
        f"resources: {len(per_job)} benchmark records across {len(jobs)} "
        f"tools; projected production cost "
        f"{tool_rows[-1]['projected_hours_3243_genomes']} core-hours of wall "
        f"time and {tool_rows[-1]['projected_gb_3243_genomes']} GB for "
        f"{PRODUCTION_N} genomes\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.params.bench, sm.params.outdir, list(sm.params.samples),
             sm.output.tsv, sm.output.per_tool)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3].split(","), sys.argv[4],
             sys.argv[5])
