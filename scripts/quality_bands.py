"""
Definitive quality-band assignment, Table 9 of the protocol.

Selection assigned a provisional band from assembly statistics alone, because
completeness and contamination are not available until CheckM2 has run. This
applies the full definition and reports every genome whose band changed.

A genome that drops out of its selected band is reported, not replaced.
Replacing it after the fact would make the pilot a quality-selected sample of a
quality-selected sample, and the resulting dispersion estimate would understate
the variance the production run will meet.
"""

import os
import sys

import yaml

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv, write_tsv, num

FIELDS = ["sample", "block", "ecological_category", "originating_study",
          "country", "assembly_level", "contig_n50", "number_of_contigs",
          "total_length", "gc_percent", "completeness", "contamination",
          "passes_broad", "passes_high_contiguity", "quality_band",
          "selected_band", "band_changed", "exclusion_reason"]


def band_of(n50, ctg, comp, cont, q):
    """Return (band, broad_pass, highc_pass, reason)."""
    missing = [k for k, v in (("contig_n50", n50), ("number_of_contigs", ctg),
                              ("completeness", comp),
                              ("contamination", cont)) if v is None]
    if missing:
        return "unknown", "NO", "NO", "missing:" + ",".join(missing)

    b, h = q["broad"], q["high_contiguity"]
    reasons = []
    broad = True
    if n50 < b["min_contig_n50"]:
        broad = False
        reasons.append(f"n50<{b['min_contig_n50']}")
    if ctg > b["max_contigs"]:
        broad = False
        reasons.append(f"contigs>{b['max_contigs']}")
    if comp < b["min_completeness"]:
        broad = False
        reasons.append(f"completeness<{b['min_completeness']}")
    if cont > b["max_contamination"]:
        broad = False
        reasons.append(f"contamination>{b['max_contamination']}")

    highc = (n50 > h["min_contig_n50"] and ctg <= h["max_contigs"]
             and comp >= h["min_completeness"] and cont <= h["max_contamination"])

    if highc:
        band = "high_contiguity"
    elif broad:
        band = "broad_only"
    else:
        band = "below_broad"
    return (band, "YES" if broad else "NO", "YES" if highc else "NO",
            ";".join(reasons))


def parse_seqkit(paths):
    out = {}
    for p in paths:
        for r in read_tsv(p):
            sample = os.path.basename(p)[:-4]
            out[sample] = {
                "contig_n50": num(r.get("N50")),
                "number_of_contigs": num(r.get("num_seqs")),
                "total_length": num(r.get("sum_len")),
                "gc_percent": num(r.get("GC(%)")),
            }
    return out


def parse_checkm2(path):
    out = {}
    for r in read_tsv(path):
        out[r.get("Name") or r.get("Bin Id") or ""] = {
            "completeness": num(r.get("Completeness")),
            "contamination": num(r.get("Contamination")),
        }
    return out


def main(seqkit_paths, checkm2_path, manifest_path, cfg_path,
         out_tsv, out_changed):
    with open(cfg_path) as fh:
        q = yaml.safe_load(fh)["quality"]

    stats = parse_seqkit(seqkit_paths)
    cm = parse_checkm2(checkm2_path)
    man = {r["assembly_accession"]: r for r in read_tsv(manifest_path)}

    rows, changed = [], []
    for sample in sorted(stats):
        s = stats[sample]
        c = cm.get(sample, {})
        m = man.get(sample, {})
        comp, cont = c.get("completeness"), c.get("contamination")
        band, broad, highc, reason = band_of(
            s["contig_n50"], s["number_of_contigs"], comp, cont, q)
        selected = m.get("quality_band", "")
        rec = {
            "sample": sample,
            "block": m.get("block", ""),
            "ecological_category": m.get("ecological_category", ""),
            "originating_study": m.get("originating_study", ""),
            "country": m.get("country", ""),
            "assembly_level": m.get("assembly_level", ""),
            "contig_n50": s["contig_n50"],
            "number_of_contigs": s["number_of_contigs"],
            "total_length": s["total_length"],
            "gc_percent": s["gc_percent"],
            "completeness": "" if comp is None else comp,
            "contamination": "" if cont is None else cont,
            "passes_broad": broad,
            "passes_high_contiguity": highc,
            "quality_band": band,
            "selected_band": selected,
            "band_changed": "YES" if (selected and band != selected) else "NO",
            "exclusion_reason": reason,
        }
        rows.append(rec)
        if rec["band_changed"] == "YES":
            changed.append(rec)

    write_tsv(out_tsv, rows, FIELDS)
    write_tsv(out_changed, changed, FIELDS)

    if not cm:
        sys.stderr.write("WARNING CheckM2 report parsed to zero genomes. Band "
                         "assignment is on assembly statistics only.\n")
    sys.stderr.write(
        f"quality bands: {len(rows)} genomes, {len(changed)} changed band "
        f"from selection\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(list(sm.input.seqkit), sm.input.checkm2, sm.input.manifest,
             sm.params.cfg, sm.output.tsv, sm.output.changed)
    else:
        main(sys.argv[1].split(","), sys.argv[2], sys.argv[3], sys.argv[4],
             sys.argv[5], sys.argv[6])
