"""
Extract integrated proviral predictions from the geNomad virus summary.

Protocol Step 5. Only predicted integrated proviral regions contribute to the
prophage-burden analysis. Predicted free viral sequences, plasmids and other
mobile elements are recorded separately rather than counted as prophages, so
the discarded rows are written to a companion file instead of being dropped
silently.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import contig_lengths, read_tsv, write_tsv, split_contig_id, num

FIELDS = ["sample", "prediction_id", "contig", "contig_length",
          "genomad_start", "genomad_end", "genomad_length", "topology",
          "virus_score", "fdr", "n_genes", "n_hallmarks",
          "marker_enrichment", "taxonomy", "boundary_class"]


def boundary_class(start, end, clen):
    """Step 3.4 classification, from the raw geNomad coordinates."""
    if clen is None or start is None:
        return "unknown"
    left = start <= 1
    right = end >= clen
    if left and right:
        return "truncated_both"
    if left:
        return "truncated_left"
    if right:
        return "truncated_right"
    return "internally_resolved"


def main(summary_path, fna_path, out_path, sample):
    lens = contig_lengths(fna_path)
    rows = read_tsv(summary_path)

    kept, dropped = [], []
    n = 0
    for r in rows:
        seq = r.get("seq_name", "")
        topo = (r.get("topology") or "").strip()
        contig, s, e = split_contig_id(seq)
        is_provirus = topo.lower() == "provirus" and s is not None

        if not is_provirus:
            dropped.append({
                "sample": sample, "seq_name": seq, "topology": topo,
                "length": r.get("length", ""),
                "virus_score": r.get("virus_score", ""),
                "taxonomy": r.get("taxonomy", ""),
                "reason": "not_an_integrated_provirus",
            })
            continue

        n += 1
        clen = lens.get(contig)
        kept.append({
            "sample": sample,
            "prediction_id": f"{sample}_gnm{n:04d}",
            "contig": contig,
            "contig_length": "" if clen is None else clen,
            "genomad_start": s,
            "genomad_end": e,
            "genomad_length": e - s + 1,
            "topology": topo,
            "virus_score": r.get("virus_score", ""),
            "fdr": r.get("fdr", ""),
            "n_genes": r.get("n_genes", ""),
            "n_hallmarks": r.get("n_hallmarks", ""),
            "marker_enrichment": r.get("marker_enrichment", ""),
            "taxonomy": r.get("taxonomy", ""),
            "boundary_class": boundary_class(s, e, clen),
        })

    write_tsv(out_path, kept, FIELDS)
    write_tsv(out_path.replace(".tsv", ".excluded.tsv"), dropped,
              ["sample", "seq_name", "topology", "length", "virus_score",
               "taxonomy", "reason"])
    sys.stderr.write(
        f"{sample}: {len(kept)} integrated proviruses, "
        f"{len(dropped)} non-provirus predictions recorded separately\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.input.summary, sm.input.fna, sm.output.tsv, sm.wildcards.sample)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])
