"""
Parse CRISPRCasFinder output into an array table and a spacer FASTA.

Protocol Step 11. For each array the following are recorded: array identifier,
contig identifier, coordinates, repeat sequence, number of spacers, individual
spacer sequences, orientation, predicted subtype, associated Cas operon and
Cas operon completeness.

Step 11.2 matters for E. coli specifically. The Type I-E CRISPR1 array lies
downstream of iap and CRISPR2 lies between ygcE and ygcF, the two are about
20 kb apart, and only CRISPR1 carries associated cas genes. CRISPR2 must
therefore not be called an orphan array merely because cas genes are not
immediately adjacent. The association window below is set at 30 kb for that
reason, and arrays beyond it are called orphaned only when no Cas locus exists
anywhere on the contig.
"""

import json
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import write_tsv

ARRAY_FIELDS = ["sample", "array_id", "contig", "contig_length", "start",
                "end", "length", "dr_consensus", "dr_length", "n_spacers",
                "orientation", "crispr_direction", "evidence_level",
                "conservation_drs", "conservation_spacers",
                "nearest_cas_type", "nearest_cas_distance",
                "cas_on_same_contig", "association_class"]

GENOME_FIELDS = ["sample", "n_arrays", "n_arrays_evidence_3_or_4",
                 "n_spacers_total", "n_cas_loci", "cas_types",
                 "n_orphan_arrays", "n_associated_arrays"]

ASSOCIATION_WINDOW = 30000


def main(samples, outdir, out_arrays, out_spacers):
    arrays, genome_rows, spacers = [], [], []

    for s in samples:
        path = f"{outdir}/crispr/{s}/result.json"
        if not os.path.exists(path):
            genome_rows.append({"sample": s, "n_arrays": 0,
                                "n_arrays_evidence_3_or_4": 0,
                                "n_spacers_total": 0, "n_cas_loci": 0,
                                "cas_types": "", "n_orphan_arrays": 0,
                                "n_associated_arrays": 0})
            continue
        with open(path) as fh:
            data = json.load(fh)

        n_arr = n_sp = n_cas = n_orph = n_assoc = n_strong = 0
        cas_types = set()
        a_idx = 0

        for seq in data.get("Sequences", []) or []:
            contig = seq.get("Id", "")
            clen = seq.get("Length", "")
            cas_loci = []
            for c in seq.get("Cas", []) or []:
                try:
                    cs, ce = int(c.get("Start", 0)), int(c.get("End", 0))
                except (TypeError, ValueError):
                    continue
                cas_loci.append((c.get("Type", ""), cs, ce,
                                 len(c.get("Genes", []) or [])))
                cas_types.add(str(c.get("Type", "")))
            n_cas += len(cas_loci)

            for cr in seq.get("Crisprs", []) or []:
                a_idx += 1
                n_arr += 1
                aid = cr.get("Name") or f"{s}_arr{a_idx:03d}"
                try:
                    st, en = int(cr.get("Start", 0)), int(cr.get("End", 0))
                except (TypeError, ValueError):
                    st = en = 0
                nsp = cr.get("Spacers", 0) or 0
                try:
                    nsp = int(nsp)
                except (TypeError, ValueError):
                    nsp = 0
                n_sp += nsp
                ev = cr.get("Evidence_Level", "")
                try:
                    if int(ev) >= 3:
                        n_strong += 1
                except (TypeError, ValueError):
                    pass

                nearest, dist = "", ""
                if cas_loci:
                    best = min(cas_loci,
                               key=lambda c: min(abs(st - c[2]), abs(c[1] - en))
                               if not (c[1] <= en and st <= c[2]) else 0)
                    overlapping = best[1] <= en and st <= best[2]
                    dist = 0 if overlapping else min(abs(st - best[2]),
                                                     abs(best[1] - en))
                    nearest = best[0]

                if not cas_loci:
                    cls = "orphan_no_cas_on_contig"
                    n_orph += 1
                elif dist != "" and dist <= ASSOCIATION_WINDOW:
                    cls = "cas_associated"
                    n_assoc += 1
                else:
                    cls = "cas_on_contig_beyond_window"
                    n_assoc += 1

                arrays.append({
                    "sample": s, "array_id": aid, "contig": contig,
                    "contig_length": clen, "start": st, "end": en,
                    "length": max(0, en - st + 1),
                    "dr_consensus": cr.get("DR_Consensus", ""),
                    "dr_length": cr.get("DR_Length", ""),
                    "n_spacers": nsp,
                    "orientation": cr.get("Potential_Orientation", ""),
                    "crispr_direction": cr.get("CRISPRDirection", ""),
                    "evidence_level": ev,
                    "conservation_drs": cr.get("Conservation_DRs", ""),
                    "conservation_spacers": cr.get("Conservation_Spacers", ""),
                    "nearest_cas_type": nearest,
                    "nearest_cas_distance": dist,
                    "cas_on_same_contig": "YES" if cas_loci else "NO",
                    "association_class": cls,
                })

                k = 0
                for reg in cr.get("Regions", []) or []:
                    if str(reg.get("Type", "")).lower() != "spacer":
                        continue
                    k += 1
                    seqstr = reg.get("Sequence", "") or ""
                    if seqstr:
                        spacers.append(
                            (f"{s}|{aid}|spacer{k:03d}", seqstr))

        genome_rows.append({
            "sample": s, "n_arrays": n_arr,
            "n_arrays_evidence_3_or_4": n_strong,
            "n_spacers_total": n_sp, "n_cas_loci": n_cas,
            "cas_types": ";".join(sorted(x for x in cas_types if x)),
            "n_orphan_arrays": n_orph, "n_associated_arrays": n_assoc,
        })

    write_tsv(out_arrays, arrays, ARRAY_FIELDS)
    write_tsv(out_arrays.replace("crispr_arrays.tsv", "crispr_by_genome.tsv"),
              genome_rows, GENOME_FIELDS)
    os.makedirs(os.path.dirname(out_spacers) or ".", exist_ok=True)
    with open(out_spacers, "w") as fh:
        for name, seqstr in spacers:
            fh.write(f">{name}\n{seqstr}\n")

    sys.stderr.write(
        f"crispr: {len(arrays)} arrays, {len(spacers)} spacers across "
        f"{len(genome_rows)} genomes\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(list(sm.params.samples), sm.params.outdir, sm.output.arrays,
             sm.output.spacers)
    else:
        main(sys.argv[1].split(","), sys.argv[2], sys.argv[3], sys.argv[4])
