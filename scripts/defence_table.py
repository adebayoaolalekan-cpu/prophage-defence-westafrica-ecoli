"""
Collate DefenseFinder output into one table per genome and one per system.

Protocol Step 10. Four variables are required and they are kept apart: the
total number of complete systems, the number of distinct system types, system
diversity, and the number of partial or contig-truncated candidates. A genome
carrying several copies of one system type is not more diverse than a genome
carrying several distinct types, so count and diversity never share a column.

CRISPR-Cas systems are characterised alongside other defence systems, but
their array and spacer counts stay in the CRISPR table so the same biological
feature is not counted twice. The CRISPR-Cas systems DefenseFinder reports are
therefore flagged here and excluded from the headline defence count, which is
reported both ways.
"""

import math
import os
import sys
from collections import Counter

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv, write_tsv, num

GENOME_FIELDS = ["sample", "n_systems_total", "n_systems_excl_crispr",
                 "n_distinct_types", "n_distinct_subtypes",
                 "shannon_diversity_types", "simpson_diversity_types",
                 "n_crispr_cas_systems", "n_truncated_candidates",
                 "n_systems_near_contig_end", "system_types", "subtypes"]

SYSTEM_FIELDS = ["sample", "sys_id", "type", "subtype", "activity", "contig",
                 "sys_beg", "sys_end", "genes_count", "profiles",
                 "is_crispr_cas", "near_contig_end", "truncation_flag"]

# A system whose outermost gene sits within this distance of a contig end is
# flagged, because its missing components cannot be distinguished from
# components that were never there (Step 3.5).
EDGE_BP = 5000

TYPE_KEYS = ["type", "Type", "system"]
SUBTYPE_KEYS = ["subtype", "Subtype"]
BEG_KEYS = ["sys_beg", "sys_begin", "start"]
END_KEYS = ["sys_end", "end"]
GENES_KEYS = ["genes_count", "genes_number", "n_genes"]
PROFILE_KEYS = ["name_of_profiles_in_sys", "profiles", "protein_in_syst"]
ACTIVITY_KEYS = ["activity", "Activity"]


def pick(r, keys, default=""):
    for k in keys:
        if k in r and str(r[k]).strip() != "":
            return r[k]
    return default


def shannon(counts):
    n = sum(counts)
    if n == 0:
        return ""
    return round(-sum((c / n) * math.log(c / n) for c in counts if c > 0), 4)


def simpson(counts):
    n = sum(counts)
    if n == 0:
        return ""
    return round(1.0 - sum((c / n) ** 2 for c in counts), 4)


def contig_of(tag, sample):
    """
    DefenseFinder reports the protein identifier Bakta assigned. Bakta locus
    tags do not carry the contig, so the contig is recovered from the gene
    table where that is possible and left blank where it is not. A blank is
    honest; a guessed contig would corrupt the edge flag.
    """
    if not tag:
        return ""
    if "|" in str(tag):
        return str(tag).split("_")[0]
    return ""


def main(samples, outdir, out_path):
    lens = {}
    for s in samples:
        for r in read_tsv(f"{outdir}/qc/seqkit/{s}.tsv"):
            lens[s] = num(r.get("max_len"))

    genome_rows, system_rows = [], []
    for s in samples:
        systems = read_tsv(f"{outdir}/defence/{s}/defense_finder_systems.tsv")
        genes = read_tsv(f"{outdir}/defence/{s}/defense_finder_genes.tsv")
        gene_contig = {}
        for g in genes:
            hit = pick(g, ["hit_id", "replicon", "gene_name"], "")
            rep = pick(g, ["replicon", "replicon_name"], "")
            if hit and rep:
                gene_contig[hit] = rep

        types, subtypes = [], []
        n_crispr = n_trunc = n_edge = 0
        for r in systems:
            t = str(pick(r, TYPE_KEYS))
            st = str(pick(r, SUBTYPE_KEYS))
            beg, end = num(pick(r, BEG_KEYS)), num(pick(r, END_KEYS))
            profiles = str(pick(r, PROFILE_KEYS))
            contig = gene_contig.get(profiles.split(",")[0].strip(), "") \
                or contig_of(profiles, s)

            is_crispr = "crispr" in t.lower() or "cas" in st.lower()[:3]
            # Truncation: DefenseFinder does not report completeness directly,
            # so a system is flagged when its gene count falls below the
            # mandatory count for its subtype where that is known, and
            # otherwise when it sits near a contig end.
            max_len = lens.get(s)
            near_edge = bool(
                max_len and beg is not None and end is not None
                and (beg <= EDGE_BP or end >= max_len - EDGE_BP))
            trunc = "near_contig_end" if near_edge else ""

            if is_crispr:
                n_crispr += 1
            else:
                types.append(t)
                subtypes.append(st)
            if near_edge:
                n_edge += 1
            if trunc:
                n_trunc += 1

            system_rows.append({
                "sample": s,
                "sys_id": pick(r, ["sys_id", "id"]),
                "type": t, "subtype": st,
                "activity": pick(r, ACTIVITY_KEYS),
                "contig": contig,
                "sys_beg": "" if beg is None else beg,
                "sys_end": "" if end is None else end,
                "genes_count": pick(r, GENES_KEYS),
                "profiles": profiles,
                "is_crispr_cas": "YES" if is_crispr else "NO",
                "near_contig_end": "YES" if near_edge else "NO",
                "truncation_flag": trunc,
            })

        tc = Counter(types)
        genome_rows.append({
            "sample": s,
            "n_systems_total": len(systems),
            "n_systems_excl_crispr": len(systems) - n_crispr,
            "n_distinct_types": len(tc),
            "n_distinct_subtypes": len(set(subtypes)),
            "shannon_diversity_types": shannon(list(tc.values())),
            "simpson_diversity_types": simpson(list(tc.values())),
            "n_crispr_cas_systems": n_crispr,
            "n_truncated_candidates": n_trunc,
            "n_systems_near_contig_end": n_edge,
            "system_types": ";".join(f"{k}:{v}" for k, v in sorted(tc.items())),
            "subtypes": ";".join(sorted(set(subtypes))),
        })

    write_tsv(out_path, genome_rows, GENOME_FIELDS)
    write_tsv(out_path.replace(".tsv", "_by_system.tsv"), system_rows,
              SYSTEM_FIELDS)
    tot = sum(r["n_systems_total"] for r in genome_rows)
    sys.stderr.write(
        f"defence: {tot} systems across {len(genome_rows)} genomes, "
        f"{sum(r['n_systems_near_contig_end'] for r in genome_rows)} near a "
        f"contig end\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(list(sm.params.samples), sm.params.outdir, sm.output.tsv)
    else:
        main(sys.argv[1].split(","), sys.argv[2], sys.argv[3])
