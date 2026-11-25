"""
Cut each geNomad prediction out of its host contig with the protocol's flank.

Protocol Step 7 and Table 10. CheckV assesses single-contig viral genomes, so
supplying a whole bacterial contig is wrong by the tool's design and supplying
an already-trimmed geNomad provirus is equally wrong, because the breakpoint
rules need host-side genes on both sides of the junction. The extraction is
therefore the prediction plus a flank.

    S = max(1, s - L)
    E = min(C, e + L)

where s and e are the geNomad coordinates, C is the contig length and L is the
flank. Where an extraction reaches the start or end of a contig the missing
flank is recorded, because a one-sided extraction gives CheckV less to work
with and that has to be visible in the output rather than inferred.

Where two predictions on one contig have overlapping flanked windows they are
processed as one shared region, and the map file records both predictions
against that region so the refinement step can split the result back out.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_fasta, read_tsv, write_tsv, wrap

MAP_FIELDS = ["region_id", "sample", "contig", "contig_length",
              "region_start", "region_end", "region_length", "flank",
              "flank_left_actual", "flank_right_actual",
              "flank_left_clipped", "flank_right_clipped",
              "prediction_id", "genomad_start", "genomad_end",
              "offset_start_in_region", "offset_end_in_region",
              "shared_region", "n_predictions_in_region"]


def build_windows(preds, clen, flank, merge):
    """Return [(start, end, [pred, ...])] in 1-based inclusive coordinates."""
    windows = []
    for p in preds:
        s, e = int(p["genomad_start"]), int(p["genomad_end"])
        S = max(1, s - flank)
        E = min(clen, e + flank)
        windows.append([S, E, [p]])
    if not merge or len(windows) < 2:
        return sorted(windows, key=lambda w: (w[0], w[1]))

    windows.sort(key=lambda w: (w[0], w[1]))
    merged = [windows[0]]
    for S, E, ps in windows[1:]:
        last = merged[-1]
        if S <= last[1]:                 # inclusive coordinates, so <= not <
            last[1] = max(last[1], E)
            last[2].extend(ps)
        else:
            merged.append([S, E, ps])
    return merged


def main(fna_path, prov_path, out_fasta, out_map, flank, merge, sample):
    flank = int(flank)
    merge = str(merge).strip().lower() in ("1", "true", "yes")

    preds = read_tsv(prov_path)
    seqs = dict(read_fasta(fna_path))

    by_contig = {}
    for p in preds:
        by_contig.setdefault(p["contig"], []).append(p)

    records, maps = [], []
    rid = 0
    for contig in sorted(by_contig):
        seq = seqs.get(contig)
        if seq is None:
            sys.stderr.write(
                f"WARNING {sample}: contig {contig} carries a prediction but "
                f"is absent from the assembly FASTA. Skipped.\n")
            continue
        clen = len(seq)
        for S, E, ps in build_windows(by_contig[contig], clen, flank, merge):
            rid += 1
            region_id = f"{sample}_reg{rid:04d}"
            records.append((region_id, seq[S - 1:E]))
            shared = len(ps) > 1
            for p in ps:
                s, e = int(p["genomad_start"]), int(p["genomad_end"])
                maps.append({
                    "region_id": region_id,
                    "sample": sample,
                    "contig": contig,
                    "contig_length": clen,
                    "region_start": S,
                    "region_end": E,
                    "region_length": E - S + 1,
                    "flank": flank,
                    "flank_left_actual": s - S,
                    "flank_right_actual": E - e,
                    "flank_left_clipped": "YES" if (s - S) < flank else "NO",
                    "flank_right_clipped": "YES" if (E - e) < flank else "NO",
                    "prediction_id": p["prediction_id"],
                    "genomad_start": s,
                    "genomad_end": e,
                    # Offsets inside the extracted region, 1-based.
                    "offset_start_in_region": s - S + 1,
                    "offset_end_in_region": e - S + 1,
                    "shared_region": "YES" if shared else "NO",
                    "n_predictions_in_region": len(ps),
                })

    os.makedirs(os.path.dirname(out_fasta) or ".", exist_ok=True)
    with open(out_fasta, "w") as fh:
        for name, s in records:
            fh.write(f">{name}\n{wrap(s)}\n")
    write_tsv(out_map, maps, MAP_FIELDS)

    n_clipped = sum(1 for m in maps
                    if m["flank_left_clipped"] == "YES"
                    or m["flank_right_clipped"] == "YES")
    sys.stderr.write(
        f"{sample} flank={flank}: {len(records)} regions from {len(preds)} "
        f"predictions, {sum(1 for m in maps if m['shared_region'] == 'YES')} "
        f"predictions in shared regions, {n_clipped} with a clipped flank\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.input.fna, sm.input.prov, sm.output.fasta, sm.output.map,
             sm.params.flank, sm.params.merge, sm.wildcards.sample)
    else:
        main(*sys.argv[1:8])
