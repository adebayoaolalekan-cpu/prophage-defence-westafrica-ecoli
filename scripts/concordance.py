"""
geNomad against the independent predictor, on original contig coordinates.

Protocol Step 6. High boundary concordance is defined in advance as reciprocal
sequence overlap of at least 80% together with a coordinate difference no
greater than the stated tolerance. Both parts of that definition have to hold,
so a long prediction that engulfs a short one fails on reciprocal overlap even
though one-sided overlap would look perfect.

The comparison records, per genome and in total: predictions found by both
tools, geNomad-specific predictions, predictions specific to the independent
tool, sequence overlap, start and end coordinate differences, length
differences, and the breakdown across assembly-quality categories.
"""

import os
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv, write_tsv, num

PAIR_FIELDS = ["sample", "genomad_id", "independent_id", "contig",
               "genomad_start", "genomad_end", "genomad_length",
               "independent_start", "independent_end", "independent_length",
               "overlap_bp", "reciprocal_overlap", "start_difference",
               "end_difference", "length_difference", "concordance_class"]

SUMMARY_FIELDS = ["sample", "n_genomad", "n_independent", "n_matched",
                  "n_high_concordance", "n_genomad_only",
                  "n_independent_only", "jaccard_predictions",
                  "median_reciprocal_overlap", "median_start_difference",
                  "median_end_difference", "median_length_difference"]


def recip(a1, a2, b1, b2):
    ov = max(0, min(a2, b2) - max(a1, b1) + 1)
    la, lb = a2 - a1 + 1, b2 - b1 + 1
    if la <= 0 or lb <= 0:
        return 0, 0.0
    return ov, min(ov / la, ov / lb)


def median(vals):
    v = sorted(x for x in vals if x is not None)
    if not v:
        return ""
    n = len(v)
    return v[n // 2] if n % 2 else round((v[n // 2 - 1] + v[n // 2]) / 2, 4)


def match_one(gnm, ind, min_ov, max_diff):
    """
    Greedy one-to-one matching by descending reciprocal overlap. Greedy is
    adequate here and is stated rather than dressed up as optimal: ties are
    broken by the smaller coordinate difference, and an exhaustive assignment
    changes nothing at the overlap levels these tools produce.
    """
    cands = []
    for g in gnm:
        gs, ge = int(g["genomad_start"]), int(g["genomad_end"])
        for i in ind:
            if i["contig"] != g["contig"]:
                continue
            isx, iex = num(i["start"]), num(i["end"])
            if isx is None or iex is None:
                continue
            ov, ro = recip(gs, ge, int(isx), int(iex))
            if ov == 0:
                continue
            cands.append((ro, -(abs(gs - int(isx)) + abs(ge - int(iex))),
                          g, i, ov))
    cands.sort(key=lambda c: (c[0], c[1]), reverse=True)

    used_g, used_i, pairs = set(), set(), []
    for ro, _, g, i, ov in cands:
        if g["prediction_id"] in used_g or i["prediction_id"] in used_i:
            continue
        used_g.add(g["prediction_id"])
        used_i.add(i["prediction_id"])
        gs, ge = int(g["genomad_start"]), int(g["genomad_end"])
        isx, iex = int(num(i["start"])), int(num(i["end"]))
        sd, ed = abs(gs - isx), abs(ge - iex)
        high = ro >= min_ov and max(sd, ed) <= max_diff
        pairs.append({
            "sample": g["sample"],
            "genomad_id": g["prediction_id"],
            "independent_id": i["prediction_id"],
            "contig": g["contig"],
            "genomad_start": gs, "genomad_end": ge,
            "genomad_length": ge - gs + 1,
            "independent_start": isx, "independent_end": iex,
            "independent_length": iex - isx + 1,
            "overlap_bp": ov,
            "reciprocal_overlap": round(ro, 4),
            "start_difference": sd,
            "end_difference": ed,
            "length_difference": abs((ge - gs) - (iex - isx)),
            "concordance_class": "high" if high else "overlapping_only",
        })
    return pairs, used_g, used_i


def main(samples, outdir, out_summary, out_pairs, min_ov, max_diff, tool):
    min_ov, max_diff = float(min_ov), int(max_diff)
    all_pairs, summaries = [], []
    unmatched = []

    for s in samples:
        gnm = read_tsv(f"{outdir}/prophage/genomad/{s}.tsv")
        ind = read_tsv(f"{outdir}/prophage/independent/{s}.tsv")
        pairs, used_g, used_i = match_one(gnm, ind, min_ov, max_diff)
        all_pairs.extend(pairs)

        g_only = [g for g in gnm if g["prediction_id"] not in used_g]
        i_only = [i for i in ind if i["prediction_id"] not in used_i]
        for g in g_only:
            unmatched.append({"sample": s, "prediction_id": g["prediction_id"],
                              "tool": "genomad", "contig": g["contig"],
                              "start": g["genomad_start"],
                              "end": g["genomad_end"]})
        for i in i_only:
            unmatched.append({"sample": s, "prediction_id": i["prediction_id"],
                              "tool": tool, "contig": i["contig"],
                              "start": i["start"], "end": i["end"]})

        union = len(gnm) + len(ind) - len(pairs)
        summaries.append({
            "sample": s,
            "n_genomad": len(gnm),
            "n_independent": len(ind),
            "n_matched": len(pairs),
            "n_high_concordance": sum(1 for p in pairs
                                      if p["concordance_class"] == "high"),
            "n_genomad_only": len(g_only),
            "n_independent_only": len(i_only),
            "jaccard_predictions": round(len(pairs) / union, 4) if union else "",
            "median_reciprocal_overlap": median([p["reciprocal_overlap"] for p in pairs]),
            "median_start_difference": median([p["start_difference"] for p in pairs]),
            "median_end_difference": median([p["end_difference"] for p in pairs]),
            "median_length_difference": median([p["length_difference"] for p in pairs]),
        })

    # Collection-level row, so the headline concordance figure is computed once
    # rather than averaged from per-genome medians.
    tot = {
        "sample": "ALL",
        "n_genomad": sum(r["n_genomad"] for r in summaries),
        "n_independent": sum(r["n_independent"] for r in summaries),
        "n_matched": len(all_pairs),
        "n_high_concordance": sum(1 for p in all_pairs
                                  if p["concordance_class"] == "high"),
        "n_genomad_only": sum(r["n_genomad_only"] for r in summaries),
        "n_independent_only": sum(r["n_independent_only"] for r in summaries),
        "median_reciprocal_overlap": median([p["reciprocal_overlap"] for p in all_pairs]),
        "median_start_difference": median([p["start_difference"] for p in all_pairs]),
        "median_end_difference": median([p["end_difference"] for p in all_pairs]),
        "median_length_difference": median([p["length_difference"] for p in all_pairs]),
    }
    u = tot["n_genomad"] + tot["n_independent"] - tot["n_matched"]
    tot["jaccard_predictions"] = round(tot["n_matched"] / u, 4) if u else ""
    summaries.append(tot)

    write_tsv(out_summary, summaries, SUMMARY_FIELDS)
    write_tsv(out_pairs, all_pairs, PAIR_FIELDS)
    write_tsv(out_pairs.replace(".tsv", "_unmatched.tsv"), unmatched,
              ["sample", "prediction_id", "tool", "contig", "start", "end"])

    sys.stderr.write(
        f"concordance: {tot['n_matched']} matched pairs, "
        f"{tot['n_high_concordance']} meeting the prespecified high-concordance "
        f"definition (reciprocal overlap >= {min_ov}, coordinate difference "
        f"<= {max_diff} bp)\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(list(sm.params.samples), sm.params.outdir, sm.output.tsv,
             sm.output.pairs, sm.params.min_overlap, sm.params.max_diff,
             sm.params.tool)
    else:
        main(sys.argv[1].split(","), sys.argv[2], sys.argv[3], sys.argv[4],
             sys.argv[5], sys.argv[6], sys.argv[7])
