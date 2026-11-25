"""
Map CheckV breakpoints back onto the original contig coordinates.

Protocol Step 7. The prophage database holds both the original geNomad
coordinates and the refined CheckV coordinates, and the refined coordinates are
used downstream. Step 3.4's boundary classification is recomputed on the
refined coordinates, since refinement can move a prediction away from a contig
end or onto one.

CheckV reports its breakpoints in the coordinate space of the sequence it was
given, which here is the flanked extraction. Converting back is a single
offset:

    contig_coordinate = region_start + region_coordinate - 1

Three cases need handling rather than hiding.

  1. CheckV calls the extraction a provirus and reports viral and host
     segments. The viral segment becomes the refined prophage.
  2. CheckV finds no host, so it reports no breakpoint. Adding a 10 kb host
     flank and then being told the whole thing is viral means refinement
     failed on that prediction. The geNomad coordinates are retained and the
     prediction is flagged.
  3. Two predictions shared one extracted region. Each viral segment is
     assigned to whichever prediction it overlaps most in region space.

Boundary stability is the coordinate movement between the 10 kb and the 5 kb
extraction. A prediction whose refined start or end shifts by more than the
concordance tolerance between the two flank lengths has an unstable boundary,
and the protocol treats stability as a condition for the high-confidence set.
"""

import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv, write_tsv, num

FIELDS = ["sample", "prediction_id", "contig", "contig_length",
          "genomad_start", "genomad_end", "genomad_length",
          "refined_start", "refined_end", "refined_length",
          "refinement_status", "region_id", "shared_region",
          "checkv_provirus", "checkv_completeness", "checkv_quality",
          "miuvig_quality", "checkv_contamination_pct",
          "host_length_in_region", "viral_length_in_region",
          "boundary_class_genomad", "boundary_class_refined",
          "start_shift", "end_shift",
          "sens_refined_start", "sens_refined_end",
          "boundary_stability", "flank_left_clipped", "flank_right_clipped",
          "virus_score", "fdr", "n_hallmarks", "taxonomy"]

STABILITY_TOLERANCE = 2000  # bp; mirrors concordance.max_coordinate_difference


def boundary_class(start, end, clen):
    if clen is None or start is None or end is None:
        return "unknown"
    left, right = start <= 1, end >= clen
    if left and right:
        return "truncated_both"
    if left:
        return "truncated_left"
    if right:
        return "truncated_right"
    return "internally_resolved"


def parse_regions(coords, types):
    """
    CheckV writes region_coords_bp as '1-5000,5001-30000' alongside
    region_types as 'host,viral'. Return [(kind, start, end)] in region space.
    """
    out = []
    if not coords or not types:
        return out
    cs = [c.strip() for c in str(coords).split(",") if c.strip()]
    ts = [t.strip().lower() for t in str(types).split(",") if t.strip()]
    for c, t in zip(cs, ts):
        if "-" not in c:
            continue
        a, b = c.split("-", 1)
        a, b = num(a), num(b)
        if a is None or b is None:
            continue
        out.append((t, int(a), int(b)))
    return out


def overlap(a1, a2, b1, b2):
    return max(0, min(a2, b2) - max(a1, b1) + 1)


def refine(prov, mapping, contam, complete):
    """Return {prediction_id: dict} of refined coordinates and CheckV fields."""
    contam_by_region = {r["contig_id"]: r for r in contam}
    complete_by_region = {r["contig_id"]: r for r in complete}
    maps_by_pred = {m["prediction_id"]: m for m in mapping}

    result = {}
    for p in prov:
        pid = p["prediction_id"]
        m = maps_by_pred.get(pid)
        gs, ge = int(p["genomad_start"]), int(p["genomad_end"])
        base = {
            "refined_start": gs, "refined_end": ge,
            "refinement_status": "no_extraction",
            "region_id": "", "shared_region": "",
            "checkv_provirus": "", "checkv_completeness": "",
            "checkv_quality": "", "miuvig_quality": "",
            "checkv_contamination_pct": "",
            "host_length_in_region": "", "viral_length_in_region": "",
            "flank_left_clipped": "", "flank_right_clipped": "",
        }
        if m is None:
            result[pid] = base
            continue

        rid = m["region_id"]
        rstart = int(m["region_start"])
        os_, oe_ = int(m["offset_start_in_region"]), int(m["offset_end_in_region"])
        base.update({
            "region_id": rid,
            "shared_region": m["shared_region"],
            "flank_left_clipped": m["flank_left_clipped"],
            "flank_right_clipped": m["flank_right_clipped"],
        })

        c = contam_by_region.get(rid)
        comp = complete_by_region.get(rid, {})
        base["checkv_completeness"] = comp.get("completeness", "")
        base["checkv_quality"] = comp.get("checkv_quality", "")
        base["miuvig_quality"] = comp.get("miuvig_quality", "")

        if c is None:
            base["refinement_status"] = "checkv_region_missing"
            result[pid] = base
            continue

        base["checkv_provirus"] = c.get("provirus", "")
        base["host_length_in_region"] = c.get("host_length", "")
        base["viral_length_in_region"] = c.get("viral_length", "")
        total = num(c.get("contig_length"))
        host = num(c.get("host_length"))
        if total and host is not None and total > 0:
            base["checkv_contamination_pct"] = round(100.0 * host / total, 2)

        viral = [r for r in parse_regions(c.get("region_coords_bp"),
                                          c.get("region_types"))
                 if r[0] == "viral"]
        if str(c.get("provirus", "")).strip().lower() not in ("yes", "true", "1") \
                or not viral:
            # Case 2. Host flank was supplied and CheckV still found no
            # junction. Keep geNomad and say so.
            base["refinement_status"] = "no_breakpoint_geNomad_retained"
            result[pid] = base
            continue

        # Case 1 and 3. Pick the viral segment that overlaps this prediction
        # most in region space.
        best, best_ov = None, 0
        for _, a, b in viral:
            ov = overlap(a, b, os_, oe_)
            if ov > best_ov:
                best, best_ov = (a, b), ov
        if best is None:
            base["refinement_status"] = "no_overlapping_viral_segment"
            result[pid] = base
            continue

        base["refined_start"] = rstart + best[0] - 1
        base["refined_end"] = rstart + best[1] - 1
        base["refinement_status"] = (
            "refined_shared_region" if m["shared_region"] == "YES"
            else "refined")
        result[pid] = base
    return result


def main(prov_path, p_contam, p_complete, p_map, s_contam, s_map,
         out_path, sample):
    prov = read_tsv(prov_path)
    primary = refine(prov, read_tsv(p_map), read_tsv(p_contam),
                     read_tsv(p_complete))
    sens = refine(prov, read_tsv(s_map), read_tsv(s_contam), [])

    rows = []
    for p in prov:
        pid = p["prediction_id"]
        pr = primary[pid]
        sr = sens.get(pid, {})
        clen = num(p.get("contig_length"))
        gs, ge = int(p["genomad_start"]), int(p["genomad_end"])
        rs, re_ = pr["refined_start"], pr["refined_end"]
        ss, se = sr.get("refined_start"), sr.get("refined_end")

        if ss is None or se is None:
            stability = "not_assessed"
            sshift = eshift = ""
        else:
            sshift, eshift = abs(rs - ss), abs(re_ - se)
            stability = ("stable"
                         if max(sshift, eshift) <= STABILITY_TOLERANCE
                         else "unstable")

        rows.append({
            "sample": sample,
            "prediction_id": pid,
            "contig": p["contig"],
            "contig_length": p.get("contig_length", ""),
            "genomad_start": gs,
            "genomad_end": ge,
            "genomad_length": ge - gs + 1,
            "refined_start": rs,
            "refined_end": re_,
            "refined_length": re_ - rs + 1,
            "refinement_status": pr["refinement_status"],
            "region_id": pr["region_id"],
            "shared_region": pr["shared_region"],
            "checkv_provirus": pr["checkv_provirus"],
            "checkv_completeness": pr["checkv_completeness"],
            "checkv_quality": pr["checkv_quality"],
            "miuvig_quality": pr["miuvig_quality"],
            "checkv_contamination_pct": pr["checkv_contamination_pct"],
            "host_length_in_region": pr["host_length_in_region"],
            "viral_length_in_region": pr["viral_length_in_region"],
            "boundary_class_genomad": p.get("boundary_class", ""),
            "boundary_class_refined": boundary_class(rs, re_, clen),
            "start_shift": sshift,
            "end_shift": eshift,
            "sens_refined_start": "" if ss is None else ss,
            "sens_refined_end": "" if se is None else se,
            "boundary_stability": stability,
            "flank_left_clipped": pr["flank_left_clipped"],
            "flank_right_clipped": pr["flank_right_clipped"],
            "virus_score": p.get("virus_score", ""),
            "fdr": p.get("fdr", ""),
            "n_hallmarks": p.get("n_hallmarks", ""),
            "taxonomy": p.get("taxonomy", ""),
        })

    write_tsv(out_path, rows, FIELDS)
    n_ref = sum(1 for r in rows if r["refinement_status"].startswith("refined"))
    n_unstable = sum(1 for r in rows if r["boundary_stability"] == "unstable")
    sys.stderr.write(
        f"{sample}: {len(rows)} predictions, {n_ref} refined, "
        f"{n_unstable} with an unstable boundary\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.input.prov, sm.input.primary_contam, sm.input.primary_complete,
             sm.input.primary_map, sm.input.sens_contam, sm.input.sens_map,
             sm.output.tsv, sm.wildcards.sample)
    else:
        main(*sys.argv[1:9])
