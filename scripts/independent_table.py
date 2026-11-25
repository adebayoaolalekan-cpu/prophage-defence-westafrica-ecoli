"""
Normalise the independent predictor's output into the geNomad table's columns.

Protocol Step 6. The comparison is made on the original bacterial contig
coordinates, not on separately extracted viral FASTA files, so whichever tool
ran has to be mapped back to contig space before the concordance step sees it.

PIDE's output layout is read from whatever the installed release writes. The
parser below tries the plausible file names and column spellings in turn and
fails with a readable message rather than silently returning zero predictions,
because a silent zero would look like a genome with no prophages.

PhiSpy writes prophage_coordinates.tsv with no header: prophage id, contig,
start, stop, then attachment-site columns.
"""

import csv
import glob
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _common import read_tsv, write_tsv, num

FIELDS = ["sample", "prediction_id", "tool", "contig", "contig_length",
          "start", "end", "length", "score"]

PIDE_CANDIDATES = ["*prophage*.tsv", "*prophage*.csv", "*predictions*.tsv",
                   "*result*.tsv", "*.prophage.tsv"]
START_KEYS = ["start", "Start", "start_pos", "prophage_start", "begin"]
END_KEYS = ["end", "End", "stop", "Stop", "end_pos", "prophage_end"]
CONTIG_KEYS = ["contig", "Contig", "seq_name", "sequence", "seqid",
               "contig_id", "scaffold", "chromosome"]
SCORE_KEYS = ["score", "Score", "probability", "confidence"]


def pick(row, keys):
    for k in keys:
        if k in row and str(row[k]).strip() != "":
            return row[k]
    return None


def parse_pide(indir):
    files = []
    for pat in PIDE_CANDIDATES:
        files.extend(glob.glob(os.path.join(indir, "**", pat), recursive=True))
    files = sorted(set(files))
    if not files:
        raise SystemExit(
            f"No PIDE output table found under {indir}. Looked for "
            f"{PIDE_CANDIDATES}. Set the pattern to match the installed "
            f"release, or run with independent_predictor: phispy.")

    out = []
    for f in files:
        delim = "," if f.endswith(".csv") else "\t"
        with open(f, newline="") as fh:
            for r in csv.DictReader(fh, delimiter=delim):
                c, s, e = pick(r, CONTIG_KEYS), pick(r, START_KEYS), pick(r, END_KEYS)
                if c is None or s is None or e is None:
                    continue
                out.append((str(c), num(s), num(e), pick(r, SCORE_KEYS) or ""))
    if not out:
        raise SystemExit(
            f"PIDE output at {indir} parsed to zero rows. Column names did "
            f"not match any of {CONTIG_KEYS}, {START_KEYS}, {END_KEYS}. "
            f"A genuine zero-prophage genome must be distinguishable from a "
            f"parse failure, so this is an error rather than an empty table.")
    return out


def parse_phispy(indir):
    path = os.path.join(indir, "prophage_coordinates.tsv")
    if not os.path.exists(path):
        # PhiSpy writes nothing when it finds nothing.
        return []
    out = []
    with open(path, newline="") as fh:
        for row in csv.reader(fh, delimiter="\t"):
            if len(row) < 4:
                continue
            out.append((row[1], num(row[2]), num(row[3]), ""))
    return out


def main(indir, seqkit_path, out_path, tool, sample):
    raw = parse_pide(indir) if tool == "pide" else parse_phispy(indir)

    # Contig lengths come from the geNomad-side table where available, so the
    # two tables agree on contig space. Falling back to blank is acceptable
    # here because the concordance step only needs coordinates.
    sibling = out_path.replace("/independent/", "/genomad/")
    lens = {r["contig"]: r["contig_length"] for r in read_tsv(sibling)}

    rows = []
    for i, (contig, s, e, score) in enumerate(sorted(raw), start=1):
        if s is None or e is None:
            continue
        if s > e:
            s, e = e, s
        # PhiSpy reports contigs under their original name; this workflow
        # prefixes them with the sample. Accept either spelling.
        if contig not in lens and f"{sample}|{contig}" in lens:
            contig = f"{sample}|{contig}"
        rows.append({
            "sample": sample,
            "prediction_id": f"{sample}_ind{i:04d}",
            "tool": tool,
            "contig": contig,
            "contig_length": lens.get(contig, ""),
            "start": s,
            "end": e,
            "length": e - s + 1,
            "score": score,
        })

    write_tsv(out_path, rows, FIELDS)
    sys.stderr.write(f"{sample}: {len(rows)} {tool} predictions\n")


if __name__ == "__main__":
    try:
        sm = snakemake  # noqa: F821
    except NameError:
        sm = None
    if sm is not None:
        main(sm.params.indir, sm.input.seqkit, sm.output.tsv, sm.params.tool,
             sm.wildcards.sample)
    else:
        main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], sys.argv[5])
