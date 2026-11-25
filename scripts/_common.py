"""Shared helpers for the pilot workflow scripts."""

import csv
import os
import re


def read_fasta(path):
    """Yield (name, sequence). Names are taken up to the first whitespace."""
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                name, chunks = line[1:].split()[0], []
            else:
                chunks.append(line)
    if name is not None:
        yield name, "".join(chunks)


def contig_lengths(path):
    return {n: len(s) for n, s in read_fasta(path)}


def read_tsv(path):
    if not os.path.exists(path):
        return []
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def write_tsv(path, rows, fields):
    d = os.path.dirname(path)
    if d:
        os.makedirs(d, exist_ok=True)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        for r in rows:
            w.writerow(r)


def num(v, default=None):
    try:
        f = float(v)
    except (TypeError, ValueError):
        return default
    return int(f) if f.is_integer() else f


def split_contig_id(seq_name):
    """
    geNomad writes a provirus as '<contig>|provirus_<start>_<end>'. Contigs in
    this workflow are already prefixed '<sample>|<original_contig>', so the
    provirus suffix is the LAST pipe-delimited field and the contig is
    everything before it.
    """
    m = re.match(r"^(?P<contig>.+)\|provirus_(?P<s>\d+)_(?P<e>\d+)$", seq_name)
    if m:
        return m.group("contig"), int(m.group("s")), int(m.group("e"))
    return seq_name, None, None


def wrap(seq, width=70):
    return "\n".join(seq[i:i + width] for i in range(0, len(seq), width))
