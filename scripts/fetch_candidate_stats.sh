#!/usr/bin/env bash
# Fetch assembly statistics for the block B candidate frame.
#
# The ENA portal does not expose contig N50 or contig count for assembly
# records, so the quality band cannot be assigned from the retrieval metadata.
# This pulls the statistics from NCBI Datasets for the candidate accessions
# only, which keeps the request small enough to run interactively.
#
# Requires the NCBI datasets command line client, v16 or later, and jq.
#   conda install -c conda-forge ncbi-datasets-cli jq
#
# Usage: bash scripts/fetch_candidate_stats.sh
set -euo pipefail

CAND=${1:-config/pilot_block_b_candidates.tsv}
OUT=${2:-resources/candidate_assembly_stats.tsv}
BATCH=${BATCH:-50}

command -v datasets >/dev/null || { echo "datasets client not found" >&2; exit 1; }
command -v jq >/dev/null || { echo "jq not found" >&2; exit 1; }

mkdir -p "$(dirname "$OUT")" tmp_stats
tail -n +2 "$CAND" | cut -f1 | sort -u > tmp_stats/accs.txt
total=$(wc -l < tmp_stats/accs.txt)
echo "resolving $total candidate accessions in batches of $BATCH" >&2

printf 'assembly_accession\tassembly_level\tcontig_n50\tnumber_of_contigs\ttotal_length\tgc_percent\tsubmitter\n' > "$OUT"

split -l "$BATCH" tmp_stats/accs.txt tmp_stats/batch_
for b in tmp_stats/batch_*; do
  echo "  $(basename "$b")" >&2
  # --inputfile takes one accession per line. The report is streamed as
  # newline-delimited JSON, one object per assembly.
  datasets summary genome accession --inputfile "$b" --as-json-lines \
    | jq -r '
        [ (.accession // ""),
          (.assembly_info.assembly_level // ""),
          (.assembly_stats.contig_n50 // ""),
          (.assembly_stats.number_of_contigs // ""),
          (.assembly_stats.total_sequence_length // ""),
          (.assembly_stats.gc_percent // ""),
          (.assembly_info.biosample.submission_date // "")
        ] | @tsv' >> "$OUT"
  sleep 1
done

got=$(( $(wc -l < "$OUT") - 1 ))
echo "wrote $got of $total records to $OUT" >&2
if [ "$got" -lt "$total" ]; then
  echo "WARNING $(( total - got )) accessions returned no record. Stage 2 drops" >&2
  echo "        them from the frame and reports the count." >&2
fi
rm -rf tmp_stats
