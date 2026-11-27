#!/usr/bin/env bash
# Build the five tool databases the pilot needs and record what was installed.
#
# Run this once, on the machine that will run the workflow, before the first
# Snakemake invocation. It is deliberately not a Snakemake rule: a database
# resolved at run time cannot be pinned, and a release that cannot be pinned
# cannot be reported in the provenance table.
#
#   bash scripts/setup_databases.sh /data/db
#
# Needs roughly 90 GB free and a few hours on a reasonable connection. Bakta's
# light database is the large one. Re-running skips anything already present.

set -euo pipefail

DB="${1:-/data/db}"
mkdir -p "$DB"
LOG="$DB/setup.log"
exec > >(tee -a "$LOG") 2>&1
echo "=== database setup started $(date -u +%Y-%m-%dT%H:%M:%SZ) into $DB ==="

need() { command -v "$1" >/dev/null 2>&1 || { echo "MISSING: $1 is not on PATH. Activate the right conda environment first."; exit 1; }; }

step() { echo; echo "--- $1 ---"; }

# geNomad
step "geNomad"
if [ -d "$DB/genomad_db" ]; then
  echo "already present, skipping"
else
  need genomad
  genomad download-database "$DB"
fi

# CheckV
step "CheckV"
if compgen -G "$DB/checkv-db-*" > /dev/null; then
  echo "already present, skipping"
else
  need checkv
  checkv download_database "$DB"
fi

# CheckM2
step "CheckM2"
if [ -f "$DB/CheckM2_database/uniref100.KO.1.dmnd" ]; then
  echo "already present, skipping"
else
  need checkm2
  checkm2 database --download --path "$DB/CheckM2_database"
fi

# Bakta
step "Bakta (light)"
if [ -d "$DB/bakta/db-light" ]; then
  echo "already present, skipping"
else
  need bakta_db
  mkdir -p "$DB/bakta"
  bakta_db download --output "$DB/bakta" --type light
fi

# DefenseFinder models
step "DefenseFinder models"
if [ -d "${HOME}/.macsyfinder/data/defense-finder-models" ]; then
  echo "already present, skipping"
else
  need defense-finder
  defense-finder update
fi

# Record what landed, so the provenance table is filled from fact rather than memory
step "recording installed versions"
VER="$DB/installed_versions.txt"
{
  echo "# recorded $(date -u +%Y-%m-%dT%H:%M:%SZ) on $(hostname)"
  for t in genomad checkv checkm2 bakta seqkit poppunk defense-finder; do
    if command -v "$t" >/dev/null 2>&1; then
      printf '%-16s %s\n' "$t" "$("$t" --version 2>&1 | head -1)"
    else
      printf '%-16s NOT INSTALLED\n' "$t"
    fi
  done
  echo
  echo "# database paths and sizes"
  for d in "$DB"/genomad_db "$DB"/checkv-db-* "$DB"/CheckM2_database "$DB"/bakta/db-light; do
    [ -e "$d" ] && printf '%-40s %s\n' "$d" "$(du -sh "$d" 2>/dev/null | cut -f1)"
  done
} > "$VER"
cat "$VER"

echo
echo "=== done. Copy the versions above into the versions block of config/config.yaml,"
echo "    and set the databases block to paths under $DB. ==="
