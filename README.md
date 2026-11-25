# Genomic pilot

Prophage Burden, Prophage-Associated Virulence Determinants and Anti-Phage
Defence Systems in West African *Escherichia coli*.

The pilot exists to close the one open item on the registration hold point.
Condition 4 of the protocol's power assessment currently holds at a rate ratio
of 2.0 on the primary inference set because the negative binomial dispersion
was assumed rather than observed. Section 10.1 cannot be re-run until a real
prophage count distribution exists. Three further questions are settled at the
same time: the PopPUNK lineage partition the clustered standard errors need,
the runtime and memory figures the production run has to be sized against, and
whether PIDE is practical at collection scale.

Nothing here is exploratory. Every parameter is the one the protocol fixes, and
changing one is an amendment rather than a configuration change.

## The selected set

96 genomes in two blocks. The selection is in `config/pilot_manifest.tsv` and
the workflow reads `config/samples.tsv`.

Block A, 48 genomes, is the inference block. It is drawn only from the two
originating studies that satisfy the eligibility rule for the human versus
poultry contrast, PRJNA293225-A and PRJNA1019343, restricted to the
high-contiguity collection and balanced 24 human against 24 poultry and 24 per
study. Dispersion estimated anywhere else would not describe the set the
primary inference runs on, which is the whole reason this block is separate.

Block B, 48 genomes, is the catalogue block. It is drawn from the rest of the
retrieval across ecological source, country and quality band, and it includes
genomes that fail the broad-quality threshold on purpose. The truncation rate
among excluded genomes is what shows whether the quality filter is earning its
place. Block B spans 30 distinct originating studies and 12 of the 13 countries
in the catalogue.

No accession appears in both blocks, and block B borrows nothing from the two
qualifying studies.

## Running it

```bash
python3 scripts/pilot_select.py --stage 1     # already run, outputs committed
bash scripts/fetch_candidate_stats.sh         # already run
python3 scripts/pilot_select.py --stage 2     # already run, outputs committed
python3 -I tests/test_chain.py                # 44 checks on the arithmetic
snakemake -n --configfile config/config.yaml  # dry run, prints the DAG
snakemake --use-conda --cores 32 --configfile config/config.yaml
```

The selection is already done and its outputs are in `config/`. Both draws use
a fixed seed, so re-running reproduces the same genomes.

### Before the first run

Two things need setting in `config/config.yaml`.

The database paths under `databases:` are absolute and the workflow does not
download any of them. A database resolved at run time cannot be pinned, and an
unpinned release cannot be reported. Build them once:

```bash
genomad download-database /data/db
checkv download_database /data/db
checkm2 database --download --path /data/db/CheckM2_database
bakta_db download --output /data/db/bakta --type light
```

The PIDE command line is deliberately left unset. Fill `pide.command` from the
installed release's own documentation; the placeholders `{input}`, `{outdir}`
and `{threads}` are substituted. Leaving the default in place makes the
workflow fail at parse time rather than run a guessed command. If PIDE proves
impractical, set `independent_predictor: phispy` and the declared substitute
runs instead. Which one ran is recorded in the provenance table, not inferred.

### What it will cost

1,356 jobs for 96 genomes. geNomad and Bakta dominate, CheckV runs twice per
genome because the 10 kb and 5 kb extractions are both needed for the boundary
stability measure, and CheckM2 and PopPUNK each run once across the whole set.
The resource report is the pilot's answer to how long the full collection will
take, so do not estimate it from this paragraph.

## What it produces

| File under `results/reports/` | Contents |
| --- | --- |
| `dispersion_estimates.json` | Negative binomial theta with confidence intervals, inference block and catalogue block separately. The number Section 10.1 is held up on. |
| `prophage_counts.tsv` | Per-genome burden, internally resolved predictions separated from truncated candidates. |
| `concordance.tsv` | geNomad against the independent predictor, at the prespecified 80% reciprocal overlap and 2 kb coordinate tolerance. |
| `defence_systems.tsv` | DefenseFinder systems and genes per genome. |
| `crispr_arrays.tsv`, `crispr_spacers.fna` | CRISPR arrays, evidence levels and spacer sequences. |
| `resource_report.tsv`, `resource_per_tool.tsv` | Wall time, CPU time and peak memory per tool per genome, with output footprints. |
| `provenance.tsv` | Software versions, database releases and every parameter actually used. |
| `pilot_summary.md` | The four pilot questions, answered against the files above. |

`results/qc/band_changes.tsv` lists every genome whose quality band changed once
CheckM2 completeness and contamination were available. A genome that drops out
of its selected band is reported, never replaced. Replacing it would make the
pilot a quality-selected sample of a quality-selected sample, and the dispersion
estimate would understate the variance the production run will meet.

## Licence and citation

The workflow, scripts and tests are MIT licensed. The selection tables and
retrieval metadata are dedicated to the public domain under CC0, since they are
derived from ENA and NCBI records that nobody owns. LICENSE and LICENSE-DATA
set out which files fall under which, and CITATION.cff carries the citation
metadata.

This pilot implements Section 10.2 of a registered protocol. Cite both:

  Olalekan, A. A. (2026). Prophages, Prophage-Encoded Virulence Genes and
  Anti-Phage Defence Systems in West African Escherichia coli
  [Preregistration]. OSF. https://doi.org/10.17605/OSF.IO/FXMEQ

## Known gaps

Four candidate accessions, GCA_980407635, GCA_980408575, GCA_980409345 and
GCA_980410785, returned no record from NCBI Datasets and were dropped from the
block B frame. All four are Liberian, so Liberia falls out of the pilot for a
mirroring reason rather than a sampling one. The catalogue holds five further
Liberian accessions that the draw did not reach. This is the completeness gap
running in the opposite direction to the one documented in Step 1.1: ENA holds
assemblies that NCBI Datasets has not yet exposed, just as NCBI held 117 that
ENA had not.

One cell, other animal at broad quality, had no candidates at all. The slot was
backfilled from another band of the same source, which holds the ecological
balance and spends the loss on the band balance. That is the right way round.
Source balance carries the catalogue description; band balance is a convenience.
