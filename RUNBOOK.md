# Running the pilot

Everything here happens after registration, which is the point. The
registration at https://doi.org/10.17605/OSF.IO/FXMEQ fixed the design and
the decision rules before any of this was run.

Budget roughly a day of wall time: a few hours building databases, a few hours
for the workflow itself on a cluster, longer on a single machine.

## 0. What you need

A Linux machine with at least 16 cores, 64 GB of memory and 200 GB of free
disk. Less will work but CheckM2 and PopPUNK are the memory-hungry steps and
will thrash below about 48 GB. A GPU is needed only if you run PIDE.

Conda or Mamba. Snakemake 8 with the Slurm executor plugin if you are on a
cluster, or Snakemake 7 with `--cluster` if that is what your site has.

Outbound network access to NCBI, Zenodo and the tool database hosts. The
workflow downloads the 96 genomes itself.

## 1. Settle the independent predictor

The workflow refuses to start until you deal with this, which is deliberate.

`config/config.yaml` has `pide.command` set to a placeholder. Fill it in from
the installed PIDE release's own documentation. The placeholders `{input}`,
`{outdir}` and `{threads}` are substituted. A guessed command line that
happens to run is worse than one that fails, because it enters the provenance
record looking correct.

If PIDE proves impractical, set `independent_predictor: phispy` and the
declared substitute runs instead. Either way the provenance table records
which predictor actually ran rather than which was intended. Switching is a
configuration change, not an amendment: both predictors are named in the
registered protocol.

## 2. Build the databases

```bash
bash scripts/setup_databases.sh /data/db
```

About 90 GB and a few hours. Re-running skips what is already there.

When it finishes it writes `/data/db/installed_versions.txt`. Copy those
versions into the `versions` block of `config/config.yaml` and point the
`databases` block at the paths under `/data/db`. Do this before the run, not
after. Versions recorded afterwards are recollections.

## 3. Dry run

```bash
snakemake -n --configfile config/config.yaml
```

Expect 1,356 jobs across the 96 genomes. If you see a different number,
something is wrong with the sample sheet and you should stop and find out
what before spending the compute.

```bash
python3 -I tests/test_chain.py
```

Forty-four checks on the coordinate arithmetic between prediction, flanked
extraction and boundary refinement. A silent error there corrupts every
downstream count without raising anything. Run it on the machine that will do
the work, not just on the machine where you read this.

## 4. Run

On a cluster, once the Slurm executor plugin is installed:

```bash
pip install snakemake-executor-plugin-slurm
snakemake --profile profiles/slurm --configfile config/config.yaml
```

Without that plugin Snakemake rejects the profile with `invalid choice:
'slurm'`, which looks like a broken profile but is only a missing dependency.
The memory and runtime figures in `profiles/slurm/config.yaml` are starting
estimates. Replace them with the observed values from the pilot before the
full collection.

On a single machine:

```bash
snakemake --use-conda --cores 16 --configfile config/config.yaml
```

The first invocation builds thirteen conda environments, which takes a while
and happens once. `--keep-going` is already set in the profile so one failing
genome does not halt the rest.

## 5. Read the output

Start with `results/reports/pilot_summary.md`, which answers the four pilot
questions directly. Then:

| File | What to do with it |
| --- | --- |
| `reports/dispersion_estimates.json` | The number Section 10.1 is held on. Take the covariate-adjusted estimate on the inference block, with its bootstrap interval, and apply the decision table in Table 17 of the protocol. |
| `reports/resource_per_tool.tsv` | Observed memory and wall time per tool per genome. Replace the estimates in `profiles/slurm/config.yaml` with these before the full collection. |
| `poppunk/lineage_assignments.tsv` | The observed lineage partition for the clustered standard errors. |
| `reports/concordance.tsv` | geNomad against the independent predictor. A low figure is a finding, not a failure. |
| `qc/band_changes.tsv` | Genomes whose quality band changed once completeness and contamination were known. Reported, never replaced. |
| `reports/provenance.tsv` | Versions, database releases and every parameter actually used. This is what goes in the manuscript. |

## 6. Re-run the power assessment

Take the observed dispersion, mean count and lineage partition and repeat the
Section 10.1 simulation. Table 17 of the protocol states in advance what each
possible outcome permits, so this step applies a rule rather than making a
decision. Record the result as an amendment to the registration.

Only then does the full collection get processed.

## If something goes wrong

A failed genome leaves its outputs missing and the run continues. Find them
with `snakemake -n` afterwards, which lists what is still to do.

Conda environment solves that fail are nearly always a channel problem. The
environment files pin channels deliberately. Do not relax them to get a solve.

CRISPRCasFinder needs `sel392v2.so` at the path in the config. It ships with
the tool but not always where the config expects.

If the run dies partway, Snakemake resumes from where it stopped. Do not
delete `results/` to start clean unless you have a reason, because the
download step is the slowest part to repeat.
