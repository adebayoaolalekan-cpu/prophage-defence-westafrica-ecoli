# ---------------------------------------------------------------------------
# Genomic pilot workflow
#
# Prophage Burden, Prophage-Associated Virulence Determinants and Anti-Phage
# Defence Systems in West African Escherichia coli
#
# The pilot has to deliver four things before the full collection is processed.
#   1. The observed prophage count distribution and the negative binomial
#      dispersion for the primary inference set, so that Section 10.1 can be
#      re-run on observed rather than assumed dispersion.
#   2. The PopPUNK lineage partition, so the clustered standard errors in the
#      power simulation use the real partition.
#   3. Runtime, peak memory and disk per genome per tool, and the workflow
#      failure rate, so the production run can be sized.
#   4. Whether PIDE is practical at collection scale, and the geNomad to PIDE
#      boundary concordance.
#
# Everything the workflow computes is written to a file. Nothing is reported
# from a log scrape.
#
# Run:
#   snakemake --use-conda --cores 32 --configfile config/config.yaml
#   snakemake --use-conda --cores 32 -n            # dry run, prints the DAG
# ---------------------------------------------------------------------------

import os
import pandas as pd

configfile: "config/config.yaml"

OUT = config["outdir"]
DB = config["databases"]
THREADS = config["threads"]

samples = pd.read_csv(config["samples"], sep="\t", dtype=str)
SAMPLES = samples["sample"].tolist()
ACC = dict(zip(samples["sample"], samples["accession"]))
BLOCK = dict(zip(samples["sample"], samples["block"]))

INDEP = config["independent_predictor"].lower()
if INDEP not in ("pide", "phispy"):
    raise WorkflowError("independent_predictor must be 'pide' or 'phispy'")
if INDEP == "pide" and config["pide"]["command"].startswith("PLEASE_SET"):
    raise WorkflowError(
        "pide.command is not set in config/config.yaml. Fill it in from the "
        "installed PIDE release's documentation, or set "
        "independent_predictor: phispy to use the declared substitute. The "
        "workflow will not run a guessed command line."
    )

wildcard_constraints:
    sample=r"[A-Za-z0-9_.\-]+",
    flank=r"primary|sensitivity",


# ---------------------------------------------------------------------------
rule all:
    input:
        f"{OUT}/reports/pilot_summary.md",
        f"{OUT}/reports/resource_report.tsv",
        f"{OUT}/reports/prophage_counts.tsv",
        f"{OUT}/reports/dispersion_estimates.json",
        f"{OUT}/reports/concordance.tsv",
        f"{OUT}/reports/defence_systems.tsv",
        f"{OUT}/reports/crispr_arrays.tsv",
        f"{OUT}/reports/provenance.tsv",
        f"{OUT}/poppunk/lineage_assignments.tsv",


# ---------------------------------------------------------------------------
# Stage 1. Genomes and quality
# ---------------------------------------------------------------------------

rule download_genome:
    """One accession per rule instance, so a single withdrawn accession fails
    one job rather than the whole download."""
    output:
        fna=f"{OUT}/genomes/{{sample}}.fna",
    params:
        acc=lambda wc: ACC[wc.sample],
    benchmark:
        f"{OUT}/benchmarks/download/{{sample}}.tsv"
    conda:
        "envs/datasets.yaml"
    retries: 3
    shell:
        r"""
        tmp=$(mktemp -d)
        trap 'rm -rf "$tmp"' EXIT
        datasets download genome accession {params.acc} \
            --include genome --filename "$tmp/g.zip" --no-progressbar
        unzip -q -o "$tmp/g.zip" -d "$tmp/unz"
        # Concatenate in case the accession carries more than one FASTA, and
        # prefix every contig with the sample identifier so coordinates stay
        # traceable to a genome after any downstream merge.
        find "$tmp/unz" -name '*.fna' -print0 \
          | xargs -0 cat \
          | awk -v s={wildcards.sample} \
              '/^>/ {{ sub(/^>/,""); print ">" s "|" $1; next }} {{ print }}' \
          > {output.fna}
        test -s {output.fna}
        """


rule seqkit_stats:
    input:
        fna=f"{OUT}/genomes/{{sample}}.fna",
    output:
        tsv=f"{OUT}/qc/seqkit/{{sample}}.tsv",
    benchmark:
        f"{OUT}/benchmarks/seqkit/{{sample}}.tsv"
    conda:
        "envs/seqkit.yaml"
    shell:
        "seqkit stats -a -T {input.fna} > {output.tsv}"


rule checkm2:
    """CheckM2 is run once over every genome. Its per-genome cost is dominated
    by the DIAMOND pass, which batches far better than it parallelises."""
    input:
        fna=expand(f"{OUT}/genomes/{{sample}}.fna", sample=SAMPLES),
    output:
        tsv=f"{OUT}/qc/checkm2/quality_report.tsv",
    params:
        indir=f"{OUT}/genomes",
        outdir=f"{OUT}/qc/checkm2",
        db=DB["checkm2"],
    threads: THREADS["checkm2"]
    benchmark:
        f"{OUT}/benchmarks/checkm2/all.tsv"
    conda:
        "envs/checkm2.yaml"
    shell:
        r"""
        export CHECKM2DB={params.db}
        checkm2 predict --threads {threads} --input {params.indir} \
            --output-directory {params.outdir} -x fna --force
        test -s {output.tsv}
        """


rule quality_bands:
    """Assign the Table 9 bands from the observed statistics. The band assigned
    at selection time used assembly statistics only; this is the definitive
    assignment, and any genome that changes band is reported rather than
    replaced."""
    input:
        seqkit=expand(f"{OUT}/qc/seqkit/{{sample}}.tsv", sample=SAMPLES),
        checkm2=f"{OUT}/qc/checkm2/quality_report.tsv",
        manifest=config["manifest"],
    output:
        tsv=f"{OUT}/qc/quality_bands.tsv",
        changed=f"{OUT}/qc/band_changes.tsv",
    params:
        cfg="config/config.yaml",
    conda:
        "envs/python.yaml"
    script:
        "scripts/quality_bands.py"


# ---------------------------------------------------------------------------
# Stage 2. Annotation
# ---------------------------------------------------------------------------

rule bakta:
    """Annotation supplies the proteins DefenseFinder consumes and the GenBank
    file PhiSpy needs if it is substituted for PIDE."""
    input:
        fna=f"{OUT}/genomes/{{sample}}.fna",
    output:
        faa=f"{OUT}/annotation/{{sample}}/{{sample}}.faa",
        gbff=f"{OUT}/annotation/{{sample}}/{{sample}}.gbff",
    params:
        outdir=lambda wc: f"{OUT}/annotation/{wc.sample}",
        db=DB["bakta"],
    threads: THREADS["bakta"]
    benchmark:
        f"{OUT}/benchmarks/bakta/{{sample}}.tsv"
    conda:
        "envs/bakta.yaml"
    shell:
        r"""
        bakta --db {params.db} --output {params.outdir} --prefix {wildcards.sample} \
            --threads {threads} --force --skip-plot {input.fna}
        """


# ---------------------------------------------------------------------------
# Stage 3. Prophage identification
# ---------------------------------------------------------------------------

rule genomad:
    input:
        fna=f"{OUT}/genomes/{{sample}}.fna",
    output:
        summary=f"{OUT}/genomad/{{sample}}/{{sample}}_summary/{{sample}}_virus_summary.tsv",
        genes=f"{OUT}/genomad/{{sample}}/{{sample}}_summary/{{sample}}_virus_genes.tsv",
    params:
        outdir=lambda wc: f"{OUT}/genomad/{wc.sample}",
        db=DB["genomad"],
        splits=config["genomad"]["splits"],
        score=lambda wc: ("" if config["genomad"]["min_score"] in (None, "")
                          else f"--min-score {config['genomad']['min_score']}"),
        extra=config["genomad"]["extra"],
    threads: THREADS["genomad"]
    benchmark:
        f"{OUT}/benchmarks/genomad/{{sample}}.tsv"
    conda:
        "envs/genomad.yaml"
    shell:
        r"""
        genomad end-to-end --threads {threads} --splits {params.splits} \
            {params.score} {params.extra} \
            {input.fna} {params.outdir} {params.db}
        """


rule provirus_table:
    """Pull the integrated proviral predictions out of the geNomad summary and
    drop everything geNomad called a free virus, a plasmid or another mobile
    element. Protocol Step 5: only integrated proviruses contribute to burden."""
    input:
        summary=f"{OUT}/genomad/{{sample}}/{{sample}}_summary/{{sample}}_virus_summary.tsv",
        seqkit=f"{OUT}/qc/seqkit/{{sample}}.tsv",
        fna=f"{OUT}/genomes/{{sample}}.fna",
    output:
        tsv=f"{OUT}/prophage/genomad/{{sample}}.tsv",
    conda:
        "envs/python.yaml"
    script:
        "scripts/provirus_table.py"


if INDEP == "pide":

    rule pide:
        input:
            fna=f"{OUT}/genomes/{{sample}}.fna",
        output:
            raw=f"{OUT}/prophage/pide_raw/{{sample}}/done.txt",
        params:
            outdir=lambda wc: f"{OUT}/prophage/pide_raw/{wc.sample}",
            command=config["pide"]["command"],
        threads: THREADS["pide"]
        resources:
            gpu=1 if config["pide"]["gpu"] else 0,
        benchmark:
            f"{OUT}/benchmarks/pide/{{sample}}.tsv"
        conda:
            "envs/pide.yaml"
        shell:
            r"""
            mkdir -p {params.outdir}
            cmd=$(printf '%s' "{params.command}" \
                  | sed -e "s#{{input}}#{input.fna}#g" \
                        -e "s#{{outdir}}#{params.outdir}#g" \
                        -e "s#{{threads}}#{threads}#g")
            echo "$cmd" > {params.outdir}/command.txt
            eval "$cmd"
            touch {output.raw}
            """

else:

    rule phispy:
        """Declared substitute, protocol Step 6 and Section 11.0."""
        input:
            gbff=f"{OUT}/annotation/{{sample}}/{{sample}}.gbff",
        output:
            raw=f"{OUT}/prophage/pide_raw/{{sample}}/done.txt",
        params:
            outdir=lambda wc: f"{OUT}/prophage/pide_raw/{wc.sample}",
            extra=config["phispy"]["extra"],
        threads: THREADS["pide"]
        benchmark:
            f"{OUT}/benchmarks/pide/{{sample}}.tsv"
        conda:
            "envs/phispy.yaml"
        shell:
            r"""
            PhiSpy.py {input.gbff} -o {params.outdir} --threads {threads} \
                {params.extra}
            touch {output.raw}
            """


rule independent_table:
    """Normalise whichever independent predictor ran into the same columns as
    the geNomad table, so the concordance step does not care which one it was."""
    input:
        raw=f"{OUT}/prophage/pide_raw/{{sample}}/done.txt",
        seqkit=f"{OUT}/qc/seqkit/{{sample}}.tsv",
    output:
        tsv=f"{OUT}/prophage/independent/{{sample}}.tsv",
    params:
        tool=INDEP,
        indir=lambda wc: f"{OUT}/prophage/pide_raw/{wc.sample}",
    conda:
        "envs/python.yaml"
    script:
        "scripts/independent_table.py"


# ---------------------------------------------------------------------------
# Stage 4. CheckV boundary refinement
# ---------------------------------------------------------------------------

rule checkv_extract:
    """Cut each geNomad prediction out of its host contig with the Table 10
    flank, clipping at contig ends and merging windows that overlap. CheckV is
    designed for single-contig viral genomes, so neither a whole bacterial
    contig nor an already-trimmed provirus is a valid input."""
    input:
        fna=f"{OUT}/genomes/{{sample}}.fna",
        prov=f"{OUT}/prophage/genomad/{{sample}}.tsv",
    output:
        fasta=f"{OUT}/checkv/extract/{{sample}}.{{flank}}.fna",
        map=f"{OUT}/checkv/extract/{{sample}}.{{flank}}.map.tsv",
    params:
        flank=lambda wc: config["checkv"][f"flank_{wc.flank}"],
        merge=config["checkv"]["merge_overlapping_windows"],
    conda:
        "envs/python.yaml"
    script:
        "scripts/checkv_extract.py"


rule checkv_contamination:
    input:
        fasta=f"{OUT}/checkv/extract/{{sample}}.{{flank}}.fna",
    output:
        contam=f"{OUT}/checkv/{{flank}}/{{sample}}/contamination.tsv",
        complete=f"{OUT}/checkv/{{flank}}/{{sample}}/completeness.tsv",
    params:
        outdir=lambda wc: f"{OUT}/checkv/{wc.flank}/{wc.sample}",
        db=DB["checkv"],
    threads: THREADS["checkv"]
    benchmark:
        f"{OUT}/benchmarks/checkv_{{flank}}/{{sample}}.tsv"
    conda:
        "envs/checkv.yaml"
    shell:
        r"""
        if [ ! -s {input.fasta} ]; then
            # No prediction in this genome. Write the headers CheckV would have
            # written so the aggregation step sees an empty result rather than a
            # missing file.
            mkdir -p {params.outdir}
            printf 'contig_id\tcontig_length\tviral_length\thost_length\tprovirus\tregion_types\tregion_lengths\tregion_coords_bp\n' > {output.contam}
            printf 'contig_id\tcontig_length\tprovirus\tproviral_length\tgene_count\tviral_genes\thost_genes\tcheckv_quality\tmiuvig_quality\tcompleteness\tcompleteness_method\tcontamination\tkmer_freq\twarnings\n' > {output.complete}
            exit 0
        fi
        checkv end_to_end {input.fasta} {params.outdir} -t {threads} -d {params.db} \
            --remove_tmp
        test -s {output.contam}
        """


rule refine_boundaries:
    """Map the CheckV breakpoints back onto the original contig coordinates and
    classify every prediction by its position relative to the contig ends
    (Step 3.4). Boundary stability is the difference between the 10 kb and 5 kb
    extractions."""
    input:
        prov=f"{OUT}/prophage/genomad/{{sample}}.tsv",
        primary_contam=f"{OUT}/checkv/primary/{{sample}}/contamination.tsv",
        primary_complete=f"{OUT}/checkv/primary/{{sample}}/completeness.tsv",
        primary_map=f"{OUT}/checkv/extract/{{sample}}.primary.map.tsv",
        sens_contam=f"{OUT}/checkv/sensitivity/{{sample}}/contamination.tsv",
        sens_map=f"{OUT}/checkv/extract/{{sample}}.sensitivity.map.tsv",
    output:
        tsv=f"{OUT}/prophage/refined/{{sample}}.tsv",
    conda:
        "envs/python.yaml"
    script:
        "scripts/refine_boundaries.py"


# ---------------------------------------------------------------------------
# Stage 5. Defence systems and CRISPR
# ---------------------------------------------------------------------------

rule defensefinder:
    """Run on the complete assembly's proteins, not on predicted viral regions
    (Step 10)."""
    input:
        faa=f"{OUT}/annotation/{{sample}}/{{sample}}.faa",
    output:
        systems=f"{OUT}/defence/{{sample}}/defense_finder_systems.tsv",
        genes=f"{OUT}/defence/{{sample}}/defense_finder_genes.tsv",
    params:
        outdir=lambda wc: f"{OUT}/defence/{wc.sample}",
    threads: THREADS["defensefinder"]
    benchmark:
        f"{OUT}/benchmarks/defensefinder/{{sample}}.tsv"
    conda:
        "envs/defensefinder.yaml"
    shell:
        r"""
        defense-finder run --out-dir {params.outdir} --workers {threads} \
            --preserve-raw {input.faa}
        test -s {output.systems}
        """


rule crisprcasfinder:
    output:
        results=f"{OUT}/crispr/{{sample}}/result.json",
    input:
        fna=f"{OUT}/genomes/{{sample}}.fna",
    params:
        outdir=lambda wc: f"{OUT}/crispr/{wc.sample}",
        so=DB["crisprcasfinder_so"],
    threads: THREADS["crisprcasfinder"]
    benchmark:
        f"{OUT}/benchmarks/crisprcasfinder/{{sample}}.tsv"
    conda:
        "envs/crisprcasfinder.yaml"
    shell:
        r"""
        rm -rf {params.outdir}
        CRISPRCasFinder.pl -in {input.fna} -out {params.outdir} \
            -so {params.so} -cas -def General -cpuM {threads} -keep -q
        test -s {output.results}
        """


rule defence_table:
    input:
        systems=expand(f"{OUT}/defence/{{sample}}/defense_finder_systems.tsv", sample=SAMPLES),
        genes=expand(f"{OUT}/defence/{{sample}}/defense_finder_genes.tsv", sample=SAMPLES),
        seqkit=expand(f"{OUT}/qc/seqkit/{{sample}}.tsv", sample=SAMPLES),
    output:
        tsv=f"{OUT}/reports/defence_systems.tsv",
    params:
        samples=SAMPLES,
        outdir=OUT,
    conda:
        "envs/python.yaml"
    script:
        "scripts/defence_table.py"


rule crispr_table:
    input:
        results=expand(f"{OUT}/crispr/{{sample}}/result.json", sample=SAMPLES),
    output:
        arrays=f"{OUT}/reports/crispr_arrays.tsv",
        spacers=f"{OUT}/reports/crispr_spacers.fna",
    params:
        samples=SAMPLES,
        outdir=OUT,
    conda:
        "envs/python.yaml"
    script:
        "scripts/crispr_table.py"


# ---------------------------------------------------------------------------
# Stage 6. Lineage partition
# ---------------------------------------------------------------------------

rule poppunk:
    """The power simulation clusters standard errors by lineage. Until this has
    run, the clustering in Section 10.1 is a generic partition rather than the
    observed one."""
    input:
        fna=expand(f"{OUT}/genomes/{{sample}}.fna", sample=SAMPLES),
    output:
        clusters=f"{OUT}/poppunk/lineage_assignments.tsv",
    params:
        db=f"{OUT}/poppunk/db",
        fit=f"{OUT}/poppunk/fit",
        model=config["poppunk"]["model"],
        max_k=config["poppunk"]["max_k"],
        outdir=f"{OUT}/poppunk",
    threads: THREADS["poppunk"]
    benchmark:
        f"{OUT}/benchmarks/poppunk/all.tsv"
    conda:
        "envs/poppunk.yaml"
    shell:
        r"""
        mkdir -p {params.outdir}
        rlist={params.outdir}/rfiles.txt
        : > "$rlist"
        for f in {input.fna}; do
            printf '%s\t%s\n' "$(basename "$f" .fna)" "$(readlink -f "$f")" >> "$rlist"
        done
        poppunk --create-db --r-files "$rlist" --output {params.db} \
            --threads {threads}
        poppunk --fit-model {params.model} --ref-db {params.db} \
            --output {params.fit} --K {params.max_k} --threads {threads}
        # The cluster file name carries the fit directory's basename.
        cp {params.fit}/$(basename {params.fit})_clusters.csv \
           {params.outdir}/clusters.csv
        python3 -c "
import csv, sys
rd = csv.DictReader(open('{params.outdir}/clusters.csv'))
w = csv.writer(open('{output.clusters}', 'w'), delimiter='\t')
w.writerow(['sample', 'lineage'])
for r in rd:
    w.writerow([r['Taxon'], r['Cluster']])
"
        test -s {output.clusters}
        """


# ---------------------------------------------------------------------------
# Stage 7. Aggregation and the numbers the protocol is waiting on
# ---------------------------------------------------------------------------

rule concordance:
    input:
        genomad=expand(f"{OUT}/prophage/genomad/{{sample}}.tsv", sample=SAMPLES),
        indep=expand(f"{OUT}/prophage/independent/{{sample}}.tsv", sample=SAMPLES),
    output:
        tsv=f"{OUT}/reports/concordance.tsv",
        pairs=f"{OUT}/reports/concordance_pairs.tsv",
    params:
        samples=SAMPLES,
        outdir=OUT,
        min_overlap=config["concordance"]["min_reciprocal_overlap"],
        max_diff=config["concordance"]["max_coordinate_difference"],
        tool=INDEP,
    conda:
        "envs/python.yaml"
    script:
        "scripts/concordance.py"


rule prophage_counts:
    """The burden table. Internally resolved predictions are the primary
    measure; truncated candidates are carried as separate columns rather than
    folded in (Step 3.4)."""
    input:
        refined=expand(f"{OUT}/prophage/refined/{{sample}}.tsv", sample=SAMPLES),
        pairs=f"{OUT}/reports/concordance_pairs.tsv",
        bands=f"{OUT}/qc/quality_bands.tsv",
        lineage=f"{OUT}/poppunk/lineage_assignments.tsv",
        manifest=config["manifest"],
    output:
        tsv=f"{OUT}/reports/prophage_counts.tsv",
    params:
        samples=SAMPLES,
        outdir=OUT,
    conda:
        "envs/python.yaml"
    script:
        "scripts/prophage_counts.py"


rule dispersion:
    """Fit the negative binomial and report theta with a confidence interval,
    on the primary inference set and on the catalogue block separately. This is
    the single number Section 10.1 is held up on."""
    input:
        counts=f"{OUT}/reports/prophage_counts.tsv",
    output:
        json=f"{OUT}/reports/dispersion_estimates.json",
        tsv=f"{OUT}/reports/dispersion_estimates.tsv",
    conda:
        "envs/stats.yaml"
    script:
        "scripts/dispersion.py"


rule resource_report:
    """Collate every benchmark file into one table: wall time, CPU time and
    peak RSS per tool per genome, plus the per-tool output footprint."""
    input:
        counts=f"{OUT}/reports/prophage_counts.tsv",
        defence=f"{OUT}/reports/defence_systems.tsv",
        crispr=f"{OUT}/reports/crispr_arrays.tsv",
        concordance=f"{OUT}/reports/concordance.tsv",
    output:
        tsv=f"{OUT}/reports/resource_report.tsv",
        per_tool=f"{OUT}/reports/resource_per_tool.tsv",
    params:
        bench=f"{OUT}/benchmarks",
        outdir=OUT,
        samples=SAMPLES,
    conda:
        "envs/python.yaml"
    script:
        "scripts/resource_report.py"


rule provenance:
    """Software versions, database releases and every parameter actually used,
    written from the config rather than transcribed by hand."""
    output:
        tsv=f"{OUT}/reports/provenance.tsv",
    params:
        cfg="config/config.yaml",
        tool=INDEP,
    conda:
        "envs/python.yaml"
    script:
        "scripts/provenance.py"


rule pilot_summary:
    input:
        counts=f"{OUT}/reports/prophage_counts.tsv",
        disp=f"{OUT}/reports/dispersion_estimates.json",
        res=f"{OUT}/reports/resource_report.tsv",
        per_tool=f"{OUT}/reports/resource_per_tool.tsv",
        conc=f"{OUT}/reports/concordance.tsv",
        defence=f"{OUT}/reports/defence_systems.tsv",
        crispr=f"{OUT}/reports/crispr_arrays.tsv",
        bands=f"{OUT}/qc/quality_bands.tsv",
        changed=f"{OUT}/qc/band_changes.tsv",
        prov=f"{OUT}/reports/provenance.tsv",
    output:
        md=f"{OUT}/reports/pilot_summary.md",
    params:
        outdir=OUT,
        n_requested=len(SAMPLES),
    conda:
        "envs/python.yaml"
    script:
        "scripts/pilot_summary.py"
