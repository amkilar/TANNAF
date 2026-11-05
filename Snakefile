"""Download reference genomes for an NCBI taxon, encode as NAF, and build taxonomy."""

import csv
import os
import sys


CONFIG_DIR = os.path.dirname(os.path.abspath(workflow.configfiles[0])) if workflow.configfiles else os.getcwd()


def config_path(path):
    path = os.path.expandvars(os.path.expanduser(path))
    return path if os.path.isabs(path) else os.path.join(CONFIG_DIR, path)


OUTDIR = config_path(config["output_dir"])
RESULTSDIR = os.path.join(OUTDIR, "results")
LOGDIR = os.path.join(OUTDIR, "logs")
BENCHDIR = os.path.join(OUTDIR, "benchmarks")
QCDIR = os.path.join(OUTDIR, "qc")
MANIFEST_OVERRIDE = config_path(config["assembly_manifest"]) if config.get("assembly_manifest") else None
COMPRESS_GENOMES = config.get("compress_genomes", True)
if isinstance(COMPRESS_GENOMES, str) and COMPRESS_GENOMES.lower() in ("true", "false"):
    COMPRESS_GENOMES = COMPRESS_GENOMES.lower() == "true"
if not isinstance(COMPRESS_GENOMES, bool):
    raise ValueError("compress_genomes must be true or false")
DEHYDRATED_THRESHOLD = int(config.get("dehydrated_threshold", 1000))
REHYDRATE_WORKERS = int(config.get("rehydrate_workers", 10))
if DEHYDRATED_THRESHOLD < 1:
    raise ValueError("dehydrated_threshold must be at least 1")
if not 1 <= REHYDRATE_WORKERS <= 30:
    raise ValueError("rehydrate_workers must be between 1 and 30")
TAXON_NAME = str(config.get("taxon") or "").strip()
TAXID = str(config.get("taxid") or "").strip()
if TAXID:
    if not TAXID.isdigit() or int(TAXID) < 1:
        raise ValueError("taxid must be a positive NCBI taxonomy ID")
    TAXON_SELECTOR = TAXID
elif TAXON_NAME:
    TAXON_SELECTOR = TAXON_NAME
else:
    raise ValueError("Set either taxon or taxid in the configuration")

ASSEMBLY_JSONL = os.path.join(RESULTSDIR, "metadata", "assemblies.jsonl")
ASSEMBLY_TSV = os.path.join(RESULTSDIR, "metadata", "assemblies.tsv")
ACCESSION_LIST = os.path.join(RESULTSDIR, "metadata", "accessions.txt")
FNA_DIR = os.path.join(RESULTSDIR, "fna")
TAXONOMY_ZIP = os.path.join(RESULTSDIR, "taxonomy", "taxonomy.zip")
TAXONOMY_RAW = os.path.join(RESULTSDIR, "taxonomy", "taxonomy_summary.tsv")
TAXONOMY_TABLE = os.path.join(RESULTSDIR, "taxonomy", "taxonomy.tsv")
HELPER = os.path.join(workflow.basedir, "scripts", "helpers.py")
PYTHON = sys.executable
NCBI_ENV = os.path.join(workflow.basedir, "envs", "ncbi-datasets.yaml")
NAF_ENV = os.path.join(workflow.basedir, "envs", "naf.yaml")


def read_accessions(path):
    with open(path, newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if not rows or "Assembly Accession" not in rows[0]:
        raise ValueError(f"Empty or invalid assembly manifest: {path}")
    accessions = [row["Assembly Accession"] for row in rows]
    if len(accessions) != len(set(accessions)):
        raise ValueError(f"Duplicate assembly accessions in: {path}")
    return accessions


if MANIFEST_OVERRIDE:
    ACCESSIONS = read_accessions(MANIFEST_OVERRIDE)

    def manifest_input(_wildcards):
        return MANIFEST_OVERRIDE

    def naf_targets(_wildcards):
        if not COMPRESS_GENOMES:
            return []
        return expand(
            os.path.join(RESULTSDIR, "naf", "{accession}.naf"),
            accession=ACCESSIONS,
        )

    def fasta_targets(_wildcards):
        return expand(
            os.path.join(FNA_DIR, "{accession}.fna"),
            accession=ACCESSIONS,
        )

    def fasta_input(wildcards):
        return os.path.join(FNA_DIR, f"{wildcards.accession}.fna")
else:
    def manifest_input(wildcards):
        return checkpoints.format_assembly_manifest.get(**wildcards).output.tsv

    def naf_targets(wildcards):
        if not COMPRESS_GENOMES:
            return []
        return expand(
            os.path.join(RESULTSDIR, "naf", "{accession}.naf"),
            accession=read_accessions(manifest_input(wildcards)),
        )

    def fasta_targets(wildcards):
        directory = checkpoints.download_genomes.get().output.directory
        return expand(
            os.path.join(directory, "{accession}.fna"),
            accession=read_accessions(manifest_input(wildcards)),
        )

    def fasta_input(wildcards):
        directory = checkpoints.download_genomes.get().output.directory
        return os.path.join(directory, f"{wildcards.accession}.fna")


rule all:
    input:
        os.path.join(OUTDIR, "run_metadata.json"),
        TAXONOMY_TABLE,
        os.path.join(QCDIR, "taxon_genomes.pass"),


if not MANIFEST_OVERRIDE:
    rule query_assemblies:
        """Query reference assemblies for the configured taxon."""
        output:
            jsonl=temp(ASSEMBLY_JSONL),
        log:
            os.path.join(LOGDIR, "query_assemblies.log"),
        benchmark:
            os.path.join(BENCHDIR, "query_assemblies.tsv")
        conda:
            NCBI_ENV
        retries: 2
        resources:
            ncbi=1
        params:
            taxon=TAXON_SELECTOR,
        shell:
            """
            datasets summary genome taxon {params.taxon:q} \
                --reference \
                --assembly-level scaffold,chromosome,complete \
                --as-json-lines > {output.jsonl:q} 2> {log:q}
            test -s {output.jsonl:q}
            """


    checkpoint format_assembly_manifest:
        """Convert the NCBI assembly report to the workflow manifest."""
        input:
            jsonl=ASSEMBLY_JSONL,
        output:
            tsv=ASSEMBLY_TSV,
        log:
            os.path.join(LOGDIR, "format_assembly_manifest.log"),
        conda:
            NCBI_ENV
        retries: 0
        shell:
            """
            dataformat tsv genome \
                --fields accession,organism-tax-id,organism-name,assminfo-level \
                < {input.jsonl:q} > {output.tsv:q} 2> {log:q}
            test "$(wc -l < {output.tsv:q})" -gt 1
            """


if MANIFEST_OVERRIDE:
    rule download_genomes:
        """Download all selected genomes in one NCBI data package."""
        input:
            manifest=manifest_input,
        output:
            accessions=ACCESSION_LIST,
            fastas=expand(os.path.join(FNA_DIR, "{accession}.fna"), accession=ACCESSIONS),
        log:
            os.path.join(LOGDIR, "download_genomes.log"),
        benchmark:
            os.path.join(BENCHDIR, "download_genomes.tsv")
        conda:
            NCBI_ENV
        retries: 2
        resources:
            ncbi=1
        params:
            helper=HELPER,
            python=PYTHON,
            output_dir=lambda wildcards, output: os.path.dirname(output.fastas[0]),
            threshold=DEHYDRATED_THRESHOLD,
            workers=REHYDRATE_WORKERS,
        shell:
            """
            {params.python:q} {params.helper:q} download-genomes \
                --manifest {input.manifest:q} \
                --accessions {output.accessions:q} \
                --output-dir {params.output_dir:q} \
                --dehydrated-threshold {params.threshold} \
                --rehydrate-workers {params.workers} > {log:q} 2>&1
            """
else:
    checkpoint download_genomes:
        """Download all selected genomes, rehydrating large collections."""
        input:
            manifest=manifest_input,
        output:
            accessions=ACCESSION_LIST,
            directory=directory(FNA_DIR),
        log:
            os.path.join(LOGDIR, "download_genomes.log"),
        benchmark:
            os.path.join(BENCHDIR, "download_genomes.tsv")
        conda:
            NCBI_ENV
        retries: 2
        resources:
            ncbi=1
        params:
            helper=HELPER,
            python=PYTHON,
            output_dir=lambda wildcards, output: output.directory,
            threshold=DEHYDRATED_THRESHOLD,
            workers=REHYDRATE_WORKERS,
        shell:
            """
            {params.python:q} {params.helper:q} download-genomes \
                --manifest {input.manifest:q} \
                --accessions {output.accessions:q} \
                --output-dir {params.output_dir:q} \
                --dehydrated-threshold {params.threshold} \
                --rehydrate-workers {params.workers} > {log:q} 2>&1
            """


rule compress_genome:
    """Encode one genomic FASTA in NAF format."""
    input:
        fasta=fasta_input,
    output:
        naf=protected(os.path.join(RESULTSDIR, "naf", "{accession}.naf")),
    log:
        os.path.join(LOGDIR, "compress_genome_{accession}.log"),
    benchmark:
        os.path.join(BENCHDIR, "compress_genome_{accession}.tsv")
    conda:
        NAF_ENV
    retries: 1
    threads: 1
    shell:
        """
        ennaf {input.fasta:q} > {output.naf:q} 2> {log:q}
        test -s {output.naf:q}
        """


rule download_taxonomy:
    """Download the NCBI taxonomy package for the taxon and all descendants."""
    output:
        package=temp(TAXONOMY_ZIP),
    log:
        os.path.join(LOGDIR, "download_taxonomy.log"),
    benchmark:
        os.path.join(BENCHDIR, "download_taxonomy.tsv")
    conda:
        NCBI_ENV
    retries: 2
    resources:
        ncbi=1
    params:
        taxon=TAXON_SELECTOR,
    shell:
        """
        datasets download taxonomy taxon {params.taxon:q} \
            --children \
            --filename {output.package:q} \
            --no-progressbar > {log:q} 2>&1
        test -s {output.package:q}
        """


rule extract_taxonomy:
    """Extract taxonomy_summary.tsv from the NCBI package."""
    input:
        package=TAXONOMY_ZIP,
    output:
        table=temp(TAXONOMY_RAW),
    log:
        os.path.join(LOGDIR, "extract_taxonomy.log"),
    conda:
        NCBI_ENV
    retries: 1
    params:
        helper=HELPER,
        python=PYTHON,
    shell:
        """
        {params.python:q} {params.helper:q} extract-taxonomy \
            --package {input.package:q} \
            --output {output.table:q} > {log:q} 2>&1
        """


rule build_taxonomy_table:
    """Join selected assemblies to their NCBI taxonomy lineages."""
    input:
        manifest=manifest_input,
        taxonomy=TAXONOMY_RAW,
    output:
        table=protected(TAXONOMY_TABLE),
    log:
        os.path.join(LOGDIR, "build_taxonomy_table.log"),
    conda:
        NCBI_ENV
    retries: 0
    params:
        helper=HELPER,
        python=PYTHON,
    shell:
        """
        {params.python:q} {params.helper:q} build-taxonomy \
            --manifest {input.manifest:q} \
            --taxonomy {input.taxonomy:q} \
            --output {output.table:q} > {log:q} 2>&1
        """


rule validate_outputs:
    """Require one FASTA and taxonomy row per assembly, plus NAF when enabled."""
    input:
        manifest=manifest_input,
        taxonomy=TAXONOMY_TABLE,
        fastas=fasta_targets,
        nafs=naf_targets,
    output:
        passed=os.path.join(QCDIR, "taxon_genomes.pass"),
    log:
        os.path.join(LOGDIR, "validate_outputs.log"),
    conda:
        NCBI_ENV
    retries: 0
    params:
        helper=HELPER,
        python=PYTHON,
        skip_naf="" if COMPRESS_GENOMES else "--skip-naf",
    shell:
        """
        {params.python:q} {params.helper:q} validate \
            --manifest {input.manifest:q} \
            --taxonomy {input.taxonomy:q} \
            --fastas {input.fastas:q} \
            --nafs {input.nafs:q} \
            {params.skip_naf} \
            --output {output.passed:q} > {log:q} 2>&1
        """


rule write_run_metadata:
    """Record the frozen inputs and tool settings for this run."""
    input:
        manifest=manifest_input,
    output:
        os.path.join(OUTDIR, "run_metadata.json"),
    log:
        os.path.join(LOGDIR, "write_run_metadata.log"),
    conda:
        NCBI_ENV
    retries: 0
    params:
        helper=HELPER,
        python=PYTHON,
        ennaf="ennaf (Conda: naf=1.3.0)",
        compress_genomes=str(COMPRESS_GENOMES).lower(),
        taxon_name=TAXON_NAME,
        taxid=TAXID,
        taxon_selector=TAXON_SELECTOR,
        config_file=workflow.configfiles[0] if workflow.configfiles else "unknown",
    shell:
        """
        {params.python:q} {params.helper:q} write-metadata \
            --manifest {input.manifest:q} \
            --ennaf {params.ennaf:q} \
            --compress-genomes {params.compress_genomes} \
            --taxon-name={params.taxon_name:q} \
            --taxid={params.taxid:q} \
            --taxon-selector {params.taxon_selector:q} \
            --config-file {params.config_file:q} \
            --output {output:q} > {log:q} 2>&1
        """
