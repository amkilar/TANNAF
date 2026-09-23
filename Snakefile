"""Download reference genomes for an NCBI taxon, encode as NAF, and build taxonomy."""

import csv
import os


OUTDIR = config["output_dir"]
RESULTSDIR = os.path.join(OUTDIR, "results")
LOGDIR = os.path.join(OUTDIR, "logs")
BENCHDIR = os.path.join(OUTDIR, "benchmarks")
QCDIR = os.path.join(OUTDIR, "qc")
ENNAF = config["ennaf"]
MANIFEST_OVERRIDE = config.get("assembly_manifest")
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
TAXONOMY_ZIP = os.path.join(RESULTSDIR, "taxonomy", "taxonomy.zip")
TAXONOMY_RAW = os.path.join(RESULTSDIR, "taxonomy", "taxonomy_summary.tsv")
TAXONOMY_TABLE = os.path.join(RESULTSDIR, "taxonomy", "taxonomy.tsv")
HELPER = os.path.join(workflow.basedir, "scripts", "helpers.py")
NCBI_ENV = os.path.join(workflow.basedir, "envs", "ncbi-datasets.yaml")


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
        return expand(
            os.path.join(RESULTSDIR, "naf", "{accession}.naf"),
            accession=ACCESSIONS,
        )
else:
    def manifest_input(wildcards):
        return checkpoints.format_assembly_manifest.get(**wildcards).output.tsv

    def naf_targets(wildcards):
        return expand(
            os.path.join(RESULTSDIR, "naf", "{accession}.naf"),
            accession=read_accessions(manifest_input(wildcards)),
        )


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


rule download_genome:
    """Download one NCBI genome data package."""
    input:
        manifest=manifest_input,
    output:
        package=temp(os.path.join(RESULTSDIR, "packages", "{accession}.zip")),
    log:
        os.path.join(LOGDIR, "download_genome_{accession}.log"),
    benchmark:
        os.path.join(BENCHDIR, "download_genome_{accession}.tsv")
    conda:
        NCBI_ENV
    retries: 2
    resources:
        ncbi=1
    shell:
        """
        datasets download genome accession {wildcards.accession:q} \
            --include genome \
            --filename {output.package:q} \
            --no-progressbar > {log:q} 2>&1
        test -s {output.package:q}
        """


rule extract_genome:
    """Extract and normalize one genomic FASTA filename."""
    input:
        package=os.path.join(RESULTSDIR, "packages", "{accession}.zip"),
    output:
        fasta=temp(os.path.join(RESULTSDIR, "fna", "{accession}.fna")),
    log:
        os.path.join(LOGDIR, "extract_genome_{accession}.log"),
    conda:
        NCBI_ENV
    retries: 1
    params:
        helper=HELPER,
    shell:
        """
        python {params.helper:q} extract-genome \
            --package {input.package:q} \
            --accession {wildcards.accession:q} \
            --output {output.fasta:q} > {log:q} 2>&1
        """


rule compress_genome:
    """Encode one genomic FASTA in NAF format."""
    input:
        fasta=os.path.join(RESULTSDIR, "fna", "{accession}.fna"),
    output:
        naf=protected(os.path.join(RESULTSDIR, "naf", "{accession}.naf")),
    log:
        os.path.join(LOGDIR, "compress_genome_{accession}.log"),
    benchmark:
        os.path.join(BENCHDIR, "compress_genome_{accession}.tsv")
    conda:
        NCBI_ENV
    retries: 1
    threads: 1
    params:
        ennaf=ENNAF,
    shell:
        """
        {params.ennaf:q} {input.fasta:q} > {output.naf:q} 2> {log:q}
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
    shell:
        """
        python {params.helper:q} extract-taxonomy \
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
    shell:
        """
        python {params.helper:q} build-taxonomy \
            --manifest {input.manifest:q} \
            --taxonomy {input.taxonomy:q} \
            --output {output.table:q} > {log:q} 2>&1
        """


rule validate_outputs:
    """Require one non-empty NAF and one taxonomy row per selected assembly."""
    input:
        manifest=manifest_input,
        taxonomy=TAXONOMY_TABLE,
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
    shell:
        """
        python {params.helper:q} validate \
            --manifest {input.manifest:q} \
            --taxonomy {input.taxonomy:q} \
            --nafs {input.nafs:q} \
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
        ennaf=ENNAF,
        taxon_name=TAXON_NAME,
        taxid=TAXID,
        taxon_selector=TAXON_SELECTOR,
        config_file=workflow.configfiles[0] if workflow.configfiles else "unknown",
    shell:
        """
        python {params.helper:q} write-metadata \
            --manifest {input.manifest:q} \
            --ennaf {params.ennaf:q} \
            --taxon-name {params.taxon_name:q} \
            --taxid {params.taxid:q} \
            --taxon-selector {params.taxon_selector:q} \
            --config-file {params.config_file:q} \
            --output {output:q} > {log:q} 2>&1
        """
