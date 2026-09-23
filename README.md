# NCBI taxon reference genomes → NAF + taxonomy

This Snakemake 9 workflow replaces the three older download, taxonomy, and compression pipelines without modifying them.

It selects NCBI **reference genomes** below a configured taxonomic group with assembly levels `scaffold`, `chromosome`, or `complete`, downloads each genome, encodes it with `ennaf`, and joins the selected assemblies to the NCBI taxonomy table.

## Configure

```bash
cd pipeline_download_naf_taxonomy/workflows/taxon_genomes
cp config_template.yaml config.yaml
```

Set `taxon` to any NCBI taxonomic name, plus `output_dir` and the absolute path to `ennaf`. Alternatively, set `taxid` to a positive NCBI taxonomy ID; when supplied, it takes precedence over `taxon`:

```yaml
taxon: null
taxid: 3055
```

Before using SLURM, also replace `FILL_ME` in `profiles/slurm/config.yaml` with the cluster account and confirm the `eddy` partition.
The launcher uses Snakemake's Conda integration; Snakemake 9.27 requires Conda 24.7.1 or newer.

## Check without running

```bash
python3 test/test_helpers.py
snakemake --snakefile Snakefile --lint --configfile test/test_config.yaml
snakemake --snakefile Snakefile --dry-run --configfile test/test_config.yaml --cores 1
snakemake --snakefile Snakefile --dag --configfile test/test_config.yaml | dot -Tpdf > test/dag.pdf
```

The test configuration uses a frozen two-row manifest, so dry-run builds the complete DAG without contacting NCBI. `/bin/cat` is only a dry-run placeholder for `ennaf`; it must not be used for a real run.

## Run on the cluster later

```bash
bash run_snakemake.sh config.yaml --dry-run
bash run_snakemake.sh config.yaml
```

## Final outputs

- `results/naf/<assembly_accession>.naf`
- `results/taxonomy/taxonomy.tsv`
- `qc/taxon_genomes.pass`
- `run_metadata.json`

Set `assembly_manifest` to a previously generated `assemblies.tsv` to reproduce the exact accession set without repeating assembly discovery. The configured taxon or taxid is still used to retrieve its taxonomy table. NCBI genome packages contain one genomic FASTA per assembly plus an assembly data report; the taxonomy package supplies `taxonomy_summary.tsv`.

Sources: [NCBI genome downloads](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/download/genome/), [genome package layout](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/data-packages/genome/), and [taxonomy metadata](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/how-tos/taxonomy/taxonomy/).
