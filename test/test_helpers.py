#!/usr/bin/env python3
"""One dependency-free smoke test for the workflow helpers."""

import importlib.util
import tempfile
from pathlib import Path
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("helpers", ROOT / "scripts" / "helpers.py")
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        accession = "GCF_000000001.1"
        genome_zip = tmp / "genome.zip"
        with ZipFile(genome_zip, "w") as archive:
            archive.writestr(
                f"ncbi_dataset/data/{accession}/{accession}_test_genomic.fna",
                ">chr1\nACGT\n",
            )
        fasta = tmp / "fna" / f"{accession}.fna"
        helpers.extract_genome(genome_zip, accession, fasta)
        assert fasta.read_text() == ">chr1\nACGT\n"

        taxonomy_zip = tmp / "taxonomy.zip"
        taxonomy_text = (
            "Taxid\tDomain/Realm name\tKingdom name\tPhylum name\tClass name\t"
            "Order name\tFamily name\tGenus name\tSpecies name\n"
            "3055\tEukaryota\tViridiplantae\tChlorophyta\tChlorophyceae\t"
            "Chlamydomonadales\tChlamydomonadaceae\tChlamydomonas\t"
            "Chlamydomonas reinhardtii\n"
        )
        with ZipFile(taxonomy_zip, "w") as archive:
            archive.writestr("ncbi_dataset/data/taxonomy_summary.tsv", taxonomy_text)
        taxonomy = tmp / "taxonomy_summary.tsv"
        helpers.extract_taxonomy(taxonomy_zip, taxonomy)

        manifest = tmp / "assemblies.tsv"
        manifest.write_text(
            "Assembly Accession\tOrganism Taxonomic ID\tOrganism Name\tAssembly Level\n"
            f"{accession}\t3055\tChlamydomonas reinhardtii\tChromosome\n"
        )
        table = tmp / "taxonomy.tsv"
        helpers.build_taxonomy(manifest, taxonomy, table)
        assert "Chlorophyta" in table.read_text()

        naf = tmp / f"{accession}.naf"
        naf.write_bytes(b"NAF-test")
        passed = tmp / "qc.pass"
        helpers.validate(manifest, table, [naf], passed)
        assert passed.read_text() == "PASS\tassemblies=1\n"

        metadata = tmp / "run_metadata.json"
        helpers.write_metadata(
            manifest,
            "/tools/ennaf",
            "Chlorophyta",
            "",
            "Chlorophyta",
            "config.yaml",
            metadata,
        )
        assert '"assembly_count": 1' in metadata.read_text()
        assert '"ncbi_taxon_selector": "Chlorophyta"' in metadata.read_text()

    print("PASS: helper smoke test")


if __name__ == "__main__":
    main()
