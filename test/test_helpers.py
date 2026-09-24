#!/usr/bin/env python3
"""One dependency-free smoke test for the workflow helpers."""

import importlib.util
import tempfile
from pathlib import Path
from unittest.mock import patch
from zipfile import ZipFile


ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("helpers", ROOT / "scripts" / "helpers.py")
helpers = importlib.util.module_from_spec(spec)
spec.loader.exec_module(helpers)


def main():
    with tempfile.TemporaryDirectory() as tmp:
        tmp = Path(tmp)
        accession = "GCF_000000001.1"
        manifest = tmp / "assemblies.tsv"
        manifest.write_text(
            "Assembly Accession\tOrganism Taxonomic ID\tOrganism Name\tAssembly Level\n"
            f"{accession}\t3055\tChlamydomonas reinhardtii\tChromosome\n"
        )

        calls = []

        def fake_datasets(command, check):
            assert check
            calls.append(command)
            if command[1:4] == ["download", "genome", "accession"]:
                package = Path(command[command.index("--filename") + 1])
                with ZipFile(package, "w") as archive:
                    if "--dehydrated" in command:
                        archive.writestr("README.md", "dehydrated")
                    else:
                        archive.writestr(
                            f"ncbi_dataset/data/{accession}/{accession}_test_genomic.fna",
                            ">chr1\nACGT\n",
                        )
            else:
                extracted = Path(command[command.index("--directory") + 1])
                genome = extracted / "ncbi_dataset" / "data" / accession / "test_genomic.fna"
                genome.parent.mkdir(parents=True)
                genome.write_text(">chr1\nACGT\n")

        with patch.object(helpers.subprocess, "run", side_effect=fake_datasets):
            direct = tmp / "direct"
            helpers.download_genomes(manifest, tmp / "direct.txt", direct, 1000, 10)
            assert (direct / f"{accession}.fna").read_text() == ">chr1\nACGT\n"
            assert "--dehydrated" not in calls[0]

            dehydrated = tmp / "dehydrated"
            helpers.download_genomes(manifest, tmp / "dehydrated.txt", dehydrated, 1, 7)
            assert (dehydrated / f"{accession}.fna").read_text() == ">chr1\nACGT\n"
            assert "--dehydrated" in calls[1]
            assert calls[2][-3:] == ["--max-workers", "7", "--no-progressbar"]

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
