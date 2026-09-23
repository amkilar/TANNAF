#!/usr/bin/env python3
"""Small, dependency-free helpers for the taxon genome workflow."""

import argparse
import csv
import json
import shutil
from datetime import datetime, timezone
from pathlib import Path
from zipfile import ZipFile


MANIFEST_COLUMNS = (
    "Assembly Accession",
    "Organism Taxonomic ID",
    "Organism Name",
    "Assembly Level",
)
OUTPUT_COLUMNS = (
    "assembly_accession",
    "organism_tax_id",
    "organism_name",
    "assembly_level",
    "domain",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
)


def log(message):
    print(f"[{datetime.now().isoformat(timespec='seconds')}] {message}", flush=True)


def rows(path):
    with open(path, newline="") as handle:
        yield from csv.DictReader(handle, delimiter="\t")


def require_columns(fieldnames, required, path):
    missing = set(required) - set(fieldnames or ())
    if missing:
        raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")


def extract_genome(package, accession, output):
    suffix = "_genomic.fna"
    prefix = f"ncbi_dataset/data/{accession}/"
    with ZipFile(package) as archive:
        matches = [n for n in archive.namelist() if n.startswith(prefix) and n.endswith(suffix)]
        if len(matches) != 1:
            raise ValueError(f"Expected one genomic FASTA for {accession}, found {len(matches)}")
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        with archive.open(matches[0]) as source, open(output, "wb") as target:
            shutil.copyfileobj(source, target)
    if Path(output).stat().st_size == 0:
        raise ValueError(f"Extracted an empty FASTA: {output}")


def extract_taxonomy(package, output):
    with ZipFile(package) as archive:
        matches = [n for n in archive.namelist() if Path(n).name == "taxonomy_summary.tsv"]
        if len(matches) != 1:
            raise ValueError(f"Expected one taxonomy_summary.tsv, found {len(matches)}")
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        with archive.open(matches[0]) as source, open(output, "wb") as target:
            shutil.copyfileobj(source, target)


def build_taxonomy(manifest, taxonomy, output):
    with open(manifest, newline="") as handle:
        manifest_reader = csv.DictReader(handle, delimiter="\t")
        require_columns(manifest_reader.fieldnames, MANIFEST_COLUMNS, manifest)
        assemblies = list(manifest_reader)
    if not assemblies:
        raise ValueError(f"Assembly manifest is empty: {manifest}")

    with open(taxonomy, newline="") as handle:
        taxonomy_reader = csv.DictReader(handle, delimiter="\t")
        require_columns(taxonomy_reader.fieldnames, ("Taxid",), taxonomy)
        lineages = {row["Taxid"]: row for row in taxonomy_reader}

    rank_columns = {
        "domain": "Domain/Realm name",
        "kingdom": "Kingdom name",
        "phylum": "Phylum name",
        "class": "Class name",
        "order": "Order name",
        "family": "Family name",
        "genus": "Genus name",
        "species": "Species name",
    }
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    with open(output, "w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_COLUMNS, delimiter="\t")
        writer.writeheader()
        for assembly in assemblies:
            taxid = assembly["Organism Taxonomic ID"]
            lineage = lineages.get(taxid, {})
            row = {
                "assembly_accession": assembly["Assembly Accession"],
                "organism_tax_id": taxid,
                "organism_name": assembly["Organism Name"],
                "assembly_level": assembly["Assembly Level"],
            }
            row.update({name: lineage.get(source) or "NA" for name, source in rank_columns.items()})
            writer.writerow(row)


def validate(manifest, taxonomy, nafs, output):
    expected = {row["Assembly Accession"] for row in rows(manifest)}
    observed_taxonomy = {row["assembly_accession"] for row in rows(taxonomy)}
    observed_nafs = {Path(path).stem for path in nafs}
    empty = [path for path in nafs if not Path(path).is_file() or Path(path).stat().st_size == 0]
    if not expected or expected != observed_taxonomy or expected != observed_nafs or empty:
        raise ValueError(
            f"Output validation failed: expected={sorted(expected)}, "
            f"taxonomy={sorted(observed_taxonomy)}, naf={sorted(observed_nafs)}, empty={empty}"
        )
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(f"PASS\tassemblies={len(expected)}\n")


def write_metadata(manifest, ennaf, taxon_name, taxid, taxon_selector, config_file, output):
    accessions = [row["Assembly Accession"] for row in rows(manifest)]
    metadata = {
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "taxon": taxon_name or None,
        "taxid": taxid or None,
        "ncbi_taxon_selector": taxon_selector,
        "reference_only": True,
        "assembly_levels": ["scaffold", "chromosome", "complete"],
        "assembly_manifest": str(Path(manifest).resolve()),
        "assembly_count": len(accessions),
        "ennaf": ennaf,
        "config_file": config_file,
    }
    Path(output).parent.mkdir(parents=True, exist_ok=True)
    Path(output).write_text(json.dumps(metadata, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    genome = subparsers.add_parser("extract-genome")
    genome.add_argument("--package", required=True)
    genome.add_argument("--accession", required=True)
    genome.add_argument("--output", required=True)

    tax = subparsers.add_parser("extract-taxonomy")
    tax.add_argument("--package", required=True)
    tax.add_argument("--output", required=True)

    build = subparsers.add_parser("build-taxonomy")
    build.add_argument("--manifest", required=True)
    build.add_argument("--taxonomy", required=True)
    build.add_argument("--output", required=True)

    check = subparsers.add_parser("validate")
    check.add_argument("--manifest", required=True)
    check.add_argument("--taxonomy", required=True)
    check.add_argument("--nafs", nargs="+", required=True)
    check.add_argument("--output", required=True)

    metadata = subparsers.add_parser("write-metadata")
    metadata.add_argument("--manifest", required=True)
    metadata.add_argument("--ennaf", required=True)
    metadata.add_argument("--taxon-name", default="")
    metadata.add_argument("--taxid", default="")
    metadata.add_argument("--taxon-selector", required=True)
    metadata.add_argument("--config-file", required=True)
    metadata.add_argument("--output", required=True)

    args = parser.parse_args()
    log(f"=== {args.command} ===")
    if args.command == "extract-genome":
        extract_genome(args.package, args.accession, args.output)
    elif args.command == "extract-taxonomy":
        extract_taxonomy(args.package, args.output)
    elif args.command == "build-taxonomy":
        build_taxonomy(args.manifest, args.taxonomy, args.output)
    elif args.command == "validate":
        validate(args.manifest, args.taxonomy, args.nafs, args.output)
    else:
        write_metadata(
            args.manifest,
            args.ennaf,
            args.taxon_name,
            args.taxid,
            args.taxon_selector,
            args.config_file,
            args.output,
        )
    log(f"=== DONE: {args.command} completed successfully ===")


if __name__ == "__main__":
    main()
