# Changelog

## 2026-09-24

- Replaced per-assembly downloads with one batch package and automatic dehydration/rehydration for large collections.
- Reduced a 100-genome run from roughly 300 SLURM jobs to roughly 100 compression jobs plus shared setup steps.

## 2026-09-23

- Generalized taxon selection to accept any NCBI taxonomic name or taxonomy ID.
- Reused the original cluster-local `$HOME/software/naf/naf/ennaf/ennaf` encoder path.
- Made the launcher use the active environment's Conda and disable its SQLite-backed sharded metadata cache for shared filesystems.
- Added the unified reference-genome download, NAF encoding, and taxonomy workflow.
