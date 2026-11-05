#!/usr/bin/env bash
set -euo pipefail
unset SLURM_MEM_PER_NODE
export CONDA_REPODATA_USE_SHARDS=false
export CONDA_CHANNEL_PRIORITY=strict

if [[ -n "${CONDA_PREFIX:-}" && -x "$CONDA_PREFIX/bin/conda" ]]; then
    export CONDA_EXE="$CONDA_PREFIX/bin/conda"
fi

if [[ $# -lt 1 ]]; then
    echo "Usage: $0 CONFIG.yaml [snakemake arguments]" >&2
    exit 2
fi

CONFIGFILE="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
shift
WORKFLOW_DIR="$(cd "$(dirname "$0")" && pwd)"

snakemake \
    --snakefile "$WORKFLOW_DIR/Snakefile" \
    --configfile "$CONFIGFILE" \
    --workflow-profile "$WORKFLOW_DIR/profiles/slurm" \
    --use-conda \
    --rerun-incomplete \
    "$@"
