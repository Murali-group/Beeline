#!/bin/bash

# Build the BEELINE Sphinx documentation.
# Output is written to docs/build/html/index.html.
# Run from the repository root or from the utils/ directory.

set -e

SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
REPO_ROOT="$(dirname "$SCRIPT_DIR")"

eval "$(conda shell.bash hook)"
conda activate BEELINE

sphinx-build -b html \
    "$REPO_ROOT/docs/source" \
    "$REPO_ROOT/docs/build/html"

echo "Documentation built: $REPO_ROOT/docs/build/html/index.html"