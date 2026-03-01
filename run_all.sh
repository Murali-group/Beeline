#!/bin/bash

set -e
BASEDIR="$(dirname "$(readlink -f "$0")")"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate BEELINE

CONFIGS=(
    "$BASEDIR/config-files/Curated/GSD.yaml"
    "$BASEDIR/config-files/Curated/HSC.yaml"
    "$BASEDIR/config-files/Curated/mCAD.yaml"
    "$BASEDIR/config-files/Curated/VSC.yaml"
    "$BASEDIR/config-files/Synthetic/dyn-LL.yaml"
    "$BASEDIR/config-files/Synthetic/dyn-TF.yaml"
)

for config in "${CONFIGS[@]}"; do
    echo "========================================"
    echo "Running: $config"
    echo "========================================"
    python "$BASEDIR/BLRunner.py" -c "$config"
done

echo "All configs complete."
