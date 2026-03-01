#!/bin/bash

set -e
BASEDIR="$(dirname "$(readlink -f "$0")")"

source "$HOME/miniconda3/etc/profile.d/conda.sh"
conda activate BEELINE

CONFIGS=(
    "$BASEDIR/config-files/Synthetic/dyn-BF.yaml"
    "$BASEDIR/config-files/Synthetic/dyn-BFC.yaml"
    "$BASEDIR/config-files/Synthetic/dyn-CY.yaml"
    "$BASEDIR/config-files/Synthetic/dyn-LI.yaml"
)

for config in "${CONFIGS[@]}"; do
    echo "========================================"
    echo "Running: $config"
    echo "========================================"
    python "$BASEDIR/BLRunner.py" -c "$config"
done

echo "All synthetic datasets complete."
