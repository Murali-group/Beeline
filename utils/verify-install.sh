#!/usr/bin/env bash
# Verify the tutorial install by exercising each of the four conda envs.
# Prints a pass/fail line per check and exits non-zero if anything failed.

set -u

PASS=0
FAIL=0
FAILED_CHECKS=()

check() {
    local name="$1"; shift
    if "$@" >/dev/null 2>&1; then
        printf "  \e[32mPASS\e[0m  %s\n" "$name"
        PASS=$((PASS + 1))
    else
        printf "  \e[31mFAIL\e[0m  %s\n" "$name"
        FAIL=$((FAIL + 1))
        FAILED_CHECKS+=("$name")
    fi
}

echo "==> conda environments"
check "conda on PATH"                command -v conda
check "BEELINE env exists"           bash -c "conda env list | grep -qE '^BEELINE\s'"
check "BEELINE-PPCOR env exists"     bash -c "conda env list | grep -qE '^BEELINE-PPCOR\s'"
check "BEELINE-LEAP env exists"      bash -c "conda env list | grep -qE '^BEELINE-LEAP\s'"
check "BEELINE-GRNBOOST2 env exists" bash -c "conda env list | grep -qE '^BEELINE-GRNBOOST2\s'"

echo "==> BEELINE pipeline imports"
check "pandas + yaml + networkx + sklearn" \
    conda run --no-capture-output -n BEELINE python -c \
        "import pandas, yaml, networkx, sklearn"
check "BLRunner.py --help" \
    conda run --no-capture-output -n BEELINE python BLRunner.py --help

echo "==> Algorithm envs"
check "BEELINE-PPCOR loads r-ppcor" \
    conda run --no-capture-output -n BEELINE-PPCOR Rscript -e 'library(ppcor)'
check "BEELINE-LEAP loads r-leap" \
    conda run --no-capture-output -n BEELINE-LEAP Rscript -e 'library(LEAP)'
check "BEELINE-GRNBOOST2 loads arboreto" \
    conda run --no-capture-output -n BEELINE-GRNBOOST2 python -c \
        "from arboreto.algo import grnboost2"

echo "==> Tutorial data"
check "Synthetic input dir present" \
    test -d inputs/Synthetic/dyn-BF/dyn-BF-500
check "All 6 networks have a ground-truth file" \
    bash -c 'for n in dyn-BF dyn-BFC dyn-CY dyn-LI dyn-LL dyn-TF; do
        test -f "inputs/Synthetic/$n/$n-500/GroundTruthNetwork.csv" || exit 1
    done'
check "tutorial.yaml parses as YAML" \
    conda run --no-capture-output -n BEELINE python -c \
        "import yaml; yaml.safe_load(open('config-files/tutorial.yaml'))"

echo
echo "------------------------------------"
printf "  %d passed, %d failed\n" "$PASS" "$FAIL"
echo "------------------------------------"

if [ "$FAIL" -gt 0 ]; then
    echo
    echo "Failures:"
    for c in "${FAILED_CHECKS[@]}"; do
        echo "  - $c"
    done
    exit 1
fi
