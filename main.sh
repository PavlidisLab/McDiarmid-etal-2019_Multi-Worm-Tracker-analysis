#!/usr/bin/env bash
set -eu

## Usage:
## ./main.sh
## This script effectively runs main.R but retains the logs in the $LOGS_DIR directory.

LOGS_DIR="logs/"
mkdir -p "${LOGS_DIR}"

echo "Running main.R"
Rscript main.R \
    2> "${LOGS_DIR}/main.sh.stderr" \
    | tee > "${LOGS_DIR}/main.sh.stdout"
