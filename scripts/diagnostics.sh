#!/bin/bash

set -euxo pipefail

export NUXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
cd "$WORKSPACE/Cactus"

export LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH:-}"
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

SIMULATION_NAME=NuXDiagnostics
PARAMETER_FILE="$NUXSPACE/diagnostics/nuX_M1_cpu.par"

time ./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    create-run "$SIMULATION_NAME" --cores 1 --ppn-used 1 --num-threads 2 \
    --parfile "$PARAMETER_FILE"

DIAGNOSTICS_OUTPUT_DIR="$(./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    get-output-dir "$SIMULATION_NAME")"

if test -n "${GITHUB_ENV:-}"; then
    echo "DIAGNOSTICS_OUTPUT_DIR=$DIAGNOSTICS_OUTPUT_DIR" >>"$GITHUB_ENV"
fi

NORM_DIR="$(find "$DIAGNOSTICS_OUTPUT_DIR" -type d -name norms -print -quit)"

# Optional diagnostic fields are temporary observability aids, not correctness
# or regression oracles.  Report whatever was produced and let the dedicated
# nuX_Tests thorn provide the numerical pass/fail decision.
if test -n "$NORM_DIR" && test -d "$NORM_DIR"; then
    find "$NORM_DIR" -maxdepth 1 -type f -name '*.tsv' -print | sort
else
    echo "No diagnostic norm directory was produced" >&2
fi
