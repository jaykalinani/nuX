#!/bin/bash

set -euxo pipefail

export NUXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
cd "$WORKSPACE/Cactus"

export LD_LIBRARY_PATH="/usr/local/lib:${LD_LIBRARY_PATH:-}"
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

time ./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    create-run NuXTest01 --cores 1 --num-threads 2 \
    --testsuite --select-tests=nuX

TEST_OUTPUT_DIR="$(./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    get-output-dir NuXTest01)/TEST/sim"

if test -n "${GITHUB_ENV:-}"; then
    echo "TEST_OUTPUT_DIR=$TEST_OUTPUT_DIR" >>"$GITHUB_ENV"
fi

cat "$TEST_OUTPUT_DIR/summary.log"
# A zero-test run also reports zero failures, so require the complete nuX
# arrangement suite to have run and passed.  This checks ordinary Cactus
# pass/fail tests and does not depend on optional diagnostic grid functions.
grep -q '^    Total available tests    -> 5$' \
    "$TEST_OUTPUT_DIR/summary.log"
grep -q '^    Number of tests passed   -> 5$' \
    "$TEST_OUTPUT_DIR/summary.log"
grep -q '^    Number failed            -> 0$' \
    "$TEST_OUTPUT_DIR/summary.log"
