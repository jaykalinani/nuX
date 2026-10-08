#!/bin/bash

set -euxo pipefail

export NUXSPACE="$PWD"
export WORKSPACE="$PWD/../workspace"
cd "$WORKSPACE/Cactus"

NUX_SCRIPTS="$NUXSPACE/scripts"
cp "$NUX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.cfg" \
    simfactory/mdb/optionlists/
cp "$NUX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.ini" \
    simfactory/mdb/machines/
cp "$NUX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.run" \
    simfactory/mdb/runscripts/
cp "$NUX_SCRIPTS/actions-$ACCELERATOR-$REAL_PRECISION.sub" \
    simfactory/mdb/submitscripts/
cp "$NUX_SCRIPTS/defs.local.ini" simfactory/etc/
cp "$NUX_SCRIPTS/nux.th" .
# Keep unit-test infrastructure out of production executables.  CI uses a
# private copy of the production ThornList and adds the test thorn there.
printf '\nnuX/nuX_Tests\n' >>nux.th

if command -v ccache >/dev/null 2>&1; then
    export CCACHE_DIR="${CCACHE_DIR:-$NUXSPACE/.ccache}"
    ccache --max-size=2G
    ccache --zero-stats || true
    sed -i -e 's/^CC = /CC = ccache /' -e 's/^CXX = /CXX = ccache /' \
        "simfactory/mdb/optionlists/actions-$ACCELERATOR-$REAL_PRECISION.cfg"
fi

set +e
time ./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    build -j "$(nproc)" sim 2>&1 | tee build.log
build_status=${PIPESTATUS[0]}
set -e

if ((build_status != 0)); then
    # GitHub does not expose public Actions logs without authentication. Keep
    # the useful compiler diagnostics visible in the job annotation instead.
    build_errors=$(
        { grep -E '(^|: )(fatal )?error:|undefined (reference|symbol)|No rule to make target|make(\[[0-9]+\])?: \*\*\*' \
            build.log || true; } | tail -n 12 | tr '\n' ' '
    )
    printf '::error title=Cactus build failed (%s)::%s\n' \
        "$ACCELERATOR" "${build_errors:-See the Build Cactus log for details.}"
    exit "$build_status"
fi

test -x exe/cactus_sim
command -v ccache >/dev/null 2>&1 && ccache --show-stats || true
