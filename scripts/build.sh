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

time ./simfactory/bin/sim \
    --machine="actions-$ACCELERATOR-$REAL_PRECISION" \
    build -j "$(nproc)" sim 2>&1 | tee build.log

test -x exe/cactus_sim
command -v ccache >/dev/null 2>&1 && ccache --show-stats || true
