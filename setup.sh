#!/bin/sh

(
cd $(dirname ${BASH_SOURCE[0]})
pushd MELA
./setup.sh "$@" || exit 1
popd
)
