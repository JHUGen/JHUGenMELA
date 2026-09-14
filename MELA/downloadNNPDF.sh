#!/bin/bash

set -euo pipefail

(
cd $(dirname $0)

if [ -f ./data/Pdfdata/NNPDF30_lo_as_0130.LHgrid ]; then
    exit 0
else
    wget --no-check-certificate http://spin.pha.jhu.edu/Generator/NNPDF30_lo_as_0130.LHgrid.tgz
    tar -zxvf ./NNPDF30_lo_as_0130.LHgrid.tgz
    rm ./NNPDF30_lo_as_0130.LHgrid.tgz
    mv ./NNPDF30_lo_as_0130.LHgrid ./data/Pdfdata/NNPDF30_lo_as_0130.LHgrid
fi
)
