#!/bin/bash
set -euo pipefail

export CTEST_OUTPUT_ON_FAILURE=1
export OMPI_MCA_rmaps_base_oversubscribe=1
export equation=elastic
export build_type=Release
export precision=double
dirname="./SeisSol/build_${equation}_${precision}_${build_type}"
mkdir -p $dirname && cd $dirname
cmake .. -DNETCDF=OFF -DMETIS=ON -DCOMMTHREAD=ON -DASAGI=OFF -DHDF5=ON \
                -DCMAKE_BUILD_TYPE=$build_type -DTESTING=ON -DLOG_LEVEL=warning \
                -DLOG_LEVEL_MASTER=info -DHOST_ARCH=$HOST_ARCH -DPRECISION=$precision \
                -DEQUATIONS=$equation -DNUMBER_OF_MECHANISMS=$mechanisms
make -j $NUM_PROC
make test
cd ../..
set +u
