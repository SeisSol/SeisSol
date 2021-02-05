#!/bin/bash

set -euo pipefail

export CUDA_PATH=$(dirname $(dirname $(readlink $(which nvcc))))
export CUDA_HOME=$CUDA_PATH
export CMAKE_PREFIX_PATH=$CMAKE_PREFIX_PATH:$CUDA_HOME
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$CUDA_HOME/lib64/stubs
export LIBRARY_PATH=$LIBRARY_PATH:$CUDA_HOME/lib64/stubs

# install gemmforge
pip3 install git+https://github.com/ravil-mobile/gemmforge.git

# prepare
cd ./tests/convergence_elastic && ./generateCubes.sh && cd ../..

for precision in double single; do
    mkdir -p ./SeisSol/_build_${precision} && cd ./SeisSol/_build_${precision}
    cp ../../tests/convergence_elastic/* .

    # compile
    cmake .. -DCMAKE_BUILD_TYPE=Release -DDEVICE_ARCH=${GPU_VENDOR} -DDEVICE_SUB_ARCH=${GPU_MODEL}  -DHOST_ARCH=${HOST} -DPRECISION=${precision}
    make -j

    # run proxy
    ./SeisSol_proxy_Release_?${GPU_VENDOR}_?_elastic 100000 50 all
    
    # run convergence test
    mpirun -n 1 ./SeisSol_Release_?nvidia_?_elastic ./parameters.par
    cd ../..
done
set +u
