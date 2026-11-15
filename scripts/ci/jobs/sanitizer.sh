#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The sanitizer jobs (.github/workflows/sanitizers.yml):
#
#   scripts/ci/jobs/sanitizer.sh build   configures and builds SeisSol with the sanitizers in build/
#   scripts/ci/jobs/sanitizer.sh test    runs its tests there, as the CI does
#
# The setup comes from the environment, as the matrix of the workflow gives it:
#
#   CC, CXX             the compilers (gcc/g++, clang/clang++ for the thread sanitizer)
#   SANITIZERS          the sanitizers, e.g. "address;undefined"
#   CXX_FLAGS           further flags for C++, e.g. -D_GLIBCXX_ASSERTIONS
#   CONFIGS             the configurations, e.g. ci-sanitize
#   DR_QUAD_RULE        the quadrature rule of the faults (dunavant, stroud)
#   TESTING_GENERATED   whether to test the generated kernels as well (ON, OFF)
#   THREADS             true to test with a communication thread and asynchronous output (thread)
#   CMAKE_PREFIX_PATH   where the dependencies are (default: /opt/dependencies)
#   BUILD_JOBS          the jobs of the build (default: 4, as on the runners of GitHub)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}"

case "${1:-}" in
build)
  mkdir -p build
  cd build
  export CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-/opt/dependencies}"
  # The verification cases run in each configuration of the build, with a limit raised for the time
  # they take with sanitizers. The snapshot references come from another environment, which this
  # one deviates from by rounding, as the CPU CI (scripts/compile/compile-variants.yaml).
  cmake .. -GNinja \
    -DCMAKE_BUILD_TYPE=Debug \
    -DCMAKE_CXX_FLAGS="${CXX_FLAGS:-}" \
    -DCMAKE_C_FLAGS_DEBUG="-O1 -g1" \
    -DCMAKE_CXX_FLAGS_DEBUG="-O1 -g1" \
    -DSANITIZERS="${SANITIZERS}" \
    -DCONFIGS="${CONFIGS}" \
    -DDR_QUAD_RULE="${DR_QUAD_RULE}" \
    -DGEMM_TOOLS_LIST=none \
    -DHOST_ARCH=hsw \
    -DNETCDF=OFF \
    -DGRAPH_PARTITIONING_LIBS=none \
    -DNUMA_AWARE_PINNING=OFF \
    -DTESTING=ON \
    -DTESTING_GENERATED="${TESTING_GENERATED}" \
    -DSEISSOL_SNAPSHOT_TOLERANCE_DOUBLE=1e-9 \
    -DSEISSOL_MINITEST_TIMEOUT=3600
  # a translation unit takes up to 4.5 GB with the sanitizers of GCC and all configurations
  ninja -j "${BUILD_JOBS:-4}"
  ;;
test)
  # (the thread sanitizer finds the runtime of OpenMP of the compiler by CXX)
  export CXX
  source .ci/sanitizers/environment.sh

  if [ "${THREADS:-false}" = true ]; then
    # A run of several ranks communicates in a thread of its own, and every run writes its output
    # from one; on the runner, with as many ranks and threads as CPUs, they share the CPUs of the
    # workers, one CPU for each rank of a node (SEISSOL_FREE_CPUS_MASK).
    cpus=$(ci_cpus)
    mask=""
    for rank in $(seq 0 $(((cpus > 4 ? cpus : 4) - 1))); do
      mask="${mask:+${mask},}$((rank % cpus))"
    done
    export SEISSOL_COMMTHREAD=1
    export SEISSOL_FREE_CPUS_MASK="${mask}"
    export ASYNC_MODE=THREAD
  fi

  cd build
  ctest -j "$(ci_cpus)" -LE convergence --output-on-failure --output-junit junit.xml
  ;;
*)
  echo "usage: $0 build|test" >&2
  exit 2
  ;;
esac
