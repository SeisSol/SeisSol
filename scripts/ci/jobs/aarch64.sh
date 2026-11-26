#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The jobs on aarch64 (aarch64 in .github/workflows/architectures.yml), on an aarch64 machine:
#
#   scripts/ci/jobs/aarch64.sh build   configures and builds SeisSol with PSpaMM in build/
#   scripts/ci/jobs/aarch64.sh test    runs its tests there
#
#   CPU                 the CPU to compile for (-mcpu), e.g. neoverse-v2
#   HOST_ARCH           the host architecture of SeisSol, e.g. neon or sve128
#   CMAKE_PREFIX_PATH   where the dependencies are (default: /opt/dependencies)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}"

case "${1:-}" in
build)
  mkdir -p build
  cd build
  export CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-/opt/dependencies}"
  # The architecture the kernels are generated for (HOST_ARCH) sets no flags of the compiler on ARM;
  # -mcpu does, for the code GCC generates and for the instructions its assembler takes. The
  # snapshot references of the verification tests come from another environment, which this one
  # deviates from by rounding, as the CPU CI (scripts/compile/compile-variants.yaml).
  cmake .. -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS="-mcpu=${CPU}" \
    -DCMAKE_CXX_FLAGS="-mcpu=${CPU}" \
    -DHOST_ARCH="${HOST_ARCH}" \
    -DGEMM_TOOLS_LIST=PSpaMM \
    -DCONFIGS=ci-cpu \
    -DNETCDF=OFF \
    -DGRAPH_PARTITIONING_LIBS=none \
    -DTESTING=ON \
    -DTESTING_GENERATED=ON \
    -DSEISSOL_SNAPSHOT_TOLERANCE_DOUBLE=1e-9
  ninja -j "$(ci_cpus)"
  ;;
test)
  # all but the cases of dynamic rupture on fault2d-periodic, which take minutes each and hold the
  # friction laws against a spectral method: with them, the job takes longer than its 150 minutes on
  # the four cores of the runner. The CPU CI runs them.
  cd build
  ctest -j "$(ci_cpus)" -LE spectral --output-on-failure --output-junit junit.xml
  ;;
*)
  echo "usage: $0 build|test" >&2
  exit 2
  ;;
esac
