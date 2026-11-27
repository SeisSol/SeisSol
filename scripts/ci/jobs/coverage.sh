#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The coverage job (.github/workflows/coverage.yml):
#
#   scripts/ci/jobs/coverage.sh build      configures and builds SeisSol for coverage in build/
#   scripts/ci/jobs/coverage.sh coverage   runs its tests and gathers the coverage (seissol-coverage)
#
#   CC, CXX             the compilers (default: gcc-15, g++-15)
#   CMAKE_PREFIX_PATH   where the dependencies are (default: /opt/dependencies)
#   INSTALL_PREFIX      where SeisSol would be installed (default: /opt/seissol)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}"

case "${1:-}" in
build)
  mkdir -p build
  cd build
  export CC="${CC:-gcc-15}"
  export CXX="${CXX:-g++-15}"
  export CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-/opt/dependencies}"
  # without NetCDF for consistency with the clang-tidy check; the configurations of the CPU tests
  # (cmake/confignames.cmake), so that the code of every material counts. The snapshot references
  # of the verification tests come from another environment, which this one deviates from by
  # rounding, as the CPU CI (scripts/compile/compile-variants.yaml). The verification cases run in
  # the first configuration only, since in this unoptimized build, the ones of all configurations
  # would take hours; the CPU CI runs them all.
  cmake .. -GNinja -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX:-/opt/seissol}" -DGEMM_TOOLS_LIST=none \
    -DNETCDF=OFF -DCOVERAGE=ON -DTESTING=ON -DTESTING_GENERATED=ON -DHOST_ARCH=hsw \
    -DCMAKE_BUILD_TYPE=Debug -DCONFIGS=ci-cpu -DGRAPH_PARTITIONING_LIBS=ptscotch \
    -DSEISSOL_SNAPSHOT_TOLERANCE_DOUBLE=1e-9 -DSEISSOL_VERIFICATION_CONFIGS=first
  ninja
  ;;
coverage)
  cd build
  ninja seissol-coverage
  ;;
*)
  echo "usage: $0 build|coverage" >&2
  exit 2
  ;;
esac
