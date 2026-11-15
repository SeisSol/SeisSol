#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The macOS job (.github/workflows/macos.yml), with the libraries of Homebrew:
#
#   scripts/ci/jobs/macos.sh build   configures and builds SeisSol with PSpaMM in build/
#   scripts/ci/jobs/macos.sh test    runs its tests there
#
#   CMAKE_PREFIX_PATH   where easi is (default: /opt/dependencies)
#   INSTALL_PREFIX      where SeisSol would be installed (default: /opt/seissol)
#
# (for the bash 3.2 of macOS as well)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}"

case "${1:-}" in
build)
  mkdir -p build
  cd build
  export CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:-/opt/dependencies}"

  # cf. https://gist.github.com/scivision/16c2ca1dc250f54d34f1a1a35596f4a0
  OpenMP_ROOT="$(brew --prefix)/opt/libomp"
  export OpenMP_ROOT
  export LDFLAGS="${LDFLAGS:-} -L${OpenMP_ROOT}/lib"
  export CFLAGS="-mcpu=native ${CFLAGS:-}"
  export CXXFLAGS="-mcpu=native ${CXXFLAGS:-}"

  cmake --version

  # without NetCDF for consistency with the clang-tidy check. The snapshot references of the
  # verification tests come from another environment, which this one deviates from by rounding, as
  # the CPU CI (scripts/compile/compile-variants.yaml). They run in the first configuration only, as
  # before the configurations were built into one executable; the CPU CI runs them in all.
  cmake .. -GNinja -DCMAKE_INSTALL_PREFIX="${INSTALL_PREFIX:-/opt/seissol}" \
    -DGEMM_TOOLS_LIST=PSpaMM \
    -DNETCDF=OFF \
    -DCOVERAGE=OFF \
    -DTESTING=ON \
    -DTESTING_GENERATED=ON \
    -DHOST_ARCH=neon \
    -DCMAKE_BUILD_TYPE=Release \
    -DCONFIGS=ci-cpu \
    -DGRAPH_PARTITIONING_LIBS=none \
    -DNUMA_AWARE_PINNING=OFF \
    -DSEISSOL_SNAPSHOT_TOLERANCE_DOUBLE=1e-9 \
    -DSEISSOL_VERIFICATION_CONFIGS=first \
    -DPython3_EXECUTABLE="$(command -v python3)"

  # one compile job per core: a translation unit of the set takes up to 1.6 GB
  ninja -j "$(ci_cpus)"
  ;;
test)
  cd build
  ctest -j --rerun-failed --output-on-failure
  ;;
*)
  echo "usage: $0 build|test" >&2
  exit 2
  ;;
esac
