#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The emulated jobs (emulated in .github/workflows/architectures.yml): SeisSol cross-compiled for
# another architecture, and its tests run with qemu-user (see .ci/emulated):
#
#   scripts/ci/jobs/emulated.sh dependencies   builds the libraries of the target into CROSS_PREFIX
#   scripts/ci/jobs/emulated.sh easi           builds easi (checked out in easi/) into it
#   scripts/ci/jobs/emulated.sh build          configures and builds SeisSol in build/
#   scripts/ci/jobs/emulated.sh test           runs its tests there
#
#   TRIPLET        the target, e.g. riscv64-linux-gnu
#   CROSS_PREFIX   the prefix of its libraries (default: /opt/cross/<TRIPLET>)
#   QEMU_CPU       the CPU that qemu emulates
#   CC, CXX        the cross compilers
#   FLAGS          the flags for C and C++, e.g. -march=...
#   HOST_ARCH      the host architecture of SeisSol
#   GEMM_TOOLS     the GEMM generators (GEMM_TOOLS_LIST)
#   CONFIGS        the configurations

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}"

CROSS_PREFIX="${CROSS_PREFIX:-/opt/cross/${TRIPLET}}"
# the cross compilation and qemu for the target
source .ci/emulated/environment.sh "${TRIPLET}" "${CROSS_PREFIX}" "${QEMU_CPU}"
toolchain="$(pwd)/.ci/emulated/toolchain.cmake"

case "${1:-}" in
dependencies)
  .ci/emulated/dependencies.sh "${TRIPLET}" "${CROSS_PREFIX}"
  ;;
easi)
  # next to the other libraries of the target, after the cache has kept those
  "${SEISSOL_CI_ROOT}/scripts/ci/build-easi.sh" easi "${CROSS_PREFIX}" \
    -DCMAKE_TOOLCHAIN_FILE="${toolchain}"
  ;;
build)
  # HOST_ARCH sets no flags of the compiler for most of these targets: the flags name the vector
  # instructions, for the code of the compiler and for its assembler, which takes those of the
  # kernels. Eigen is the one of the build machine, whose headers are the same for every target.
  # The verification cases take longer than natively, and the snapshot references come from another
  # environment, which this one deviates from by rounding, as the CPU CI
  # (scripts/compile/compile-variants.yaml).
  cmake -S . -B build -GNinja \
    -DCMAKE_TOOLCHAIN_FILE="${toolchain}" \
    -DEigen3_DIR=/usr/share/eigen3/cmake \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_C_FLAGS="${FLAGS}" \
    -DCMAKE_CXX_FLAGS="${FLAGS}" \
    -DHOST_ARCH="${HOST_ARCH}" \
    -DGEMM_TOOLS_LIST="${GEMM_TOOLS}" \
    -DCONFIGS="${CONFIGS}" \
    -DNETCDF=OFF \
    -DGRAPH_PARTITIONING_LIBS=none \
    -DNUMA_AWARE_PINNING=OFF \
    -DTESTING=ON \
    -DTESTING_GENERATED=ON \
    -DSEISSOL_VERIFICATION_CONFIGS=first \
    -DSEISSOL_MINITEST_TIMEOUT=1800 \
    -DSEISSOL_SNAPSHOT_TOLERANCE_DOUBLE=1e-9
  # one compile job per core, as for the native builds of the set
  cmake --build build --parallel "$(ci_cpus)"
  ;;
test)
  cd build
  ctest -j "$(ci_cpus)" -LE convergence --output-on-failure --output-junit junit.xml
  ;;
*)
  echo "usage: $0 dependencies|easi|build|test" >&2
  exit 2
  ;;
esac
