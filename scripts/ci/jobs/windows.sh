#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The Windows job (.github/workflows/windows.yml), in the UCRT64 environment of MSYS2, in the
# folder that holds the checkouts of SeisSol (seissol/) and easi (easi/):
#
#   seissol/scripts/ci/jobs/windows.sh hdf5    builds HDF5 with MS-MPI into DEPS/hdf5
#   seissol/scripts/ci/jobs/windows.sh easi    builds easi into DEPS/easi
#   seissol/scripts/ci/jobs/windows.sh build   configures and builds SeisSol in build/
#   seissol/scripts/ci/jobs/windows.sh test    runs its tests there
#
#   DEPS           where the dependencies go (a path in Windows form, e.g. D:/a/_temp/dependencies)
#   PYTHON         the Python of the code generators and the tests, with numpy, scipy, h5py and
#                  PSpaMM (in Windows form)
#   MPI_HINTS      the arguments that point CMake to the MS-MPI of MSYS2
#   MSMPI_BIN      the folder of mpiexec (of the MS-MPI runtime)
#   HDF5_VERSION   the version of HDF5 to build

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"

case "${1:-}" in
hdf5)
  # MSYS2 only packages a serial HDF5; SeisSol needs a parallel one
  git clone --depth 1 --branch "${HDF5_VERSION}" https://github.com/HDFGroup/hdf5.git hdf5
  # (MPI_HINTS are several arguments)
  # shellcheck disable=SC2086
  cmake -S hdf5 -B hdf5-build -GNinja \
    -DCMAKE_BUILD_TYPE=Release \
    -DCMAKE_INSTALL_PREFIX="${DEPS}/hdf5" \
    ${MPI_HINTS} \
    -DHDF5_ENABLE_PARALLEL=ON \
    -DHDF5_ENABLE_ZLIB_SUPPORT=ON \
    -DBUILD_SHARED_LIBS=ON \
    -DBUILD_STATIC_LIBS=OFF \
    -DBUILD_TESTING=OFF \
    -DHDF5_BUILD_EXAMPLES=OFF \
    -DHDF5_BUILD_TOOLS=OFF \
    -DHDF5_BUILD_UTILS=OFF \
    -DHDF5_BUILD_FORTRAN=OFF \
    -DHDF5_BUILD_CPP_LIB=OFF \
    -DHDF5_BUILD_HL_LIB=ON
  ninja -C hdf5-build install
  ;;
easi)
  EASI_BUILD_DIR=easi-build "${SEISSOL_CI_ROOT}/scripts/ci/build-easi.sh" easi "${DEPS}/easi"
  ;;
build)
  # the tests are run once after linking them (to collect the test cases); hence, the HDF5 DLLs and
  # mpiexec are already needed here. pspamm-generator lies next to the Python
  PATH="$(cygpath -u "${DEPS}/hdf5/bin"):${PATH}:$(cygpath -u "${MSMPI_BIN}"):$(dirname "$(cygpath -u "${PYTHON}")")/Scripts"
  export PATH

  cmake --version

  # without NetCDF, as on macOS. The snapshot references of the verification tests come from
  # another environment, which this one deviates from by rounding; they run in the first
  # configuration only. No NUMA-aware pinning (it is Linux-only), and mpiexec for the MPI tests, by
  # default
  # shellcheck disable=SC2086
  cmake -S seissol -B build -GNinja \
    -DCMAKE_PREFIX_PATH="${DEPS}/hdf5;${DEPS}/easi" \
    ${MPI_HINTS} \
    -DGEMM_TOOLS_LIST=PSpaMM \
    -DNETCDF=OFF \
    -DCOVERAGE=OFF \
    -DTESTING=ON \
    -DTESTING_GENERATED=ON \
    -DHOST_ARCH=hsw \
    -DCMAKE_BUILD_TYPE=Release \
    -DCONFIGS=ci-cpu \
    -DGRAPH_PARTITIONING_LIBS=none \
    -DSEISSOL_SNAPSHOT_TOLERANCE_DOUBLE=1e-9 \
    -DSEISSOL_VERIFICATION_CONFIGS=first \
    -DPython3_EXECUTABLE="${PYTHON}"

  # one compile job per core: a translation unit of the set takes up to 1.6 GB
  ninja -C build -j "$(ci_cpus)"
  ;;
test)
  PATH="$(cygpath -u "${DEPS}/hdf5/bin"):${PATH}:$(cygpath -u "${MSMPI_BIN}")"
  export PATH
  cd build
  ctest -j --rerun-failed --output-on-failure
  ;;
*)
  echo "usage: $0 hdf5|easi|build|test" >&2
  exit 2
  ;;
esac
