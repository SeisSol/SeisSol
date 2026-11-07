#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Builds the libraries that SeisSol needs (Open MPI, zlib, HDF5, BLAS and LAPACK, yaml-cpp) for the
# target of a GNU cross compiler, into a prefix:
#
#   .ci/emulated/dependencies.sh <triplet> <prefix>
#
# with CC, CXX and FC the compilers, if they are not <triplet>-gcc, <triplet>-g++ and
# <triplet>-gfortran. The sources are those of the archive of Ubuntu, the same for every target;
# Ubuntu keeps the packages of most other architectures in an archive of its own, and has none for
# LoongArch.
#
# Their build systems run the test programs they compile: the target has to run on QEMU,
# registered with binfmt_misc, and with the environment of .ci/emulated/environment.sh.

set -euo pipefail

if [ $# -ne 2 ]; then
  echo "usage: $0 <triplet> <prefix>" >&2
  exit 2
fi

triplet=$1
prefix=$2
export CC=${CC:-${triplet}-gcc}
export CXX=${CXX:-${triplet}-g++}
export FC=${FC:-${triplet}-gfortran}

here="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
pool=http://archive.ubuntu.com/ubuntu/pool
work=$(mktemp -d)
trap 'rm -rf "${work}"' EXIT

# <url> <sha256>: download a source archive and unpack it into the working directory; the checksum
# is what makes the download over HTTP safe
fetch() {
  local archive="${work}/${1##*/}"
  curl --fail --silent --show-error --location --retry 3 --output "${archive}" "$1"
  echo "$2  ${archive}" | sha256sum --check --quiet
  tar -xf "${archive}" -C "${work}"
}

mkdir -p "${prefix}"

fetch "${pool}/universe/o/openmpi/openmpi_5.0.8.orig.tar.xz" \
  c10cfb9e4939c227e1f8faa0a1b5829da58f04379c2fb8a55af4aa0b217978d2
(
  cd "${work}/openmpi-5.0.8"
  # The libraries of the tarball, not those of the build machine that the cross compiler finds in
  # /usr/include as well. MPI-IO is that of OMPIO: ROMIO, the other implementation, compiles its
  # part of the shared libmpi without -fPIC, which the linker rejects on AArch64.
  ./configure --host="${triplet}" --prefix="${prefix}" \
    --disable-mpi-fortran --disable-oshmem --disable-io-romio \
    --with-hwloc=internal --with-libevent=internal --with-pmix=internal --with-prrte=internal
  make -j "$(nproc)"
  make install
)

# for the compression of HDF5, which the writers of SeisSol use
fetch "${pool}/main/z/zlib/zlib_1.3.dfsg.orig.tar.xz" \
  5eea0322c1c21c75cad3b607ac1c43ff5c71e014b8ac4a34300b5e2b80d02e70
(
  cd "${work}/zlib-1.3.dfsg"
  # CHOST names the archiver of the target
  CHOST="${triplet}" ./configure --prefix="${prefix}"
  make -j "$(nproc)"
  make install
)

fetch "${pool}/universe/h/hdf5/hdf5_1.14.6+repack.orig.tar.gz" \
  c18bc47d07f76a5c2cfb70773bfe4152857e48ca4c96e7aef2684c96fa0e6a7c
(
  cd "${work}/hdf5-1.14.6+repack"
  # with the zlib of the target, and without szip, whose header of the build machine the cross
  # compiler would find
  ./configure --host="${triplet}" --prefix="${prefix}" CC="${prefix}/bin/mpicc" \
    --enable-parallel --disable-fortran --disable-cxx --disable-tests --disable-tools \
    --disable-hl --with-zlib="${prefix}" --with-szlib=no
  make -j "$(nproc)"
  make install
)

# the reference BLAS and LAPACK, which SeisSol takes for poroelastic materials
fetch "${pool}/main/l/lapack/lapack_3.12.0.orig.tar.gz" \
  eac9570f8e0ad6f30ce4b963f4f033f0f643e7c3912fc9ee6cd99120675ad48b
SEISSOL_CROSS_TRIPLET="${triplet}" SEISSOL_CROSS_PREFIX="${prefix}" \
  cmake -S "${work}/lapack-3.12.0" -B "${work}/lapack-build" -GNinja \
  -DCMAKE_TOOLCHAIN_FILE="${here}/toolchain.cmake" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${prefix}" -DBUILD_SHARED_LIBS=ON -DBUILD_TESTING=OFF \
  -DCBLAS=OFF -DLAPACKE=OFF
cmake --build "${work}/lapack-build"
cmake --install "${work}/lapack-build"

fetch "${pool}/main/y/yaml-cpp/yaml-cpp_0.8.0+dfsg.orig.tar.xz" \
  32c85a5f36633a9005b20261b35fccb0d8ec69e9ef16c67de7d5d007944e3709
SEISSOL_CROSS_TRIPLET="${triplet}" SEISSOL_CROSS_PREFIX="${prefix}" \
  cmake -S "${work}"/jbeder-yaml-cpp-* -B "${work}/yaml-cpp-build" -GNinja \
  -DCMAKE_TOOLCHAIN_FILE="${here}/toolchain.cmake" -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_INSTALL_PREFIX="${prefix}" -DYAML_BUILD_SHARED_LIBS=ON \
  -DYAML_CPP_BUILD_TESTS=OFF -DYAML_CPP_BUILD_TOOLS=OFF
cmake --build "${work}/yaml-cpp-build"
cmake --install "${work}/yaml-cpp-build"
