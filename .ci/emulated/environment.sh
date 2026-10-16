# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The environment that builds SeisSol for another architecture and runs it on QEMU. In bash, with
# the triplet of the target, the prefix of its libraries (.ci/emulated/dependencies.sh) and the CPU
# that QEMU emulates, in the syntax of its option -cpu:
#
#   source .ci/emulated/environment.sh riscv64-linux-gnu /opt/cross/riscv64-linux-gnu \
#     rv64,v=true,vlen=128,vext_spec=v1.0
#
# A program of the target starts on QEMU through binfmt_misc, as the qemu-user-static of Ubuntu
# registers it; QEMU reads its options from the environment, which every process it starts
# inherits, also the ranks that mpirun starts.

export SEISSOL_CROSS_TRIPLET="$1"
export SEISSOL_CROSS_PREFIX="$2"

# the dynamic loader and the C and C++ runtime libraries of the cross compiler
export QEMU_LD_PREFIX="/usr/${SEISSOL_CROSS_TRIPLET}"
export QEMU_CPU="$3"

# mpirun and the compiler wrappers of the MPI library of the target
export PATH="${SEISSOL_CROSS_PREFIX}/bin:${PATH}"
export LD_LIBRARY_PATH="${SEISSOL_CROSS_PREFIX}/lib:${LD_LIBRARY_PATH:-}"
# PMIx of the target has no compression library to compress large messages with, which it would
# point out in every run of mpirun
export PMIX_MCA_pcompress_base_silence_warning=1
