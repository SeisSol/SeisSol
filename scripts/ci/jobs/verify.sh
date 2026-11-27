#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The verification of an installed SeisSol against the precomputed solutions (seissol-run-test in
# .github/workflows/build-seissol-cpu.yml and build-seissol-gpu.yml), with
# scripts/validate/verify.py, in the current directory:
#
#   scripts/ci/jobs/verify.sh
#
#   INSTALL_DIR      the installation of SeisSol (default: install)
#   PRECOMPUTED_DIR  the checkout of seissol/precomputed-seissol (default: precomputed)
#   MPI              the MPI of the build: openmpi or another one (default: openmpi)
#   RANKS, THREADS   the ranks and the threads per rank
#   COMMTHREAD       whether to run a communication thread (true, false)
#   REPORT           the report of verify.py (default: report.txt)

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"

install_dir="$(realpath "${INSTALL_DIR:-install}")"

mpi_options=""
if [ "${MPI:-openmpi}" = openmpi ]; then
  mpi_options="--oversubscribe --allow-run-as-root -bind-to none"
fi

export PATH="${install_dir}/bin:${PATH}"
export LD_LIBRARY_PATH="${install_dir}/lib:${install_dir}/lib64:${LD_LIBRARY_PATH:-}"
chmod +x "${install_dir}"/bin/*
# the executables verify.py chooses from
ls -l "${install_dir}/bin"

export OMP_NUM_THREADS="${THREADS}"
export SEISSOL_COMMTHREAD="${COMMTHREAD:-false}"
export SEISSOL_MINI_NUM_ELEMENTS=100

python3 "${SEISSOL_CI_ROOT}/scripts/validate/verify.py" \
  --precomputed-dir "${PRECOMPUTED_DIR:-precomputed}" \
  --binary-dir "${install_dir}" \
  --all \
  --mpi-cmd mpirun \
  --mpi-args "${mpi_options}" \
  --ranks "${RANKS}" \
  --threads "${THREADS}" \
  --color always \
  --report "${REPORT:-report.txt}"
