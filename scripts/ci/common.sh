# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# Sourced by the scripts of the CI jobs (scripts/ci/jobs): the environment MPI runs in on the
# runners and in the containers of the CI, with more ranks than cores and as root; the root of the
# repository in SEISSOL_CI_ROOT; and ci_cpus, the number of CPUs.

export OMPI_MCA_rmaps_base_oversubscribe=1
export PRTE_MCA_rmaps_base_oversubscribe=1
export PRTE_MCA_rmaps_default_mapping_policy=:oversubscribe
export OMPI_ALLOW_RUN_AS_ROOT=1
export OMPI_ALLOW_RUN_AS_ROOT_CONFIRM=1

SEISSOL_CI_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)"
export SEISSOL_CI_ROOT

ci_cpus() {
  if command -v nproc > /dev/null; then
    nproc
  else
    sysctl -n hw.ncpu
  fi
}
