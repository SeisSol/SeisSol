#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The builds and tests of seissol-build-test (.github/workflows/build-seissol-cpu.yml and
# build-seissol-gpu.yml): scripts/compile/compile_variants.py with the given arguments, in the
# environment that MPI needs on the runners and in the containers of the CI. E.g. the job of GCC:
#
#   scripts/ci/jobs/build-variants.sh --preset ci-cpu --toolchain gcc-13 \
#     --set defines.WARNINGS_AS_ERRORS=ON --install-prefix /opt/seissol --max-parallel 1 \
#     --report build-report.json
#
# The build of each variant goes to build/<variant>/ (with its log in build.log, and the results of
# its tests in ctest-results.xml).

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}"

python3 scripts/compile/compile_variants.py "$@"
