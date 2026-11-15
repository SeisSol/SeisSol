#!/usr/bin/env bash
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

# The tests of the Python code (.github/workflows/pytest-python.yml), with the coverage of the code
# generator and the scripts; the results also go to tests/python/pytest-results.xml (JUnit), the
# coverage to tests/python/coverage.xml. Further arguments go to pytest, e.g. -x or -k <name>.
#
#   scripts/ci/jobs/pytest.sh [<pytest arguments>...]
#
#   PYTHON   the Python to run pytest with (default: python3); it needs pytest, pytest-cov and
#            hypothesis, and the packages the workflow installs

set -euo pipefail
source "$(dirname "${BASH_SOURCE[0]}")/../common.sh"
cd "${SEISSOL_CI_ROOT}/tests/python"

"${PYTHON:-python3}" -m pytest -v --tb=short \
  --cov=../../codegen \
  --cov=../../scripts/validate \
  --cov=../../scripts/perf \
  --cov=../../scripts/benchmark \
  --cov=../../.ci \
  --cov-report=xml \
  --cov-report=term-missing \
  --cov-branch \
  --junitxml=pytest-results.xml \
  "$@"
