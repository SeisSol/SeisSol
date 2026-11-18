# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The package seissol of the repository, without installing it."""

import sys
from pathlib import Path

PYTHON = Path(__file__).resolve().parents[3] / "python"
sys.path.insert(0, str(PYTHON))
