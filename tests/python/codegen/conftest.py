# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Resolve SeisSol's codegen tree as a sys.path entry.

This mirrors how `codegen/generate.py` is actually invoked from CMake
(with cwd inside `codegen/`), which is why we add `codegen/` itself —
its submodules are imported as `kernels.*`, not `codegen.kernels.*`.
"""

import sys
from pathlib import Path

# tests/python/codegen/conftest.py -> repo root is three parents up
REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "codegen"))
