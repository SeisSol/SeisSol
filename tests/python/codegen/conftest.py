# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Resolve SeisSol's codegen tree and the yateto submodule as sys.path entries.

This mirrors how `codegen/generate.py` is actually invoked from CMake
(with cwd inside `codegen/`), which is why we add `codegen/` itself —
its submodules are imported as `kernels.*`, not `codegen.kernels.*`.
CMake hands yateto to it through PYTHONPATH; so does this file, for the
tests that run generate.py in a process of its own.
"""

import os
import sys
from pathlib import Path

# tests/python/codegen/conftest.py -> repo root is three parents up
REPO_ROOT = Path(__file__).resolve().parents[3]
YATETO_ROOT = REPO_ROOT / "submodules" / "yateto"
sys.path.insert(0, str(REPO_ROOT / "codegen"))
sys.path.insert(1, str(YATETO_ROOT))

# (ahead of the Python path that the environment gives, as in CMake)
os.environ["PYTHONPATH"] = os.pathsep.join(
    path for path in (str(YATETO_ROOT), os.environ.get("PYTHONPATH")) if path
)
