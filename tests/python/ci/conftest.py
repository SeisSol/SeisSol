# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""Load the `.ci/` scripts (they're not a package; hyphen-free names)."""

import importlib.util
import sys
from pathlib import Path

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
_CI = SEISSOL_ROOT / ".ci"


def _load(name):
    spec = importlib.util.spec_from_file_location(name, _CI / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod  # allow re-import in other test files
    spec.loader.exec_module(mod)
    return mod
