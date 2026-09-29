# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "postprocessing" / "validation"))


def pytest_addoption(parser):
    parser.addoption(
        "--vtkhdf", action="store", default=None, help="a VTKHDF file to check"
    )


@pytest.fixture
def vtkhdf_path(request):
    path = request.config.getoption("--vtkhdf")
    if path is None:
        pytest.skip("no --vtkhdf given")
    return path
