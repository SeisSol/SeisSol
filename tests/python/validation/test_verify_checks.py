# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""A failed check of a compare script (such as meshcompare's global-id verdict) fails the
comparison in verify.py, also where per-quantity thresholds decide the verdict."""

import sys
from pathlib import Path

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SEISSOL_ROOT / "scripts" / "validate"))

import validation_report  # noqa: E402
import verify  # noqa: E402

TPV_DATA = {"surface": {"enabled": True, "quantities": {"v1": 1e-6}}}


def summary(tmp_path, checks):
    path = tmp_path / "report.json"
    validation_report.write_report_json(
        str(path), "surface", 1e-6, True, {"v1": 1e-9}, checks=checks
    )
    return verify._read_summary(path)


def evaluate(tmp_path, checks):
    result = verify._apply_summary(
        verify.CompareResult("surface", False, 1e-6), summary(tmp_path, checks)
    )
    verify._evaluate_quantities(result, TPV_DATA, "double", None)
    return result


def test_passed_check_passes(tmp_path):
    assert evaluate(tmp_path, {"global-id": True}).passed


def test_failed_check_fails_despite_thresholds(tmp_path):
    result = evaluate(tmp_path, {"global-id": False})
    assert not result.passed
    assert [name for name, _, _ in result.failures] == ["check:global-id"]


def test_reports_without_checks_still_read(tmp_path):
    assert "checks" not in summary(tmp_path, None)
    assert evaluate(tmp_path, None).passed
