# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause

"""A case-level xfail of cases.json expects the run of the case to fail as well, as where
the libraries of an environment are too old for the case; without one, a failed run
fails."""

import sys
from pathlib import Path

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(SEISSOL_ROOT / "scripts" / "validate"))

import verify  # noqa: E402

TOO_OLD = "easi is too old"


def failed(xfail, labels=()):
    case = verify.CaseSpec(
        key="loh1-ic/double",
        config="elastic-o6-f64",
        params_path=Path("loh1-ic/parameters.par"),
        case_dir=Path("loh1-ic"),
        prefix="tpv",
        xfail=xfail,
    )
    compares = verify.failed_run(case, 133, set(labels))
    return verify.CaseResult(case, True, 133, 1.0, compares)


def test_failed_run_fails():
    result = failed([])
    assert not result.passed
    assert verify._state_of(result) == "RC=133"


def test_xfail_expects_the_run_to_fail():
    result = failed([{"when": {}, "reason": TOO_OLD}])
    assert result.passed
    assert verify._state_of(result) == "PASS"
    assert [(c.category, c.xfail_outcome, c.xfail_reason) for c in result.compares] == [
        ("run", "xfail", TOO_OLD)
    ]


def test_xfail_holds_only_where_its_rule_matches():
    rule = [{"when": {"label": "gpu"}, "reason": TOO_OLD}]
    assert not failed(rule).passed
    assert failed(rule, labels=["gpu"]).passed
