# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""Tests for the smoke-test driver in scripts/validate/smoke.py."""

from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
DRIVER = REPO_ROOT / "scripts" / "validate" / "smoke.py"

_spec = importlib.util.spec_from_file_location("seissol_smoke", DRIVER)
assert _spec is not None and _spec.loader is not None
smoke = importlib.util.module_from_spec(_spec)
sys.modules["seissol_smoke"] = smoke
_spec.loader.exec_module(smoke)


VALID_FIELDS = {
    "name": "all",
    "cycle-source": "tsc",
    "time": 0.0123,
    "cycles": 42.0,
    "gflop-libxsmm": 0.0,
    "gflop-pspamm": 0.0,
    "gflop-libxsmm-pspamm": 0.0,
    "gflop-nz": 1.5,
    "gflop-hw": 2.1,
    "gib": 0.4,
    "gib-kernel": 0.3,
    "gflopcycle-nz": 1.0,
    "gflopcycle-hw": 1.4,
    "gibcycle": 0.1,
    "gibcycle-kernel": 0.08,
    "gflops-nz": 122.0,
    "gflops-hw": 170.0,
    "gibs": 32.0,
    "gibs-kernel": 24.0,
}


def test_field_table_matches_the_reference_document():
    """The table and a known-good document have to agree, or every run fails."""
    assert set(smoke.PROXY_FIELDS) == set(VALID_FIELDS)


def test_valid_document_passes():
    smoke.check_proxy_fields(dict(VALID_FIELDS), "all")


@pytest.mark.parametrize("field", ["gflopcycle-nz", "gibcycle", "time"])
def test_non_finite_value_is_rejected(field):
    document = dict(VALID_FIELDS)
    document[field] = float("inf")
    with pytest.raises(smoke.CheckFailed, match="not finite"):
        smoke.check_proxy_fields(document, "all")


def test_zero_time_is_rejected():
    document = dict(VALID_FIELDS)
    document["time"] = 0.0
    with pytest.raises(smoke.CheckFailed, match="positive"):
        smoke.check_proxy_fields(document, "all")


def test_negative_counter_is_rejected():
    document = dict(VALID_FIELDS)
    document["gib"] = -1.0
    with pytest.raises(smoke.CheckFailed, match="non-negative"):
        smoke.check_proxy_fields(document, "all")


def test_missing_field_is_reported_by_name():
    document = dict(VALID_FIELDS)
    del document["cycles"]
    with pytest.raises(smoke.CheckFailed, match="missing field"):
        smoke.check_proxy_fields(document, "all")


def test_unexpected_field_is_reported_by_name():
    document = dict(VALID_FIELDS)
    document["gflop-newthing"] = 1.0
    with pytest.raises(smoke.CheckFailed, match="unexpected field"):
        smoke.check_proxy_fields(document, "all")


def test_kernel_name_has_to_match_what_was_requested():
    with pytest.raises(smoke.CheckFailed, match="name"):
        smoke.check_proxy_fields(dict(VALID_FIELDS), "ader")


def test_hardware_flops_below_nonzero_flops_is_rejected():
    document = dict(VALID_FIELDS)
    document["gflop-nz"] = 9.9
    with pytest.raises(smoke.CheckFailed, match="gflop-hw"):
        smoke.check_proxy_fields(document, "all")


def test_bare_inf_is_diagnosed_as_a_writer_problem():
    text = '{"name":"all","gflopcycle-nz":inf}'
    with pytest.raises(smoke.CheckFailed, match="non-finite number reached the writer"):
        smoke.parse_proxy_json(text)


def test_json_infinity_literal_is_also_rejected():
    text = '{"name":"all","gflopcycle-nz":Infinity}'
    with pytest.raises(smoke.CheckFailed, match="not valid JSON"):
        smoke.parse_proxy_json(text)


def test_preamble_on_stdout_is_reported():
    stdout = 'Cycles via __rdtsc()\n{"name":"all"}\n'
    with pytest.raises(smoke.CheckFailed, match="non-JSON lines"):
        smoke.extract_json_object(stdout)


def test_plain_output_is_not_mistaken_for_json():
    with pytest.raises(smoke.CheckFailed, match="not a JSON object"):
        smoke.extract_json_object("=== PERFORMANCE SUMMARY ===\ntime : 1\n")


def test_signal_death_counts_as_failure():
    smoke.check_termination(-6, expect_failure=True)
    with pytest.raises(smoke.CheckFailed, match="SIGABRT"):
        smoke.check_termination(-6, expect_failure=False)


def test_clean_nonzero_exit_counts_as_failure():
    smoke.check_termination(255, expect_failure=True)
    with pytest.raises(smoke.CheckFailed, match="exit code 255"):
        smoke.check_termination(255, expect_failure=False)


def test_unexpected_success_is_rejected():
    with pytest.raises(smoke.CheckFailed, match="expected an unsuccessful exit"):
        smoke.check_termination(0, expect_failure=True)


def test_unknown_cycle_source_is_rejected():
    document = dict(VALID_FIELDS)
    document["cycle-source"] = "rdpmc"
    with pytest.raises(smoke.CheckFailed, match="cycle-source"):
        smoke.check_proxy_fields(document, "all")


def test_absent_counter_may_not_report_ticks():
    document = dict(VALID_FIELDS)
    document["cycle-source"] = "none"
    document["cycles"] = 17.0
    with pytest.raises(smoke.CheckFailed, match="cycles is non-zero"):
        smoke.check_proxy_fields(document, "all")


def test_absent_counter_with_zero_ticks_is_fine():
    document = dict(VALID_FIELDS)
    document["cycle-source"] = "none"
    document["cycles"] = 0.0
    document["gflopcycle-nz"] = 0.0
    document["gflopcycle-hw"] = 0.0
    document["gibcycle"] = 0.0
    document["gibcycle-kernel"] = 0.0
    smoke.check_proxy_fields(document, "all")


def test_combined_gflop_has_to_be_the_sum():
    document = dict(VALID_FIELDS)
    document["gflop-libxsmm"] = 1.0
    document["gflop-pspamm"] = 2.0
    document["gflop-libxsmm-pspamm"] = 5.0
    with pytest.raises(smoke.CheckFailed, match="is not the sum of"):
        smoke.check_proxy_fields(document, "all")
