# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""What is read back from runs, and the scaling of a series."""

import pytest
from strongscaling import collect, report


def test_cpu_sets():
    assert collect.cpu_set("0,1,4-6") == {0, 1, 4, 5, 6}
    assert collect.cpu_set("none") is None
    assert collect.cpu_set("no-info") is None


def pinning_row(host, rank, workers, comm="no-info", numa="0"):
    return {
        "hostname": host,
        "rank": str(rank),
        "workermask": workers,
        "workernuma": numa,
        "commthread_mask": comm,
    }


def test_pinning_problems():
    fine = [pinning_row("a", 0, "0,1", "2"), pinning_row("a", 1, "3,4", "5")]
    assert collect.pinning_problems(fine, 2) == []
    shared = [pinning_row("a", 0, "0,1"), pinning_row("a", 1, "1,2")]
    assert collect.pinning_problems(shared, 2) == ["a: ranks 0 and 1 share 1 CPU(s)"]
    # on different nodes, the same CPUs are fine
    assert (
        collect.pinning_problems(
            [pinning_row("a", 0, "0,1"), pinning_row("b", 1, "0,1")], 2
        )
        == []
    )
    comm = [pinning_row("a", 0, "0,1", "1")]
    assert "communication thread" in collect.pinning_problems(comm, 2)[0]
    few = [pinning_row("a", 0, "0")]
    assert "1 CPUs for 2 threads" in collect.pinning_problems(few, 2)[0]
    spread = [pinning_row("a", 0, "0,1", numa="0,1")]
    assert "NUMA" in collect.pinning_problems(spread, 2)[0]


def test_steady_state():
    epochs = [
        {"epoch": "1", "timeFrom": "0", "timeTo": "0.5", "wallMax": "3"},
        {"epoch": "2", "timeFrom": "0.5", "timeTo": "1", "wallMax": "1"},
        {"epoch": "3", "timeFrom": "1", "timeTo": "1.5", "wallMax": "1"},
    ]
    assert collect.steady_state(epochs, 1) == (2, 2.0, 1.0)
    assert collect.steady_state(epochs, 3) == (0, None, None)


def row(nodes, time, status="completed", repetition=1, updates=50000.0):
    return {
        "run": f"n{nodes}-r{repetition}",
        "build": "b",
        "variant": "v",
        "repetition": float(repetition),
        "nodes": float(nodes),
        "ranks": float(nodes),
        "status": status,
        "wallTime": time,
        "computeTime": time,
        "cellUpdates": updates,
        "endTime": 1.0,
        "steadyWallTime": time / 2,
        "steadySimulatedTime": 0.5,
    }


def test_series_scaling():
    rows = [row(1, 100.0), row(1, 110.0, repetition=2), row(2, 52.5), row(4, 40.0)]
    (series,) = report.make_series(rows, "total")
    points = {point.nodes: point for point in series.points}
    assert points[1].time == pytest.approx(105.0)
    assert points[2].speedup == pytest.approx(2.0)
    assert points[2].efficiency == pytest.approx(1.0)
    assert points[4].efficiency == pytest.approx(105.0 / 40.0 / 4)
    table = report.scaling_rows([series])
    assert table[1]["cellUpdatesPerSecondPerNode"] == pytest.approx(50000 / 52.5 / 2)


def test_steady_metric():
    assert report.metric_value(row(1, 100.0), "steady") == pytest.approx(100.0)
    assert report.metric_value(row(1, 100.0), "compute") == 100.0


def test_result_checks():
    rows = [row(1, 100.0), row(1, 130.0, repetition=2), row(2, 50.0, updates=60000.0)]
    rows.append(dict(row(4, 0.0, status="killed"), wallTime=None))
    series = report.make_series(rows, "total")
    findings = report.checks(rows, series)
    assert any("n4-r1: killed" in finding for finding in findings)
    assert any("cell updates differ" in finding for finding in findings)
    assert any("spread" in finding for finding in findings)
