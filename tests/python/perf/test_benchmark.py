# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""Tests for the proxy benchmark in scripts/perf/benchmark.py."""

from __future__ import annotations

import importlib.util
import json
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
SCRIPT = REPO_ROOT / "scripts" / "perf" / "benchmark.py"

_spec = importlib.util.spec_from_file_location("seissol_benchmark", SCRIPT)
assert _spec is not None and _spec.loader is not None
benchmark = importlib.util.module_from_spec(_spec)
sys.modules["seissol_benchmark"] = benchmark
_spec.loader.exec_module(benchmark)


HELP = textwrap.dedent(
    """\
    The SeisSol proxy is used to benchmark the kernels used in the SeisSol earthquake
    simulation software. This version of SeisSol proxy was built for the following
    configurations:
    elastic-linearck-o6-f64-stroud
      material: elastic
      order: 6
    viscoelastic-linearckanelastic-m3-o6-f32-stroud
      material: viscoelastic
    elastic-linearck-o3-f64-stroud-s8
      simulations: 8
    """
)


def fake_proxy(tmp_path: Path, seconds: float = 0.5) -> Path:
    """A stand-in for the proxy: its help, and a JSON document in the format of the proxy."""
    script = tmp_path / "proxy.py"
    script.write_text(
        textwrap.dedent(
            f"""\
            import json, sys
            if "--help" in sys.argv:
                print({HELP!r})
                sys.exit(0)
            cells, steps, kernel = sys.argv[1:4]
            print("Running on host")
            print(json.dumps({{"name": kernel, "time": {seconds}, "gflops-hw": 12.0,
                               "gflops-nz": 9.0, "gibs": 3.0}}))
            """
        ),
        encoding="utf-8",
    )
    return script


def report(times: dict[str, list[float]], **environment) -> dict:
    base = {
        "cpu": "Test CPU",
        "cpus": 8,
        "omp_num_threads": "8",
        "omp_places": "cores",
        "omp_proc_bind": "close",
    }
    base.update(environment)
    return {
        "schema": benchmark.REPORT_SCHEMA,
        "environment": base,
        "benchmarks": {
            name: {"ns_per_cell_step": values, "gflops_hw": [1.0] * len(values)}
            for name, values in times.items()
        },
    }


def test_spread_ignores_a_single_outlier():
    steady = [100.0, 100.5, 99.5, 100.2, 99.8]
    assert benchmark.robust_spread(steady) < 0.01
    assert benchmark.robust_spread(steady + [180.0]) < 0.01


def test_spread_of_a_single_value_is_zero():
    assert benchmark.robust_spread([42.0]) == 0.0


def test_configs_are_read_from_the_help(tmp_path):
    proxy = [sys.executable, str(fake_proxy(tmp_path))]
    assert benchmark.built_configs(proxy) == [
        "elastic-linearck-o6-f64-stroud",
        "viscoelastic-linearckanelastic-m3-o6-f32-stroud",
        "elastic-linearck-o3-f64-stroud-s8",
    ]


def test_run_writes_every_sample(tmp_path):
    proxy = fake_proxy(tmp_path, seconds=0.25)
    launcher = tmp_path / "launch.sh"
    launcher.write_text(f'#!/bin/sh\nexec {sys.executable} "$@"\n', encoding="utf-8")
    launcher.chmod(0o755)
    output = tmp_path / "report.json"
    entries = tmp_path / "entries.json"
    code = benchmark.main(
        [
            "run",
            "--proxy",
            str(proxy),
            "--launcher",
            str(launcher),
            "--config",
            "elastic-linearck-o6-f64-stroud",
            "--kernel",
            "all",
            "--cells",
            "1000",
            "--timesteps",
            "5",
            "--repetitions",
            "3",
            "--warmup",
            "0",
            "--output",
            str(output),
            "--github-benchmark",
            str(entries),
        ]
    )
    assert code == 0
    written = json.loads(output.read_text())
    data = written["benchmarks"]["elastic-linearck-o6-f64-stroud/all"]
    # 0.25 s for 1000 cells and 5 steps
    assert data["ns_per_cell_step"] == pytest.approx([50000.0] * 3)
    assert written["parameters"]["repetitions"] == 3
    assert json.loads(entries.read_text())[0]["unit"] == "ns/(cell step)"


def test_a_proxy_in_the_working_directory_is_found(tmp_path, monkeypatch):
    """./proxyseissol is no name to look up in the PATH."""
    script = fake_proxy(tmp_path)
    proxy = tmp_path / "proxyseissol"
    proxy.write_text(
        f"#!{sys.executable}\n" + script.read_text(encoding="utf-8"), encoding="utf-8"
    )
    proxy.chmod(0o755)
    monkeypatch.chdir(tmp_path)
    output = tmp_path / "report.json"
    arguments = ["run", "--proxy", "./proxyseissol", "--kernel", "all"]
    arguments += ["--repetitions", "1", "--warmup", "0", "--output", str(output)]
    assert benchmark.main(arguments) == 0
    assert len(json.loads(output.read_text())["benchmarks"]) == 3


def test_a_failing_proxy_is_reported(tmp_path):
    proxy = tmp_path / "proxy.py"
    proxy.write_text("import sys\nsys.exit(3)\n", encoding="utf-8")
    with pytest.raises(benchmark.BenchmarkError, match="exit code 3"):
        benchmark.run_proxy([sys.executable, str(proxy)], "c", "all", 1, 1)


def test_change_within_the_noise_is_no_regression():
    noisy = [100.0, 104.0, 96.0, 102.0, 98.0]
    rows = benchmark.compare_reports(
        report({"a/all": noisy}),
        report({"a/all": [value * 1.06 for value in noisy]}),
        threshold=0.05,
        noise_factor=3.0,
    )
    assert rows[0]["verdict"] == "unchanged"


def test_change_beyond_threshold_and_noise_is_a_regression():
    steady = [100.0, 100.2, 99.8, 100.1, 99.9]
    rows = benchmark.compare_reports(
        report({"a/all": steady}),
        report({"a/all": [value * 1.08 for value in steady]}),
        threshold=0.05,
        noise_factor=3.0,
    )
    assert rows[0]["verdict"] == "slower"
    assert rows[0]["change"] == pytest.approx(0.08)


def test_speedup_and_unmatched_benchmarks():
    steady = [100.0, 100.2, 99.8]
    rows = benchmark.compare_reports(
        report({"a/all": steady, "b/all": steady}),
        report({"a/all": [value * 0.8 for value in steady], "c/all": steady}),
        threshold=0.05,
        noise_factor=3.0,
    )
    verdicts = {row["name"]: row["verdict"] for row in rows}
    assert verdicts == {"a/all": "faster", "b/all": "gone", "c/all": "new"}


def test_reports_of_different_machines_are_not_compared(tmp_path):
    baseline = tmp_path / "baseline.json"
    current = tmp_path / "current.json"
    baseline.write_text(json.dumps(report({"a/all": [1.0]}, cpu="CPU A")))
    current.write_text(json.dumps(report({"a/all": [2.0]}, cpu="CPU B")))
    summary = tmp_path / "summary.md"
    code = benchmark.main(
        [
            "compare",
            "--baseline",
            str(baseline),
            "--current",
            str(current),
            "--summary",
            str(summary),
            "--fail-on-regression",
        ]
    )
    assert code == 0
    assert "CPU A -> CPU B" in summary.read_text()


def test_regression_fails_when_asked_to(tmp_path):
    baseline = tmp_path / "baseline.json"
    current = tmp_path / "current.json"
    baseline.write_text(json.dumps(report({"a/all": [100.0, 100.1, 99.9]})))
    current.write_text(json.dumps(report({"a/all": [130.0, 130.1, 129.9]})))
    arguments = ["compare", "--baseline", str(baseline), "--current", str(current)]
    assert benchmark.main(arguments) == 0
    assert benchmark.main(arguments + ["--fail-on-regression"]) == 1


def test_noise_mode_prints_the_spread(tmp_path, capsys):
    path = tmp_path / "report.json"
    path.write_text(json.dumps(report({"a/all": [100.0, 102.0, 98.0]})))
    assert benchmark.main(["noise", str(path)]) == 0
    assert "largest spread" in capsys.readouterr().out
