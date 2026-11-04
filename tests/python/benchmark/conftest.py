# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Fixtures of the tests of scripts/benchmark: a fake SeisSol, launcher and Slurm.

The fake SeisSol answers --check with a report, and a run with the files a run
of SeisSol leaves: its summary, epochs, ranks and pinning, with times that
shrink with the ranks. The fake launcher tells it the number of ranks.
"""

import json
import os
import stat
import sys
import textwrap
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[3]
sys.path.insert(0, str(REPO_ROOT / "scripts" / "benchmark"))

FAKE_SEISSOL = r"""#!/usr/bin/env python3
import json
import os
import sys
from pathlib import Path

args = sys.argv[1:]
sets = dict(a[len("--set="):].split("=", 1) for a in args if a.startswith("--set="))
check = [a for a in args if a == "--check" or a.startswith("--check=")]
if os.environ.get("FAKE_SIGNAL"):
    import signal

    os.kill(os.getpid(), signal.SIGILL)
if check:
    report = check[0].split("=", 1)[1] if "=" in check[0] else None
    data = {
        "schema": "seissol-check/1",
        "status": "passed",
        "problems": [],
        "mesh": {"file": "mesh.puml.h5", "format": "puml", "cells": 1000, "vertices": 300},
        "output": {"enabled": [], "writable": True},
        "unknownFields": [],
    }
    if os.environ.get("FAKE_CHECK_FAILS"):
        data["status"] = "failed"
        data["problems"] = ["The file x.yaml that equations.materialfilename names does not exist."]
    if report:
        Path(report).write_text(json.dumps(data))
    sys.exit(1 if os.environ.get("FAKE_CHECK_FAILS") else 0)

ranks = int(os.environ.get("FAKE_RANKS", "1"))
marker = os.environ.get("FAKE_FAIL_ONCE")
if marker and not Path(marker).exists():
    Path(marker).write_text("failed once")
    sys.exit(3)
prefix = Path(sets["output.outputfile"])
prefix.parent.mkdir(parents=True, exist_ok=True)
run_id = os.environ["SEISSOL_RUN_ID"]
total = 100.0 / ranks + 2.0
summary = {
    "schema": "seissol-run-summary/1",
    "runId": run_id,
    "status": "completed",
    "parallel": {"ranks": ranks, "nodes": ranks},
    "mesh": {"cells": {"sum": 1000}},
    "simulation": {
        "endTime": 1.0,
        "currentTime": 1.0,
        "cellUpdates": 50000,
        "duplicateCellUpdates": 100 * (ranks - 1),
    },
    "time": {
        "init": {"total": {"max": 1.0}},
        "simulation": {
            "total": {"max": total},
            "compute": {"max": 0.9 * total},
            "io": {"max": 0.01},
            "outsideKernels": {"max": 0.1 * total, "mean": 0.05 * total},
        },
    },
    "kernels": {"total": {"max": 0.8 * total, "mean": 0.7 * total}, "loadImbalance": 0.125},
    "performance": {
        "cellUpdatesPerSecond": 50000 / total,
        "hardwareFlopPerSecond": 1e9,
        "nonzeroFlopPerSecond": 5e8,
    },
    "memory": {"atEnd": {"hostPeak": {"max": 1.0e9 / ranks}}},
    "build": {"version": "fake", "commit": "0" * 40},
}
Path(f"{prefix}-summary.json").write_text(json.dumps(summary))
epochs = ["epoch,timeFrom,timeTo,wallMin,wallMean,wallMax,runId"]
for epoch in range(4):
    wall = total / 4 * (1.5 if epoch == 0 else 1.0)
    epochs.append(f"{epoch + 1},{epoch / 4},{(epoch + 1) / 4},{wall},{wall},{wall},{run_id}")
Path(f"{prefix}-epochs.csv").write_text("\n".join(epochs) + "\n")
lines = ["rank,host,localRank,cells,simulationTotal,runId"]
pinning = ["hostname,device,rank,localRank,workermask,workernuma,commthread_mask,commthread_numa,nproc,runId"]
for rank in range(ranks):
    lines.append(f'{rank},"node{rank}",0,{1000 // ranks},{total},{run_id}')
    pinning.append(f'"node{rank}","",{rank},0,"0,1,2","0","3","0",{ranks},{run_id}')
Path(f"{prefix}-ranks.csv").write_text("\n".join(lines) + "\n")
Path(f"{prefix}-threadPinning.csv").write_text("\n".join(pinning) + "\n")
print("fake SeisSol ran on", ranks, "ranks")
"""

FAKE_LAUNCHER = r"""#!/usr/bin/env bash
# fake-launch -np N executable arguments...
if [ "$1" = "-np" ]; then
  export FAKE_RANKS="$2"
  shift 2
fi
exec "$@"
"""

FAKE_SBATCH = r"""#!/usr/bin/env bash
# records the submission; answers with a job ID as sbatch --parsable does
dir="$(dirname "$0")"
count=$(( $(cat "$dir/count" 2>/dev/null || echo 0) + 1 ))
echo "$count" > "$dir/count"
echo "$@" >> "$dir/submissions"
echo "$((1000 + count))"
"""

FAKE_SQUEUE = r"""#!/usr/bin/env bash
# the queue: the lines of the file queue, as ID|STATE
dir="$(dirname "$0")"
cat "$dir/queue" 2>/dev/null
"""

FAKE_SACCT = r"""#!/usr/bin/env bash
dir="$(dirname "$0")"
cat "$dir/accounting" 2>/dev/null
"""


def write_executable(path: Path, text: str) -> Path:
    path.write_text(text, encoding="utf-8")
    path.chmod(path.stat().st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)
    return path


@pytest.fixture
def fakes(tmp_path):
    """The fake SeisSol and launcher, in a directory of their own."""
    folder = tmp_path / "fakes"
    folder.mkdir()
    return {
        "seissol": write_executable(folder / "seissol", FAKE_SEISSOL),
        "launcher": write_executable(folder / "fake-launch", FAKE_LAUNCHER),
    }


@pytest.fixture
def fake_slurm(tmp_path, monkeypatch):
    """sbatch, squeue and sacct on the PATH, which record and answer from files."""
    folder = tmp_path / "slurm"
    folder.mkdir()
    write_executable(folder / "sbatch", FAKE_SBATCH)
    write_executable(folder / "squeue", FAKE_SQUEUE)
    write_executable(folder / "sacct", FAKE_SACCT)
    write_executable(
        folder / "scancel",
        '#!/usr/bin/env bash\necho "$@" >> "$(dirname "$0")/cancelled"\n',
    )
    monkeypatch.setenv("PATH", f"{folder}{os.pathsep}{os.environ['PATH']}")
    return folder


@pytest.fixture
def campaign_file(tmp_path, fakes):
    """A campaign of the fake SeisSol on the local machine, as JSON."""

    def make(**changes):
        case = tmp_path / "case"
        case.mkdir(exist_ok=True)
        (case / "parameters.par").write_text("&equations\n/\n")
        config = {
            "version": 1,
            "name": "fake",
            "machine": {
                "scheduler": "local",
                "node": {"cores": 4},
                "layout": {
                    "ranks_per_node": 1,
                    "cpus_per_rank": 4,
                    "threads_per_rank": 3,
                },
                "launcher": f"{fakes['launcher']} -np {{ranks}}",
            },
            "case": {
                "parameter_file": "case/parameters.par",
                "set": {"abortcriteria.endtime": 1.0},
            },
            "builds": [{"name": "fake", "binary": str(fakes["seissol"])}],
            "variants": [
                {"name": "a"},
                {"name": "b", "env": {"SEISSOL_COMMTHREAD": 0}},
            ],
            "scaling": {"nodes": [1, 2, 4], "repetitions": 2},
            "job": {"time_limit": "00:10:00"},
        }
        for key, value in changes.items():
            config[key] = value
        path = tmp_path / "campaign.json"
        path.write_text(json.dumps(config, indent=2))
        return path

    return make


def dedent(text: str) -> str:
    return textwrap.dedent(text).lstrip()
