# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""A machine profile drafted from a job script that is known to run SeisSol well.

The #SBATCH options give the layout of a node and the scheduler settings, the
exports the environment, the module loads and other set-up lines the setup,
and the line that starts SeisSol the launcher and the wrapper. The counts of
nodes and ranks become placeholders, which every run fills in. What a job
script does not say, the cores and the memory of a node, the draft guesses;
it is meant to be read and corrected.
"""

from __future__ import annotations

import json
import re
import shlex
from pathlib import Path
from typing import Any

#: #SBATCH options that the campaign sets per job
MANAGED_OPTIONS = {
    "nodes",
    "ntasks",
    "time",
    "job-name",
    "output",
    "error",
    "chdir",
    "mail-type",
    "mail-user",
    "no-requeue",
    "requeue",
}
SHORT_OPTIONS = {
    "N": "nodes",
    "n": "ntasks",
    "c": "cpus-per-task",
    "p": "partition",
    "A": "account",
    "q": "qos",
    "C": "constraint",
    "t": "time",
    "J": "job-name",
    "o": "output",
    "e": "error",
    "D": "chdir",
}
SETUP_COMMANDS = (
    "module",
    "ml",
    "source",
    ".",
    "ulimit",
    "unset",
    "spack",
    "conda",
    "cd",
)
LAUNCHERS = ("srun", "mpirun", "mpiexec", "mpiexec.hydra", "aprun", "jsrun", "flux")


def sbatch_options(line: str) -> list[tuple[str, str | None]]:
    """The options of an #SBATCH line, by their long names."""
    words = shlex.split(line.split("#SBATCH", 1)[1], comments=False, posix=True)
    options = []
    index = 0
    while index < len(words):
        word = words[index]
        if word.startswith("#"):
            break
        if word.startswith("--"):
            name, _, value = word[2:].partition("=")
            options.append((name, value if _ else None))
        elif word.startswith("-") and len(word) >= 2:
            name = SHORT_OPTIONS.get(word[1], word[1])
            value = word[2:] or None
            if (
                value is None
                and index + 1 < len(words)
                and not words[index + 1].startswith("-")
            ):
                index += 1
                value = words[index]
            options.append((name, value))
        index += 1
    return options


def _count(value: str | None) -> int | None:
    match = re.match(r"\d+", value or "")
    return int(match.group(0)) if match else None


def launch_line(lines: list[str]) -> tuple[str, str | None] | None:
    """The command that starts SeisSol (up to the wrapper), and the wrapper."""
    for line in lines:
        stripped = line.strip()
        if not stripped or stripped.startswith("#"):
            continue
        words = shlex.split(stripped, posix=False)
        # leading assignments, e.g. VAR=1 srun ...
        start = 0
        while start < len(words) and re.fullmatch(r"[A-Za-z_]\w*=.*", words[start]):
            start += 1
        if start >= len(words) or Path(words[start]).name not in LAUNCHERS:
            continue
        candidates = [
            i
            for i in range(start + 1, len(words))
            if "seissol" in Path(words[i]).name.lower()
        ]
        if not candidates:
            continue
        binary = candidates[-1]
        wrapper = None
        end = binary
        before = words[binary - 1] if binary - 1 > start else ""
        if before and not before.startswith("-") and "=" not in before:
            wrapper = before
            end = binary - 1
        return " ".join(words[start:end]), wrapper
    return None


def placeholders(launcher: str) -> str:
    """The launcher with its counts of nodes and ranks replaced by the placeholders of a run."""
    words = launcher.split()
    program = Path(words[0]).name if words else ""
    replacements = [
        (r"--nodes[= ]\d+", "--nodes={nodes}"),
        (r"(?<!\S)-N ?\d+", "--nodes={nodes}"),
        (r"--ntasks-per-node[= ]\d+", "--ntasks-per-node={ranks_per_node}"),
        (r"--ntasks[= ]\d+", "--ntasks={ranks}"),
        (r"--cpus-per-task[= ]\d+", "--cpus-per-task={cpus_per_rank}"),
        (r"(?<!\S)-c ?\d+", "--cpus-per-task={cpus_per_rank}"),
    ]
    if program == "srun":
        replacements.append((r"(?<!\S)-n ?\d+", "--ntasks={ranks}"))
    else:
        replacements.append((r"(?<!\S)-np? ?\d+", "-n {ranks}"))
    for pattern, replacement in replacements:
        launcher = re.sub(pattern, replacement, launcher)
    words = launcher.split()
    if program == "srun":
        # explicit, as the variants of a job may lay out their ranks differently
        missing = [
            f"{option}={value}"
            for option, value in (
                ("--nodes", "{nodes}"),
                ("--ntasks", "{ranks}"),
                ("--ntasks-per-node", "{ranks_per_node}"),
            )
            if not re.search(rf"{option}=", launcher)
        ]
        words[1:1] = missing
    elif (
        program in ("mpirun", "mpiexec", "mpiexec.hydra") and "{ranks}" not in launcher
    ):
        words.insert(1, "-n {ranks}")
    return " ".join(words)


def draft_profile(script: Path, name: str | None = None) -> dict[str, Any]:
    """A machine profile from the job script @p script."""
    lines = script.read_text(encoding="utf-8").splitlines()
    slurm: dict[str, Any] = {"exclusive": False, "extra": []}
    layout: dict[str, Any] = {}
    node: dict[str, Any] = {}
    env: dict[str, str] = {}
    setup: list[str] = []
    for line in lines:
        stripped = line.strip()
        if stripped.startswith("#SBATCH"):
            for option, value in sbatch_options(stripped):
                if option in MANAGED_OPTIONS:
                    continue
                if option == "ntasks-per-node":
                    layout["ranks_per_node"] = _count(value)
                elif option == "cpus-per-task":
                    layout["cpus_per_rank"] = _count(value)
                elif option in ("gpus-per-node", "gres"):
                    gpus = _count((value or "").split(":")[-1])
                    if gpus:
                        node["gpus"] = gpus
                elif option in ("partition", "account", "qos", "constraint"):
                    slurm[option] = value
                elif option == "exclusive":
                    slurm["exclusive"] = True
                else:
                    slurm["extra"].append(
                        f"--{option}" + (f"={value}" if value is not None else "")
                    )
            continue
        if not stripped or stripped.startswith("#"):
            continue
        export = re.fullmatch(r"export\s+([A-Za-z_]\w*)=(.*)", stripped)
        if export:
            variable, value = export.group(1), export.group(2).strip()
            if len(value) >= 2 and value[0] == value[-1] and value[0] in "'\"":
                value = value[1:-1]
            if variable == "OMP_NUM_THREADS" and value.isdigit():
                layout["threads_per_rank"] = int(value)
            elif "$" not in value:
                env[variable] = value
            else:
                setup.append(stripped)
            continue
        first = stripped.split()[0]
        if first in SETUP_COMMANDS or re.fullmatch(r"[A-Za-z_]\w*=\S*", stripped):
            setup.append(stripped)
    launch = launch_line(lines)
    ranks = layout.get("ranks_per_node") or 1
    cpus = layout.get("cpus_per_rank")
    threads = layout.get("threads_per_rank")
    if cpus is None:
        cpus = (threads + 1) if threads else 1
    if node.get("gpus"):
        layout["gpus_per_rank"] = max(1, node["gpus"] // ranks)
    profile: dict[str, Any] = {
        "name": name or script.stem,
        "description": f"drafted from {script.name}; check the node before using it",
        "scheduler": "slurm",
        "node": dict(
            {"cores": cpus * ranks, "threads_per_core": 1},
            **({"gpus": node["gpus"]} if node.get("gpus") else {}),
        ),
        "layout": {
            "ranks_per_node": ranks,
            "cpus_per_rank": cpus,
            "threads_per_rank": threads or max(1, cpus - 1),
            **(
                {"gpus_per_rank": layout["gpus_per_rank"]}
                if "gpus_per_rank" in layout
                else {}
            ),
        },
        "slurm": {
            key: value for key, value in slurm.items() if value not in (None, [], "")
        },
        "setup": setup,
        "env": env,
    }
    if launch is not None:
        launcher, wrapper = launch
        profile["launcher"] = placeholders(launcher)
        if wrapper:
            profile["wrapper"] = wrapper
    return profile


def render(profile: dict[str, Any], source: Path) -> str:
    header = (
        f"# A machine profile drafted from {source}.\n"
        "# node.cores is the CPUs the job script gives a node; correct it to the cores of\n"
        "# a node (and set threads_per_core and memory), and check the layout.\n"
    )
    try:
        import yaml

        return header + yaml.safe_dump(
            profile, sort_keys=False, default_flow_style=False
        )
    except ImportError:
        return json.dumps(profile, indent=2) + "\n"
