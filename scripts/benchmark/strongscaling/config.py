# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The campaign file and the machine profiles: reading and validating them.

Both are YAML (with PyYAML) or JSON. Every mapping is checked against the keys
it may have, so that a misspelt key is an error rather than a setting that is
silently ignored.
"""

from __future__ import annotations

import copy
import json
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

MACHINE_DIRECTORY = Path(__file__).resolve().parent.parent / "machines"
SCHEDULERS = ("slurm", "local")


class ConfigError(Exception):
    """A campaign file or a machine profile that cannot be used, and why."""


def load_document(path: Path) -> dict[str, Any]:
    """The mapping a YAML or JSON file holds."""
    try:
        text = path.read_text(encoding="utf-8")
    except OSError as error:
        raise ConfigError(f"cannot read {path}: {error.strerror}") from error
    if path.suffix.lower() in (".yaml", ".yml"):
        try:
            import yaml
        except ImportError as error:
            raise ConfigError(
                f"{path} is YAML, but PyYAML is not installed; install it "
                "(pip install pyyaml) or write the file as JSON"
            ) from error
        try:
            data = yaml.safe_load(text)
        except yaml.YAMLError as error:
            raise ConfigError(f"{path} is not valid YAML: {error}") from error
    elif path.suffix.lower() == ".json":
        try:
            data = json.loads(text)
        except json.JSONDecodeError as error:
            raise ConfigError(f"{path} is not valid JSON: {error}") from error
    else:
        raise ConfigError(f"{path}: use a .yaml, .yml or .json file")
    if not isinstance(data, dict):
        raise ConfigError(f"{path}: the top level has to be a mapping")
    return data


def deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    """A copy of base with overlay merged in; mappings merge, all else is replaced."""
    merged = copy.deepcopy(base)
    for key, value in overlay.items():
        if isinstance(value, dict) and isinstance(merged.get(key), dict):
            merged[key] = deep_merge(merged[key], value)
        else:
            merged[key] = copy.deepcopy(value)
    return merged


def _check_keys(data: Any, allowed: set[str], where: str) -> dict[str, Any]:
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"{where}: expected a mapping, got {type(data).__name__}")
    unknown = sorted(set(data) - allowed)
    if unknown:
        raise ConfigError(
            f"{where}: unknown key(s) {', '.join(unknown)}; "
            f"known are {', '.join(sorted(allowed))}"
        )
    return data


def _positive_int(value: Any, where: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int) or value < 1:
        raise ConfigError(f"{where}: expected a positive integer, got {value!r}")
    return value


def _optional_positive_int(value: Any, where: str) -> int | None:
    return None if value is None else _positive_int(value, where)


def _number(value: Any, where: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        raise ConfigError(f"{where}: expected a number, got {value!r}")
    return float(value)


def as_text(value: Any) -> str:
    """A setting as the text it is passed on as: to the environment or to --set."""
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (list, tuple)):
        return " ".join(as_text(item) for item in value)
    return str(value)


def _text_mapping(data: Any, where: str) -> dict[str, str]:
    if data is None:
        return {}
    if not isinstance(data, dict):
        raise ConfigError(f"{where}: expected a mapping of names to values")
    return {str(key): as_text(value) for key, value in data.items()}


def _overrides(data: Any, where: str) -> dict[str, str]:
    overrides = _text_mapping(data, where)
    for name in overrides:
        if not re.fullmatch(r"[A-Za-z_][\w]*\.[A-Za-z_][\w]*", name):
            raise ConfigError(
                f"{where}: {name!r} is not a parameter as section.field, "
                "e.g. abortcriteria.endtime"
            )
    return {name.lower(): value for name, value in overrides.items()}


def parse_duration(value: Any, where: str) -> int:
    """Seconds of a duration given as minutes, or as [D-]HH:MM[:SS] as Slurm takes it."""
    if isinstance(value, bool):
        raise ConfigError(f"{where}: expected a duration, got {value!r}")
    if isinstance(value, (int, float)):
        if value <= 0:
            raise ConfigError(f"{where}: a duration has to be positive")
        return int(round(value * 60))
    match = re.fullmatch(r"(?:(\d+)-)?(\d+):(\d{2})(?::(\d{2}))?", str(value).strip())
    if match is None:
        raise ConfigError(
            f"{where}: {value!r} is no duration; give minutes or [D-]HH:MM[:SS]"
        )
    days, hours, minutes, seconds = (int(part or 0) for part in match.groups())
    total = ((days * 24 + hours) * 60 + minutes) * 60 + seconds
    if total <= 0:
        raise ConfigError(f"{where}: a duration has to be positive")
    return total


def format_duration(seconds: int) -> str:
    """A duration as Slurm takes it: [D-]HH:MM:SS."""
    days, rest = divmod(int(seconds), 86400)
    hours, rest = divmod(rest, 3600)
    minutes, secs = divmod(rest, 60)
    text = f"{hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{days}-{text}" if days else text


def parse_size(value: Any, where: str) -> float:
    """Bytes of a size given as a number of bytes or with a unit, e.g. 96G or 512 GiB."""
    if isinstance(value, bool):
        raise ConfigError(f"{where}: expected a size, got {value!r}")
    if isinstance(value, (int, float)):
        return float(value)
    match = re.fullmatch(r"\s*([0-9.]+)\s*([kKmMgGtT]?)(i?[bB])?\s*", str(value))
    if match is None:
        raise ConfigError(f"{where}: {value!r} is no size; e.g. 96G or 512GiB")
    factor = {"": 1, "k": 2**10, "m": 2**20, "g": 2**30, "t": 2**40}[
        match.group(2).lower()
    ]
    return float(match.group(1)) * factor


# ---------------------------------------------------------------------------
# Machine profiles
# ---------------------------------------------------------------------------

MACHINE_KEYS = {
    "name",
    "description",
    "scheduler",
    "node",
    "layout",
    "slurm",
    "launcher",
    "wrapper",
    "setup",
    "env",
    "python",
}
NODE_KEYS = {"cores", "threads_per_core", "gpus", "memory", "numa_domains"}
LAYOUT_KEYS = {"ranks_per_node", "threads_per_rank", "cpus_per_rank", "gpus_per_rank"}
SLURM_KEYS = {
    "partition",
    "account",
    "qos",
    "constraint",
    "exclusive",
    "cpus_per_task",
    "extra",
}


@dataclass
class Layout:
    """How the ranks of a run sit on a node."""

    ranks_per_node: int
    #: the OpenMP workers of a rank, OMP_NUM_THREADS
    threads_per_rank: int
    #: the CPUs a rank holds (the workers and the communication thread), --cpus-per-task
    cpus_per_rank: int
    gpus_per_rank: int = 0

    def to_dict(self) -> dict[str, int]:
        return {
            "ranks_per_node": self.ranks_per_node,
            "threads_per_rank": self.threads_per_rank,
            "cpus_per_rank": self.cpus_per_rank,
            "gpus_per_rank": self.gpus_per_rank,
        }


def partition_for(slurm: dict[str, Any], nodes: int) -> str | None:
    """The partition of a job on @p nodes nodes: the first whose max_nodes it fits."""
    partition = slurm.get("partition")
    if not isinstance(partition, list):
        return partition
    for rule in partition:
        if nodes <= rule["max_nodes"]:
            return rule["partition"]
    raise ConfigError(
        f"no partition of the machine takes a job of {nodes} nodes; the largest takes "
        f"{partition[-1]['max_nodes']}"
    )


@dataclass
class Machine:
    name: str
    scheduler: str
    cores_per_node: int
    threads_per_core: int
    gpus_per_node: int
    memory_per_node: float | None
    numa_domains: int | None
    layout: Layout
    #: the command a run is started with, before the wrapper and the executable;
    #: {nodes}, {ranks}, {ranks_per_node}, {cpus_per_rank}, {threads_per_rank} and
    #: {gpus_per_rank} stand for the values of the run
    launcher: str
    #: a script started per rank in front of the executable, e.g. for its pinning
    wrapper: str | None
    #: shell lines that set up the environment of a job, e.g. module loads
    setup: list[str]
    env: dict[str, str]
    slurm: dict[str, Any]
    #: how the jobs run this tool, for the report job
    python: str
    description: str = ""
    source: str = ""

    def to_dict(self) -> dict[str, Any]:
        return {
            "name": self.name,
            "description": self.description,
            "scheduler": self.scheduler,
            "node": {
                "cores": self.cores_per_node,
                "threads_per_core": self.threads_per_core,
                "gpus": self.gpus_per_node,
                "memory": self.memory_per_node,
                "numa_domains": self.numa_domains,
            },
            "layout": self.layout.to_dict(),
            "slurm": self.slurm,
            "launcher": self.launcher,
            "wrapper": self.wrapper,
            "setup": self.setup,
            "env": self.env,
            "python": self.python,
            "source": self.source,
        }


DEFAULT_LAUNCHERS = {
    "slurm": "srun --nodes={nodes} --ntasks={ranks} --ntasks-per-node={ranks_per_node} "
    "--cpus-per-task={cpus_per_rank}",
    "local": "mpirun -np {ranks} --oversubscribe",
}


def _resolve_relative(text: str | None, base: Path) -> str | None:
    """A path relative to the file that names it; anything else is left alone."""
    if text is None:
        return None
    path = Path(text)
    if path.is_absolute() or "/" not in text and not (base / text).exists():
        # a command found on the PATH, e.g. seissol-launch
        return text
    return str((base / path).resolve())


def find_machine_profile(name: str, base: Path) -> Path:
    """The file of the profile @p name: a path, or a profile of machines/."""
    candidate = Path(name)
    if candidate.suffix:
        path = candidate if candidate.is_absolute() else base / candidate
        if not path.exists():
            raise ConfigError(f"the machine profile {path} does not exist")
        return path
    for suffix in (".yaml", ".yml", ".json"):
        path = MACHINE_DIRECTORY / f"{name}{suffix}"
        if path.exists():
            return path
    known = sorted(p.stem for p in MACHINE_DIRECTORY.glob("*.*"))
    raise ConfigError(
        f"there is no machine profile {name!r}; known are {', '.join(known)}, "
        "or give the path of a profile"
    )


def parse_machine(data: dict[str, Any], source: Path) -> Machine:
    where = f"machine profile {source}"
    data = _check_keys(data, MACHINE_KEYS, where)
    scheduler = data.get("scheduler", "slurm")
    if scheduler not in SCHEDULERS:
        raise ConfigError(
            f"{where}: scheduler {scheduler!r} is not supported; "
            f"supported are {', '.join(SCHEDULERS)}"
        )
    node = _check_keys(data.get("node"), NODE_KEYS, f"{where}, node")
    if "cores" not in node:
        raise ConfigError(f"{where}: node.cores, the cores of a node, is missing")
    cores = _positive_int(node["cores"], f"{where}, node.cores")
    threads_per_core = _positive_int(
        node.get("threads_per_core", 1), f"{where}, node.threads_per_core"
    )
    gpus = node.get("gpus", 0) or 0
    if gpus:
        gpus = _positive_int(gpus, f"{where}, node.gpus")
    memory = node.get("memory")
    memory_bytes = (
        None if memory is None else parse_size(memory, f"{where}, node.memory")
    )
    numa = _optional_positive_int(
        node.get("numa_domains"), f"{where}, node.numa_domains"
    )

    layout_data = _check_keys(data.get("layout"), LAYOUT_KEYS, f"{where}, layout")
    ranks_per_node = _positive_int(
        layout_data.get("ranks_per_node", 1), f"{where}, layout.ranks_per_node"
    )
    cpus_default = max(1, cores * threads_per_core // ranks_per_node)
    cpus_per_rank = _positive_int(
        layout_data.get("cpus_per_rank", cpus_default), f"{where}, layout.cpus_per_rank"
    )
    # one CPU of a rank is left to its communication thread
    threads_per_rank = _positive_int(
        layout_data.get("threads_per_rank", max(1, cpus_per_rank - 1)),
        f"{where}, layout.threads_per_rank",
    )
    gpus_per_rank = layout_data.get("gpus_per_rank", 0) or 0
    if gpus_per_rank:
        gpus_per_rank = _positive_int(gpus_per_rank, f"{where}, layout.gpus_per_rank")
    layout = Layout(ranks_per_node, threads_per_rank, cpus_per_rank, gpus_per_rank)

    slurm = _check_keys(data.get("slurm"), SLURM_KEYS, f"{where}, slurm")
    slurm = dict(slurm)
    slurm.setdefault("exclusive", True)
    # whether a job asks for the CPUs of a rank; where Slurm counts the hyperthreads of a node as
    # its CPUs and where it counts its cores differs between sites
    slurm.setdefault("cpus_per_task", True)
    partition = slurm.get("partition")
    if isinstance(partition, list):
        rules = []
        for index, rule in enumerate(partition):
            rule = _check_keys(
                rule, {"max_nodes", "partition"}, f"{where}, slurm.partition[{index}]"
            )
            if "partition" not in rule or "max_nodes" not in rule:
                raise ConfigError(
                    f"{where}, slurm.partition[{index}]: give partition and max_nodes"
                )
            rules.append(
                {
                    "max_nodes": _positive_int(
                        rule["max_nodes"],
                        f"{where}, slurm.partition[{index}].max_nodes",
                    ),
                    "partition": str(rule["partition"]),
                }
            )
        slurm["partition"] = sorted(rules, key=lambda rule: rule["max_nodes"])
    elif partition is not None and not isinstance(partition, str):
        raise ConfigError(
            f"{where}, slurm.partition: a partition, or a list of "
            "{max_nodes: N, partition: name} by the size of a job"
        )
    extra = slurm.get("extra") or []
    if not isinstance(extra, list) or not all(isinstance(e, str) for e in extra):
        raise ConfigError(f"{where}, slurm.extra: expected a list of sbatch options")
    slurm["extra"] = list(extra)

    setup = data.get("setup") or []
    if isinstance(setup, str):
        setup = setup.splitlines()
    if not isinstance(setup, list) or not all(isinstance(line, str) for line in setup):
        raise ConfigError(f"{where}, setup: expected shell lines")

    launcher = data.get("launcher") or DEFAULT_LAUNCHERS[scheduler]
    if not isinstance(launcher, str):
        raise ConfigError(f"{where}, launcher: expected a command line")

    return Machine(
        name=str(data.get("name", source.stem)),
        description=str(data.get("description", "")),
        scheduler=scheduler,
        cores_per_node=cores,
        threads_per_core=threads_per_core,
        gpus_per_node=gpus,
        memory_per_node=memory_bytes,
        numa_domains=numa,
        layout=layout,
        launcher=launcher,
        wrapper=_resolve_relative(data.get("wrapper"), source.parent),
        setup=list(setup),
        env=_text_mapping(data.get("env"), f"{where}, env"),
        slurm=slurm,
        python=str(data.get("python", "python3")),
        source=str(source),
    )


def load_machine(spec: Any, base: Path) -> Machine:
    """The machine a campaign names: a profile, possibly with settings of its own."""
    if isinstance(spec, str):
        path = find_machine_profile(spec, base)
        return parse_machine(load_document(path), path)
    if isinstance(spec, dict):
        spec = dict(spec)
        profile = spec.pop("profile", None)
        if profile is None:
            return parse_machine(spec, base / "campaign")
        path = find_machine_profile(str(profile), base)
        return parse_machine(deep_merge(load_document(path), spec), path)
    raise ConfigError("machine: name a profile, or give a mapping with 'profile'")


# ---------------------------------------------------------------------------
# Campaigns
# ---------------------------------------------------------------------------

CAMPAIGN_KEYS = {
    "version",
    "name",
    "directory",
    "machine",
    "case",
    "builds",
    "variants",
    "scaling",
    "job",
}
CASE_KEYS = {"parameter_file", "set", "epoch_interval", "cells", "env"}
BUILD_KEYS = {"name", "binary", "variants", "preset", "toolchain", "set", "select"}
VARIANT_KEYS = {
    "name",
    "env",
    "set",
    "ranks_per_node",
    "threads_per_rank",
    "cpus_per_rank",
    "launcher_args",
}
SCALING_KEYS = {
    "nodes",
    "repetitions",
    "min_cells_per_rank",
    "max_cells_per_node",
    "order",
}
NODE_RANGE_KEYS = {"min", "max", "factor"}
JOB_KEYS = {
    "time_limit",
    "run_time_limit",
    "max_attempts",
    "report_job",
    "check_in_job",
    "warmup_epochs",
}


@dataclass
class Case:
    parameter_file: Path
    #: --set overrides of every run, section.field -> value
    overrides: dict[str, str]
    epoch_interval: float | None
    #: the cells of the mesh, if known without checking it
    cells: int | None
    env: dict[str, str]


@dataclass
class Build:
    name: str
    #: an executable built already
    binary: Path | None = None
    #: a variant-set file of compile_variants.py, which builds the executable
    variants_file: Path | None = None
    preset: str | None = None
    toolchain: str | None = None
    #: --set options of compile_variants.py
    overrides: list[str] = field(default_factory=list)
    #: the executable to take of a build of several
    select: str | None = None


@dataclass
class Variant:
    name: str
    env: dict[str, str] = field(default_factory=dict)
    overrides: dict[str, str] = field(default_factory=dict)
    ranks_per_node: int | None = None
    threads_per_rank: int | None = None
    cpus_per_rank: int | None = None
    launcher_args: list[str] = field(default_factory=list)

    def layout(self, machine: Machine) -> Layout:
        """The layout of a run of this variant on @p machine."""
        base = machine.layout
        ranks = self.ranks_per_node or base.ranks_per_node
        if self.cpus_per_rank is not None:
            cpus = self.cpus_per_rank
        elif self.ranks_per_node is not None:
            cpus = max(1, base.cpus_per_rank * base.ranks_per_node // ranks)
        else:
            cpus = base.cpus_per_rank
        if self.threads_per_rank is not None:
            threads = self.threads_per_rank
        elif self.cpus_per_rank is not None or self.ranks_per_node is not None:
            threads = max(1, cpus - (base.cpus_per_rank - base.threads_per_rank))
        else:
            threads = base.threads_per_rank
        return Layout(ranks, threads, cpus, base.gpus_per_rank)


@dataclass
class Scaling:
    nodes: list[int] | None
    node_range: tuple[int | None, int, int] | None
    repetitions: int
    min_cells_per_rank: float | None
    max_cells_per_node: float | None
    order: str


@dataclass
class JobSettings:
    time_limit: int
    run_time_limit: int | None
    max_attempts: int
    report_job: bool
    check_in_job: bool
    warmup_epochs: int


@dataclass
class Campaign:
    name: str
    path: Path
    directory: Path
    machine: Machine
    case: Case
    builds: list[Build]
    variants: list[Variant]
    scaling: Scaling
    job: JobSettings
    raw: dict[str, Any]


def _name(value: Any, where: str) -> str:
    text = str(value) if value is not None else ""
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.+-]*", text):
        raise ConfigError(
            f"{where}: {value!r} is no name; use letters, digits and _.+- only"
        )
    return text


def _parse_case(data: Any, base: Path) -> Case:
    data = _check_keys(data, CASE_KEYS, "case")
    if "parameter_file" not in data:
        raise ConfigError(
            "case.parameter_file, the parameter file of the runs, is missing"
        )
    parameter_file = Path(str(data["parameter_file"]))
    if not parameter_file.is_absolute():
        parameter_file = (base / parameter_file).resolve()
    overrides = _overrides(data.get("set"), "case.set")
    if "output.outputfile" in overrides:
        raise ConfigError(
            "case.set: output.outputfile is set per run, in the run's directory"
        )
    epoch_interval = data.get("epoch_interval")
    if epoch_interval is not None:
        epoch_interval = _number(epoch_interval, "case.epoch_interval")
        if epoch_interval <= 0:
            raise ConfigError("case.epoch_interval has to be positive")
    cells = data.get("cells")
    if cells is not None:
        cells = _positive_int(cells, "case.cells")
    return Case(
        parameter_file=parameter_file,
        overrides=overrides,
        epoch_interval=epoch_interval,
        cells=cells,
        env=_text_mapping(data.get("env"), "case.env"),
    )


def _parse_builds(data: Any, base: Path) -> list[Build]:
    if data is None:
        raise ConfigError("builds, the executables to compare, is missing")
    if not isinstance(data, list) or not data:
        raise ConfigError("builds: expected a list of builds")
    builds = []
    for index, entry in enumerate(data):
        where = f"builds[{index}]"
        entry = _check_keys(entry, BUILD_KEYS, where)
        name = _name(entry.get("name", f"build{index}"), f"{where}.name")
        binary = entry.get("binary")
        variants = entry.get("variants")
        if (binary is None) == (variants is None):
            raise ConfigError(
                f"{where}: give either binary, an executable built already, or "
                "variants, a variant-set file of compile_variants.py"
            )
        overrides = entry.get("set") or []
        if isinstance(overrides, dict):
            overrides = [f"{key}={as_text(value)}" for key, value in overrides.items()]
        if not isinstance(overrides, list):
            raise ConfigError(f"{where}.set: expected KEY=VALUE entries")
        build = Build(
            name=name,
            preset=entry.get("preset"),
            toolchain=entry.get("toolchain"),
            overrides=[str(item) for item in overrides],
            select=entry.get("select"),
        )
        if binary is not None:
            path = Path(str(binary)).expanduser()
            build.binary = path if path.is_absolute() else (base / path).resolve()
        else:
            path = Path(str(variants)).expanduser()
            build.variants_file = (
                path if path.is_absolute() else (base / path).resolve()
            )
        builds.append(build)
    names = [build.name for build in builds]
    if len(set(names)) != len(names):
        raise ConfigError("builds: two builds have the same name")
    return builds


def _parse_variants(data: Any) -> list[Variant]:
    if data is None:
        return [Variant(name="default")]
    if not isinstance(data, list) or not data:
        raise ConfigError("variants: expected a list of variants")
    variants = []
    for index, entry in enumerate(data):
        where = f"variants[{index}]"
        entry = _check_keys(entry, VARIANT_KEYS, where)
        overrides = _overrides(entry.get("set"), f"{where}.set")
        if "output.outputfile" in overrides:
            raise ConfigError(f"{where}.set: output.outputfile is set per run")
        launcher_args = entry.get("launcher_args") or []
        if isinstance(launcher_args, str):
            launcher_args = launcher_args.split()
        variants.append(
            Variant(
                name=_name(entry.get("name", f"variant{index}"), f"{where}.name"),
                env=_text_mapping(entry.get("env"), f"{where}.env"),
                overrides=overrides,
                ranks_per_node=_optional_positive_int(
                    entry.get("ranks_per_node"), f"{where}.ranks_per_node"
                ),
                threads_per_rank=_optional_positive_int(
                    entry.get("threads_per_rank"), f"{where}.threads_per_rank"
                ),
                cpus_per_rank=_optional_positive_int(
                    entry.get("cpus_per_rank"), f"{where}.cpus_per_rank"
                ),
                launcher_args=[str(arg) for arg in launcher_args],
            )
        )
    names = [variant.name for variant in variants]
    if len(set(names)) != len(names):
        raise ConfigError("variants: two variants have the same name")
    return variants


def _parse_scaling(data: Any) -> Scaling:
    data = _check_keys(data, SCALING_KEYS, "scaling")
    nodes = data.get("nodes")
    explicit = None
    node_range = None
    if isinstance(nodes, list):
        explicit = sorted({_positive_int(n, "scaling.nodes") for n in nodes})
        if not explicit:
            raise ConfigError("scaling.nodes: give at least one node count")
    elif isinstance(nodes, dict):
        nodes = _check_keys(nodes, NODE_RANGE_KEYS, "scaling.nodes")
        if "max" not in nodes:
            raise ConfigError("scaling.nodes.max, the largest node count, is missing")
        low = nodes.get("min")
        low = None if low in (None, "auto") else _positive_int(low, "scaling.nodes.min")
        high = _positive_int(nodes["max"], "scaling.nodes.max")
        factor = _positive_int(nodes.get("factor", 2), "scaling.nodes.factor")
        if factor < 2:
            raise ConfigError("scaling.nodes.factor has to be 2 or more")
        if low is not None and low > high:
            raise ConfigError("scaling.nodes: min is larger than max")
        node_range = (low, high, factor)
    else:
        raise ConfigError(
            "scaling.nodes: give a list of node counts, or a range "
            "{min: 1, max: 64, factor: 2}"
        )
    order = data.get("order", "interleaved")
    if order not in ("interleaved", "blocked"):
        raise ConfigError("scaling.order: interleaved or blocked")
    min_cells = data.get("min_cells_per_rank")
    max_cells = data.get("max_cells_per_node")
    return Scaling(
        nodes=explicit,
        node_range=node_range,
        repetitions=_positive_int(data.get("repetitions", 1), "scaling.repetitions"),
        min_cells_per_rank=(
            None
            if min_cells is None
            else _number(min_cells, "scaling.min_cells_per_rank")
        ),
        max_cells_per_node=(
            None
            if max_cells is None
            else _number(max_cells, "scaling.max_cells_per_node")
        ),
        order=order,
    )


def _parse_job(data: Any) -> JobSettings:
    data = _check_keys(data, JOB_KEYS, "job")
    if "time_limit" not in data:
        raise ConfigError("job.time_limit, the time limit of a job, is missing")
    run_limit = data.get("run_time_limit")
    return JobSettings(
        time_limit=parse_duration(data["time_limit"], "job.time_limit"),
        run_time_limit=(
            None
            if run_limit is None
            else parse_duration(run_limit, "job.run_time_limit")
        ),
        max_attempts=_positive_int(data.get("max_attempts", 3), "job.max_attempts"),
        report_job=bool(data.get("report_job", False)),
        check_in_job=bool(data.get("check_in_job", True)),
        warmup_epochs=int(data.get("warmup_epochs", 1)),
    )


def parse_campaign(data: dict[str, Any], path: Path) -> Campaign:
    data = _check_keys(data, CAMPAIGN_KEYS, f"campaign {path}")
    version = data.get("version", 1)
    if version != 1:
        raise ConfigError(f"campaign {path}: version {version!r} is not supported")
    base = path.resolve().parent
    if "name" not in data:
        raise ConfigError("name, the name of the campaign, is missing")
    name = _name(data["name"], "name")
    directory = Path(str(data.get("directory", name))).expanduser()
    if not directory.is_absolute():
        directory = (base / directory).resolve()
    if "machine" not in data:
        raise ConfigError("machine, the machine profile, is missing")
    machine = load_machine(data["machine"], base)
    campaign = Campaign(
        name=name,
        path=path.resolve(),
        directory=directory,
        machine=machine,
        case=_parse_case(data.get("case"), base),
        builds=_parse_builds(data.get("builds"), base),
        variants=_parse_variants(data.get("variants")),
        scaling=_parse_scaling(data.get("scaling")),
        job=_parse_job(data.get("job")),
        raw=data,
    )
    return campaign


def load_campaign(path: Path) -> Campaign:
    """The campaign of a campaign file, or of the copy a campaign directory keeps."""
    if path.is_dir():
        copy_path = path / "campaign.json"
        if not copy_path.exists():
            raise ConfigError(
                f"{path} is a directory, but no campaign directory: it has no campaign.json"
            )
        stored = load_document(copy_path)
        return parse_campaign(stored["config"], Path(stored["path"]))
    return parse_campaign(load_document(path), path)
