#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""Build several SeisSol variants from a single declarative variant-set file.

SeisSol bakes ``EQUATIONS``, ``ORDER`` and ``PRECISION`` (and friends) into the
generated kernels at configure time, so one build tree produces exactly one
variant.  The binaries encode the variant in their name
(``seissol-<eq>[-<mech>]-o<order>-f{32,64}[-s<multisim>][-<suffix>]``) and are
meant to be installed side by side.  This tool is a *thin orchestrator*: it
expands a variant matrix, configures/builds/installs each variant in its own
build tree, and reports the result.  Dependency tracking stays with CMake/Ninja.

Compilers live in a named ``toolchains:`` registry and are selected like any
other matrix axis (``toolchain: [gcc-13, clang-19]``) or pinned from the command
line with ``--toolchain``.  The latter is what CI uses: the compiler there is
tied to the container image, so the workflow matrix owns it and the variant-set
file stays compiler-agnostic.

After a successful build each variant can run its own CTest suite
(``--test``/``testing:``), including the smoke tests registered by SeisSol's
CMake.  Per-test results are read back from CTest's JUnit output so that failing
test names reach the JSON report and the step summary.

Scope: core (declarative matrix, per-variant build dir, side-by-side install,
fail isolation) + operational (cross-variant parallelism with a core budget,
--dry-run, JSON report, --step-summary, --rerun-failed, predicted-binary check).

Deliberately *out* of scope for now (ask if wanted): ccache launcher wiring,
spec-hash incremental skipping, cluster module-load hooks.
"""

from __future__ import annotations

import argparse
import concurrent.futures
import copy
import hashlib
import itertools
import json
import os
import shlex
import shutil
import subprocess
import sys
import threading
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, NoReturn

REPORT_SCHEMA = "seissol-build-report/2"

# Mapping: canonical variant field -> CMake cache variable.
# Only fields present on a variant are passed to CMake.
FIELD_TO_CMAKE: dict[str, str] = {
    "equations": "EQUATIONS",
    "order": "ORDER",
    "precision": "PRECISION",
    "mechanisms": "NUMBER_OF_MECHANISMS",
    "multisim": "NUMBER_OF_FUSED_SIMULATIONS",
    "host_arch": "HOST_ARCH",
    "device_backend": "DEVICE_BACKEND",
    "device_arch": "DEVICE_ARCH",
    "device_codegen": "DEVICE_CODEGEN",
    "build_type": "CMAKE_BUILD_TYPE",
    "custom_suffix": "CUSTOM_BINARY_SUFFIX",
}

# Fields that identify a variant (must resolve to something).
REQUIRED_FIELDS = ("equations", "order", "precision")

# Settings that cascade global-defaults -> preset-defaults -> toolchain -> variant.
META_FIELDS = (
    "source_dir",
    "build_root",
    "install_prefix",
    "generator",
    "jobs",
    "toolchain",
)

# Path-valued meta fields in which {placeholder} expansion is performed.
TEMPLATED_META_FIELDS = ("build_root", "install_prefix")

# Defaults for the per-variant `testing:` block.
TESTING_DEFAULTS: dict[str, Any] = {
    "enabled": False,
    "build_only": False,
    "jobs": 0,
    "labels": [],
    "exclude_labels": [],
    "ctest_args": [],
    "timeout": None,
}

# CMake cache entries implied by `testing.enabled`. Anything the config sets
# explicitly under `defines:` wins over these.
TESTING_DEFINES: dict[str, str] = {
    "TESTING": "ON",
    "TESTING_GENERATED": "ON",
}

# ---------------------------------------------------------------------------
# Console helpers
# ---------------------------------------------------------------------------

_PRINT_LOCK = threading.Lock()
_USE_COLOR = sys.stdout.isatty() and os.environ.get("NO_COLOR") is None


def _c(code: str, text: str) -> str:
    return f"\033[{code}m{text}\033[0m" if _USE_COLOR else text


def green(t: str) -> str:
    return _c("32", t)


def red(t: str) -> str:
    return _c("31", t)


def yellow(t: str) -> str:
    return _c("33", t)


def dim(t: str) -> str:
    return _c("2", t)


def log(msg: str) -> None:
    with _PRINT_LOCK:
        print(msg, flush=True)


def warn(msg: str) -> None:
    log(yellow(f"warning: {msg}"))


def die(msg: str) -> NoReturn:
    print(red(f"error: {msg}"), file=sys.stderr, flush=True)
    raise SystemExit(2)


# ---------------------------------------------------------------------------
# Config loading
# ---------------------------------------------------------------------------


def load_config(path: Path) -> dict[str, Any]:
    text = path.read_text(encoding="utf-8")
    suffix = path.suffix.lower()
    if suffix in (".yaml", ".yml"):
        try:
            import yaml  # type: ignore
        except ImportError:
            die(
                f"{path} is YAML but PyYAML is not installed. "
                "Install it (pip install pyyaml) or use a .json config."
            )
        data = yaml.safe_load(text)
    elif suffix == ".json":
        data = json.loads(text)
    else:
        die(f"unsupported config extension {suffix!r} (use .yaml/.yml/.json)")
    if not isinstance(data, dict):
        die(f"{path}: top level must be a mapping")
    # Path settings are anchored to this directory, see _anchor().
    data["_config_dir"] = str(path.resolve().parent)
    return data


def _deep_merge(base: dict[str, Any], overlay: dict[str, Any]) -> dict[str, Any]:
    """Merge overlay onto a copy of base. Dict values are merged one level deep
    (used for ``env`` and ``defines``); everything else is replaced."""
    out = copy.deepcopy(base)
    for key, val in overlay.items():
        if isinstance(val, dict) and isinstance(out.get(key), dict):
            merged = copy.deepcopy(out[key])
            merged.update(val)
            out[key] = merged
        else:
            out[key] = copy.deepcopy(val)
    return out


# ---------------------------------------------------------------------------
# Toolchains
# ---------------------------------------------------------------------------


def toolchain_layer(config: dict[str, Any], name: str) -> dict[str, Any]:
    """Return the settings bundle a toolchain name stands for.

    The bundle is merged like any other defaults layer, so a toolchain may carry
    ``env``, ``defines`` or plain variant fields. An empty bundle is legal and
    means "take CC/CXX/FC from the ambient environment".
    """
    registry = config.get("toolchains") or {}
    if not isinstance(registry, dict):
        die("'toolchains' must be a mapping of name -> settings")
    if name not in registry:
        available = ", ".join(sorted(registry)) or "(none defined)"
        die(f"unknown toolchain {name!r}. Available: {available}")
    entry = registry[name]
    if entry is None:
        return {}
    if not isinstance(entry, dict):
        die(f"toolchain {name!r} must be a mapping")
    return copy.deepcopy(entry)


def _expand_placeholders(text: str, values: dict[str, Any]) -> str:
    """Substitute ``{field}`` references in a path-valued setting.

    Unset fields collapse to the empty string, and the resulting empty path
    segments are folded away so that ``install/{toolchain}`` stays usable when
    no toolchain is selected.
    """
    out = text
    for key, val in values.items():
        out = out.replace("{" + key + "}", "" if val is None else str(val))
    while "//" in out:
        out = out.replace("//", "/")
    return out.rstrip("/") or "."


# ---------------------------------------------------------------------------
# Matrix expansion (GitHub-Actions-style semantics)
# ---------------------------------------------------------------------------


def _axis_values(values: Any) -> list[dict[str, Any]]:
    """Normalise one matrix axis into a list of partial-variant dicts.

    A scalar axis ``order: [4, 6]`` under key ``order`` yields
    ``[{order: 4}, {order: 6}]``.  An axis whose values are dicts (e.g. the
    GitHub-style ``equation: [{type: elastic, multisim: 8}]``) yields those
    dicts verbatim, with ``type`` renamed to ``equations``.
    """
    if not isinstance(values, list):
        die(f"matrix axis must be a list, got {type(values).__name__}")
    out: list[dict[str, Any]] = []
    for v in values:
        if isinstance(v, dict):
            out.append(_normalise_partial(v))
        else:
            out.append({"__scalar__": v})
    return out


def _normalise_partial(part: dict[str, Any]) -> dict[str, Any]:
    """Rename convenience keys to canonical field names."""
    part = dict(part)
    if "type" in part and "equations" not in part:
        part["equations"] = part.pop("type")
    if "equation" in part and "equations" not in part:
        part["equations"] = part.pop("equation")
    if "arch" in part and "host_arch" not in part:
        part["host_arch"] = part.pop("arch")
    return part


def expand_matrix(matrix: dict[str, Any]) -> list[dict[str, Any]]:
    """Cross-product of axes, then apply include/exclude (GH semantics)."""
    matrix = dict(matrix)
    include = matrix.pop("include", []) or []
    exclude = matrix.pop("exclude", []) or []

    axis_names = list(matrix.keys())
    axis_lists = [_axis_values(matrix[name]) for name in axis_names]

    combos: list[dict[str, Any]] = []
    for choice in itertools.product(*axis_lists) if axis_lists else [()]:
        entry: dict[str, Any] = {}
        for name, part in zip(axis_names, choice):
            if "__scalar__" in part:
                entry[name] = part["__scalar__"]
            else:
                # dict-valued axis (e.g. `equation`): merge its keys in.
                entry.update(part)
        combos.append(entry)

    # exclude: drop any combo that matches ALL keys of an exclude entry.
    def matches(entry: dict[str, Any], pattern: dict[str, Any]) -> bool:
        return all(entry.get(k) == v for k, v in _normalise_partial(pattern).items())

    combos = [e for e in combos if not any(matches(e, ex) for ex in exclude)]

    # include: append extra explicit variants.
    for inc in include:
        combos.append(_normalise_partial(inc))

    return combos


def preset_variants(preset: dict[str, Any]) -> list[dict[str, Any]]:
    if "matrix" in preset:
        return expand_matrix(preset["matrix"])
    if "variants" in preset:
        return [_normalise_partial(v) for v in preset["variants"]]
    die("preset must contain either 'matrix' or 'variants'")


# ---------------------------------------------------------------------------
# Variant model
# ---------------------------------------------------------------------------


@dataclass
class Variant:
    spec: dict[str, Any]  # canonical field -> value (equations, order, ...)
    meta: dict[str, Any]  # source_dir, build_root, install_prefix, generator, jobs
    env: dict[str, str]
    defines: dict[str, Any]
    testing: dict[str, Any] = field(default_factory=lambda: dict(TESTING_DEFAULTS))
    slug: str = ""  # build-dir name, unique across the selection

    @property
    def toolchain(self) -> str | None:
        return self.meta.get("toolchain")

    @property
    def runs_tests(self) -> bool:
        return bool(self.testing.get("enabled")) and not self.testing.get("build_only")

    @property
    def name_suffix(self) -> str:
        """Reproduces the NEW_BINARY_NAMING NAME_SUFFIX from SeisSol's CMake."""
        s = self.spec
        eq = str(s["equations"])
        # viscoelastic2 shares the binary name 'viscoelastic'
        base = "viscoelastic" if eq == "viscoelastic2" else eq
        suffix = base
        mech = int(s.get("mechanisms", 0) or 0)
        if mech > 0:
            suffix += f"-{mech}"
        suffix += f"-o{int(s['order'])}"
        suffix += "-f32" if str(s["precision"]) == "single" else "-f64"
        multisim = int(s.get("multisim", 1) or 1)
        if multisim > 1:
            suffix += f"-s{multisim}"
        custom = str(s.get("custom_suffix", "") or "")
        if custom:
            suffix += f"-{custom}"
        return suffix

    @property
    def binary_name(self) -> str:
        return f"seissol-{self.name_suffix}"

    def base_slug(self) -> str:
        """Human-readable build-dir slug (includes arch/build_type so that
        arch-only differences do not collide, unlike the binary name)."""
        parts = [self.name_suffix]
        if self.spec.get("device_backend") and self.spec["device_backend"] != "none":
            parts.append(str(self.spec["device_backend"]))
            if self.spec.get("device_arch"):
                parts.append(str(self.spec["device_arch"]))
        elif self.spec.get("host_arch"):
            parts.append(str(self.spec["host_arch"]))
        if self.toolchain:
            parts.append(str(self.toolchain))
        bt = str(self.meta.get("build_type", "")).strip()
        if bt and bt.lower() != "release":
            parts.append(bt)
        return "-".join(p for p in parts if p).replace("/", "_")

    def spec_hash(self) -> str:
        payload = json.dumps(
            {
                "spec": self.spec,
                "env": self.env,
                "defines": self.defines,
                "toolchain": self.toolchain,
                "generator": self.meta.get("generator"),
                "build_type": self.meta.get("build_type"),
            },
            sort_keys=True,
        )
        # Not a security primitive: this only disambiguates build-dir names.
        return hashlib.sha1(payload.encode(), usedforsecurity=False).hexdigest()[:6]

    @property
    def build_dir(self) -> Path:
        return Path(self.meta["build_root"]) / self.slug

    @property
    def install_prefix(self) -> Path:
        return Path(self.meta["install_prefix"])

    @property
    def binary_path(self) -> Path:
        return self.install_prefix / "bin" / self.binary_name

    def cmake_configure_cmd(self) -> list[str]:
        cmd = [
            "cmake",
            "-S",
            str(self.meta["source_dir"]),
            "-B",
            str(self.build_dir),
            "-G",
            str(self.meta["generator"]),
            f"-DCMAKE_INSTALL_PREFIX={self.install_prefix}",
        ]
        for field_name, cmake_var in FIELD_TO_CMAKE.items():
            if field_name in self.spec and self.spec[field_name] is not None:
                cmd.append(f"-D{cmake_var}={self.spec[field_name]}")
        for key, val in sorted(self.defines.items()):
            cmd.append(f"-D{key}={val}")
        return cmd

    def cmake_build_cmd(self, jobs: int) -> list[str]:
        return [
            "cmake",
            "--build",
            str(self.build_dir),
            "--target",
            "install",
            "-j",
            str(jobs),
        ]

    @property
    def junit_path(self) -> Path:
        return self.build_dir / "ctest-results.xml"

    def ctest_cmd(self, jobs: int) -> list[str]:
        cmd = [
            "ctest",
            "--test-dir",
            str(self.build_dir),
            "--output-on-failure",
            "--output-junit",
            str(self.junit_path.name),
            "-j",
            str(jobs),
        ]
        for label in self.testing.get("labels") or []:
            cmd += ["-L", str(label)]
        for label in self.testing.get("exclude_labels") or []:
            cmd += ["-LE", str(label)]
        timeout = self.testing.get("timeout")
        if timeout:
            cmd += ["--timeout", str(int(timeout))]
        cmd += [str(a) for a in (self.testing.get("ctest_args") or [])]
        return cmd

    def effective_env(self) -> dict[str, str]:
        env = os.environ.copy()
        env.update({k: str(v) for k, v in self.env.items()})
        return env


def _anchor(path_value: str, base: Path) -> str:
    """Resolve a path setting against the variant-set file's directory.

    Anchoring to the file rather than to the working directory means the same
    invocation works from anywhere in the tree. Absolute paths pass through, so
    a caller can point the install prefix at a fixed location.
    """
    candidate = Path(path_value)
    return str(candidate if candidate.is_absolute() else (base / candidate))


def build_variants(config: dict[str, Any], preset_name: str | None) -> list[Variant]:
    global_defaults = config.get("defaults", {}) or {}

    presets = config.get("presets")
    if presets:
        if preset_name is None:
            if "default" in presets:
                preset_name = "default"
            else:
                die(
                    "config has presets; choose one with --preset. Available: "
                    + ", ".join(sorted(presets))
                )
        if preset_name not in presets:
            die(
                f"unknown preset {preset_name!r}. Available: "
                + ", ".join(sorted(presets))
            )
        preset = presets[preset_name]
        preset_defaults = preset.get("defaults", {}) or {}
        raw = preset_variants(preset)
    else:
        # top-level implicit preset
        preset_defaults = {}
        raw = preset_variants(config)

    merged_defaults = _deep_merge(global_defaults, preset_defaults)

    known_keys = set(FIELD_TO_CMAKE) | set(META_FIELDS) | {"env", "defines", "testing"}

    variants: list[Variant] = []
    for raw_entry in raw:
        entry = _normalise_partial(raw_entry)

        # Cascade: global defaults -> preset defaults -> toolchain -> variant.
        # The toolchain sits below the variant so that a single `include:` leg
        # can still override individual settings the toolchain brings along.
        toolchain = entry.get("toolchain", merged_defaults.get("toolchain"))
        layer = merged_defaults
        if toolchain:
            layer = _deep_merge(layer, toolchain_layer(config, str(toolchain)))
        resolved = _deep_merge(layer, entry)
        resolved.setdefault("build_type", "Release")

        env = {k: str(v) for k, v in (resolved.pop("env", {}) or {}).items()}
        defines = resolved.pop("defines", {}) or {}
        testing = _deep_merge(TESTING_DEFAULTS, resolved.pop("testing", {}) or {})
        if testing["enabled"]:
            for key, val in TESTING_DEFINES.items():
                defines.setdefault(key, val)
            if testing["build_only"]:
                defines.setdefault("TESTING_BUILD_ONLY", "ON")

        meta = {
            "source_dir": resolved.get("source_dir", "."),
            "build_root": resolved.get("build_root", "build"),
            "install_prefix": resolved.get("install_prefix", "install"),
            "generator": resolved.get("generator", "Ninja"),
            "jobs": int(resolved.get("jobs", 0) or 0),
            "build_type": resolved.get("build_type", "Release"),
            "toolchain": str(toolchain) if toolchain else None,
            "preset": preset_name or "default",
            "config_dir": config.get("_config_dir", "."),
        }
        # The unexpanded forms are kept so that a variant can be re-resolved
        # against a different toolchain without losing its placeholders.
        meta["templates"] = {k: str(meta[k]) for k in TEMPLATED_META_FIELDS}

        # spec = every resolved key that maps to a CMake variable
        # (this correctly picks up host_arch/device_* set in defaults).
        spec = {
            k: v for k, v in resolved.items() if k in FIELD_TO_CMAKE and v is not None
        }

        unknown = set(resolved) - known_keys
        if unknown:
            warn(
                f"variant {raw_entry!r}: ignoring unknown key(s) "
                + ", ".join(sorted(unknown))
            )

        for req in REQUIRED_FIELDS:
            if spec.get(req) in (None, ""):
                die(
                    f"variant {raw_entry!r} is missing required field {req!r} "
                    "(set it on the variant or in defaults)"
                )
        variant = Variant(
            spec=spec, meta=meta, env=env, defines=defines, testing=testing
        )
        _expand_paths(variant)
        variants.append(variant)

    _assign_slugs(variants)
    _check_binary_collisions(variants)
    _domain_warnings(variants)
    return variants


def _expand_paths(variant: Variant) -> None:
    """Resolve {placeholder} references in the variant's path settings."""
    placeholders: dict[str, Any] = {
        "preset": variant.meta.get("preset"),
        "toolchain": variant.meta.get("toolchain"),
        "build_type": variant.meta.get("build_type"),
        **{k: variant.spec.get(k) for k in FIELD_TO_CMAKE},
    }
    base = Path(variant.meta["config_dir"])
    for key, template in variant.meta["templates"].items():
        variant.meta[key] = _anchor(_expand_placeholders(template, placeholders), base)
    variant.meta["source_dir"] = _anchor(str(variant.meta["source_dir"]), base)


def apply_overrides(variants: list[Variant], assignments: list[str]) -> None:
    """Apply `--set key=value` to every selected variant.

    Accepts a variant field (`device_arch=sm_80`), or a nested setting via a
    dotted prefix (`env.CC=clang-19`, `defines.ASAGI=OFF`). Meant for settings
    that are orthogonal to the variant set itself, such as the device backend
    a CI job happens to target.
    """
    for assignment in assignments:
        key, sep, value = assignment.partition("=")
        if not sep:
            die(f"--set expects key=value, got {assignment!r}")
        key = key.strip()
        for v in variants:
            if key.startswith("env."):
                v.env[key[4:]] = value
            elif key.startswith("defines."):
                v.defines[key[8:]] = value
            elif key in FIELD_TO_CMAKE:
                v.spec[key] = value
            elif key in META_FIELDS:
                v.meta[key] = value
            else:
                die(
                    f"--set: unknown key {key!r}. Use a variant field, "
                    "env.NAME or defines.NAME."
                )
    if assignments:
        _assign_slugs(variants)
        _check_binary_collisions(variants)


def apply_toolchains(
    config: dict[str, Any], variants: list[Variant], names: list[str]
) -> list[Variant]:
    """Re-resolve every variant against each named toolchain.

    This is the command-line override: the selection is crossed with the given
    toolchains, whatever the variant-set file said. CI uses it with a single
    name, because there the compiler is fixed by the container image.
    """
    out: list[Variant] = []
    for name in names:
        layer = toolchain_layer(config, name)
        layer_env = {k: str(v) for k, v in (layer.get("env") or {}).items()}
        for base in variants:
            v = copy.deepcopy(base)
            v.meta["toolchain"] = name
            v.env = _deep_merge(v.env, layer_env)
            v.defines = _deep_merge(v.defines, layer.get("defines") or {})
            for key, val in layer.items():
                if key in FIELD_TO_CMAKE and val is not None:
                    v.spec[key] = val
            _expand_paths(v)
            out.append(v)

    _assign_slugs(out)
    _check_binary_collisions(out)
    return out


def _assign_slugs(variants: list[Variant]) -> None:
    """Give every variant a unique, mostly-readable build-dir slug."""
    counts: dict[str, int] = {}
    for v in variants:
        counts[v.base_slug()] = counts.get(v.base_slug(), 0) + 1
    for v in variants:
        base = v.base_slug()
        v.slug = base if counts[base] == 1 else f"{base}-{v.spec_hash()}"


def _check_binary_collisions(variants: list[Variant]) -> None:
    """The binary name does NOT encode arch, so two variants differing only in
    arch would overwrite each other in a shared install prefix. Warn loudly."""
    seen: dict[tuple[str, str], Variant] = {}
    for v in variants:
        key = (str(v.install_prefix), v.binary_name)
        if key in seen:
            other = seen[key]
            hint = "Use distinct install_prefix values or a custom_suffix."
            if other.toolchain != v.toolchain:
                # A custom_suffix would work too, but it lands in the binary
                # name and thus breaks the lookup verify.py does via cases.json.
                hint = "Set install_prefix to something like 'install/{toolchain}'."
            warn(
                f"binary-name collision: {v.binary_name!r} produced by both "
                f"build dir {other.slug!r} and {v.slug!r} into the same prefix "
                f"({v.install_prefix}). The second install will overwrite the "
                f"first. {hint}"
            )
        else:
            seen[key] = v


def _domain_warnings(variants: list[Variant]) -> None:
    for v in variants:
        eq = str(v.spec.get("equations"))
        mech = int(v.spec.get("mechanisms", 0) or 0)
        if eq in ("viscoelastic", "viscoelastic2") and mech <= 0:
            warn(
                f"variant {v.slug!r}: {eq} usually needs NUMBER_OF_MECHANISMS > 0 "
                "(set 'mechanisms')."
            )


# ---------------------------------------------------------------------------
# Execution
# ---------------------------------------------------------------------------


@dataclass
class Result:
    variant: Variant
    status: str = "pending"  # ok | failed | skipped
    phase_failed: str | None = None  # configure | build | test | None
    returncode: int = 0
    configure_seconds: float = 0.0
    build_seconds: float = 0.0
    test_seconds: float = 0.0
    total_seconds: float = 0.0
    log_path: Path | None = None
    binary_found: bool = False
    test_status: str | None = None  # ok | failed | skipped | None (not requested)
    tests_total: int = 0
    tests_failed: int = 0
    failed_tests: list[str] = field(default_factory=list)


def _tail(path: Path, n: int = 40) -> str:
    try:
        lines = path.read_text(errors="replace").splitlines()
    except OSError:
        return ""
    return "\n".join(lines[-n:])


def _run_phase(cmd: list[str], env: dict[str, str], log_handle, stream: bool) -> int:
    log_handle.write("$ " + shlex.join(cmd) + "\n")
    log_handle.flush()
    if stream:
        proc = subprocess.Popen(
            cmd,
            env=env,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            text=True,
            bufsize=1,
        )
        assert proc.stdout is not None
        for line in proc.stdout:
            log_handle.write(line)
            sys.stdout.write(line)
        log_handle.flush()
        return proc.wait()
    proc = subprocess.run(
        cmd, env=env, stdout=log_handle, stderr=subprocess.STDOUT, text=True
    )
    return proc.returncode


def read_junit(path: Path) -> tuple[int, list[str]]:
    """Read CTest's JUnit output. Returns (total, names of failing tests)."""
    try:
        root = ET.parse(path).getroot()
    except (OSError, ET.ParseError):
        return 0, []
    cases = root.iter("testcase")
    total = 0
    failed: list[str] = []
    for case in cases:
        total += 1
        if case.find("failure") is not None or case.find("error") is not None:
            failed.append(case.get("name") or "<unnamed>")
    return total, failed


def test_one(variant: Variant, jobs: int, stream: bool) -> Result:
    """Run the variant's CTest suite against an existing build tree."""
    res = Result(variant=variant)
    if not variant.build_dir.exists():
        res.status, res.test_status = "skipped", "skipped"
        return res

    log_path = variant.build_dir / "test.log"
    res.log_path = log_path
    env = variant.effective_env()
    variant.junit_path.unlink(missing_ok=True)

    t0 = time.monotonic()
    with open(log_path, "w", encoding="utf-8") as fh:
        rc = _run_phase(variant.ctest_cmd(jobs), env, fh, stream)
    res.test_seconds = time.monotonic() - t0
    res.total_seconds = res.test_seconds
    res.tests_total, res.failed_tests = read_junit(variant.junit_path)
    res.tests_failed = len(res.failed_tests)
    res.binary_found = variant.binary_path.exists()

    if rc != 0:
        res.status, res.phase_failed, res.returncode = "failed", "test", rc
        res.test_status = "failed"
    else:
        res.status, res.test_status = "ok", "ok"
    return res


def build_one(
    variant: Variant, jobs: int, clean: bool, stream: bool, test_jobs: int = 0
) -> Result:
    res = Result(variant=variant)
    build_dir = variant.build_dir
    if clean and build_dir.exists():
        shutil.rmtree(build_dir)
    build_dir.mkdir(parents=True, exist_ok=True)
    variant.install_prefix.mkdir(parents=True, exist_ok=True)

    log_path = build_dir / "build.log"
    res.log_path = log_path
    env = variant.effective_env()
    start = time.monotonic()

    with open(log_path, "w", encoding="utf-8") as fh:
        # configure
        t0 = time.monotonic()
        rc = _run_phase(variant.cmake_configure_cmd(), env, fh, stream)
        res.configure_seconds = time.monotonic() - t0
        if rc != 0:
            res.status, res.phase_failed, res.returncode = "failed", "configure", rc
            res.total_seconds = time.monotonic() - start
            return res
        # build + install
        t0 = time.monotonic()
        rc = _run_phase(variant.cmake_build_cmd(jobs), env, fh, stream)
        res.build_seconds = time.monotonic() - t0
        if rc != 0:
            res.status, res.phase_failed, res.returncode = "failed", "build", rc
            res.total_seconds = time.monotonic() - start
            return res

        res.binary_found = variant.binary_path.exists()

        # unit + smoke tests, run inside the build tree: the test binaries are
        # not part of the install set.
        if variant.runs_tests:
            variant.junit_path.unlink(missing_ok=True)
            t0 = time.monotonic()
            rc = _run_phase(variant.ctest_cmd(test_jobs or jobs), env, fh, stream)
            res.test_seconds = time.monotonic() - t0
            res.tests_total, res.failed_tests = read_junit(variant.junit_path)
            res.tests_failed = len(res.failed_tests)
            if rc != 0:
                res.status, res.phase_failed, res.returncode = "failed", "test", rc
                res.test_status = "failed"
                res.total_seconds = time.monotonic() - start
                return res
            res.test_status = "ok"

    res.total_seconds = time.monotonic() - start
    res.status = "ok"
    return res


def compute_jobs(meta_jobs: int, cli_jobs: int, max_parallel: int) -> int:
    """Per-build -j so that max_parallel builds fit the core budget."""
    if cli_jobs > 0:
        return cli_jobs
    if meta_jobs > 0:
        return meta_jobs
    cores = os.cpu_count() or 1
    return max(1, cores // max(1, max_parallel))


def run_builds(
    variants: list[Variant],
    max_parallel: int,
    cli_jobs: int,
    clean: bool,
    fail_fast: bool,
    verbose: bool,
    test_jobs: int = 0,
    test_only: bool = False,
) -> list[Result]:
    results: dict[str, Result] = {}
    stop = threading.Event()
    total = len(variants)
    done = 0
    done_lock = threading.Lock()
    stream = verbose and max_parallel == 1

    def worker(idx: int, variant: Variant) -> Result:
        nonlocal done
        if stop.is_set():
            r = Result(variant=variant, status="skipped")
            return r
        jobs = compute_jobs(variant.meta["jobs"], cli_jobs, max_parallel)
        tjobs = compute_jobs(variant.testing.get("jobs", 0), test_jobs, max_parallel)
        if test_only:
            log(dim(f"[{idx}/{total}] testing {variant.slug}  (-j{tjobs})"))
            r = test_one(variant, tjobs, stream)
        else:
            log(dim(f"[{idx}/{total}] building {variant.slug}  (-j{jobs})"))
            r = build_one(variant, jobs, clean, stream, tjobs)
        with done_lock:
            done += 1
            marker = green("ok") if r.status == "ok" else red(r.status)
            extra = ""
            if r.status == "ok" and not r.binary_found:
                extra = yellow(f"  [binary {variant.binary_name} not found!]")
            if r.test_status is not None and r.tests_total:
                passed = r.tests_total - r.tests_failed
                tally = f"  [{passed}/{r.tests_total} tests]"
                extra += green(tally) if r.tests_failed == 0 else red(tally)
            log(
                f"[{done}/{total}] {marker:>16} {variant.slug} "
                f"{dim(f'({r.total_seconds:.0f}s)')}{extra}"
            )
            if r.status == "failed":
                log(
                    red(
                        f"    failed in {r.phase_failed} (rc={r.returncode}); "
                        f"log: {r.log_path}"
                    )
                )
                for name in r.failed_tests[:20]:
                    log(red(f"      test failed: {name}"))
                if len(r.failed_tests) > 20:
                    log(red(f"      ... and {len(r.failed_tests) - 20} more"))
                tail = _tail(r.log_path) if r.log_path else ""
                if tail:
                    log(dim("    --- last log lines ---"))
                    for line in tail.splitlines():
                        log(dim("    " + line))
                if fail_fast:
                    stop.set()
        return r

    with concurrent.futures.ThreadPoolExecutor(max_workers=max_parallel) as ex:
        futures = {ex.submit(worker, i + 1, v): v for i, v in enumerate(variants)}
        for fut in concurrent.futures.as_completed(futures):
            r = fut.result()
            results[r.variant.slug] = r

    # keep input order in the returned list
    return [results[v.slug] for v in variants]


# ---------------------------------------------------------------------------
# Reporting
# ---------------------------------------------------------------------------


def result_to_dict(r: Result) -> dict[str, Any]:
    v = r.variant
    return {
        "id": v.slug,
        "spec": v.spec,
        "toolchain": v.toolchain,
        "env": v.env,
        "defines": v.defines,
        "cmake_configure": shlex.join(v.cmake_configure_cmd()),
        "build_dir": str(v.build_dir),
        "install_prefix": str(v.install_prefix),
        "binary_name": v.binary_name,
        "binary_path": str(v.binary_path),
        "status": r.status,
        "phase_failed": r.phase_failed,
        "returncode": r.returncode,
        "configure_seconds": round(r.configure_seconds, 2),
        "build_seconds": round(r.build_seconds, 2),
        "test_seconds": round(r.test_seconds, 2),
        "total_seconds": round(r.total_seconds, 2),
        "log": str(r.log_path) if r.log_path else None,
        "binary_found": r.binary_found,
        "testing": {
            "status": r.test_status,
            "total": r.tests_total,
            "failed": r.tests_failed,
            "failed_tests": r.failed_tests,
            "junit": str(v.junit_path) if r.test_status else None,
        },
    }


def write_report(
    path: Path, results: list[Result], preset: str | None, wall: float
) -> None:
    n_ok = sum(r.status == "ok" for r in results)
    n_failed = sum(r.status == "failed" for r in results)
    n_skipped = sum(r.status == "skipped" for r in results)
    prefix = results[0].variant.install_prefix if results else Path("install")
    doc = {
        "schema": REPORT_SCHEMA,
        "generated": datetime.now(timezone.utc).isoformat(),
        "preset": preset,
        "install_prefix": str(prefix),
        "total_seconds": round(wall, 2),
        "summary": {
            "total": len(results),
            "ok": n_ok,
            "failed": n_failed,
            "skipped": n_skipped,
        },
        "variants": [result_to_dict(r) for r in results],
    }
    path.write_text(json.dumps(doc, indent=2), encoding="utf-8")
    log(f"wrote JSON report: {path}")


def write_step_summary(
    path: Path, results: list[Result], preset: str | None, wall: float
) -> None:
    n_ok = sum(r.status == "ok" for r in results)
    n_failed = sum(r.status == "failed" for r in results)
    toolchains = sorted({r.variant.toolchain for r in results if r.variant.toolchain})
    heading = f"## SeisSol build — {preset or 'default'}"
    if toolchains:
        heading += " (" + ", ".join(toolchains) + ")"
    lines = [
        heading,
        "",
        f"**{n_ok}/{len(results)} ok**, {n_failed} failed " f"· {wall:.0f}s wall",
        "",
        "| Variant | Status | Binary | Tests | Time |",
        "| --- | --- | --- | --- | --- |",
    ]
    for r in results:
        icon = {"ok": "✅", "failed": "❌", "skipped": "⏭️"}.get(r.status, "❔")
        binary = r.variant.binary_name
        if r.status == "ok" and not r.binary_found:
            binary += " ⚠️ missing"
        note = f" ({r.phase_failed})" if r.status == "failed" else ""
        if r.test_status is None:
            tests = "—"
        elif r.test_status == "skipped":
            tests = "skipped"
        else:
            tests = f"{r.tests_total - r.tests_failed}/{r.tests_total}"
        lines.append(
            f"| `{r.variant.slug}` | {icon} {r.status}{note} | "
            f"`{binary}` | {tests} | {r.total_seconds:.0f}s |"
        )

    failing = [(r.variant.slug, r.failed_tests) for r in results if r.failed_tests]
    if failing:
        lines += ["", "<details><summary>Failing tests</summary>", ""]
        for slug, names in failing:
            lines.append(f"**{slug}**")
            lines += [f"- `{n}`" for n in names[:50]]
            if len(names) > 50:
                lines.append(f"- ... and {len(names) - 50} more")
            lines.append("")
        lines.append("</details>")

    # append (GitHub step summaries are appended across steps)
    with open(path, "a", encoding="utf-8") as fh:
        fh.write("\n".join(lines) + "\n")
    log(f"wrote step summary: {path}")


def print_summary(results: list[Result]) -> None:
    n_ok = sum(r.status == "ok" for r in results)
    n_failed = sum(r.status == "failed" for r in results)
    n_skipped = sum(r.status == "skipped" for r in results)
    log("")
    log("summary:")
    for r in results:
        icon = {
            "ok": green("ok"),
            "failed": red("failed"),
            "skipped": yellow("skipped"),
        }.get(r.status, r.status)
        tests = ""
        if r.test_status is not None and r.tests_total:
            tests = dim(f"  {r.tests_total - r.tests_failed}/{r.tests_total} tests")
        log(
            f"  {icon:>16}  {r.variant.slug}  "
            f"{dim(f'{r.total_seconds:.0f}s')}{tests}"
        )
    tally = f"{n_ok} ok, {n_failed} failed"
    if n_skipped:
        tally += f", {n_skipped} skipped"
    tests_total = sum(r.tests_total for r in results)
    if tests_total:
        tests_failed = sum(r.tests_failed for r in results)
        tally += f" · {tests_total - tests_failed}/{tests_total} tests passed"
    log("")
    log((green if n_failed == 0 else red)(tally))


# ---------------------------------------------------------------------------
# Listing / dry-run
# ---------------------------------------------------------------------------


def list_variants(variants: list[Variant]) -> None:
    log(f"{len(variants)} variant(s):\n")
    for v in variants:
        log(green(v.slug))
        log(f"  binary : {v.binary_name}")
        log(f"  build  : {v.build_dir}")
        log(f"  install: {v.install_prefix}")
        spec = "  ".join(f"{k}={v.spec[k]}" for k in sorted(v.spec))
        log(f"  spec   : {spec}")
        if v.toolchain:
            log(f"  chain  : {v.toolchain}")
        if v.env:
            log(f"  env    : {' '.join(f'{k}={val}' for k, val in v.env.items())}")
        if v.testing.get("enabled"):
            mode = "build-only" if v.testing.get("build_only") else "build + ctest"
            log(f"  tests  : {mode}")
        log("")


def list_toolchains(config: dict[str, Any]) -> None:
    registry = config.get("toolchains") or {}
    if not registry:
        log("no toolchains defined")
        return
    log(f"{len(registry)} toolchain(s):\n")
    for name in sorted(registry):
        layer = toolchain_layer(config, name)
        log(green(name))
        env = layer.get("env") or {}
        log(
            "  env    : "
            + (" ".join(f"{k}={v}" for k, v in env.items()) or "(ambient)")
        )
        defines = layer.get("defines") or {}
        if defines:
            log("  defines: " + " ".join(f"{k}={v}" for k, v in defines.items()))
        log("")


def dry_run(
    variants: list[Variant], max_parallel: int, cli_jobs: int, test_jobs: int
) -> None:
    log(f"dry-run: {len(variants)} variant(s), max_parallel={max_parallel}\n")
    for v in variants:
        jobs = compute_jobs(v.meta["jobs"], cli_jobs, max_parallel)
        tjobs = compute_jobs(v.testing.get("jobs", 0), test_jobs, max_parallel)
        log(green(f"# {v.slug}  ->  {v.binary_name}"))
        if v.env:
            log(dim("  env: " + " ".join(f"{k}={val}" for k, val in v.env.items())))
        log("  " + shlex.join(v.cmake_configure_cmd()))
        log("  " + shlex.join(v.cmake_build_cmd(jobs)))
        if v.runs_tests:
            log("  " + shlex.join(v.ctest_cmd(tjobs)))
        log("")


# ---------------------------------------------------------------------------
# --rerun-failed
# ---------------------------------------------------------------------------


def filter_failed(variants: list[Variant], report_path: Path) -> list[Variant]:
    doc = json.loads(report_path.read_text(encoding="utf-8"))
    failed_ids = {v["id"] for v in doc.get("variants", []) if v.get("status") != "ok"}
    kept = [v for v in variants if v.slug in failed_ids]
    if not kept:
        log(green("nothing to rerun: all variants in the report succeeded."))
    else:
        log(f"rerun-failed: {len(kept)} of {len(variants)} variant(s)")
    return kept


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Build several SeisSol variants from a variant-set file.",
    )
    p.add_argument(
        "-c",
        "--config",
        type=Path,
        default=Path(__file__).resolve().parent / "compile-variants.yaml",
        help="variant-set file (.yaml/.yml/.json)",
    )
    p.add_argument(
        "-p",
        "--preset",
        default=None,
        help="which preset to build (if the config defines presets)",
    )
    p.add_argument(
        "-t",
        "--toolchain",
        action="append",
        default=None,
        metavar="NAME",
        help="build with this toolchain from the config's 'toolchains' block, "
        "overriding whatever the preset selects. Repeatable: the selection is "
        "crossed with every name given.",
    )
    p.add_argument(
        "--list-toolchains",
        action="store_true",
        help="list the toolchains the config defines and exit",
    )
    p.add_argument(
        "-j",
        "--jobs",
        type=int,
        default=0,
        help="compiler jobs per build (0 = auto from core budget)",
    )
    p.add_argument(
        "-P",
        "--max-parallel",
        type=int,
        default=1,
        help="how many variant builds to run concurrently (default 1)",
    )
    p.add_argument(
        "--test",
        dest="test",
        action="store_true",
        default=None,
        help="run each variant's CTest suite after building it "
        "(overrides 'testing.enabled' in the config)",
    )
    p.add_argument(
        "--no-test",
        dest="test",
        action="store_false",
        help="never run CTest, whatever the config says",
    )
    p.add_argument(
        "--test-only",
        action="store_true",
        help="skip configure/build and only run CTest in the existing build trees",
    )
    p.add_argument(
        "--test-jobs",
        type=int,
        default=0,
        help="parallel tests per variant (0 = auto from core budget)",
    )
    p.add_argument(
        "--ctest-arg",
        action="append",
        default=[],
        dest="ctest_args",
        metavar="ARG",
        help="extra argument passed to ctest verbatim (repeatable)",
    )
    p.add_argument(
        "--set",
        action="append",
        default=[],
        dest="overrides",
        metavar="KEY=VALUE",
        help="override a variant field (device_arch=sm_80), or a nested "
        "setting (env.CC=clang-19, defines.ASAGI=OFF). Repeatable.",
    )
    p.add_argument(
        "--install-prefix",
        type=Path,
        default=None,
        help="install every selected variant here, overriding the config",
    )
    p.add_argument(
        "--build-root",
        type=Path,
        default=None,
        help="place the per-variant build trees here, overriding the config",
    )
    p.add_argument(
        "--clean", action="store_true", help="remove each build dir before configuring"
    )
    p.add_argument(
        "--fail-fast",
        action="store_true",
        help="stop scheduling new builds after the first failure",
    )
    p.add_argument(
        "--dry-run", action="store_true", help="print the exact cmake commands and exit"
    )
    p.add_argument(
        "--list", action="store_true", help="list the selected variants and exit"
    )
    p.add_argument(
        "--report",
        type=Path,
        default=Path("build-report.json"),
        help="write a machine-readable JSON report here",
    )
    p.add_argument(
        "--step-summary",
        type=Path,
        default=None,
        help="write a Markdown summary here " "(falls back to $GITHUB_STEP_SUMMARY)",
    )
    p.add_argument(
        "--rerun-failed",
        type=Path,
        default=None,
        metavar="REPORT",
        help="only rebuild variants that were not 'ok' in a prior report",
    )
    p.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        help="stream build output live (only with --max-parallel 1)",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.config.exists():
        die(f"config not found: {args.config}")
    if args.max_parallel < 1:
        die("--max-parallel must be >= 1")

    config = load_config(args.config)

    if args.list_toolchains:
        list_toolchains(config)
        return 0

    variants = build_variants(config, args.preset)
    if args.toolchain:
        variants = apply_toolchains(config, variants, args.toolchain)
    apply_overrides(variants, args.overrides)
    if not variants:
        die("no variants selected")

    for v in variants:
        if args.install_prefix is not None:
            v.meta["install_prefix"] = str(args.install_prefix)
        if args.build_root is not None:
            v.meta["build_root"] = str(args.build_root)
    if args.install_prefix is not None:
        _check_binary_collisions(variants)

    # The CLI wins over the config for the test phase; --test-only implies it.
    for v in variants:
        if args.test_only:
            v.testing["enabled"] = True
            v.testing["build_only"] = False
        elif args.test is not None:
            v.testing["enabled"] = args.test
        if args.ctest_args:
            v.testing["ctest_args"] = list(v.testing.get("ctest_args") or []) + list(
                args.ctest_args
            )

    if args.rerun_failed is not None:
        if not args.rerun_failed.exists():
            die(f"report not found: {args.rerun_failed}")
        variants = filter_failed(variants, args.rerun_failed)
        if not variants:
            return 0

    if args.list:
        list_variants(variants)
        return 0
    if args.dry_run:
        dry_run(variants, args.max_parallel, args.jobs, args.test_jobs)
        return 0

    if args.verbose and args.max_parallel > 1:
        warn("--verbose is ignored with --max-parallel > 1 (logs go to files)")

    verb = "testing" if args.test_only else "building"
    log(f"{verb} {len(variants)} variant(s), " f"max_parallel={args.max_parallel}")
    wall_start = time.monotonic()
    results = run_builds(
        variants,
        max_parallel=args.max_parallel,
        cli_jobs=args.jobs,
        clean=args.clean,
        fail_fast=args.fail_fast,
        verbose=args.verbose,
        test_jobs=args.test_jobs,
        test_only=args.test_only,
    )
    wall = time.monotonic() - wall_start

    print_summary(results)
    write_report(args.report, results, args.preset, wall)

    step_summary = args.step_summary
    if step_summary is None and os.environ.get("GITHUB_STEP_SUMMARY"):
        step_summary = Path(os.environ["GITHUB_STEP_SUMMARY"])
    if step_summary is not None:
        write_step_summary(step_summary, results, args.preset, wall)

    return 1 if any(r.status == "failed" for r in results) else 0


if __name__ == "__main__":
    raise SystemExit(main())
