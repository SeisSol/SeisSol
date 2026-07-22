#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""Run SeisSol against precomputed reference solutions and verify the outputs.

Reads ``cases.json`` (in --precomputed-dir) and per-case ``tpv-data.json`` files,
dispatches the appropriate SeisSol binary from --binary-dir per case, then calls
the existing ``compare-{mesh,receivers,energies}.py`` scripts next to this file.

Binary naming (pre-PR #1421): one binary per (equation, precision, fused-count)
combination, named ``seissol-<config>`` where ``<config>`` is read directly
from cases.json (e.g. ``elastic-o6-f64``, ``viscoelastic-3-o6-f64-s8``). Once
PR #1421 lands, all configs collapse to a single ``seissol`` binary; this script
falls back to that name automatically.

Examples
--------
Run one case (both precisions are picked up when no precision is given)::

    verify.py --precomputed-dir /path/to/precomputed-seissol \
              --binary-dir /opt/seissol/bin --case tpv5

Run everything in cases.json::

    verify.py --precomputed-dir . --binary-dir /opt/seissol/bin --all

Re-record the reference for one case against the currently-built SeisSol::

    verify.py --precomputed-dir . --binary-dir /opt/seissol/bin \
              --case tpv5/double --record --seissol-commit $(git -C ../SeisSol rev-parse HEAD)

Re-record everything (used by the recompute-all GH workflow)::

    verify.py --precomputed-dir . --binary-dir /opt/seissol/bin \
              --all --record --seissol-commit $SEISSOL_COMMIT --note "Q1 2026 refresh"
"""

from __future__ import annotations

import argparse
import dataclasses
import datetime as _dt
import json
import math
import os
import re
import shlex
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterable, Optional, Sequence

SCRIPT_DIR = Path(__file__).resolve().parent
COMPARE_MESH = SCRIPT_DIR / "compare-mesh.py"
COMPARE_RECEIVERS = SCRIPT_DIR / "compare-receivers.py"
COMPARE_ENERGIES = SCRIPT_DIR / "compare-energies.py"

DEFAULT_EPSILON = 0.05
OUTPUT_SUBDIR = "output"  # relative to the case directory (fresh run output)
# Placeholder threshold written into an auto-generated <prefix>-data.json. A
# recording has zero self-error, so a meaningful epsilon can't be derived here --
# these defaults exist to be reviewed and tuned by hand (or via --suggest-epsilons).
SKELETON_EPSILON = 0.01


# ---------------------------------------------------------------------------
# ANSI color
# ---------------------------------------------------------------------------


class _Ansi:
    RESET = "\033[0m"
    BOLD = "\033[1m"
    DIM = "\033[2m"
    RED = "\033[31m"
    GREEN = "\033[32m"
    YELLOW = "\033[33m"
    BLUE = "\033[34m"
    MAGENTA = "\033[35m"
    CYAN = "\033[36m"


_color_enabled = False


def setup_color(mode: str) -> None:
    """Decide whether to emit ANSI escapes. Respects ``NO_COLOR`` and ``TERM=dumb``."""
    global _color_enabled
    if mode == "always":
        _color_enabled = True
    elif mode == "never":
        _color_enabled = False
    else:  # "auto"
        _color_enabled = (
            sys.stdout.isatty()
            and os.environ.get("NO_COLOR") is None
            and os.environ.get("TERM", "") != "dumb"
        )


def c(text: str, *codes: str) -> str:
    """Wrap ``text`` in the given ANSI codes if color is enabled."""
    if not _color_enabled or not codes:
        return text
    return "".join(codes) + text + _Ansi.RESET


# Convenience accessors for common roles. They take the *already-padded* text
# so the visible width is preserved in tables.
def _ok(s: str) -> str:    return c(s, _Ansi.GREEN, _Ansi.BOLD)
def _bad(s: str) -> str:   return c(s, _Ansi.RED, _Ansi.BOLD)
def _warn(s: str) -> str:  return c(s, _Ansi.YELLOW)
def _info(s: str) -> str:  return c(s, _Ansi.CYAN)
def _dim(s: str) -> str:   return c(s, _Ansi.DIM)


# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclasses.dataclass
class CaseSpec:
    """One ``cases.json`` entry, keyed by an opaque ``key``.

    A case declares its own output via ``output`` (directory + prefix) and/or a
    ``reference`` (another case's key) to compare against; at least one must be
    present. With only ``output`` the case owns and self-checks its reference;
    with only ``reference`` it borrows another case's output and inherits that
    case's prefix; with both, it compares against ``reference`` but keeps its own
    output directory. Fresh SeisSol output always goes to ``workdir``; the
    ``output`` directory holds the committed reference (and is the record target).
    """

    key: str            # opaque identifier from cases.json (e.g. "tpv5/double")
    config: str         # binary selector, e.g. "elastic-o6-f64"
    params_path: Path   # absolute path to the parameter file
    case_dir: Path      # parent of params_path (SeisSol run cwd)
    prefix: str         # output file prefix, e.g. "tpv"
    output_dir: Optional[Path] = None  # own committed output dir (record target)
    reference: str = ""  # another case's key; "" => compare against own output
    # Resolved compare target, filled in by load_cases (own output_dir, or the
    # referenced case's output_dir):
    reference_dir: Optional[Path] = None
    # Case-level conditional rules (each: {"when": {...}, "reason": ...}); apply
    # to the whole case / all its categories. Category-level rules live in the
    # per-category <prefix>-data.json instead.
    xfail: list = dataclasses.field(default_factory=list)
    skip: list = dataclasses.field(default_factory=list)

    @property
    def is_derived(self) -> bool:
        """True when the compare target is another case's output."""
        return bool(self.reference) and self.reference != self.key

    @property
    def workdir(self) -> Path:
        """Where SeisSol writes fresh output (the parameter-file-driven location)."""
        return self.case_dir / OUTPUT_SUBDIR


@dataclasses.dataclass
class CompareResult:
    category: str
    passed: bool
    epsilon: float                       # the category threshold that was applied
    detail: str = ""
    achieved: Optional[float] = None     # worst error actually measured
    quantities: dict[str, float] = dataclasses.field(default_factory=dict)
    has_summary: bool = False            # a report JSON was produced
    # per-quantity threshold violations: (quantity, achieved, threshold)
    failures: list[tuple[str, float, float]] = dataclasses.field(default_factory=list)
    # conditional expected-failure state (see the xfail/skip machinery)
    xfail_outcome: str = ""   # "" | "xfail" (failed as expected) | "xpass" (passed unexpectedly)
    xfail_reason: str = ""
    skipped: bool = False     # category-level conditional skip (not compared)
    skip_reason: str = ""


def _compare_is_failure(cr: "CompareResult", xfail_strict: bool) -> bool:
    """Whether a comparison counts as a real failure for the run's verdict.

    A conditionally-skipped category never fails; an expected failure (xfail)
    never fails; an unexpected pass (xpass) fails only under ``--xfail-strict``.
    """
    if cr.skipped:
        return False
    if cr.xfail_outcome == "xfail":
        return False
    if cr.xfail_outcome == "xpass":
        return xfail_strict
    return not cr.passed


@dataclasses.dataclass
class CaseResult:
    case: CaseSpec
    ran: bool
    run_returncode: int
    duration_s: float
    compares: list[CompareResult]
    recorded: bool = False
    skipped: bool = False   # not evaluated (no binary / derived record / conditional skip)
    note: str = ""          # why it was skipped, shown in the report
    xfail_strict: bool = False  # if set, an unexpected pass (xpass) fails the run

    @property
    def passed(self) -> bool:
        if self.skipped:
            return True
        if self.ran and self.run_returncode != 0:
            return False
        if self.recorded:
            return True
        if not self.compares:
            return False
        return not any(_compare_is_failure(c, self.xfail_strict) for c in self.compares)


# ---------------------------------------------------------------------------
# Loading: cases.json and tpv-data.json
# ---------------------------------------------------------------------------


def load_cases(precomputed_dir: Path) -> dict[str, CaseSpec]:
    cases_file = precomputed_dir / "cases.json"
    if not cases_file.is_file():
        raise FileNotFoundError(f"cases.json not found in {precomputed_dir}")
    raw = json.loads(cases_file.read_text())

    cases: dict[str, CaseSpec] = {}
    for key, entry in raw.items():
        if "output" not in entry and "reference" not in entry:
            raise ValueError(
                f"Case '{key}' has neither 'output' nor 'reference'; one is required"
            )
        params_path = (precomputed_dir / entry["parameters"]).resolve()
        prefix, output_dir = "", None
        out = entry.get("output")
        if isinstance(out, str):        # lenient: a bare string is the prefix
            prefix = out
        elif isinstance(out, dict):
            prefix = out.get("prefix", "")
            if out.get("directory"):
                output_dir = (precomputed_dir / out["directory"]).resolve()
        cases[key] = CaseSpec(
            key=key,
            config=entry["config"],
            params_path=params_path,
            case_dir=params_path.parent,
            prefix=prefix,
            output_dir=output_dir,
            reference=entry.get("reference", ""),
            xfail=entry.get("xfail", []),
            skip=entry.get("skip", []),
        )

    # Second pass: resolve each case's compare target and inherit the prefix from
    # the referenced case when this case declared no output of its own.
    for case in cases.values():
        if case.reference:
            if case.reference not in cases:
                raise ValueError(
                    f"Case '{case.key}' references unknown case '{case.reference}'"
                )
            ref = cases[case.reference]
            if ref.output_dir is None:
                raise ValueError(
                    f"Case '{case.key}' references '{case.reference}', which has no "
                    f"output directory to compare against"
                )
            case.reference_dir = ref.output_dir
            if not case.prefix:
                case.prefix = ref.prefix
        else:
            if case.output_dir is None:
                raise ValueError(
                    f"Case '{case.key}' has no 'reference' and no output directory"
                )
            case.reference_dir = case.output_dir
        if not case.prefix:
            raise ValueError(
                f"Case '{case.key}': no output prefix (set output.prefix, or "
                f"reference a case that has one)"
            )
    return cases


def load_tpv_data(case: CaseSpec) -> dict:
    """Per-case validation config; missing file => empty (all categories off)."""
    f = case.reference_dir / f"{case.prefix}-data.json"
    if not f.is_file():
        return {}
    return json.loads(f.read_text())


def _enabled(tpv_data: dict, category: str) -> bool:
    return bool(tpv_data.get(category, {}).get("enabled", False))


def _epsilon(
    tpv_data: dict, category: str, override: Optional[float], variant: str
) -> float:
    """Resolve the threshold for one category.

    ``epsilon`` in tpv-data.json may be either a plain number (applies
    everywhere) or an object keyed by variant label with a ``default`` fallback::

        "receiver": {"enabled": true, "epsilon": {"default": 0.02, "single": 0.15}}

    The dict is resolved by the *verified entry's own* variant label, else
    ``default``. ``default`` therefore carries the reference-vs-reference
    tolerance, and each derived variant that needs a looser (or tighter) bound
    adds its own key. Keys are opaque -- they are whatever the repo names its
    precomputed subtrees (precisions, orders, ...), never interpreted here.

    Note: keying by the *reference* label would be ambiguous, since a reference
    compared against itself and a derived variant compared against it both point
    at the same tree; the discriminator is the entry being verified.
    """
    if override is not None:
        return override
    raw = tpv_data.get(category, {}).get("epsilon", DEFAULT_EPSILON)
    if isinstance(raw, dict):
        return float(raw.get(variant, raw.get("default", DEFAULT_EPSILON)))
    return float(raw)


def _has_quantity_thresholds(tpv_data: dict, category: str) -> bool:
    """True if this category defines any per-quantity override."""
    return bool(tpv_data.get(category, {}).get("quantities"))


def _quantity_epsilon(
    tpv_data: dict, category: str, quantity: str, variant: str,
    override: Optional[float],
) -> float:
    """Threshold for a single quantity, falling back to the category threshold.

    ``quantities`` in tpv-data.json maps an *exact* quantity name (as it appears
    in the report JSON, e.g. ``receiver:v1``) to either a scalar or a per-variant
    dict, so the per-quantity and per-variant axes compose::

        "receiver": {
          "epsilon":    {"default": 0.02, "single": 0.15},
          "quantities": {"receiver:SRs": {"default": 0.05, "single": 0.30},
                         "receiver:v1":  0.01}
        }

    A ``--epsilon`` override on the command line still wins over everything.
    """
    if override is not None:
        return override
    q = tpv_data.get(category, {}).get("quantities", {}).get(quantity)
    if q is None:
        return _epsilon(tpv_data, category, None, variant)
    if isinstance(q, dict):
        return float(q.get(variant, q.get("default", DEFAULT_EPSILON)))
    return float(q)


def _evaluate_quantities(
    result: "CompareResult", tpv_data: dict, variant: str,
    override: Optional[float],
) -> None:
    """Re-decide one comparison's verdict from its per-quantity achieved errors.

    Only called when the category defines a ``quantities`` map and a report JSON
    was produced. Reduces exactly to the category-epsilon verdict when no
    quantity has an override, and generalises to per-quantity thresholds
    otherwise. Warns about configured quantity names that never appeared in the
    output (exact-match; likely a typo or a rename).
    """
    configured = set(tpv_data.get(result.category, {}).get("quantities", {}))
    seen = set(result.quantities)
    for missing in sorted(configured - seen):
        print(f"{_warn('[warn]')} {result.category}: per-quantity threshold for "
              f"'{missing}' set but that quantity was not found in the output",
              file=sys.stderr)

    failures: list[tuple[str, float, float]] = []
    for q, achieved in result.quantities.items():
        thr = _quantity_epsilon(tpv_data, result.category, q, variant, override)
        # ``not <=`` so NaN and +inf both count as violations.
        if not (achieved <= thr):
            failures.append((q, float(achieved), thr))
    result.failures = failures
    result.passed = not failures


# ---------------------------------------------------------------------------
# Binary resolution
# ---------------------------------------------------------------------------


def resolve_binary(binary_dir: Path, config: str) -> Path:
    """Locate a SeisSol executable for the given config. Returns the absolute
    path to it.

    Pre-PR #1421 layout: ``binary_dir/seissol-<config>``.

    Post-PR #1421 layout: a single ``binary_dir/seissol`` (the config string
    becomes informational); this is tried last as a fallback for any variant.
    Each prefix is searched in both ``binary_dir`` and ``binary_dir/bin``.
    """
    candidates: list[Path] = []
    candidates.append(binary_dir / f"seissol-{config}")
    candidates.append(binary_dir / "bin" / f"seissol-{config}")
    candidates.append(binary_dir / "seissol")
    candidates.append(binary_dir / "bin" / "seissol")

    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return cand.absolute()
    raise FileNotFoundError(
        f"No SeisSol binary found for config '{config}'. "
        f"Tried:\n  " + "\n  ".join(str(c) for c in candidates)
    )


# ---------------------------------------------------------------------------
# Running SeisSol
# ---------------------------------------------------------------------------


def run_seissol(
    case: CaseSpec,
    binary: Path,
    mpi_cmd: list[str],
    ranks: int,
    threads: int,
    clean_output: bool = True,
) -> tuple[int, float]:
    """Run SeisSol for one case, returning ``(returncode, duration_seconds)``.

    The simulation is executed with ``cwd = case.case_dir``, which is also where
    the SeisSol parameter file expects to find relative paths (mesh, yaml, ...).
    Output is written to ``case.case_dir/output/`` per parameters.par convention.
    """
    if clean_output and case.workdir.exists():
        shutil.rmtree(case.workdir)
    case.workdir.mkdir(parents=True, exist_ok=True)

    env = os.environ.copy()
    env["OMP_NUM_THREADS"] = str(threads)

    cmd = list(mpi_cmd) + ["-np", str(ranks), str(binary), str(case.params_path)]
    print(f"{_info('[run]')} {case.key}: (cd {case.case_dir}) {shlex.join(cmd)}", flush=True)

    t0 = _dt.datetime.now()
    proc = subprocess.run(cmd, cwd=case.case_dir, env=env)
    dt = (_dt.datetime.now() - t0).total_seconds()
    return proc.returncode, dt


# ---------------------------------------------------------------------------
# Comparison dispatch
# ---------------------------------------------------------------------------


def _read_summary(path: Path) -> Optional[dict]:
    """Read a compare-*.py ``--report-json`` file, decoding inf/nan sentinels."""
    if not path.is_file():
        return None
    try:
        payload = json.loads(path.read_text())
    except (OSError, json.JSONDecodeError):
        return None

    def _decode(v):
        if isinstance(v, str):
            return {"inf": float("inf"), "-inf": float("-inf"),
                    "nan": float("nan")}.get(v, v)
        return v

    payload["max_error"] = _decode(payload.get("max_error"))
    payload["quantities"] = {
        k: _decode(v) for k, v in payload.get("quantities", {}).items()
    }
    return payload


def _run_compare(
    script: Path, args: Sequence[str], cwd: Path
) -> tuple[int, Optional[dict]]:
    """Run a compare-*.py script, returning ``(returncode, summary)``.

    A ``--report-json`` file is requested transparently so we can show the
    achieved error for every comparison -- even the ones that pass -- without
    parsing stdout. The summary is ``None`` if the script produced no report
    (e.g. an older compare-*.py that doesn't understand the flag).
    """
    with tempfile.NamedTemporaryFile(
        prefix="seissol-cmp-", suffix=".json", delete=False
    ) as tmp:
        report_path = Path(tmp.name)
    try:
        full_args = [*args, "--report-json", str(report_path)]
        cmd = [sys.executable, str(script), *full_args]
        print(f"{c('[cmp]', _Ansi.YELLOW)} {script.name} {' '.join(args)}", flush=True)
        rc = subprocess.run(cmd, cwd=cwd).returncode
        return rc, _read_summary(report_path)
    finally:
        report_path.unlink(missing_ok=True)


def _apply_summary(result: CompareResult, summary: Optional[dict]) -> CompareResult:
    """Copy the achieved error / per-quantity map from a summary into a result."""
    if summary is not None:
        result.has_summary = True
        result.achieved = summary.get("max_error")
        result.quantities = dict(summary.get("quantities", {}))
    return result


def _expect(path: Path, what: str) -> Optional[str]:
    if not path.is_file():
        return f"{what} missing: {path}"
    return None


def compare_mesh(case: CaseSpec, suffix: str, category: str, epsilon: float) -> CompareResult:
    """Volume/fault/surface XDMF comparison via compare-mesh.py."""
    out_xdmf = case.workdir / f"{case.prefix}{suffix}.xdmf"
    ref_xdmf = case.reference_dir / f"{case.prefix}{suffix}.xdmf"
    for path, what in [(out_xdmf, "output"), (ref_xdmf, "reference")]:
        msg = _expect(path, what)
        if msg:
            return CompareResult(category, False, epsilon, msg)
    rc, summary = _run_compare(
        COMPARE_MESH,
        [str(out_xdmf), str(ref_xdmf), "--epsilon", str(epsilon),
         "--category", category],
        cwd=case.case_dir,
    )
    return _apply_summary(CompareResult(category, rc == 0, epsilon), summary)


def compare_receivers(case: CaseSpec, epsilon: float) -> CompareResult:
    rc, summary = _run_compare(
        COMPARE_RECEIVERS,
        [str(case.workdir), str(case.reference_dir),
         "--prefix", case.prefix, "--epsilon", str(epsilon)],
        cwd=case.case_dir,
    )
    return _apply_summary(CompareResult("receiver", rc == 0, epsilon), summary)


def compare_energies(case: CaseSpec, epsilon: float) -> CompareResult:
    out_csv = case.workdir / f"{case.prefix}-energy.csv"
    ref_csv = case.reference_dir / f"{case.prefix}-energy.csv"
    for path, what in [(out_csv, "output"), (ref_csv, "reference")]:
        msg = _expect(path, what)
        if msg:
            return CompareResult("energy", False, epsilon, msg)
    rc, summary = _run_compare(
        COMPARE_ENERGIES,
        [str(out_csv), str(ref_csv), "--epsilon", str(epsilon)],
        cwd=case.case_dir,
    )
    return _apply_summary(CompareResult("energy", rc == 0, epsilon), summary)


def _rule_matches(when: dict, labels: set, case_key: str, config: str) -> bool:
    """A ``when`` spec matches when ALL of its keys match (AND within a rule).

    Keys: ``label`` (string or list -> any listed label is active), ``case``
    (exact key or a "key/" prefix), ``config`` (regex searched in the config).
    An empty ``when`` matches unconditionally.
    """
    if "label" in when:
        wanted = when["label"]
        wanted = [wanted] if isinstance(wanted, str) else list(wanted)
        if not (set(wanted) & labels):
            return False
    if "case" in when:
        c = when["case"]
        if not (case_key == c or case_key.startswith(c + "/")):
            return False
    if "config" in when:
        if not re.search(str(when["config"]), config):
            return False
    return True


def _first_matching_reason(rules, labels: set, case_key: str, config: str):
    """Reason of the first rule whose ``when`` matches (OR across rules), else None."""
    for rule in rules or []:
        if _rule_matches(rule.get("when", {}), labels, case_key, config):
            return rule.get("reason", "")
    return None


def _apply_xfail(cr: CompareResult, case: CaseSpec, tpv_data: dict,
                 category: str, labels: set) -> None:
    """Tag a comparison as xfail/xpass if a category- or case-level rule matches."""
    rules = list(tpv_data.get(category, {}).get("xfail", [])) + list(case.xfail)
    reason = _first_matching_reason(rules, labels, case.key, case.config)
    if reason is not None:
        cr.xfail_reason = reason
        cr.xfail_outcome = "xfail" if not cr.passed else "xpass"


def _category_skip_reason(case: CaseSpec, tpv_data: dict, category: str,
                          labels: set):
    """Reason if this category is conditionally skipped (category-level rules only)."""
    return _first_matching_reason(
        tpv_data.get(category, {}).get("skip", []), labels, case.key, case.config)


def verify_case(
    case: CaseSpec, tpv_data: dict, epsilon_override: Optional[float],
    labels: Optional[set] = None,
) -> list[CompareResult]:
    """Run all comparisons enabled in <prefix>-data.json for one case.

    Applies, per category: a conditional skip (rule match -> not compared),
    per-quantity thresholds, and an xfail/xpass tag (rule match -> a failure is
    expected). ``labels`` is the active environment label set (from ``--label``).
    """
    labels = labels or set()
    results: list[CompareResult] = []
    v = case.key

    def run(category: str, compare_fn) -> None:
        if not _enabled(tpv_data, category):
            return
        eps = _epsilon(tpv_data, category, epsilon_override, v)
        skip_reason = _category_skip_reason(case, tpv_data, category, labels)
        if skip_reason is not None:
            cr = CompareResult(category, True, eps)
            cr.skipped = True
            cr.skip_reason = skip_reason
            results.append(cr)
            return
        cr = compare_fn(eps)
        if cr.has_summary and _has_quantity_thresholds(tpv_data, category):
            _evaluate_quantities(cr, tpv_data, v, epsilon_override)
        _apply_xfail(cr, case, tpv_data, category, labels)
        results.append(cr)

    run("volume", lambda e: compare_mesh(case, "", "volume", e))
    run("fault", lambda e: compare_mesh(case, "-fault", "fault", e))
    run("surface", lambda e: compare_mesh(case, "-surface", "surface", e))
    run("energy", lambda e: compare_energies(case, e))
    run("receiver", lambda e: compare_receivers(case, e))
    return results


# ---------------------------------------------------------------------------
# Recording (--record)
# ---------------------------------------------------------------------------


def _list_quantities(script: Path, args: Sequence[str], cwd: Path) -> list[str]:
    """Run a compare-*.py in ``--list-quantities`` mode; return its key list.

    Returns ``[]`` on any error so a failed enumeration never breaks recording.
    """
    cmd = [sys.executable, str(script), *args, "--list-quantities"]
    try:
        proc = subprocess.run(cmd, cwd=cwd, capture_output=True, text=True)
        if proc.returncode != 0:
            print(f"{_warn('[warn]')} {script.name} --list-quantities failed: "
                  f"{proc.stderr.strip()}", file=sys.stderr)
            return []
        return list(json.loads(proc.stdout))
    except (OSError, json.JSONDecodeError) as e:
        print(f"{_warn('[warn]')} could not list quantities via {script.name}: {e}",
              file=sys.stderr)
        return []


def generate_data_skeleton(case: CaseSpec, target: Path) -> None:
    """Write a starter ``<prefix>-data.json`` derived from the recorded output.

    Enabled categories are detected from which output files are present; the
    per-quantity names come from the compare scripts' ``--list-quantities`` mode
    (so the keys match a real comparison exactly). Every threshold is the
    placeholder ``SKELETON_EPSILON`` -- see that constant. Only called when no
    data file already exists.
    """
    p = case.prefix
    data: dict = {}

    # Volume / fault / surface: one XDMF file each; quantities are its data fields.
    for category, suffix in (("volume", ""), ("fault", "-fault"), ("surface", "-surface")):
        xdmf = target / f"{p}{suffix}.xdmf"
        if xdmf.is_file():
            quants = _list_quantities(COMPARE_MESH, [str(xdmf)], cwd=case.case_dir)
            data[category] = _skeleton_category(quants)

    # Energy: a single CSV.
    energy_csv = target / f"{p}-energy.csv"
    if energy_csv.is_file():
        quants = _list_quantities(COMPARE_ENERGIES, [str(energy_csv)], cwd=case.case_dir)
        data["energy"] = _skeleton_category(quants)

    # Receivers: any receiver or faultreceiver .dat files.
    has_receivers = any(target.glob(f"{p}-receiver-*.dat")) or \
        any(target.glob(f"{p}-faultreceiver-*.dat"))
    if has_receivers:
        quants = _list_quantities(
            COMPARE_RECEIVERS, [str(target), "--prefix", p], cwd=case.case_dir)
        data["receiver"] = _skeleton_category(quants)

    out = target / f"{p}-data.json"
    out.write_text(json.dumps(data, indent=2) + "\n")
    cats = ", ".join(data) or "none"
    print(f"{c('[record]', _Ansi.MAGENTA)} {case.key}: generated {out.name} "
          f"(categories: {cats}; review the placeholder epsilons)")


def _skeleton_category(quantities: list[str]) -> dict:
    """A skeleton category block: enabled, a placeholder epsilon, and per-quantity
    placeholders for every enumerated name (prune/tune as needed)."""
    return {
        "enabled": True,
        "epsilon": SKELETON_EPSILON,
        "quantities": {q: SKELETON_EPSILON for q in quantities},
    }


def record_case(case: CaseSpec) -> int:
    """Overwrite this case's own ``output`` directory with the fresh output.

    Only cases that own an output directory are recorded; derived cases (which
    compare against another case's output) are skipped by the caller. Writing
    always targets the case's own ``output_dir`` -- never a borrowed reference.

    The per-case ``<prefix>-data.json`` (sitting in the output directory) is
    preserved; everything else matching the output prefix is replaced.
    """
    if case.output_dir is None:
        raise RuntimeError(f"Cannot record {case.key}: it declares no output directory")
    if not case.workdir.is_dir():
        raise RuntimeError(f"Cannot record: {case.workdir} does not exist")
    target = case.output_dir
    target.mkdir(parents=True, exist_ok=True)

    preserved = target / f"{case.prefix}-data.json"
    keep = preserved.read_bytes() if preserved.is_file() else None

    # Wipe and copy: deletes refs that no longer exist in the new output (e.g.
    # if you renamed a receiver). The preserved <prefix>-data.json is restored after.
    for f in target.iterdir():
        if f.is_file():
            f.unlink()
        if f.is_dir():
            shutil.rmtree(f)

    n = 0
    m = 0
    for src in case.workdir.iterdir():
        if src.name.startswith(case.prefix):
            if src.is_file():
                shutil.copy2(src, target / src.name)
                n += 1
            if src.is_dir():
                shutil.copytree(src, target / src.name)
                m += 1
    print(f"{c('[record]', _Ansi.MAGENTA)} {case.key}: wrote {n} files and {m} directories to {target}")

    if keep is not None:
        preserved.write_bytes(keep)          # keep a hand-authored data file as-is
    else:
        generate_data_skeleton(case, target)  # none existed -> start one automatically
    return n


def update_references_json(
    precomputed_dir: Path,
    case: CaseSpec,
    seissol_commit: str,
    note: str = "",
) -> None:
    """Append/update an entry in ``references.json`` (one entry per case name).

    Structure::

        {
          "tpv5": {
            "commit": "<sha>",
            "date":   "2026-05-21",
            "configs": ["elastic-o6-f32", "elastic-o6-f64"],
            "note":   "..."
          },
          ...
        }
    """
    f = precomputed_dir / "references.json"
    db = json.loads(f.read_text()) if f.is_file() else {}
    existing = db.get(case.key, {})
    configs = sorted(set(existing.get("configs", []) + [case.config]))
    db[case.key] = {
        "commit": seissol_commit,
        "date": _dt.date.today().isoformat(),
        "configs": configs,
        "note": note or existing.get("note", ""),
    }
    f.write_text(json.dumps(db, indent=2, sort_keys=True) + "\n")


# ---------------------------------------------------------------------------
# Selection and reporting
# ---------------------------------------------------------------------------


def select_cases(
    cases: dict[str, CaseSpec], requested: Iterable[str], all_: bool,
) -> list[CaseSpec]:
    if all_:
        return list(cases.values())
    selected: list[CaseSpec] = []
    seen: set[str] = set()
    for r in requested:
        # Match an exact key, or every key under the "<r>/" prefix (so "tpv5"
        # selects "tpv5/double", "tpv5/single", ... without a formal name split).
        matched = [c for c in cases.values()
                   if c.key == r or c.key.startswith(r + "/")]
        if not matched:
            raise KeyError(f"No case matches '{r}'")
        for c in matched:
            if c.key not in seen:
                selected.append(c)
                seen.add(c.key)
    return selected


def load_failed_keys(path: Path) -> list[str]:
    """Case keys that did not pass in a previous ``--report-json`` file."""
    data = json.loads(path.read_text())
    return [c["case"] for c in data.get("cases", []) if not c.get("passed", True)]


def list_cases(selected: list[CaseSpec], epsilon_override: Optional[float]) -> None:
    """Dry-run: show each selected case's reference and per-category epsilons."""
    for case in selected:
        tag = f"   (derived <- {case.reference})" if case.is_derived else ""
        print(f"{c(case.key, _Ansi.BOLD)}   config={case.config}   "
              f"ref={case.reference or 'self'}{tag}")
        tpv = load_tpv_data(case)
        enabled = [cat for cat in _CANONICAL_CATEGORIES if _enabled(tpv, cat)]
        if not enabled:
            print(f"    {_dim('(no categories enabled / no tpv-data.json)')}")
            continue
        for cat in enabled:
            eps = _epsilon(tpv, cat, epsilon_override, case.key)
            qmap = tpv.get(cat, {}).get("quantities", {})
            extra = f"   (+{len(qmap)} per-quantity)" if qmap else ""
            print(f"    {cat:<9} \u03b5={eps:g}{extra}")


def _state_of(r: "CaseResult") -> str:
    """One-word status for a case: PASS / FAIL / SKIP / REC / RC=<n>."""
    if r.skipped:
        return "SKIP"
    if r.recorded:
        return "REC"
    if r.ran and r.run_returncode != 0:
        return f"RC={r.run_returncode}"
    if not r.compares:
        return "FAIL"
    return "FAIL" if any(_compare_is_failure(c, r.xfail_strict) for c in r.compares) else "PASS"


_CANONICAL_CATEGORIES = ["volume", "fault", "surface", "energy", "receiver"]


def _category_columns(results: list["CaseResult"]) -> list[str]:
    """Category columns present across the run, canonical order first."""
    present = {cr.category for r in results for cr in r.compares}
    return ([c for c in _CANONICAL_CATEGORIES if c in present]
            + sorted(present - set(_CANONICAL_CATEGORIES)))


def _offenders(cr: "CompareResult") -> list[tuple[str, float, float]]:
    """Failing (quantity, achieved, threshold), worst first.

    Uses the per-quantity verdict when present; otherwise derives the offenders
    from the achieved map against the category epsilon, so a category-level
    failure lists the same way. Empty if there is no summary (hard failure).
    """
    if cr.failures:
        lst = list(cr.failures)
    elif cr.quantities:
        lst = [(q, e, cr.epsilon) for q, e in cr.quantities.items()
               if not (e <= cr.epsilon)]
    else:
        return []

    def sev(t: tuple[str, float, float]) -> float:
        _, a, thr = t
        return math.inf if (not math.isfinite(a) or thr <= 0) else a / thr

    return sorted(lst, key=sev, reverse=True)


def write_report(results: list[CaseResult], report_path: Optional[Path],
                 near_miss: float = 0.5, drift: Optional[list[tuple]] = None,
                 drift_factor: float = 2.0) -> None:
    """Render the run report.

    Two parts: a per-category grid (one row per case, one column per enabled
    category, cell = worst achieved error, ``✗`` on failure, ``—`` when the
    category is disabled for that case), followed by a ``Failures:`` block that
    breaks every failing case down per offending quantity, each shown against
    *its own* epsilon. The grid carries the full at-a-glance picture; the block
    is the readable "what actually broke" view.
    """
    MAX_QTY_ROWS = 6  # per category in the failure block, before "(+N more)"

    cols = _category_columns(results)

    key_w = max([len("case")] + [len(r.case.key) for r in results])
    CASE_W = min(max(key_w, 12), 44)
    STATE_W, TIME_W, COL_W = 6, 8, 10

    def num(x: Optional[float], prec: int) -> str:
        if x is None:
            return "n/a"
        if math.isinf(x):
            return "inf"
        if math.isnan(x):
            return "nan"
        return f"{x:.{prec}e}"

    state_of = _state_of        # module-level; shared with the JSON report
    offenders = _offenders      # module-level; shared with the Markdown summary

    # Pre-compute the failure block rows (color-independent) so widths line up
    # globally across every failing case.
    fail_blocks: list[tuple[str, list[tuple[str, str, Optional[float], Optional[float]]]]] = []
    for r in results:
        st = state_of(r)
        if st in ("PASS", "SKIP", "REC"):
            continue
        if st.startswith("RC="):
            fail_blocks.append((r.case.key,
                                [("run", f"exited with code {r.run_returncode}", None, None)]))
            continue
        cmap = {cr.category: cr for cr in r.compares}
        rows: list[tuple[str, str, Optional[float], Optional[float]]] = []
        for cat in cols:
            cr = cmap.get(cat)
            if cr is None or cr.passed or cr.skipped or cr.xfail_outcome == "xfail":
                continue  # skipped/expected-failure categories aren't real failures
            offs = offenders(cr)
            if offs:
                for q, a, thr in offs[:MAX_QTY_ROWS]:
                    rows.append((cat, q, a, thr))
                if len(offs) > MAX_QTY_ROWS:
                    rows.append((cat, f"(+{len(offs) - MAX_QTY_ROWS} more)", None, None))
            else:
                rows.append((cat, cr.detail or "failed", None, None))
        fail_blocks.append((r.case.key, rows))

    all_rows = [row for _, rows in fail_blocks for row in rows]
    cat_w = max([0] + [len(cat) for cat, _, _, _ in all_rows])
    qty_w = max([0] + [len(q) for _, q, _, _ in all_rows])
    num_w = max([0] + [len(num(v, 2)) for _, _, a, thr in all_rows
                       for v in (a, thr) if v is not None])

    drift = drift or []
    d_cat_w = max([0] + [len(cat) for _, cat, _, _, _, _ in drift])
    d_qty_w = max([0] + [len(q) for _, _, q, _, _, _ in drift])
    d_num_w = max([0] + [len(num(v, 2)) for _, _, _, b, cu, _ in drift
                         for v in (b, cu)])

    buckets = {"PASS": 0, "FAIL": 0, "SKIP": 0, "REC": 0}
    for r in results:
        st = state_of(r)
        buckets["FAIL" if st.startswith("RC=") else st] += 1

    # Per-category conditional outcomes: xfail (expected fail), xpass (unexpected
    # pass), skip (conditionally not compared) -- listed with their reasons.
    xfail_rows: list[tuple[str, str, str, str]] = []  # (case, category, kind, reason)
    for r in results:
        for cr in r.compares:
            if cr.skipped:
                xfail_rows.append((r.case.key, cr.category, "skip", cr.skip_reason))
            elif cr.xfail_outcome:
                xfail_rows.append((r.case.key, cr.category, cr.xfail_outcome, cr.xfail_reason))
    x_cat_w = max([0] + [len(cat) for _, cat, _, _ in xfail_rows])
    x_kind_w = max([0] + [len(k) for _, _, k, _ in xfail_rows])
    n_xfail = sum(1 for _, _, k, _ in xfail_rows if k == "xfail")
    n_xpass = sum(1 for _, _, k, _ in xfail_rows if k == "xpass")
    n_cskip = sum(1 for _, _, k, _ in xfail_rows if k == "skip")

    def render(colorize: bool) -> str:
        ok_ = (lambda s: _ok(s)) if colorize else (lambda s: s)
        bad_ = (lambda s: _bad(s)) if colorize else (lambda s: s)
        dim_ = (lambda s: _dim(s)) if colorize else (lambda s: s)
        bold_ = (lambda s: c(s, _Ansi.BOLD)) if colorize else (lambda s: s)
        mag_ = (lambda s: c(s, _Ansi.MAGENTA, _Ansi.BOLD)) if colorize else (lambda s: s)

        warn_ = (lambda s: _warn(s)) if colorize else (lambda s: s)

        def cell_render(cr: CompareResult) -> str:
            if cr.skipped:
                return dim_(f"{'skip':>{COL_W}} ")   # conditional category skip
            # Near-miss: passing, but the achieved max sits in [near_miss·ε, ε).
            ratio = 0.0
            if cr.achieved is not None and math.isfinite(cr.achieved) and cr.epsilon > 0:
                ratio = cr.achieved / cr.epsilon
            near = cr.passed and near_miss > 0 and near_miss <= ratio < 1.0
            # 'err' for a hard failure (no achieved value, e.g. missing file);
            # the reason is spelled out in the Failures block below.
            disp = "err" if cr.achieved is None else num(cr.achieved, 1)
            # Expected-failure states get their own low-key markers: 'x' for an
            # xfail (dim), '!' for an unexpected pass (yellow, wants attention).
            if cr.xfail_outcome == "xfail":
                return dim_(f"{disp:>{COL_W}}x")
            if cr.xfail_outcome == "xpass":
                return warn_(f"{disp:>{COL_W}}!")
            # Number right-aligned to the column edge (like the header); the mark
            # lives in its own trailing column so it never shifts the digits.
            mark = "✗" if not cr.passed else ("~" if near else " ")
            cell = f"{disp:>{COL_W}}{mark}"
            if not cr.passed:
                return bad_(cell)
            return warn_(cell) if near else ok_(cell)

        def state_cell(st: str) -> str:
            pad = f"{st:<{STATE_W}}"
            if st == "PASS":
                return ok_(pad)
            if st == "FAIL" or st.startswith("RC="):
                return bad_(pad)
            if st == "REC":
                return mag_(pad)
            if st == "SKIP":
                return dim_(pad)
            return pad

        def time_cell(r: CaseResult) -> str:
            # No run happened for skipped entries, so time carries no meaning.
            text = "—" if r.skipped else f"{r.duration_s:.1f}s"
            return dim_(f"{text:>{TIME_W}}")

        head = (f"{'case':<{CASE_W}} {'state':<{STATE_W}} {'time':>{TIME_W}}"
                + "".join(f"{cat:>{COL_W}} " for cat in cols))
        lines = [bold_(head), "─" * len(head)]

        for r in results:
            st = state_of(r)
            row = f"{r.case.key:<{CASE_W}} {state_cell(st)} {time_cell(r)}"
            if not r.compares:
                note = r.note or ("recorded" if r.recorded else "")
                lines.append(row + ("  " + dim_(note) if note else ""))
                continue
            cmap = {cr.category: cr for cr in r.compares}
            for cat in cols:
                cr = cmap.get(cat)
                if cr is None:
                    row += dim_(f"{'—':>{COL_W}} ")
                    continue
                row += cell_render(cr)
            lines.append(row.rstrip())

        lines.append("─" * len(head))
        evaluated = buckets["PASS"] + buckets["FAIL"]
        total_time = sum(r.duration_s for r in results)
        parts = [f"{buckets['PASS']}/{evaluated} passed"] if evaluated else ["no tests run"]
        if buckets["SKIP"]:
            parts.append(f"{buckets['SKIP']} skipped")
        if buckets["REC"]:
            parts.append(f"{buckets['REC']} recorded")
        if n_xfail:
            parts.append(f"{n_xfail} xfailed")
        if n_xpass:
            parts.append(f"{n_xpass} xpassed")
        if n_cskip:
            parts.append(f"{n_cskip} cat-skipped")
        parts.append(f"{total_time:.1f}s")
        summary = " · ".join(parts)
        lines.append(ok_(summary) if buckets["FAIL"] == 0 else bad_(summary))

        if fail_blocks:
            lines.append("")
            lines.append(bold_("Failures:"))
            for key, rows in fail_blocks:
                lines.append(f"  {key}")
                last_cat = None
                for cat, q, a, thr in rows:
                    catcell = f"{cat if cat != last_cat else '':<{cat_w}}"
                    last_cat = cat
                    if a is None and thr is None:
                        lines.append(f"    {catcell}  {dim_(q)}")
                    else:
                        line = (f"    {catcell}  {q:<{qty_w}}  "
                                f"{bad_(f'{num(a, 2):>{num_w}}')}  >  "
                                f"{dim_(f'{num(thr, 2):>{num_w}}')}")
                        lines.append(line)

        if xfail_rows:
            lines.append("")
            lines.append(bold_("Expected failures / skips:"))
            kind_color = {"xfail": dim_, "xpass": warn_, "skip": dim_}
            last_ck = None
            for ck, cat, kind, reason in xfail_rows:
                if ck != last_ck:
                    lines.append(f"  {ck}")
                    last_ck = ck
                paint = kind_color.get(kind, dim_)
                lines.append(f"    {cat:<{x_cat_w}}  {paint(f'{kind:<{x_kind_w}}')}  "
                             f"{dim_(reason)}")

        if drift:
            last_ck = None
            for ck, cat, q, base, cur, ratio in drift:
                if ck != last_ck:
                    lines.append(f"  {ck}")
                    last_ck = ck
                rs = "\u221e" if not math.isfinite(ratio) else f"{ratio:.1f}\u00d7"
                lines.append(
                    f"    {cat:<{d_cat_w}}  {q:<{d_qty_w}}  "
                    f"{dim_(f'{num(base, 2):>{d_num_w}}')} \u2192 "
                    f"{warn_(f'{num(cur, 2):>{d_num_w}}')}  ({rs})")
        return "\n".join(lines)

    print("\n" + render(colorize=True))
    if report_path:
        report_path.write_text(render(colorize=False) + "\n")


def _json_num(x: Optional[float]):
    """JSON-safe number: pass ints/floats through, encode inf/nan as strings."""
    if x is None:
        return None
    x = float(x)
    if math.isinf(x):
        return "inf" if x > 0 else "-inf"
    if math.isnan(x):
        return "nan"
    return x


def build_run_report(results: list[CaseResult],
                     drift: Optional[list[tuple]] = None) -> dict:
    """Assemble the whole run as a JSON-serialisable dict (schema 1).

    Carries, per case: identity (case key, config, reference, prefix,
    output_dir), state, wall time, and every comparison's category threshold,
    achieved max error, the full per-quantity achieved map, and any per-quantity
    threshold violations. inf/nan are encoded as the strings "inf"/"-inf"/"nan"
    so the file stays valid JSON (mirrors the compare-*.py report format).
    """
    buckets = {"PASS": 0, "FAIL": 0, "SKIP": 0, "REC": 0}
    for r in results:
        st = _state_of(r)
        buckets["FAIL" if st.startswith("RC=") else st] += 1

    cases = []
    for r in results:
        comparisons = [
            {
                "category": cr.category,
                "passed": cr.passed,
                "epsilon": cr.epsilon,
                "achieved": _json_num(cr.achieved),
                "quantities": {k: _json_num(v) for k, v in cr.quantities.items()},
                "failures": [
                    {"quantity": q, "achieved": _json_num(a), "threshold": _json_num(t)}
                    for (q, a, t) in cr.failures
                ],
                "xfail": cr.xfail_outcome or None,
                "xfail_reason": cr.xfail_reason or None,
                "skipped": cr.skipped,
                "skip_reason": cr.skip_reason or None,
            }
            for cr in r.compares
        ]
        cases.append({
            "case": r.case.key,
            "reference": r.case.reference,
            "config": r.case.config,
            "prefix": r.case.prefix,
            "output_dir": str(r.case.output_dir) if r.case.output_dir else None,
            "state": _state_of(r),
            "passed": r.passed,
            "ran": r.ran,
            "returncode": r.run_returncode,
            "time_s": round(r.duration_s, 3),
            "note": r.note,
            "comparisons": comparisons,
        })

    return {
        "schema": 1,
        "generated": _dt.datetime.now(_dt.timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "summary": {
            "passed": buckets["PASS"],
            "failed": buckets["FAIL"],
            "skipped": buckets["SKIP"],
            "recorded": buckets["REC"],
            "evaluated": buckets["PASS"] + buckets["FAIL"],
            "total_time_s": round(sum(r.duration_s for r in results), 3),
        },
        "cases": cases,
        "drift": [
            {"case": ck, "category": cat, "quantity": q,
             "baseline": _json_num(base), "current": _json_num(cur),
             "ratio": _json_num(ratio)}
            for (ck, cat, q, base, cur, ratio) in (drift or [])
        ],
    }


def write_run_report_json(results: list[CaseResult], path: Path,
                          drift: Optional[list[tuple]] = None) -> None:
    path.write_text(json.dumps(build_run_report(results, drift), indent=2) + "\n")


def render_markdown(results: list[CaseResult], near_miss: float = 0.5,
                    drift: Optional[list[tuple]] = None) -> str:
    """GitHub-flavored Markdown version of the report (grid + failures)."""
    cols = _category_columns(results)
    buckets = {"PASS": 0, "FAIL": 0, "SKIP": 0, "REC": 0}
    for r in results:
        st = _state_of(r)
        buckets["FAIL" if st.startswith("RC=") else st] += 1
    evaluated = buckets["PASS"] + buckets["FAIL"]
    total_time = sum(r.duration_s for r in results)

    def mdnum(x: Optional[float]) -> str:
        if x is None:
            return "—"
        if math.isinf(x):
            return "inf"
        if math.isnan(x):
            return "nan"
        return f"{x:.2e}"

    emoji = {"PASS": "✅", "FAIL": "❌", "SKIP": "⏭️", "REC": "💾"}

    parts = [f"{buckets['PASS']}/{evaluated} passed"] if evaluated else ["no tests run"]
    if buckets["SKIP"]:
        parts.append(f"{buckets['SKIP']} skipped")
    if buckets["REC"]:
        parts.append(f"{buckets['REC']} recorded")
    _nx = sum(1 for r in results for cr in r.compares if cr.xfail_outcome == "xfail")
    _np = sum(1 for r in results for cr in r.compares if cr.xfail_outcome == "xpass")
    _ns = sum(1 for r in results for cr in r.compares if cr.skipped)
    if _nx:
        parts.append(f"{_nx} xfailed")
    if _np:
        parts.append(f"{_np} xpassed")
    if _ns:
        parts.append(f"{_ns} cat-skipped")
    parts.append(f"{total_time:.1f}s")
    out = [f"## SeisSol verification — {' · '.join(parts)}", ""]

    header = "| case | state | time | " + " | ".join(cols) + " |"
    out.append(header)
    out.append("|" + "---|" * (3 + len(cols)))
    for r in results:
        st = _state_of(r)
        em = emoji.get("FAIL" if st.startswith("RC=") else st, "")
        state_txt = f"{em} {st}"
        tcell = "—" if r.skipped else f"{r.duration_s:.1f}s"
        if not r.compares:
            note = r.note or ("recorded" if r.recorded else "")
            if note:
                state_txt += f" — {note}"
            out.append(f"| {r.case.key} | {state_txt} | {tcell} | "
                       + " | ".join("" for _ in cols) + " |")
            continue
        cmap = {cr.category: cr for cr in r.compares}
        cells = []
        for cat in cols:
            cr = cmap.get(cat)
            if cr is None:
                cells.append("—")
                continue
            if cr.skipped:
                cells.append("skip")
                continue
            val = "err" if cr.achieved is None else mdnum(cr.achieved)
            ratio = (cr.achieved / cr.epsilon
                     if cr.achieved is not None and math.isfinite(cr.achieved)
                     and cr.epsilon > 0 else 0.0)
            if cr.xfail_outcome == "xfail":
                cells.append(f"~~{val}~~")           # expected failure
            elif cr.xfail_outcome == "xpass":
                cells.append(f"{val} ⚠️")            # unexpected pass
            elif not cr.passed:
                cells.append(f"**{val}**")           # failure
            elif near_miss > 0 and near_miss <= ratio < 1.0:
                cells.append(f"_{val}_")             # near-miss
            else:
                cells.append(val)
        out.append(f"| {r.case.key} | {state_txt} | {tcell} | " + " | ".join(cells) + " |")

    fails = [r for r in results if not r.passed and not r.skipped]
    if fails:
        out += ["", "<details><summary>Failures</summary>", ""]
        for r in fails:
            out.append(f"**{r.case.key}**")
            if r.ran and r.run_returncode != 0:
                out.append(f"- run exited with code {r.run_returncode}")
                out.append("")
                continue
            cmap = {cr.category: cr for cr in r.compares}
            for cat in _category_columns([r]):
                cr = cmap.get(cat)
                if cr is None or cr.passed or cr.skipped or cr.xfail_outcome == "xfail":
                    continue
                offs = _offenders(cr)
                if offs:
                    for q, a, thr in offs:
                        out.append(f"- `{cat}` — `{q}` {mdnum(a)} > {mdnum(thr)}")
                else:
                    out.append(f"- `{cat}` — {cr.detail or 'failed'}")
            out.append("")
        out.append("</details>")

    xfail_rows = [(r.case.key, cr.category,
                   "skip" if cr.skipped else cr.xfail_outcome,
                   cr.skip_reason if cr.skipped else cr.xfail_reason)
                  for r in results for cr in r.compares
                  if cr.skipped or cr.xfail_outcome]
    if xfail_rows:
        out += ["", "<details><summary>Expected failures / skips</summary>", ""]
        last = None
        for ck, cat, kind, reason in xfail_rows:
            if ck != last:
                out.append(f"**{ck}**")
                last = ck
            out.append(f"- `{cat}` — {kind}: {reason}")
        out.append("")
        out.append("</details>")

    if drift:
        out += ["", "<details><summary>Drift vs baseline</summary>", ""]
        last = None
        for ck, cat, q, base, cur, ratio in drift:
            if ck != last:
                out.append(f"**{ck}**")
                last = ck
            rs = "∞" if not math.isfinite(ratio) else f"{ratio:.1f}×"
            out.append(f"- `{cat}` — `{q}` {mdnum(base)} → {mdnum(cur)} ({rs})")
        out.append("")
        out.append("</details>")

    return "\n".join(out) + "\n"


def write_step_summary(results: list[CaseResult], mode: str,
                       near_miss: float = 0.5,
                       drift: Optional[list[tuple]] = None) -> None:
    """Append the Markdown report to ``$GITHUB_STEP_SUMMARY`` when appropriate.

    mode ``auto`` writes only inside GitHub Actions (the env var is set),
    ``always`` writes whenever the env var is set (and warns if it is not),
    ``never`` disables it. The target path is whatever ``GITHUB_STEP_SUMMARY``
    points at, so it can be redirected for local testing.
    """
    if mode == "never":
        return
    target = os.environ.get("GITHUB_STEP_SUMMARY")
    if not target:
        if mode == "always":
            print(f"{_warn('[warn]')} --step-summary=always but "
                  f"GITHUB_STEP_SUMMARY is not set; skipping", file=sys.stderr)
        return
    with open(target, "a") as fh:  # GitHub step summaries are append-only
        fh.write(render_markdown(results, near_miss, drift))


# --- drift detection --------------------------------------------------

DRIFT_FLOOR = 1e-9  # ignore drift in the numerical noise floor


def _decode_num(x):
    """Decode a JSON number that may be an "inf"/"-inf"/"nan" sentinel string."""
    if isinstance(x, str):
        return {"inf": math.inf, "-inf": -math.inf, "nan": math.nan}.get(x, math.nan)
    return x


def load_baseline_achieved(path: Path) -> dict:
    """Map (case, category, quantity) -> achieved error from a prior report JSON."""
    data = json.loads(path.read_text())
    out: dict = {}
    for case in data.get("cases", []):
        ck = case.get("case", "")
        for cmp_ in case.get("comparisons", []):
            cat = cmp_.get("category", "")
            for q, v in cmp_.get("quantities", {}).items():
                out[(ck, cat, q)] = _decode_num(v)
    return out


def compute_drift(results: list[CaseResult], baseline: dict,
                  factor: float) -> list[tuple]:
    """Find (case, category, quantity) whose achieved error grew >= factor.

    Returns ``(case, category, quantity, baseline, current, ratio)`` worst first.
    Only meaningful, above-floor values count; a baseline of ~0 (perfect match)
    is skipped since any ratio there is meaningless.
    """
    rows = []
    for r in results:
        for cr in r.compares:
            for q, cur in cr.quantities.items():
                base = baseline.get((r.case.key, cr.category, q))
                if base is None or not math.isfinite(base) or base <= DRIFT_FLOOR:
                    continue
                if math.isfinite(cur) and cur < DRIFT_FLOOR:
                    continue
                ratio = math.inf if not math.isfinite(cur) else cur / base
                if ratio >= factor:
                    rows.append((r.case.key, cr.category, q, base, cur, ratio))
    rows.sort(key=lambda t: t[5], reverse=True)
    return rows


# --- epsilon suggestion (#4) -----------------------------------------------

def _ceil_1sig(x: float) -> float:
    """Round x UP to one significant figure (clean, never below the input)."""
    if not math.isfinite(x) or x <= 0:
        return x
    step = 10.0 ** math.floor(math.log10(x))
    return math.ceil(x / step) * step


def suggest_epsilons(results: list[CaseResult], factor: float) -> dict:
    """Per case, suggest a category epsilon = ceil_1sig(max achieved * factor).

    Uses only comparisons that produced a summary; categories whose achieved
    error is non-finite are skipped (they can't seed a threshold). Returns
    ``{case_key: {category: epsilon}}`` ready to drop into tpv-data.json.
    """
    out: dict = {}
    for r in results:
        cats: dict = {}
        for cr in r.compares:
            if not cr.has_summary or cr.achieved is None or not math.isfinite(cr.achieved):
                continue
            cats[cr.category] = _ceil_1sig(cr.achieved * factor)
        if cats:
            out[r.case.key] = cats
    return out


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description=__doc__.split("\n", 1)[0] if __doc__ else None,
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    p.add_argument("--precomputed-dir", type=Path, default=Path.cwd(),
                   help="precomputed-seissol checkout (default: cwd).")
    p.add_argument("--binary-dir", type=Path, default=None,
                   help="Directory holding seissol-<config> binaries "
                        "(not required with --list).")
    sel = p.add_argument_group("Case selection")
    sel.add_argument("--case", action="append", default=[],
                     help="Case key to run: an exact key like 'tpv5/double', or "
                          "a prefix like 'tpv5' to select every 'tpv5/...' key. "
                          "Repeatable.")
    sel.add_argument("--all", action="store_true",
                     help="Run every entry in cases.json.")
    sel.add_argument("--rerun-failed", type=Path, default=None, metavar="REPORT_JSON",
                     help="Re-run only the cases that failed in a previous "
                          "--report-json file (overrides --case/--all).")
    sel.add_argument("--list", action="store_true",
                     help="List the selected cases with their resolved reference "
                          "and per-category epsilons, then exit without running.")

    run = p.add_argument_group("Run options")
    run.add_argument("--mpi-cmd", default="mpirun",
                     help="MPI launcher (default: mpirun).")
    run.add_argument("--mpi-args", default="",
                     help="Extra MPI launcher args, shell-split "
                          "(e.g. \"--oversubscribe --allow-run-as-root\").")
    run.add_argument("--ranks", type=int, default=1)
    run.add_argument("--threads", type=int, default=1)
    run.add_argument("--no-run", action="store_true",
                     help="Skip the simulation; verify the existing output/.")
    run.add_argument("--no-clean", action="store_true",
                     help="Don't wipe output/ before running.")

    cmp_ = p.add_argument_group("Verification options")
    cmp_.add_argument("--epsilon", type=float, default=None,
                      help="Override comparison epsilon for ALL categories.")
    cmp_.add_argument("--label", action="append", default=[], metavar="NAME",
                      help="Active environment label for xfail/skip rules (e.g. "
                           "'gpu', 'cpu', 'cuda'). Repeatable.")
    cmp_.add_argument("--xfail-strict", action="store_true",
                      help="Treat an unexpected pass (xpass) as a failure "
                           "(default: warn only).")

    rec = p.add_argument_group("Re-record options (--record)")
    rec.add_argument("--record", action="store_true",
                     help="Overwrite the precomputed reference with the new output.")
    rec.add_argument("--seissol-commit", default="",
                     help="SeisSol commit hash, stored in references.json.")
    rec.add_argument("--note", default="",
                     help="Free-form note stored in references.json.")

    out = p.add_argument_group("Output")
    out.add_argument("--report", type=Path, default=None,
                     help="Also write a plain-text summary report to this file.")
    out.add_argument("--report-json", type=Path, default=None,
                     help="Also write the full run report as JSON to this file "
                          "(per-case state, time, achieved epsilon per category "
                          "and per quantity, and per-quantity failures).")
    out.add_argument("--color", choices=("auto", "never", "always"), default="auto",
                     help="Colorize stdout. 'auto' = TTY+!NO_COLOR. (default: auto)")
    out.add_argument("--step-summary", choices=("auto", "always", "never"),
                     default="auto",
                     help="Append a Markdown report to $GITHUB_STEP_SUMMARY. "
                          "'auto' = only when that variable is set (i.e. in "
                          "GitHub Actions). (default: auto)")
    out.add_argument("--near-miss", type=float, default=0.5, metavar="FRAC",
                     help="Highlight passing categories whose achieved error is "
                          "at least FRAC of their epsilon (0 disables). "
                          "(default: 0.5)")
    out.add_argument("--fail-fast", action="store_true",
                     help="Abort after the first failing case.")

    ana = p.add_argument_group("Analysis (post-run)")
    ana.add_argument("--drift-baseline", type=Path, default=None, metavar="REPORT_JSON",
                     help="Compare achieved errors against a previous "
                          "--report-json and warn about per-quantity regressions "
                          "(even when still passing).")
    ana.add_argument("--drift-factor", type=float, default=2.0, metavar="X",
                     help="Flag a quantity when its achieved error grew by at "
                          "least this factor vs the baseline. (default: 2.0)")
    ana.add_argument("--suggest-epsilons", action="store_true",
                     help="After the run, print suggested per-category epsilons "
                          "(max achieved * factor, rounded up) as a tpv-data.json "
                          "fragment.")
    ana.add_argument("--suggest-factor", type=float, default=3.0, metavar="X",
                     help="Safety factor for --suggest-epsilons. (default: 3.0)")
    ana.add_argument("--suggest-out", type=Path, default=None, metavar="PATH",
                     help="Write the epsilon suggestions here instead of stdout.")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    setup_color(args.color)

    if not (args.case or args.all or args.rerun_failed):
        print("error: pass --case ..., --all, or --rerun-failed", file=sys.stderr)
        return 2
    if args.record and args.no_run:
        print("error: --record requires a fresh run; remove --no-run", file=sys.stderr)
        return 2

    precomputed_dir = args.precomputed_dir.resolve()
    cases = load_cases(precomputed_dir)

    if args.rerun_failed:
        keys = load_failed_keys(args.rerun_failed)
        selected = []
        for k in keys:
            if k in cases:
                selected.append(cases[k])
            else:
                print(f"{_warn('[warn]')} {k}: in report but not in cases.json; "
                      f"skipping", file=sys.stderr)
        if not selected:
            print("No failed cases to re-run.", file=sys.stderr)
            return 0
    else:
        selected = select_cases(cases, args.case, args.all)
    if not selected:
        print("No cases selected.", file=sys.stderr)
        return 2

    if args.list:
        list_cases(selected, args.epsilon)
        return 0

    if args.binary_dir is None:
        print("error: --binary-dir is required (except with --list)", file=sys.stderr)
        return 2

    mpi_cmd = [args.mpi_cmd, *shlex.split(args.mpi_args)]
    labels = set(args.label)
    results: list[CaseResult] = []

    for case in selected:
        print(f"\n{c('=== ' + case.key + '  config=' + case.config + ' ===', _Ansi.CYAN, _Ansi.BOLD)}")

        # Case-level conditional skip (cases.json "skip"): matches the active
        # labels / this case's key / config -> not run, recorded, or verified.
        skip_reason = _first_matching_reason(case.skip, labels, case.key, case.config)
        if skip_reason is not None:
            print(f"{_warn('[skip]')} {case.key}: {skip_reason}")
            results.append(CaseResult(case, False, 0, 0.0, [], skipped=True,
                                      note=f"skip: {skip_reason}",
                                      xfail_strict=args.xfail_strict))
            continue

        # A derived entry owns no reference, so there is nothing to record for
        # it: its reference belongs to another variant and is recorded under that
        # variant's key. Skip it entirely in record mode.
        if args.record and case.is_derived:
            note = f"derived <- {case.reference}"
            print(f"{_warn('[skip]')} {case.key}: {note}; nothing to record")
            results.append(CaseResult(case, False, 0, 0.0, [], skipped=True, note=note))
            continue

        try:
            binary = resolve_binary(args.binary_dir, case.config)
        except FileNotFoundError as e:
            print(f"{_warn('[skip]')} {case.key}: {e}", file=sys.stderr)
            results.append(CaseResult(case, False, -1, 0.0, [], skipped=True,
                                      note="no binary"))
            if args.fail_fast:
                break
            continue

        ran, rc, duration = False, 0, 0.0
        if not args.no_run:
            rc, duration = run_seissol(
                case, binary, mpi_cmd, args.ranks, args.threads,
                clean_output=not args.no_clean,
            )
            ran = True
            if rc != 0:
                print(f"{_bad('[fail]')} {case.key}: SeisSol exited with {rc}", file=sys.stderr)
                results.append(CaseResult(case, True, rc, duration, []))
                if args.fail_fast:
                    break
                continue

        if args.record:
            try:
                record_case(case)
                if args.seissol_commit:
                    update_references_json(
                        precomputed_dir, case, args.seissol_commit, args.note
                    )
                results.append(CaseResult(case, ran, rc, duration, [], recorded=True))
            except Exception as e:
                print(f"{_bad('[fail]')} {case.key}: record failed: {e}", file=sys.stderr)
                results.append(CaseResult(
                    case, ran, rc, duration,
                    [CompareResult("record", False, 0.0, str(e))],
                ))
                if args.fail_fast:
                    break
            continue

        tpv_data = load_tpv_data(case)
        compares = verify_case(case, tpv_data, args.epsilon, labels)
        cr = CaseResult(case, ran, rc, duration, compares,
                        xfail_strict=args.xfail_strict)
        results.append(cr)
        if not cr.passed and args.fail_fast:
            break

    drift_rows: list[tuple] = []
    if args.drift_baseline:
        baseline = load_baseline_achieved(args.drift_baseline)
        drift_rows = compute_drift(results, baseline, args.drift_factor)
        if drift_rows:
            print(f"{_warn('[drift]')} {len(drift_rows)} quantity(ies) grew "
                  f"\u2265{args.drift_factor:g}\u00d7 vs baseline "
                  f"(see the Drift section)", file=sys.stderr)

    write_report(results, args.report, near_miss=args.near_miss,
                 drift=drift_rows, drift_factor=args.drift_factor)
    if args.report_json:
        write_run_report_json(results, args.report_json, drift=drift_rows)
    write_step_summary(results, args.step_summary, near_miss=args.near_miss,
                       drift=drift_rows)

    if args.suggest_epsilons:
        suggestions = suggest_epsilons(results, args.suggest_factor)
        text = json.dumps(suggestions, indent=2) + "\n"
        if any(not r.passed for r in results):
            print(f"{_warn('[warn]')} suggesting epsilons from a run with "
                  f"failures; review before use", file=sys.stderr)
        if args.suggest_out:
            args.suggest_out.write_text(text)
            print(f"Wrote epsilon suggestions to {args.suggest_out}")
        else:
            print(f"\nSuggested epsilons (max achieved \u00d7 "
                  f"{args.suggest_factor:g}, rounded up):")
            print(text)

    return 0 if all(r.passed for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
