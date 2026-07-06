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
OUTPUT_PREFIX = "tpv"     # SeisSol output basename used across all cases
OUTPUT_SUBDIR = "output"  # relative to the case directory


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
    """One ``(case, precision)`` entry from ``cases.json``."""

    name: str          # e.g. "tpv5"
    precision: str     # own variant label / 2nd key component (conventionally a precision)
    config: str        # e.g. "elastic-o6-f64-s8"
    params_path: Path  # absolute path to parameters.par
    case_dir: Path     # parent of parameters.par
    output_prefix: str # value of the ``output`` field, kept for diagnostics
    reference: str = ""  # label of the precomputed tree to compare against; empty => own

    def __post_init__(self) -> None:
        # ``reference`` is an opaque label naming a sibling ``precomputed/<label>/``
        # tree to compare against. It is NOT tied to precision: it can point at a
        # different precision, a different order, or any other recorded variant.
        # When it differs from this entry's own label, the entry is *derived*: it
        # is run and compared against another variant's reference, and is skipped
        # during ``--record`` (its reference is owned elsewhere).
        if not self.reference:
            self.reference = self.precision

    @property
    def key(self) -> str:
        return f"{self.name}/{self.precision}"

    @property
    def is_derived(self) -> bool:
        """True if this entry borrows another variant's reference tree."""
        return self.reference != self.precision

    @property
    def precomputed_dir(self) -> Path:
        """Reference tree to compare *against* (may belong to another variant)."""
        return self.case_dir / "precomputed" / self.reference

    @property
    def own_precomputed_dir(self) -> Path:
        """This entry's own reference tree -- the target when recording."""
        return self.case_dir / "precomputed" / self.precision

    @property
    def workdir(self) -> Path:
        """Where SeisSol writes its output."""
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


@dataclasses.dataclass
class CaseResult:
    case: CaseSpec
    ran: bool
    run_returncode: int
    duration_s: float
    compares: list[CompareResult]
    recorded: bool = False
    skipped: bool = False   # not evaluated (no binary / derived-precision record)
    note: str = ""          # why it was skipped, shown in the report

    @property
    def passed(self) -> bool:
        if self.skipped:
            return True
        if self.ran and self.run_returncode != 0:
            return False
        if self.recorded:
            return True
        return bool(self.compares) and all(c.passed for c in self.compares)


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
        if "/" not in key:
            raise ValueError(
                f"Malformed key '{key}' in cases.json (expected 'name/precision')"
            )
        name, precision = key.split("/", 1)
        params_path = (precomputed_dir / entry["parameters"]).resolve()
        cases[key] = CaseSpec(
            name=name,
            precision=precision,
            config=entry["config"],
            params_path=params_path,
            case_dir=params_path.parent,
            output_prefix=entry["output"],
            reference=entry.get("reference", ""),
        )
    return cases


def load_tpv_data(case: CaseSpec) -> dict:
    """Per-case validation config; missing file => empty (all categories off)."""
    f = case.precomputed_dir / "tpv-data.json"
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
    out_xdmf = case.workdir / f"{OUTPUT_PREFIX}{suffix}.xdmf"
    ref_xdmf = case.precomputed_dir / f"{OUTPUT_PREFIX}{suffix}.xdmf"
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
        [str(case.workdir), str(case.precomputed_dir),
         "--prefix", OUTPUT_PREFIX, "--epsilon", str(epsilon)],
        cwd=case.case_dir,
    )
    return _apply_summary(CompareResult("receiver", rc == 0, epsilon), summary)


def compare_energies(case: CaseSpec, epsilon: float) -> CompareResult:
    out_csv = case.workdir / f"{OUTPUT_PREFIX}-energy.csv"
    ref_csv = case.precomputed_dir / f"{OUTPUT_PREFIX}-energy.csv"
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


def verify_case(
    case: CaseSpec, tpv_data: dict, epsilon_override: Optional[float]
) -> list[CompareResult]:
    """Run all comparisons enabled in tpv-data.json for one case."""
    results: list[CompareResult] = []
    v = case.precision
    if _enabled(tpv_data, "volume"):
        results.append(compare_mesh(case, "", "volume",
                                    _epsilon(tpv_data, "volume", epsilon_override, v)))
    if _enabled(tpv_data, "fault"):
        results.append(compare_mesh(case, "-fault", "fault",
                                    _epsilon(tpv_data, "fault", epsilon_override, v)))
    if _enabled(tpv_data, "surface"):
        results.append(compare_mesh(case, "-surface", "surface",
                                    _epsilon(tpv_data, "surface", epsilon_override, v)))
    if _enabled(tpv_data, "energy"):
        results.append(compare_energies(case,
                                        _epsilon(tpv_data, "energy", epsilon_override, v)))
    if _enabled(tpv_data, "receiver"):
        results.append(compare_receivers(
            case, _epsilon(tpv_data, "receiver", epsilon_override, v),
        ))

    # Per-quantity thresholds (opt-in): when a category defines a "quantities"
    # map and a report JSON was produced, verify.py -- not the compare script's
    # exit code -- decides that category's verdict from the achieved per-quantity
    # errors. With no "quantities" map this pass is skipped, so behaviour is
    # unchanged. When the summary is missing (hard failure) the rc-based verdict
    # stands.
    for r in results:
        if r.has_summary and _has_quantity_thresholds(tpv_data, r.category):
            _evaluate_quantities(r, tpv_data, v, epsilon_override)
    return results


# ---------------------------------------------------------------------------
# Recording (--record)
# ---------------------------------------------------------------------------


def record_case(case: CaseSpec) -> int:
    """Overwrite this precision's own ``precomputed/<precision>/`` with fresh output.

    Only entries that own a reference are recorded; derived entries (whose
    ``reference`` label points elsewhere) are skipped by the caller. Writing
    always targets ``own_precomputed_dir`` -- never the borrowed reference tree.

    The per-case ``tpv-data.json`` (sitting inside the precomputed dir) is
    preserved; everything else matching the output prefix is replaced.
    """
    if not case.workdir.is_dir():
        raise RuntimeError(f"Cannot record: {case.workdir} does not exist")
    target = case.own_precomputed_dir
    target.mkdir(parents=True, exist_ok=True)

    preserved = target / "tpv-data.json"
    keep = preserved.read_bytes() if preserved.is_file() else None

    # Wipe and copy: deletes refs that no longer exist in the new output (e.g.
    # if you renamed a receiver). The preserved tpv-data.json is restored after.
    for f in target.iterdir():
        if f.is_file():
            f.unlink()

    n = 0
    for src in case.workdir.iterdir():
        if src.is_file() and src.name.startswith(OUTPUT_PREFIX):
            shutil.copy2(src, target / src.name)
            n += 1
    if keep is not None:
        preserved.write_bytes(keep)

    print(f"{c('[record]', _Ansi.MAGENTA)} {case.key}: wrote {n} files to {target}")
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
    existing = db.get(case.name, {})
    configs = sorted(set(existing.get("configs", []) + [case.config]))
    db[case.name] = {
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
        if "/" in r:
            if r not in cases:
                raise KeyError(f"Unknown case key: {r}")
            if r not in seen:
                selected.append(cases[r])
                seen.add(r)
        else:
            matched = [c for c in cases.values() if c.name == r]
            if not matched:
                raise KeyError(f"Unknown case name: {r}")
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
        tag = "   (derived, records nothing)" if case.is_derived else ""
        print(f"{c(case.key, _Ansi.BOLD)}   config={case.config}   "
              f"ref={case.reference}{tag}")
        tpv = load_tpv_data(case)
        enabled = [cat for cat in _CANONICAL_CATEGORIES if _enabled(tpv, cat)]
        if not enabled:
            print(f"    {_dim('(no categories enabled / no tpv-data.json)')}")
            continue
        for cat in enabled:
            eps = _epsilon(tpv, cat, epsilon_override, case.precision)
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
    ok = sum(cr.passed for cr in r.compares)
    return "PASS" if r.compares and ok == len(r.compares) else "FAIL"


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
                 near_miss: float = 0.5) -> None:
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
            if cr is None or cr.passed:
                continue
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

    buckets = {"PASS": 0, "FAIL": 0, "SKIP": 0, "REC": 0}
    for r in results:
        st = state_of(r)
        buckets["FAIL" if st.startswith("RC=") else st] += 1

    def render(colorize: bool) -> str:
        ok_ = (lambda s: _ok(s)) if colorize else (lambda s: s)
        bad_ = (lambda s: _bad(s)) if colorize else (lambda s: s)
        dim_ = (lambda s: _dim(s)) if colorize else (lambda s: s)
        bold_ = (lambda s: c(s, _Ansi.BOLD)) if colorize else (lambda s: s)
        mag_ = (lambda s: c(s, _Ansi.MAGENTA, _Ansi.BOLD)) if colorize else (lambda s: s)

        warn_ = (lambda s: _warn(s)) if colorize else (lambda s: s)

        def cell_render(cr: CompareResult) -> str:
            # Near-miss: passing, but the achieved max sits in [near_miss·ε, ε).
            ratio = 0.0
            if cr.achieved is not None and math.isfinite(cr.achieved) and cr.epsilon > 0:
                ratio = cr.achieved / cr.epsilon
            near = cr.passed and near_miss > 0 and near_miss <= ratio < 1.0
            mark = "✗" if not cr.passed else ("~" if near else " ")
            cell = f"{num(cr.achieved, 1):>{COL_W - 1}}{mark}"
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
                + "".join(f"{cat:>{COL_W}}" for cat in cols))
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
                    row += dim_(f"{'—':>{COL_W}}")
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


def build_run_report(results: list[CaseResult]) -> dict:
    """Assemble the whole run as a JSON-serialisable dict (schema 1).

    Carries, per case: identity (name/precision/reference/config), state, wall
    time, and every comparison's category threshold, achieved max error, the
    full per-quantity achieved map, and any per-quantity threshold violations.
    inf/nan are encoded as the strings "inf"/"-inf"/"nan" so the file stays valid
    JSON (mirrors the compare-*.py report format).
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
            }
            for cr in r.compares
        ]
        cases.append({
            "case": r.case.key,
            "name": r.case.name,
            "precision": r.case.precision,
            "reference": r.case.reference,
            "config": r.case.config,
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
    }


def write_run_report_json(results: list[CaseResult], path: Path) -> None:
    path.write_text(json.dumps(build_run_report(results), indent=2) + "\n")


def render_markdown(results: list[CaseResult], near_miss: float = 0.5) -> str:
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
            val = mdnum(cr.achieved)
            ratio = (cr.achieved / cr.epsilon
                     if cr.achieved is not None and math.isfinite(cr.achieved)
                     and cr.epsilon > 0 else 0.0)
            if not cr.passed:
                cells.append(f"**{val}**")          # failure
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
                if cr is None or cr.passed:
                    continue
                offs = _offenders(cr)
                if offs:
                    for q, a, thr in offs:
                        out.append(f"- `{cat}` — `{q}` {mdnum(a)} > {mdnum(thr)}")
                else:
                    out.append(f"- `{cat}` — {cr.detail or 'failed'}")
            out.append("")
        out.append("</details>")

    return "\n".join(out) + "\n"


def write_step_summary(results: list[CaseResult], mode: str,
                       near_miss: float = 0.5) -> None:
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
        fh.write(render_markdown(results, near_miss))


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
                     help="Case to run: 'tpv5' for both precisions, "
                          "'tpv5/double' for one. Repeatable.")
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
    results: list[CaseResult] = []

    for case in selected:
        print(f"\n{c('=== ' + case.key + '  config=' + case.config + ' ===', _Ansi.CYAN, _Ansi.BOLD)}")

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
        compares = verify_case(case, tpv_data, args.epsilon)
        cr = CaseResult(case, ran, rc, duration, compares)
        results.append(cr)
        if not cr.passed and args.fail_fast:
            break

    write_report(results, args.report, near_miss=args.near_miss)
    if args.report_json:
        write_run_report_json(results, args.report_json)
    write_step_summary(results, args.step_summary, near_miss=args.near_miss)
    return 0 if all(r.passed for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
