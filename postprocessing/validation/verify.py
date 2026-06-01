#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
"""Run SeisSol against precomputed reference solutions and verify the outputs.

Reads ``cases.json`` (in --precomputed-dir) and per-case ``tpv-data.json`` files,
dispatches the appropriate SeisSol binary from --binary-dir per case, then calls
the existing ``compare-{mesh,receivers,energies}.py`` scripts next to this file.

Binary naming (pre-PR #1421): one binary per (equation, precision, fused-count)
combination, named ``seissol-cpu-<config>`` where ``<config>`` is read directly
from cases.json (e.g. ``elastic-p6-f64``, ``viscoelastic-3-p6-f64-s8``). Once
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
import os
import shlex
import shutil
import subprocess
import sys
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
    precision: str     # "double" or "single"
    config: str        # e.g. "elastic-p6-f64-s8"
    params_path: Path  # absolute path to parameters.par
    case_dir: Path     # parent of parameters.par
    output_prefix: str # value of the ``output`` field, kept for diagnostics

    @property
    def key(self) -> str:
        return f"{self.name}/{self.precision}"

    @property
    def precomputed_dir(self) -> Path:
        return self.case_dir / "precomputed" / self.precision

    @property
    def workdir(self) -> Path:
        """Where SeisSol writes its output."""
        return self.case_dir / OUTPUT_SUBDIR


@dataclasses.dataclass
class CompareResult:
    category: str
    passed: bool
    epsilon: float
    detail: str = ""


@dataclasses.dataclass
class CaseResult:
    case: CaseSpec
    ran: bool
    run_returncode: int
    duration_s: float
    compares: list[CompareResult]
    recorded: bool = False

    @property
    def passed(self) -> bool:
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


def _epsilon(tpv_data: dict, category: str, override: Optional[float]) -> float:
    if override is not None:
        return override
    return float(tpv_data.get(category, {}).get("epsilon", DEFAULT_EPSILON))


# ---------------------------------------------------------------------------
# Binary resolution
# ---------------------------------------------------------------------------


VARIANT_ORDER = {
    "auto": ("cpu", "gpu"),
    "cpu": ("cpu",),
    "gpu": ("gpu",),
}


def resolve_binary(binary_dir: Path, config: str, variant: str = "auto") -> Path:
    """Locate a SeisSol executable for the given config.

    Pre-PR #1421 layout: ``binary_dir/seissol-{cpu,gpu}-<config>``. Which prefix
    is preferred depends on ``variant``:

    * ``auto`` (default): try ``cpu`` first, then ``gpu``. Predictable on a
      mixed install; pass ``--variant gpu`` explicitly when CPU binaries also
      happen to be present on a GPU host.
    * ``cpu`` / ``gpu``: only that prefix is searched.

    Post-PR #1421 layout: a single ``binary_dir/seissol`` (the config string
    becomes informational); this is tried last as a fallback for any variant.
    Each prefix is searched in both ``binary_dir`` and ``binary_dir/bin``.
    """
    prefixes = VARIANT_ORDER[variant]
    candidates: list[Path] = []
    for pfx in prefixes:
        candidates.append(binary_dir / f"seissol-{pfx}-{config}")
        candidates.append(binary_dir / "bin" / f"seissol-{pfx}-{config}")
    candidates.append(binary_dir / "seissol")
    candidates.append(binary_dir / "bin" / "seissol")

    for cand in candidates:
        if cand.is_file() and os.access(cand, os.X_OK):
            return cand
    raise FileNotFoundError(
        f"No SeisSol binary found for config '{config}' (variant={variant}). "
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


def _call_compare(script: Path, args: Sequence[str], cwd: Path) -> int:
    cmd = [sys.executable, str(script), *args]
    print(f"{c('[cmp]', _Ansi.YELLOW)} {script.name} {' '.join(args)}", flush=True)
    return subprocess.run(cmd, cwd=cwd).returncode


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
    rc = _call_compare(
        COMPARE_MESH,
        [str(out_xdmf), str(ref_xdmf), "--epsilon", str(epsilon)],
        cwd=case.case_dir,
    )
    return CompareResult(category, rc == 0, epsilon)


def compare_receivers(case: CaseSpec, epsilon: float) -> CompareResult:
    rc = _call_compare(
        COMPARE_RECEIVERS,
        [str(case.workdir), str(case.precomputed_dir),
         "--prefix", OUTPUT_PREFIX, "--epsilon", str(epsilon)],
        cwd=case.case_dir,
    )
    return CompareResult("receiver", rc == 0, epsilon)


def compare_energies(case: CaseSpec, epsilon: float) -> CompareResult:
    out_csv = case.workdir / f"{OUTPUT_PREFIX}-energy.csv"
    ref_csv = case.precomputed_dir / f"{OUTPUT_PREFIX}-energy.csv"
    for path, what in [(out_csv, "output"), (ref_csv, "reference")]:
        msg = _expect(path, what)
        if msg:
            return CompareResult("energy", False, epsilon, msg)
    rc = _call_compare(
        COMPARE_ENERGIES,
        [str(out_csv), str(ref_csv), "--epsilon", str(epsilon)],
        cwd=case.case_dir,
    )
    return CompareResult("energy", rc == 0, epsilon)


def verify_case(
    case: CaseSpec, tpv_data: dict, epsilon_override: Optional[float]
) -> list[CompareResult]:
    """Run all comparisons enabled in tpv-data.json for one case."""
    results: list[CompareResult] = []
    if _enabled(tpv_data, "volume"):
        results.append(compare_mesh(case, "", "volume",
                                    _epsilon(tpv_data, "volume", epsilon_override)))
    if _enabled(tpv_data, "fault"):
        results.append(compare_mesh(case, "-fault", "fault",
                                    _epsilon(tpv_data, "fault", epsilon_override)))
    if _enabled(tpv_data, "surface"):
        results.append(compare_mesh(case, "-surface", "surface",
                                    _epsilon(tpv_data, "surface", epsilon_override)))
    if _enabled(tpv_data, "energy"):
        results.append(compare_energies(case,
                                        _epsilon(tpv_data, "energy", epsilon_override)))
    if _enabled(tpv_data, "receiver"):
        results.append(compare_receivers(
            case, _epsilon(tpv_data, "receiver", epsilon_override),
        ))
    return results


# ---------------------------------------------------------------------------
# Recording (--record)
# ---------------------------------------------------------------------------


def record_case(case: CaseSpec) -> int:
    """Overwrite ``precomputed/<precision>/`` with the fresh output.

    The per-case ``tpv-data.json`` (sitting inside the precomputed dir) is
    preserved; everything else matching the output prefix is replaced.
    """
    if not case.workdir.is_dir():
        raise RuntimeError(f"Cannot record: {case.workdir} does not exist")
    case.precomputed_dir.mkdir(parents=True, exist_ok=True)

    preserved = case.precomputed_dir / "tpv-data.json"
    keep = preserved.read_bytes() if preserved.is_file() else None

    # Wipe and copy: deletes refs that no longer exist in the new output (e.g.
    # if you renamed a receiver). The preserved tpv-data.json is restored after.
    for f in case.precomputed_dir.iterdir():
        if f.is_file():
            f.unlink()

    n = 0
    for src in case.workdir.iterdir():
        if src.is_file() and src.name.startswith(OUTPUT_PREFIX):
            shutil.copy2(src, case.precomputed_dir / src.name)
            n += 1
    if keep is not None:
        preserved.write_bytes(keep)

    print(f"{c('[record]', _Ansi.MAGENTA)} {case.key}: wrote {n} files to {case.precomputed_dir}")
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
            "configs": ["elastic-p6-f32", "elastic-p6-f64"],
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


def write_report(results: list[CaseResult], report_path: Optional[Path]) -> None:
    # Compute plain and colored versions of each row in lockstep so the file
    # gets clean ASCII and stdout gets ANSI. Padding is done on plain text.
    CASE_W, STATE_W, TIME_W, CMP_W = 32, 6, 8, 10

    def state_color(state: str) -> str:
        if state == "PASS":
            return _ok(state.rjust(STATE_W))
        if state == "REC":
            return c(state.rjust(STATE_W), _Ansi.MAGENTA, _Ansi.BOLD)
        if state == "SKIP":
            return _warn(state.rjust(STATE_W))
        # FAIL or RC=N
        return _bad(state.rjust(STATE_W))

    def cmp_color(state: str, comp: str) -> str:
        padded = comp.rjust(CMP_W)
        if comp == "-":
            return _dim(padded)
        if state == "PASS":
            return _ok(padded)
        if state == "FAIL":
            return _bad(padded)
        return padded

    header = (f"{'case':<{CASE_W}} {'state':>{STATE_W}} "
              f"{'time':>{TIME_W}} {'compares':>{CMP_W}}  detail")
    sep = "-" * (CASE_W + 1 + STATE_W + 1 + TIME_W + 1 + CMP_W + 2 + 30)

    plain_rows = [header, sep]
    color_rows = [c(header, _Ansi.BOLD), sep]

    for r in results:
        if not r.ran:
            state, comp = "SKIP", "-"
        elif r.run_returncode != 0:
            state, comp = f"RC={r.run_returncode}", "-"
        elif r.recorded:
            state, comp = "REC", "-"
        else:
            ok = sum(c_.passed for c_ in r.compares)
            state = "PASS" if ok == len(r.compares) and r.compares else "FAIL"
            comp = f"{ok}/{len(r.compares)}"

        detail = ""
        if not r.recorded:
            fails = [c_ for c_ in r.compares if not c_.passed]
            if fails:
                detail = ", ".join(f"{c_.category} ε={c_.epsilon}" for c_ in fails)

        time_str = f"{r.duration_s:>{TIME_W - 1}.1f}s"
        plain_rows.append(
            f"{r.case.key:<{CASE_W}} {state:>{STATE_W}} {time_str} "
            f"{comp:>{CMP_W}}  {detail}"
        )
        color_rows.append(
            f"{r.case.key:<{CASE_W}} {state_color(state)} {time_str} "
            f"{cmp_color(state, comp)}  {_dim(detail) if detail else ''}"
        )

    n_pass = sum(r.passed for r in results)
    n_total = len(results)
    plain_total = f"Total: {n_pass}/{n_total} passed"
    if n_pass == n_total:
        color_total = _ok(plain_total)
    else:
        color_total = _bad(f"Total: {n_pass}/{n_total} passed ({n_total - n_pass} failed)")
        plain_total = f"Total: {n_pass}/{n_total} passed ({n_total - n_pass} failed)"

    plain_rows.append(sep)
    plain_rows.append(plain_total)
    color_rows.append(sep)
    color_rows.append(color_total)

    print("\n" + "\n".join(color_rows))
    if report_path:
        report_path.write_text("\n".join(plain_rows) + "\n")


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
    p.add_argument("--binary-dir", type=Path, required=True,
                   help="Directory holding seissol-{cpu,gpu}-<config> binaries.")
    p.add_argument("--variant", choices=("auto", "cpu", "gpu"), default="auto",
                   help="Binary variant to prefer: 'auto' tries cpu then gpu, "
                        "'cpu'/'gpu' restrict to that prefix. (default: auto)")
    sel = p.add_argument_group("Case selection")
    sel.add_argument("--case", action="append", default=[],
                     help="Case to run: 'tpv5' for both precisions, "
                          "'tpv5/double' for one. Repeatable.")
    sel.add_argument("--all", action="store_true",
                     help="Run every entry in cases.json.")

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
    out.add_argument("--color", choices=("auto", "never", "always"), default="auto",
                     help="Colorize stdout. 'auto' = TTY+!NO_COLOR. (default: auto)")
    out.add_argument("--fail-fast", action="store_true",
                     help="Abort after the first failing case.")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    setup_color(args.color)

    if not args.case and not args.all:
        print("error: pass --case ... or --all", file=sys.stderr)
        return 2
    if args.record and args.no_run:
        print("error: --record requires a fresh run; remove --no-run", file=sys.stderr)
        return 2

    precomputed_dir = args.precomputed_dir.resolve()
    cases = load_cases(precomputed_dir)
    selected = select_cases(cases, args.case, args.all)
    if not selected:
        print("No cases selected.", file=sys.stderr)
        return 2

    mpi_cmd = [args.mpi_cmd, *shlex.split(args.mpi_args)]
    results: list[CaseResult] = []

    for case in selected:
        print(f"\n{c('=== ' + case.key + '  config=' + case.config + ' ===', _Ansi.CYAN, _Ansi.BOLD)}")

        try:
            binary = resolve_binary(args.binary_dir, case.config, args.variant)
        except FileNotFoundError as e:
            print(f"{_warn('[skip]')} {case.key}: {e}", file=sys.stderr)
            results.append(CaseResult(case, False, -1, 0.0, []))
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

    write_report(results, args.report)
    return 0 if all(r.passed for r in results) else 1


if __name__ == "__main__":
    raise SystemExit(main())
