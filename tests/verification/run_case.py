#!/usr/bin/env python3
# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Run one verification case and judge the result.

A case is a base parameter file plus a few overrides, run on a generated mini
mesh in its own working directory. This wrapper exists because three things
cannot be expressed in ``add_test`` alone:

* the build may not contain the configuration a case needs. The case is still
  registered so that it stays visible, and is reported as skipped rather than
  as a failure or as nothing at all.
* the parameter file has no include or override mechanism, so the overrides
  have to be applied to a copy of it here.
* a run that finishes with return code zero has not necessarily computed
  anything. The checks below insist that the output is finite and that the
  case actually did what it claims to test.

Every run leaves a small JSON record behind, which the coverage command turns
into the answer to "was every declared case executed by some build".
"""

import argparse
import csv
import hashlib
import json
import math
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

# the return code CTest is told to read as "skipped"
SKIP_RETURN_CODE = 77

SECTION = re.compile(r"^\s*&(\w+)")
ASSIGNMENT = re.compile(r"^(\s*)([A-Za-z_]\w*)(\s*=\s*)(.*)$")
COMPARISON = re.compile(r"^([\w.-]+)\s*(>=|<=|!=|=)\s*(.*)$")


# --------------------------------------------------------------------------
# parameter files
# --------------------------------------------------------------------------


class ParameterFile:
    """A minimal editor for the ``&section key = value /`` parameter format.

    Comments and untouched lines are preserved verbatim, because the base file
    is meant to stay readable and reviewable. This is a stopgap: an override
    mechanism in SeisSol itself would remove the need to rewrite the file at
    all, and would also let a test state its configuration in one place
    instead of in a file plus a command line.
    """

    def __init__(self, text):
        # an ordered list of items, each either a raw line outside any section
        # or a [name, lines] section, so that rendering reproduces the input
        self.items = []
        current = None
        for line in text.splitlines():
            match = SECTION.match(line)
            if match is not None:
                current = [match.group(1), []]
                self.items.append(current)
            elif current is not None and line.strip() == "/":
                current = None
            elif current is not None:
                current[1].append(line)
            else:
                self.items.append(line)

    @property
    def sections(self):
        return [item for item in self.items if not isinstance(item, str)]

    def set(self, section, key, value):
        for name, lines in self.sections:
            if name.lower() != section.lower():
                continue
            for index, line in enumerate(lines):
                match = ASSIGNMENT.match(line)
                if match is not None and match.group(2).lower() == key.lower():
                    # keep the comment, it usually documents the units
                    comment = ""
                    position = match.group(4).find("!")
                    if position >= 0:
                        comment = " " + match.group(4)[position:]
                    lines[index] = (
                        match.group(1)
                        + match.group(2)
                        + match.group(3)
                        + value
                        + comment
                    )
                    return
            lines.append(f"{key} = {value}")
            return
        self.items.append([section, [f"{key} = {value}"]])

    def get(self, section, key):
        for name, lines in self.sections:
            if name.lower() != section.lower():
                continue
            for line in lines:
                match = ASSIGNMENT.match(line)
                if match is not None and match.group(2).lower() == key.lower():
                    value = match.group(4)
                    position = value.find("!")
                    if position >= 0:
                        value = value[:position]
                    return value.strip().strip("'\"")
        return None

    def render(self):
        parts = []
        for item in self.items:
            if isinstance(item, str):
                parts.append(item)
            else:
                parts.append(f"&{item[0]}")
                parts.extend(item[1])
                parts.append("/")
        return "\n".join(parts) + "\n"


def apply_overrides(parameters, overrides):
    for override in overrides:
        target, _, value = override.partition("=")
        if not value:
            raise SystemExit(
                f"malformed override {override!r}, expected Section.Key=Value"
            )
        section, _, key = target.rpartition(".")
        if not section:
            raise SystemExit(f"override {override!r} does not name a section")
        parameters.set(section, key, value)


# --------------------------------------------------------------------------
# capabilities
# --------------------------------------------------------------------------


def unmet_requirements(capabilities, requirements):
    """Return the requirements this build does not satisfy.

    A requirement is ``key=a|b``, ``key!=a``, ``key>=n`` or ``key<=n``. An
    unknown key is a failure of the test declaration, not of the build, so it
    raises instead of quietly skipping: a typo that skips everything would
    leave the suite green and empty.
    """
    unmet = []
    for requirement in requirements:
        match = COMPARISON.match(requirement)
        if match is None:
            raise SystemExit(f"malformed requirement {requirement!r}")
        key, operator, wanted = match.groups()
        if key not in capabilities:
            raise SystemExit(
                f"requirement {requirement!r} names an unknown capability; "
                f"known: {', '.join(sorted(capabilities))}"
            )
        have = capabilities[key]
        if operator in ("=", "!="):
            alternatives = [item.strip().lower() for item in wanted.split("|")]
            matched = str(have).lower() in alternatives
            if matched != (operator == "="):
                unmet.append(f"{key} is {have}, required {operator} {wanted}")
        else:
            try:
                left, right = float(have), float(wanted)
            except ValueError as error:
                raise SystemExit(
                    f"requirement {requirement!r} compares non-numbers"
                ) from error
            if (operator == ">=" and left < right) or (
                operator == "<=" and left > right
            ):
                unmet.append(f"{key} is {have}, required {operator} {wanted}")
    return unmet


# --------------------------------------------------------------------------
# output checks
# --------------------------------------------------------------------------


def read_receiver(path):
    """Read a receiver file into column names and rows of floats."""
    names = []
    rows = []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        stripped = line.strip()
        if not stripped:
            continue
        if stripped.upper().startswith("VARIABLES"):
            names = re.findall(r'"([^"]*)"', stripped)
        elif stripped.startswith("#") or stripped.upper().startswith("TITLE"):
            continue
        else:
            rows.append([float(item) for item in stripped.split()])
    return names, rows


def check_outputs(work, prefix, require_finite, activity, expect_zero=False):
    """Check that the run produced finite output, and that it did something.

    Both the volume receivers and the on-fault receivers are read. With
    ``expect_zero`` the demand is the opposite of activity: every value the
    volume receivers record, except the time, has to be exactly zero, as it
    must be for a zero initial state without sources. The on-fault receivers
    are left out of that demand, since friction coefficients and state
    variables are not zero on a fault at rest.
    """
    problems = []
    directory = work / Path(prefix).parent
    name = Path(prefix).name
    volume = sorted(directory.glob(f"{name}-receiver-*.dat"))
    fault = sorted(directory.glob(f"{name}-faultreceiver-*.dat"))
    if not volume:
        return ["no receiver output was written"]
    receivers = volume + fault

    extrema = {}
    for path in receivers:
        try:
            names, rows = read_receiver(path)
        except ValueError as error:
            problems.append(f"{path.name}: unreadable ({error})")
            continue
        if not rows:
            problems.append(f"{path.name}: no samples")
            continue
        if expect_zero and path in volume:
            nonzero = [
                (abs(value), column)
                for row in rows
                for column, value in enumerate(row[1:], start=1)
                if value != 0.0
            ]
            if nonzero:
                size, column = max(nonzero)
                name = names[column] if column < len(names) else str(column)
                problems.append(
                    f"{path.name}: {len(nonzero)} nonzero values where zero is "
                    f"required, largest |{name}| = {size:.3g}"
                )
        for row in rows:
            for column, value in enumerate(row):
                if require_finite and not math.isfinite(value):
                    problems.append(f"{path.name}: non-finite value in column {column}")
                    break
                name = names[column] if column < len(names) else str(column)
                extrema[name] = max(extrema.get(name, 0.0), abs(value))

    for demand in activity:
        # "SRs|SRd:1e-6" is satisfied by either column; which fault direction
        # carries the slip depends on the orientation convention, which the
        # check has no business assuming
        names, _, threshold = demand.partition(":")
        threshold = float(threshold)
        candidates = names.split("|")
        seen = [extrema[candidate] for candidate in candidates if candidate in extrema]
        if not seen:
            problems.append(
                f"activity check wants {names!r}, which the receivers do not "
                f"contain (have: {', '.join(sorted(extrema))})"
            )
        elif max(seen) <= threshold:
            # a case that runs cleanly but never excites anything tests nothing
            problems.append(
                f"max |{names}| is {max(seen):g}, expected more than {threshold:g}"
            )
    return problems


def read_analysis(path):
    """Read ``<prefix>-analysis.csv``, which SeisSol writes at the end of a run.

    The file holds one row per quantity and norm. SeisSol writes it through its
    CSV table, which quotes text, so the norm comes in quotes, and writes every
    number in its shortest exact form.

    A quantity the wave does not excite, e.g. a stress component off the plane
    of a plane wave along an axis, has no relative error: SeisSol divides by the
    zero norm of its analytical solution and writes inf or nan. Its relative
    norms are left out; its deviation from zero stays in the absolute ones.
    """
    norms = {}
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.reader(stream))
    for row in rows[1:]:
        fields = [field.strip() for field in row]
        if len(fields) != 3:
            continue
        variable, norm, value = fields
        norms.setdefault(norm, {})[variable] = float(value)
    for norm in [name for name in norms if name.endswith("_rel")]:
        absolute = norms.get(norm.removesuffix("_rel"), {})
        norms[norm] = {
            variable: value
            for variable, value in norms[norm].items()
            if math.isfinite(value)
            or not math.isfinite(absolute.get(variable, math.nan))
        }
    return norms


def check_analysis(work, prefix, thresholds_file, key, record):
    """Judge the analytical error, or record it when no threshold is set."""
    path = work / f"{prefix}-analysis.csv"
    if not path.exists():
        return [
            "no analysis output; the initial condition has no analytical solution"
        ], {}
    norms = read_analysis(path)
    observed = {norm: max(values.values()) for norm, values in norms.items() if values}

    table = json.loads(Path(thresholds_file).read_text(encoding="utf-8"))
    wanted = table.get("entries", {}).get(key)
    if wanted is None:
        print(
            f"no thresholds for {key}; observed "
            + ", ".join(f"{n}={v:.6g}" for n, v in sorted(observed.items()))
        )
        return [], observed

    problems = []
    for norm, limit in sorted(wanted.items()):
        seen = observed.get(norm)
        if seen is None:
            problems.append(f"{norm} is not in the analysis output")
        elif limit is None:
            # recorded, not judged: a threshold nobody measured is a threshold
            # that either passes everything or fails on the first compiler
            print(f"{key} {norm} = {seen:.6g} (no threshold set)")
        elif seen > limit:
            problems.append(f"{norm} is {seen:.6g}, above the threshold {limit:.6g}")
    if record:
        print("record:", json.dumps({key: observed}))
    return problems, observed


def compare_receivers(reference, candidate, prefix, tolerance, floor):
    """Compare two runs' receiver output.

    With a tolerance of zero the files are compared byte for byte, which is the
    only way to state bit identity about an ASCII output. A mismatch is then
    quantified anyway, because a difference of 1e-16 and one of 1e-3 call for
    completely different investigations.
    """
    problems = []
    worst = {}

    def receivers(root):
        directory = root / Path(prefix).parent
        name = Path(prefix).name
        found = {}
        for kind in ("receiver", "faultreceiver"):
            for path in directory.glob(f"{name}-{kind}-*.dat"):
                found[path.name] = path
        return found

    left, right = receivers(reference), receivers(candidate)
    if not left:
        return ["the reference run wrote no receiver output"], worst
    missing = sorted(set(left) - set(right))
    extra = sorted(set(right) - set(left))
    if missing:
        problems.append(f"missing in the candidate: {', '.join(missing)}")
    if extra:
        problems.append(f"only in the candidate: {', '.join(extra)}")

    identical = True
    for name in sorted(set(left) & set(right)):
        if left[name].read_bytes() == right[name].read_bytes():
            continue
        identical = False
        names, rows_left = read_receiver(left[name])
        _, rows_right = read_receiver(right[name])
        if len(rows_left) != len(rows_right):
            problems.append(
                f"{name}: {len(rows_left)} samples against {len(rows_right)}"
            )
            continue
        for row_left, row_right in zip(rows_left, rows_right):
            for column, (a, b) in enumerate(zip(row_left, row_right)):
                label = names[column] if column < len(names) else str(column)
                relative = abs(a - b) / max(abs(a), abs(b), floor)
                worst[label] = max(worst.get(label, 0.0), relative)

    if identical:
        print("receiver output is byte for byte identical")
    else:
        summary = ", ".join(
            f"{label}={value:.3g}" for label, value in sorted(worst.items())
        )
        print(f"largest relative differences: {summary}")
        exceeded = {label: value for label, value in worst.items() if value > tolerance}
        if exceeded:
            problems.append(
                "relative difference above the tolerance "
                f"{tolerance:.3g}: "
                + ", ".join(
                    f"{label}={value:.3g}" for label, value in sorted(exceeded.items())
                )
            )
    return problems, worst


# --------------------------------------------------------------------------
# running
# --------------------------------------------------------------------------


# SeisSol's own report of a parameter it does not understand. Turned into a
# failure here: a misspelt key otherwise runs with the default and the case
# verifies something other than what it claims.
UNKNOWN_PARAMETER = re.compile(
    r"The field\s+(\S+)\s+in\s+(\S*)\s*was given in the parameter file, "
    r"but is unknown to SeisSol"
)
TIMESTEP_LINE = re.compile(r"Minimum timestep[^:]*:\s*([-+0-9.eE]+)\s*(\S*?)s\b")
SI_PREFIXES = {
    "": 1.0,
    "m": 1e-3,
    "\u00b5": 1e-6,
    "u": 1e-6,
    "n": 1e-9,
    "p": 1e-12,
}


def check_log(log, timestep):
    """Check what only the log can tell: the parameters and the time step."""
    problems = []

    unknown = sorted(
        {
            f"{section}.{field}" if section else field
            for field, section in UNKNOWN_PARAMETER.findall(log)
        }
    )
    if unknown:
        problems.append(f"parameters unknown to SeisSol: {', '.join(unknown)}")

    if timestep is not None:
        matches = TIMESTEP_LINE.findall(log)
        if not matches:
            problems.append(
                "the log does not report the time step, so the requested width "
                "cannot be confirmed"
            )
        else:
            # the last report is the effective one; with a wiggle factor
            # SeisSol reports the step before and after applying it
            value, prefix = matches[-1]
            if prefix not in SI_PREFIXES:
                problems.append(f"cannot interpret the time step unit {prefix!r}s")
            else:
                effective = float(value) * SI_PREFIXES[prefix]
                # the value is printed with four decimals after the prefix,
                # which bounds how closely it can be compared
                if abs(effective - timestep) > 1e-5 * timestep:
                    problems.append(
                        f"the time step is {effective:.6g} s, not the requested "
                        f"{timestep:.6g} s: FixTimeStep only caps the step, and the "
                        "CFL limit of this mesh, material and order is below it"
                    )
    return problems


def record(path, payload):
    if path is None:
        return
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


def _command_run(args):
    if args.expect_zero and args.activity:
        raise SystemExit("--expect-zero and --activity contradict each other")
    capabilities = json.loads(Path(args.capabilities).read_text(encoding="utf-8"))
    result = {
        "name": args.name,
        "requirements": args.require,
        "capabilities": capabilities,
    }

    unmet = unmet_requirements(capabilities, args.require)
    if unmet:
        result["status"] = "skipped"
        result["reason"] = unmet
        record(args.result, result)
        print(f"skipped: {'; '.join(unmet)}")
        return SKIP_RETURN_CODE

    work = Path(args.work_dir)
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True)

    case_dir = Path(args.case_dir)
    parameters = ParameterFile((case_dir / args.parameters).read_text(encoding="utf-8"))

    # file names and overrides may depend on the build, e.g. the material file
    # differs per equation system; a file this suite does not yet provide for
    # the configuration is a gap in the suite, reported as a skip
    substitutions = {key: str(value) for key, value in capabilities.items()}
    files = [name.format(**substitutions) for name in args.file]
    overrides = [entry.format(**substitutions) for entry in args.set]
    missing = [name for name in files if not (case_dir / name).exists()]
    if missing:
        result["status"] = "skipped"
        result["reason"] = [f"the suite has no {name} yet" for name in missing]
        record(args.result, result)
        print("skipped: " + "; ".join(result["reason"]))
        return SKIP_RETURN_CODE
    for name in files:
        shutil.copy(case_dir / name, work / Path(name).name)

    # what the run was computed from, so that a snapshot comparison can tell a
    # change of the inputs from a change of the numerics. The absolute mesh
    # path is left out, since it differs between machines; the mesh is
    # identified by its own content hash instead.
    inputs = {
        "parameters": hashlib.sha256(
            (case_dir / args.parameters).read_bytes()
        ).hexdigest(),
        "overrides": overrides,
        "files": {
            name: hashlib.sha256((case_dir / name).read_bytes()).hexdigest()
            for name in files
        },
        "timestep": args.timestep,
        "steps": args.steps,
    }
    result["inputs-id"] = hashlib.sha256(
        json.dumps(inputs, sort_keys=True).encode("utf-8")
    ).hexdigest()[:16]
    if args.mesh_manifest:
        manifest = json.loads(Path(args.mesh_manifest).read_text(encoding="utf-8"))
        result["mesh-id"] = manifest["mesh-id"]

    # an absolute mesh path keeps the working directory free of copies
    parameters.set("MeshNml", "MeshFile", f"'{Path(args.mesh).resolve()}'")
    # the contiguous distribution from the file order is the only decomposition
    # that is reproducible across rank counts and partitioner versions
    parameters.set("MeshNml", "PartitioningLib", "'none'")
    if args.timestep is not None:
        # FixTimeStep caps the step; below the CFL limit of every cell it is
        # the step, which makes runs comparable across meshes, materials and
        # orders. check_log confirms from the log that the cap is the active
        # constraint. A power of two keeps N * dt, EndTime / dt and the
        # accumulated time exact, so SeisSol computes exactly N steps of
        # exactly dt: with a decimal width, the step count ceil(EndTime / dt)
        # comes out as N + 1 for about one N in fourteen.
        mantissa, _ = math.frexp(args.timestep)
        if args.timestep <= 0 or mantissa != 0.5:
            raise SystemExit(
                f"the time step {args.timestep!r} is not a power of two; "
                f"use e.g. {2.0 ** round(math.log2(abs(args.timestep) or 1.0))!r}"
            )
        parameters.set("Discretization", "FixTimeStep", repr(args.timestep))
        # one receiver sample per step, at exactly the step times
        parameters.set("Output", "pickdt", repr(args.timestep))
        if args.steps is not None:
            parameters.set("AbortCriteria", "EndTime", repr(args.steps * args.timestep))
    apply_overrides(parameters, overrides)

    prefix = parameters.get("Output", "OutputFile") or "output/out"
    (work / Path(prefix).parent).mkdir(parents=True, exist_ok=True)
    (work / "parameters.par").write_text(parameters.render(), encoding="utf-8")

    command = []
    if args.ranks > 1:
        if not args.mpiexec:
            raise SystemExit(
                "more than one rank requested but no launcher was configured"
            )
        command += [args.mpiexec, args.numproc_flag, str(args.ranks)]
        command += args.mpiexec_flags
    command += [args.binary, "parameters.par"]

    environment = dict(os.environ)
    environment["OMP_NUM_THREADS"] = str(args.threads)
    environment.setdefault("OMP_PLACES", "cores")
    environment.setdefault("OMP_PROC_BIND", "close")

    print("running:", " ".join(command))
    completed = subprocess.run(
        command,
        cwd=work,
        env=environment,
        check=False,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
    )
    (work / "seissol.log").write_bytes(completed.stdout)
    log = completed.stdout.decode("utf-8", errors="replace")
    sys.stdout.write(log)

    result["status"] = "run"
    result["command"] = command
    result["returncode"] = completed.returncode
    problems = []
    if args.expect_failure is not None:
        # a configuration that SeisSol refuses on purpose has to be refused
        # cleanly, with the reason in the log, rather than crash or run
        if completed.returncode == 0:
            problems.append("the run succeeded, but a refusal was expected")
        elif re.search(args.expect_failure, log) is None:
            problems.append(
                f"the run failed with exit code {completed.returncode}, but the "
                f"log does not contain {args.expect_failure!r}"
            )
        result["problems"] = problems
        record(args.result, result)
        for problem in problems:
            print(f"error: {problem}", file=sys.stderr)
        return 1 if problems else 0
    if completed.returncode != 0:
        problems.append(f"exit code {completed.returncode}")
    else:
        problems += check_log(log, args.timestep)
        problems += check_outputs(
            work, prefix, args.require_finite, args.activity, args.expect_zero
        )

    if not problems and args.thresholds:
        key = "{}/o{}/{}".format(
            capabilities["equations"], capabilities["order"], capabilities["precision"]
        )
        analysis_problems, observed = check_analysis(
            work, prefix, args.thresholds, key, args.record
        )
        problems += analysis_problems
        result["analysis"] = observed
        result["configuration"] = key

    result["problems"] = problems
    record(args.result, result)
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    return 1 if problems else 0


def _command_compare(args):
    """Compare two runs, propagating a skip when either of them was skipped."""
    result = {
        "name": args.name,
        "reference": args.reference,
        "candidate": args.candidate,
    }

    inputs = {}
    for role, path in (
        ("reference", args.reference_result),
        ("candidate", args.candidate_result),
    ):
        if path is None or not Path(path).exists():
            result["status"] = "skipped"
            result["reason"] = [f"the {role} run left no result behind"]
            record(args.result, result)
            print(f"skipped: the {role} run did not report")
            return SKIP_RETURN_CODE
        inputs[role] = json.loads(Path(path).read_text(encoding="utf-8"))

    skipped = [role for role, payload in inputs.items() if payload["status"] != "run"]
    if skipped:
        # a comparison whose inputs this build cannot produce is not a failure
        reason = [f"the {role} run was skipped" for role in skipped]
        result["status"] = "skipped"
        result["reason"] = reason
        record(args.result, result)
        print("skipped: " + "; ".join(reason))
        return SKIP_RETURN_CODE

    problems, worst = compare_receivers(
        Path(args.reference),
        Path(args.candidate),
        args.prefix,
        args.tolerance,
        args.floor,
    )
    result["status"] = "run"
    result["tolerance"] = args.tolerance
    result["worst-relative"] = worst
    result["problems"] = problems
    record(args.result, result)
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    return 1 if problems else 0


# --------------------------------------------------------------------------
# snapshots
# --------------------------------------------------------------------------


def _fingerprint_receiver(path):
    """Condense a receiver file into what a later run can be held against.

    The hash says whether anything changed at all; the statistics say by how
    much, and stay comparable across compilers, where the hash cannot. The
    norm is summed with fsum, which is exact and independent of order.
    """
    names, rows = read_receiver(path)
    count = len(rows)
    entry = {
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "samples": count,
        "end-time": rows[-1][0] if rows else None,
        "columns": {},
    }
    for column in range(1, len(rows[0]) if rows else 0):
        name = names[column] if column < len(names) else str(column)
        if name in entry["columns"]:
            name = f"{name}#{column}"
        values = [row[column] for row in rows]
        entry["columns"][name] = {
            "l2": math.sqrt(math.fsum(value * value for value in values)),
            "max": max(abs(value) for value in values),
            "final": values[-1],
            "quarters": [values[(k * (count - 1)) // 4] for k in (1, 2, 3)],
        }
    return entry


def _fingerprint_analysis(path):
    norms = read_analysis(path)
    return {
        "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "norms": {
            f"{norm}/{variable}": value
            for norm, values in norms.items()
            for variable, value in values.items()
        },
    }


def fingerprint(work, prefix):
    directory = work / Path(prefix).parent
    name = Path(prefix).name
    outputs = {}
    for kind in ("receiver", "faultreceiver"):
        for path in sorted(directory.glob(f"{name}-{kind}-*.dat")):
            outputs[path.name] = _fingerprint_receiver(path)
    analysis = directory / f"{name}-analysis.csv"
    if analysis.exists():
        outputs[analysis.name] = _fingerprint_analysis(analysis)
    return outputs


def _relative(a, b, scale):
    if a == b:
        return 0.0
    return abs(a - b) / scale if scale > 0 else math.inf


def compare_fingerprints(reference, current):
    """Compare two sets of fingerprints.

    Differences are measured relative to the amplitude of the signal they
    belong to, the column's maximum over time, so that a value near a zero
    crossing of a large signal does not count as a large relative change.
    Returns the structural problems, the largest deviation per quantity, and
    the files that are identical to the byte.
    """
    problems = []
    worst = {}
    identical = []
    for name in sorted(set(reference) - set(current)):
        problems.append(f"{name} is in the reference but was not written")
    for name in sorted(set(current) - set(reference)):
        problems.append(f"{name} was written but is not in the reference")

    for name in sorted(set(reference) & set(current)):
        old, new = reference[name], current[name]
        if old["sha256"] == new["sha256"]:
            identical.append(name)
            continue
        if "norms" in old:
            for norm in sorted(set(old["norms"]) | set(new["norms"])):
                if norm not in old["norms"] or norm not in new["norms"]:
                    problems.append(f"{name}: {norm} is only in one of the two")
                    continue
                a, b = old["norms"][norm], new["norms"][norm]
                label = f"analysis {norm.split('/')[0]}"
                deviation = _relative(a, b, max(abs(a), abs(b)))
                worst[label] = max(worst.get(label, 0.0), deviation)
            continue
        if (old["samples"], old["end-time"]) != (new["samples"], new["end-time"]):
            problems.append(
                f"{name}: {new['samples']} samples up to t = {new['end-time']}, "
                f"the reference has {old['samples']} up to t = {old['end-time']}"
            )
            continue
        for column in sorted(set(old["columns"]) | set(new["columns"])):
            if column not in old["columns"] or column not in new["columns"]:
                problems.append(f"{name}: column {column} is only in one of the two")
                continue
            a, b = old["columns"][column], new["columns"][column]
            scale = max(a["max"], b["max"])
            deviations = [_relative(a["l2"], b["l2"], max(a["l2"], b["l2"]))]
            deviations += [_relative(a[key], b[key], scale) for key in ("max", "final")]
            deviations += [
                _relative(x, y, scale) for x, y in zip(a["quarters"], b["quarters"])
            ]
            worst[column] = max(worst.get(column, 0.0), max(deviations))
    return problems, worst, identical


def source_revision(source_dir, ignore):
    """The commit of the source tree, marked dirty for uncommitted changes.

    Changes below ``ignore`` do not count: updating references writes there,
    and would otherwise mark every reference after the first as dirty.
    """
    if not source_dir:
        return "unknown"
    try:
        commit = subprocess.run(
            ["git", "-C", source_dir, "rev-parse", "--short=12", "HEAD"],
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
        pathspec = ["--", "."]
        relative = os.path.relpath(ignore, source_dir)
        if not relative.startswith(".."):
            pathspec.append(f":(exclude){relative}")
        status = subprocess.run(
            ["git", "-C", source_dir, "status", "--porcelain", "--untracked-files=no"]
            + pathspec,
            capture_output=True,
            text=True,
            check=True,
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        return "unknown"
    return commit + ("+dirty" if status else "")


def _describe(worst, limit=6):
    ranked = sorted(worst.items(), key=lambda item: -item[1])[:limit]
    return ", ".join(f"{name} {value:.3g}" for name, value in ranked)


def _command_snapshot(args):
    """Hold a run against its stored reference, or store it as the reference."""
    update = args.update or bool(os.environ.get("SEISSOL_UPDATE_SNAPSHOTS"))
    capabilities = json.loads(Path(args.capabilities).read_text(encoding="utf-8"))
    environment = json.loads(Path(args.environment).read_text(encoding="utf-8"))
    key = args.key.format(**{k: str(v) for k, v in capabilities.items()})
    reference_path = Path(args.reference_dir) / args.name / f"{key}.json"
    result = {"name": f"snapshot-{args.name}", "configuration": key}

    case_path = Path(args.case_result)
    if not case_path.exists():
        result["status"] = "skipped"
        result["reason"] = ["the case left no result behind"]
        record(args.result, result)
        print("skipped: the case did not report")
        return SKIP_RETURN_CODE
    if args.only_configurations and args.configuration not in args.only_configurations:
        # references are kept for a small set of configurations, so that one
        # build refreshes all of them after a deliberate change of the numerics
        result["status"] = "skipped"
        result["reason"] = [f"{args.configuration} keeps no snapshot references"]
        record(args.result, result)
        print(f"skipped: {args.configuration} keeps no snapshot references")
        return SKIP_RETURN_CODE
    case = json.loads(case_path.read_text(encoding="utf-8"))
    if case["status"] != "run":
        result["status"] = "skipped"
        result["reason"] = [f"the case was skipped: {'; '.join(case['reason'])}"]
        record(args.result, result)
        print(result["reason"][0])
        return SKIP_RETURN_CODE
    if case.get("problems"):
        # a failed run is not a snapshot of anything, and never a reference
        print(
            "error: the case itself failed; there is nothing to compare",
            file=sys.stderr,
        )
        return 1

    environment["commit"] = source_revision(args.source_dir, args.reference_dir)
    current = {
        "case": args.name,
        "configuration": key,
        "mesh-id": case.get("mesh-id"),
        "inputs-id": case.get("inputs-id"),
        "environment": environment,
        "outputs": fingerprint(Path(args.work), args.prefix),
    }

    reference = None
    if reference_path.exists():
        reference = json.loads(reference_path.read_text(encoding="utf-8"))

    if update:
        if reference is None:
            print(f"new reference {reference_path}")
        else:
            problems, worst, identical = compare_fingerprints(
                reference["outputs"], current["outputs"]
            )
            if not problems and not worst:
                print(f"reference unchanged: {reference_path}")
                return 0
            changed = []
            for field in ("mesh-id", "inputs-id"):
                if reference.get(field) != current[field]:
                    changed.append(field.split("-")[0])
            print(f"replacing reference {reference_path}")
            print(f"  it was computed at {reference['environment'].get('commit')}")
            if changed:
                print(f"  the {' and '.join(changed)} changed since")
            for problem in problems:
                print(f"  {problem}")
            if worst:
                print(f"  largest deviations accepted: {_describe(worst)}")
        reference_path.parent.mkdir(parents=True, exist_ok=True)
        reference_path.write_text(
            json.dumps(current, indent=1, sort_keys=True) + "\n", encoding="utf-8"
        )
        return 0

    if reference is None:
        result["status"] = "skipped"
        result["reason"] = [f"no reference for configuration {key}"]
        record(args.result, result)
        print(
            f"skipped: no reference for configuration {key}; to create it, run "
            "with SEISSOL_UPDATE_SNAPSHOTS=1 and commit the result deliberately"
        )
        return SKIP_RETURN_CODE

    result["status"] = "run"
    problems = []
    # a reference computed from other inputs cannot say anything about the
    # numerics; this is reported as what it is instead of as a deviation
    if reference.get("mesh-id") != current["mesh-id"]:
        problems.append(
            f"the mesh changed (reference {reference.get('mesh-id')}, now "
            f"{current['mesh-id']}); the reference is stale and has to be "
            "regenerated deliberately"
        )
    if reference.get("inputs-id") != current["inputs-id"]:
        problems.append(
            "the parameters, case files, step width or step count changed since "
            "the reference was computed; it is stale and has to be regenerated "
            "deliberately"
        )
    if problems:
        result["problems"] = problems
        record(args.result, result)
        for problem in problems:
            print(f"error: {problem}", file=sys.stderr)
        return 1

    structure, worst, identical = compare_fingerprints(
        reference["outputs"], current["outputs"]
    )
    tolerance = (
        args.tolerance_single
        if capabilities.get("precision") == "single"
        else args.tolerance_double
    )
    problems += structure
    exceeded = {name: value for name, value in worst.items() if value > tolerance}
    if exceeded:
        problems.append(
            f"deviation from the reference above the tolerance {tolerance:.3g}: "
            + _describe(exceeded)
        )

    reference_environment = dict(reference.get("environment", {}))
    reference_environment.pop("commit", None)
    here = dict(environment)
    here.pop("commit", None)
    if exceeded and reference_environment != here:
        differences = sorted(
            key
            for key in set(reference_environment) | set(here)
            if reference_environment.get(key) != here.get(key)
        )
        problems.append(
            "the reference comes from another environment ("
            + ", ".join(
                f"{key}: {reference_environment.get(key)} -> {here.get(key)}"
                for key in differences
            )
            + "); bit identity cannot be expected there, and the tolerance for "
            "this environment is to be set from the deviation measured here"
        )

    total = len(current["outputs"])
    if len(identical) == total and not structure:
        print(f"identical to the reference in all {total} outputs")
    elif not problems:
        print(
            f"{len(identical)} of {total} outputs identical, the rest within "
            f"tolerance: {_describe(worst)}"
        )
    result["worst"] = worst
    result["identical"] = len(identical)
    result["problems"] = problems
    record(args.result, result)
    for problem in problems:
        print(f"error: {problem}", file=sys.stderr)
    return 1 if problems else 0


def _command_coverage(args):
    """Check that every declared case was run somewhere, not just skipped."""
    declared = json.loads(Path(args.declared).read_text(encoding="utf-8"))
    results = {}
    for path in sorted(Path(args.results).glob("*.json")):
        payload = json.loads(path.read_text(encoding="utf-8"))
        results[payload["name"]] = payload

    ran = [name for name, payload in results.items() if payload["status"] == "run"]
    skipped = {
        name: payload["reason"]
        for name, payload in results.items()
        if payload["status"] == "skipped"
    }
    absent = [case["name"] for case in declared["cases"] if case["name"] not in results]

    print(f"declared {len(declared['cases'])}, ran {len(ran)}, skipped {len(skipped)}")
    for name, reason in sorted(skipped.items()):
        print(f"  skipped {name}: {'; '.join(reason)}")
    for name in absent:
        print(f"  no result {name}")

    if args.history:
        history = Path(args.history)
        seen = (
            json.loads(history.read_text(encoding="utf-8")) if history.exists() else {}
        )
        for name in ran:
            seen[name] = True
        history.write_text(json.dumps(seen, indent=2, sort_keys=True) + "\n", "utf-8")
        never = [
            case["name"] for case in declared["cases"] if not seen.get(case["name"])
        ]
        if never:
            # a case no build has ever run is a hole in the matrix, and it looks
            # exactly like a passing suite unless it is called out
            print("never run by any build:", file=sys.stderr)
            for name in sorted(never):
                print(f"  {name}", file=sys.stderr)
            return 1

    if absent and args.require_results:
        return 1
    return 0


def main(argv=None):
    parser = argparse.ArgumentParser(description="run a SeisSol verification case")
    commands = parser.add_subparsers(dest="command", required=True)

    run = commands.add_parser("run", help="run a single case")
    run.add_argument("--name", required=True)
    run.add_argument("--binary", required=True)
    run.add_argument("--capabilities", required=True)
    run.add_argument("--case-dir", required=True)
    run.add_argument("--work-dir", required=True)
    run.add_argument("--mesh", required=True)
    run.add_argument("--mesh-manifest", help="the generator's manifest of the mesh")
    run.add_argument("--parameters", required=True)
    run.add_argument("--result")
    run.add_argument("--ranks", type=int, default=1)
    run.add_argument("--threads", type=int, default=1)
    run.add_argument("--steps", type=int)
    run.add_argument("--timestep", type=float)
    run.add_argument("--mpiexec", default="")
    # these are passed as --flag=value because their values start with a dash,
    # which argparse would otherwise read as the next option
    run.add_argument("--numproc-flag", default="-n")
    run.add_argument(
        "--mpiexec-flag", action="append", default=[], dest="mpiexec_flags"
    )
    run.add_argument("--set", action="append", default=[], metavar="SECTION.KEY=VALUE")
    run.add_argument("--file", action="append", default=[], metavar="NAME")
    run.add_argument("--require", action="append", default=[], metavar="KEY=VALUE")
    run.add_argument(
        "--activity",
        action="append",
        default=[],
        metavar="COLUMN:THRESHOLD",
        help="insist that the run actually excited something",
    )
    run.add_argument(
        "--thresholds", help="judge the analytical error against this table"
    )
    run.add_argument(
        "--record",
        action="store_true",
        default=bool(os.environ.get("SEISSOL_RECORD_THRESHOLDS")),
        help="print the observed errors in a form that can be pasted into the table",
    )
    run.add_argument(
        "--expect-failure",
        metavar="REGEX",
        help="require the run to fail with a log line matching REGEX",
    )
    run.add_argument(
        "--expect-zero",
        action="store_true",
        help="require every recorded value to be exactly zero",
    )
    run.add_argument("--no-finite-check", dest="require_finite", action="store_false")
    run.set_defaults(require_finite=True, func=_command_run)

    compare = commands.add_parser("compare", help="compare two runs against each other")
    compare.add_argument("--name", required=True)
    compare.add_argument(
        "--reference", required=True, help="working directory of the reference"
    )
    compare.add_argument(
        "--candidate", required=True, help="working directory of the candidate"
    )
    compare.add_argument("--reference-result", required=True)
    compare.add_argument("--candidate-result", required=True)
    compare.add_argument("--prefix", default="output/mini")
    compare.add_argument(
        "--tolerance",
        type=float,
        default=0.0,
        help="largest permitted relative difference; zero demands identical files",
    )
    compare.add_argument(
        "--floor",
        type=float,
        default=1e-30,
        help="denominator floor, so that two near-zero values do not look far apart",
    )
    compare.add_argument("--result")
    compare.set_defaults(func=_command_compare)

    snapshot = commands.add_parser(
        "snapshot", help="compare a run with its stored reference, or store it"
    )
    snapshot.add_argument("--name", required=True)
    snapshot.add_argument("--case-result", required=True)
    snapshot.add_argument("--work", required=True)
    snapshot.add_argument("--prefix", default="output/mini")
    snapshot.add_argument("--reference-dir", required=True)
    snapshot.add_argument("--key", required=True, help="configuration key template")
    snapshot.add_argument("--capabilities", required=True)
    snapshot.add_argument("--environment", required=True)
    snapshot.add_argument("--source-dir", help="for recording the commit")
    snapshot.add_argument("--configuration", default="", help="what this build is")
    snapshot.add_argument(
        "--only-configurations",
        nargs="*",
        default=[],
        help="build configurations that keep references; empty means all of them",
    )
    snapshot.add_argument("--tolerance-double", type=float, default=0.0)
    snapshot.add_argument("--tolerance-single", type=float, default=0.0)
    snapshot.add_argument(
        "--update",
        action="store_true",
        help="store the run as the reference; SEISSOL_UPDATE_SNAPSHOTS=1 does the same",
    )
    snapshot.add_argument("--result")
    snapshot.set_defaults(func=_command_snapshot)

    coverage = commands.add_parser(
        "coverage", help="account for declared versus run cases"
    )
    coverage.add_argument("--declared", required=True)
    coverage.add_argument("--results", required=True)
    coverage.add_argument(
        "--history", help="accumulate across builds and fail on holes"
    )
    coverage.add_argument("--require-results", action="store_true")
    coverage.set_defaults(func=_command_coverage)

    args = parser.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
