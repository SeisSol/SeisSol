#!/usr/bin/env python3

# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Compare a series of runs of one dynamic rupture benchmark with each other, with a
reference run and with the solutions of other codes from the SCEC/USGS code
verification website (CVWS).

The runs come in a manifest (YAML or JSON), in the order of the series, e.g. one run per
SeisSol version:

    scenario: tpv104              # sets the region of the rupture time comparison
    reference: head-stroud        # default: the last run of the series "main"
    end_time: 12.0                # optional: compare up to this time only
    scec: scec/tpv104             # optional: the uploads of other codes, see below
    runs:
      - id: v0.9.0
        output: runs/v0.9.0/output/tpv     # the OutputFile of the run
        commit: eadae61c2                 # optional, shown in the report
        date: 2022-08-29                  # optional, shown in the report
      - id: head-stroud
        output: runs/head-stroud/output/tpv
      - id: head-dunavant
        output: runs/head-dunavant/output/tpv
        series: dunavant                  # a run outside the main series

Paths are relative to the manifest. Each run is compared with the run before it in its
series (the jump), with the reference and with every SCEC upload. The script reads what
SeisSol writes, in the formats of all versions since 2021: the fault output (Xdmf), the
on-fault and off-fault receivers (text files or HDF5) and the energy output (CSV).

The SCEC directory holds the uploads in the CVWS file format, one subdirectory per code:
the on-fault stations (faultst<strike>dp<depth>), the off-fault stations
(body<normal>st<strike>dp<depth>), each named so in the file name or described in a
"# location" line, and the rupture contour data, in a file whose name contains "cplot"
(columns: along strike [m], down dip [m], rupture time [s]).

The metrics are modelled on those of the CVWS (Barall & Harris, SRL 86(1), 2015), which
normalize the RMS difference of two time series to 0..200, after the relative time shift
that fits best, and take the RMS difference of the rupture times in ms:
- two time series: the RMS of their difference over the mean of their RMS, in percent
  (0: identical; 200: opposite, or one of them zero). Against an upload of another code,
  after the relative time shift that fits best ("nrms_shifted", with the shift as
  "tshift"); between two runs of SeisSol also without a shift ("nrms"), since there a
  shift is a change in itself;
- rupture time: the RMS difference in ms over the area that ruptured in both.
In addition: the bias (the mean difference, which a difference at the scale of the
elements leaves alone), the largest difference and the area mismatch of the rupture
time, the final slip and the peak slip rate on the fault, and the seismic moment.

A run may name further runs to compare with ("compare_with: [id, ...]"), e.g. the same
commit in single precision against double precision.

The signs of the on-fault and off-fault quantities of other codes need not be those of
SeisSol: per quantity, the script takes the sign under which the reference run agrees
better with the uploads (the report states it); "scec_signs" in the manifest fixes it.

Writes metrics.csv (every value), summary.json, report.md and the plots to --out.
"""

from __future__ import annotations

import argparse
import glob
import importlib.util
import json
import math
import os
import re
import sys
from dataclasses import dataclass, field

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)


def _load_sibling(name: str, filename: str):
    """A script next to this one whose file name is no module name."""
    spec = importlib.util.spec_from_file_location(name, os.path.join(HERE, filename))
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


compare_receivers = _load_sibling("compare_receivers", "compare-receivers.py")

# The region of the rupture time comparison, along strike and down dip in m. For
# TPV101-105, the velocity-weakening patch, on which the CVWS asks for rupture times;
# for TPV5, the fault.
SCENARIOS = {
    name: (-15e3, 15e3, 0.0, 15e3)
    for name in ("tpv5", "tpv101", "tpv102", "tpv103", "tpv104", "tpv105")
}

# the CVWS fields of an on-fault station: the SeisSol quantity and its factor to CVWS
# units (MPa for stresses and pressure)
FAULT_FIELDS = {
    "h-slip": ("Sls", 1.0),
    "v-slip": ("Sld", 1.0),
    "h-slip-rate": ("SRs", 1.0),
    "v-slip-rate": ("SRd", 1.0),
    "h-shear-stress": ("Ts0", 1e-6),
    "v-shear-stress": ("Td0", 1e-6),
    "n-stress": ("Pn0", 1e-6),
    "psi": ("StV", 1.0),
    "pressure": ("P_f", 1e-6),
    "temperature": ("Tmp", 1.0),
}
BODY_FIELDS = ("h-disp", "h-vel", "v-disp", "v-vel", "n-disp", "n-vel")
# the vectors the CVWS compares as a whole
VECTORS = {
    "2d-slip": ("h-slip", "v-slip"),
    "2d-rate": ("h-slip-rate", "v-slip-rate"),
    "2d-stress": ("h-shear-stress", "v-shear-stress"),
    "vel": ("h-vel", "v-vel", "n-vel"),
    "moment-rate": ("moment-rate",),
}
QUANTITIES = {
    "fault": (
        "2d-slip",
        "2d-rate",
        "2d-stress",
        "n-stress",
        "psi",
        "pressure",
        "temperature",
    ),
    "body": ("h-vel", "v-vel", "n-vel", "vel", "h-disp", "v-disp", "n-disp"),
    "moment": ("moment-rate",),
}
# the quantity whose fit determines the time shift of a station
PRIMARY = {"fault": "2d-rate", "body": "vel", "moment": "moment-rate"}

# a rupture time at or above this marks a point that did not rupture
NO_RUPTURE = 1e8


# --- geometry -------------------------------------------------------------------------


def parse_axis(spec: str) -> np.ndarray:
    """The unit vector of an axis given as "+x", "-z", "y"."""
    match = re.fullmatch(r"\s*([+-]?)\s*([xyz])\s*", str(spec))
    if match is None:
        raise ValueError(f"not an axis: {spec!r} (expected e.g. +x, -z)")
    vector = np.zeros(3)
    vector["xyz".index(match.group(2))] = -1.0 if match.group(1) == "-" else 1.0
    return vector


@dataclass
class Axes:
    """The directions of the benchmark frame in the frame of the mesh."""

    strike: np.ndarray
    normal: np.ndarray
    down: np.ndarray

    @classmethod
    def from_spec(cls, spec: dict | None) -> "Axes":
        spec = spec or {}
        return cls(
            parse_axis(spec.get("strike", "+x")),
            parse_axis(spec.get("normal", "+y")),
            parse_axis(spec.get("down", "-z")),
        )

    def on_fault(self, points: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        """Along strike and down dip, in m."""
        points = np.atleast_2d(np.asarray(points, dtype=float))
        return points @ self.strike, points @ self.down

    def off_fault(self, points: np.ndarray) -> tuple[np.ndarray, ...]:
        """Normal to the fault, along strike and down dip, in m."""
        points = np.atleast_2d(np.asarray(points, dtype=float))
        return points @ self.normal, points @ self.strike, points @ self.down


def _hectometres(value: float) -> str:
    count = int(round(value / 100.0))
    return f"-{-count:03d}" if count < 0 else f"{count:03d}"


def fault_station(strike: float, depth: float) -> str:
    """The CVWS name of an on-fault station, from its position in m."""
    return f"faultst{_hectometres(strike)}dp{_hectometres(depth)}"


def body_station(normal: float, strike: float, depth: float) -> str:
    """The CVWS name of an off-fault station, from its position in m."""
    return f"body{_hectometres(normal)}st{_hectometres(strike)}dp{_hectometres(depth)}"


_BODY_NAME = re.compile(r"body(-?\d+)st(-?\d+)dp(-?\d+)")
_FAULT_NAME = re.compile(r"faultst(-?\d+)dp(-?\d+)")


def station_from_name(text: str) -> tuple[str, str] | None:
    """The kind ("fault" or "body") and canonical name of a station named in a text."""
    match = _BODY_NAME.search(text)
    if match:
        return "body", body_station(*(100.0 * int(g) for g in match.groups()))
    match = _FAULT_NAME.search(text)
    if match:
        return "fault", fault_station(*(100.0 * int(g) for g in match.groups()))
    return None


def station_from_location(header: list[str]) -> tuple[str, str] | None:
    """The station a "# location" line of a CVWS file describes, if it is readable."""
    text = " ".join(line.lstrip("#") for line in header if "location" in line.lower())
    if not text:
        return None

    def kilometres(pattern):
        match = re.search(
            r"(-?\d+(?:\.\d*)?)\s*km\s*(?:" + pattern + ")", text, re.IGNORECASE
        )
        return 1000.0 * float(match.group(1)) if match else None

    strike = kilometres(r"along[- ]?strike|strike")
    depth = kilometres(r"down[- ]?dip|depth|deep|dip")
    normal = kilometres(r"off[- ]?fault|normal|from (?:the )?fault")
    if strike is None:
        return None
    if normal is None or re.search(r"\bon[- ]fault\b", text, re.IGNORECASE):
        return ("fault", fault_station(strike, depth)) if depth is not None else None
    return "body", body_station(normal, strike, depth or 0.0)


# --- what a solution provides ---------------------------------------------------------


@dataclass
class Station:
    """The time series of one station, in CVWS fields and units."""

    name: str
    time: np.ndarray
    fields: dict


@dataclass
class Solution:
    """What a run of SeisSol, or the upload of another code, provides."""

    id: str
    label: str
    series: str = "main"
    meta: dict = field(default_factory=dict)
    scec: bool = False
    # time up to which the fault data holds; inf for an upload
    fault_time: float = math.inf
    # rupture time at points of the fault: along strike, down dip, time
    rupture: tuple | None = None
    # values at the cells of the fault output: along strike, down dip, {name: values}
    maps: tuple | None = None
    stations: dict = field(default_factory=lambda: {"fault": {}, "body": {}})
    moment: Station | None = None
    warnings: list = field(default_factory=list)
    # ids of further runs to compare this one with
    compare_with: list = field(default_factory=list)

    @property
    def has_data(self) -> bool:
        return (
            self.rupture is not None
            or self.maps is not None
            or any(self.stations.values())
            or self.moment is not None
        )


def _increasing(time: np.ndarray) -> np.ndarray:
    """The indices of the samples in increasing time, finite and without repetition."""
    time = np.asarray(time, dtype=float)
    finite = np.flatnonzero(np.isfinite(time))
    _, first = np.unique(time[finite], return_index=True)
    return finite[first]


def read_text_coordinates(path: str) -> tuple | None:
    """The position a text receiver file states in its "# x1", "# x2", "# x3" lines.

    From 16 to 25 June 2025, the fault receivers put the numbers of the receiver and
    the simulation before them ("# 1,1 x1").
    """
    point = {}
    with open(path) as handle:
        for number, line in enumerate(handle):
            match = re.match(r"\s*#\s*(?:\d+(?:,\d+)*\s+)?x([123])\s+(\S+)", line)
            if match:
                point[int(match.group(1))] = float(
                    match.group(2).replace("D", "E").replace("d", "e")
                )
            if len(point) == 3 or number > 32:
                break
    return tuple(point[i] for i in (1, 2, 3)) if len(point) == 3 else None


_HDF5_LAYOUT = {
    "receiver": ("receivers", "PointId"),
    "faultreceiver": ("faultreceivers", "ReceiverId"),
}


def receiver_coordinates(directory: str, prefix: str, file_type: str) -> dict:
    """The position of every receiver of one kind, by the number of its text file."""
    path = compare_receivers.hdf5_receiver_file(directory, prefix, file_type)
    if path is not None:
        import h5py

        group_name, id_name = _HDF5_LAYOUT[file_type]
        with h5py.File(path, "r") as handle:
            if group_name in handle:
                group = handle[group_name]
                numbers = group[id_name][:] + 1
                points = group["Coordinates"][:]
            else:
                # the off-fault receivers of April to September 2026, which list all
                # points in the order of the input (see compare_receivers)
                points = handle["Coordinates"][:]
                numbers = np.arange(1, len(points) + 1)
        coordinates = {}
        for number, point in zip(numbers, points):
            coordinates.setdefault(int(number), tuple(float(x) for x in point))
        return coordinates
    coordinates = {}
    for number in compare_receivers.find_all_receivers(directory, prefix, file_type):
        files = sorted(
            glob.glob(
                os.path.join(directory, f"{prefix}-{file_type}-{number:05d}*.dat")
            )
        )
        point = read_text_coordinates(files[0]) if files else None
        if point is not None:
            coordinates[int(number)] = point
    return coordinates


def _unique_name(name: str, taken: dict) -> str:
    candidate, count = name, 1
    while candidate in taken:
        count += 1
        candidate = f"{name}#{count}"
    return candidate


def read_fault_receivers(directory: str, prefix: str, axes: Axes) -> dict:
    receivers = compare_receivers.load_receivers(directory, prefix, "faultreceiver")
    coordinates = receiver_coordinates(directory, prefix, "faultreceiver")
    stations = {}
    for number, frame in sorted(receivers.items()):
        if number not in coordinates or len(frame) < 2:
            continue
        # the Fortran writer of the receivers named the absolute slip ASI
        frame = frame.rename(columns={"ASI": "ASl"})
        strike, depth = axes.on_fault(coordinates[number])
        name = _unique_name(fault_station(strike[0], depth[0]), stations)
        keep = _increasing(frame["Time"].to_numpy())
        fields = {
            key: frame[quantity].to_numpy(dtype=float)[keep] * factor
            for key, (quantity, factor) in FAULT_FIELDS.items()
            if quantity in frame.columns
        }
        stations[name] = Station(
            name, frame["Time"].to_numpy(dtype=float)[keep], fields
        )
    return stations


def read_receivers(directory: str, prefix: str, axes: Axes) -> dict:
    from scipy.integrate import cumulative_trapezoid

    receivers = compare_receivers.load_receivers(directory, prefix, "receiver")
    coordinates = receiver_coordinates(directory, prefix, "receiver")
    stations = {}
    for number, frame in sorted(receivers.items()):
        if number not in coordinates or len(frame) < 2:
            continue
        if not {"v1", "v2", "v3"} <= set(frame.columns):
            continue
        normal, strike, depth = axes.off_fault(coordinates[number])
        name = _unique_name(body_station(normal[0], strike[0], depth[0]), stations)
        keep = _increasing(frame["Time"].to_numpy())
        time = frame["Time"].to_numpy(dtype=float)[keep]
        velocity = frame[["v1", "v2", "v3"]].to_numpy(dtype=float)[keep]
        fields = {
            "h-vel": velocity @ axes.strike,
            "n-vel": velocity @ axes.normal,
            "v-vel": -(velocity @ axes.down),
        }
        for component in ("h", "n", "v"):
            fields[f"{component}-disp"] = cumulative_trapezoid(
                fields[f"{component}-vel"], time, initial=0.0
            )
        stations[name] = Station(name, time, fields)
    return stations


def energy_table(frame: pd.DataFrame) -> pd.DataFrame:
    """The energy output with a column per quantity, at the times it was written."""
    if "variable" in frame.columns:
        index = ["time"]
        if "simulation_index" in frame.columns:
            index.append("simulation_index")
        frame = frame.pivot_table(
            index=index, columns="variable", values="measurement", aggfunc="first"
        ).reset_index()
    if "simulation_index" in frame.columns:
        frame = frame[frame["simulation_index"] == frame["simulation_index"].min()]
    return frame


def read_energy_csv(path: str) -> pd.DataFrame:
    """The energy output as written.

    From 12 to 25 March 2025, SeisSol wrote the momentum rows with an empty field
    too many ("0.0,momentumX,,1.5"); such a row is read without its empty fields.
    """
    with open(path) as handle:
        width = len(handle.readline().split(","))

    def repair(fields: list[str]) -> list[str] | None:
        kept = [field for field in fields if field.strip()]
        return kept if len(kept) == width else None

    return pd.read_csv(path, engine="python", on_bad_lines=repair)


def read_moment(directory: str, prefix: str) -> Station | None:
    """The moment rate, from the seismic moment of the energy output."""
    path = os.path.join(directory, f"{prefix}-energy.csv")
    if not os.path.isfile(path):
        return None
    frame = energy_table(read_energy_csv(path))
    if "seismic_moment" not in frame.columns:
        return None
    time = frame["time"].to_numpy(dtype=float)
    moment = frame["seismic_moment"].to_numpy(dtype=float)
    keep = _increasing(time)
    time, moment = time[keep], moment[keep]
    rate = np.gradient(moment, time) if len(time) > 1 else np.zeros_like(moment)
    return Station("moment", time, {"moment": moment, "moment-rate": rate})


def read_fault_output(path: str, axes: Axes) -> tuple[float, tuple]:
    """The last snapshot of the fault output: its time and the values at the cells."""
    import seissolxdmf

    xdmf = seissolxdmf.seissolxdmf(path)
    available = set(xdmf.ReadAvailableDataFields())
    step = xdmf.ReadNdt() - 1
    try:
        time = float(np.asarray(xdmf.ReadTimes(), dtype=float)[step])
    except Exception:
        time = math.inf
    geometry = np.asarray(xdmf.ReadGeometry(), dtype=float)
    connect = np.asarray(xdmf.ReadConnect(), dtype=np.int64)
    strike, depth = axes.on_fault(geometry[connect].mean(axis=1))

    def read(name):
        return np.asarray(xdmf.ReadData(name, step), dtype=float).ravel()

    values = {}
    if "RT" in available:
        values["rt"] = read("RT")
    if "ASl" in available:
        values["slip"] = read("ASl")
    elif {"Sls", "Sld"} <= available:
        values["slip"] = np.hypot(read("Sls"), read("Sld"))
    if "PSR" in available:
        values["psr"] = read("PSR")
    return time, (strike, depth, values)


def split_output(output: str, prefix: str | None, base: str) -> tuple[str, str]:
    """The directory and prefix of an OutputFile given as a prefix or as a directory."""
    path = output if os.path.isabs(output) else os.path.join(base, output)
    if os.path.isdir(path):
        return path, prefix or "tpv"
    return os.path.dirname(path), prefix or os.path.basename(path)


def load_run(entry: dict, axes: Axes, base: str) -> Solution:
    run_id = str(entry["id"])
    compare_with = entry.get("compare_with") or []
    solution = Solution(
        id=run_id,
        label=str(entry.get("label", run_id)),
        series=str(entry.get("series", "main")),
        meta={
            key: str(entry[key])
            for key in ("commit", "date", "note", "config")
            if key in entry
        },
        compare_with=[str(other) for other in np.atleast_1d(compare_with)],
    )
    directory, prefix = split_output(str(entry["output"]), entry.get("prefix"), base)
    if not os.path.isdir(directory):
        solution.warnings.append(f"no output directory {directory}")
        return solution

    fault_xdmf = os.path.join(directory, f"{prefix}-fault.xdmf")
    if os.path.isfile(fault_xdmf):
        try:
            solution.fault_time, solution.maps = read_fault_output(fault_xdmf, axes)
            strike, depth, values = solution.maps
            if "rt" in values:
                solution.rupture = (strike, depth, values["rt"])
            else:
                solution.warnings.append("the fault output holds no rupture time (RT)")
        except Exception as error:
            solution.warnings.append(f"cannot read {fault_xdmf}: {error}")
    else:
        solution.warnings.append(f"no fault output {fault_xdmf}")

    readers = (
        ("fault", read_fault_receivers, "on-fault receivers"),
        ("body", read_receivers, "off-fault receivers"),
    )
    for kind, reader, what in readers:
        try:
            solution.stations[kind] = reader(directory, prefix, axes)
        except Exception as error:
            solution.warnings.append(f"cannot read the {what}: {error}")
        if not solution.stations[kind]:
            solution.warnings.append(f"no {what}")
    try:
        solution.moment = read_moment(directory, prefix)
    except Exception as error:
        solution.warnings.append(f"cannot read the energy output: {error}")
    return solution


# --- the uploads to the CVWS ----------------------------------------------------------


def _number(token: str) -> float:
    return float(token.replace("D", "E").replace("d", "e"))


def read_cvws_table(path: str) -> tuple[list, np.ndarray, list]:
    """The field names, the rows and the comment lines of a file in the CVWS format.

    Lines starting with "#" are comments; the first other line that is not all numbers
    names the fields; the lines of numbers are the rows.
    """
    names, rows, header = None, [], []
    with open(path, errors="replace") as handle:
        for line in handle:
            stripped = line.strip()
            if not stripped:
                continue
            if stripped.startswith("#"):
                header.append(stripped)
                continue
            tokens = stripped.replace(",", " ").split()
            try:
                rows.append([_number(token) for token in tokens])
            except ValueError:
                if names is None and not rows:
                    names = [token.lower() for token in tokens]
    width = len(names) if names else max((len(row) for row in rows), default=0)
    rows = [row for row in rows if len(row) == width]
    return names or [], np.array(rows, dtype=float).reshape(-1, width), header


def _upload_code(path: str, root: str, token: str) -> str:
    relative = os.path.relpath(path, root)
    parts = relative.split(os.sep)
    if len(parts) > 1:
        return parts[0]
    stem = os.path.basename(path)
    index = stem.lower().find(token)
    code = stem[:index].strip("_-. ") if index > 0 else ""
    return code or "scec"


def load_scec(root: str) -> tuple[list[Solution], list[str]]:
    """The uploads of the codes under a directory, as solutions, and what was skipped."""
    codes: dict[str, Solution] = {}
    skipped = []
    for path in sorted(glob.glob(os.path.join(root, "**", "*"), recursive=True)):
        if not os.path.isfile(path):
            continue
        base = os.path.basename(path).lower()
        try:
            names, rows, header = read_cvws_table(path)
        except (OSError, UnicodeDecodeError) as error:
            skipped.append(f"{path}: {error}")
            continue
        if rows.size == 0:
            skipped.append(f"{path}: no data")
            continue
        comments = " ".join(header).lower()
        is_contour = (
            "cplot" in base
            or names[:3] == ["j", "k", "t"]
            or ("rupture time" in comments and rows.shape[1] == 3)
        )
        station = None if is_contour else station_from_name(base)
        station = station or (None if is_contour else station_from_location(header))
        if is_contour:
            code = _upload_code(path, root, "cplot")
        elif station is not None:
            code = _upload_code(path, root, station[0])
        else:
            skipped.append(f"{path}: neither a station nor rupture contour data")
            continue
        solution = codes.setdefault(code, Solution(id=code, label=code, scec=True))
        if is_contour:
            along, down, time = rows[:, 0], rows[:, 1], rows[:, 2]
            if np.nanmax(np.abs(np.r_[along, down])) < 100.0:
                along, down = 1000.0 * along, 1000.0 * down
                solution.warnings.append("rupture contour coordinates read as km")
            if np.nanmedian(down) < 0.0:
                down = -down
                solution.warnings.append("rupture contour depths read as negative")
            solution.rupture = (along, down, time)
            continue
        kind, name = station
        if not names or len(names) != rows.shape[1]:
            skipped.append(f"{path}: no line naming the fields")
            continue
        time_column = next(
            (names.index(name_) for name_ in ("t", "time") if name_ in names), 0
        )
        keep = _increasing(rows[:, time_column])
        fields = {
            name_: rows[keep, i]
            for i, name_ in enumerate(names)
            if i != time_column and name_ in (*FAULT_FIELDS, *BODY_FIELDS)
        }
        solution.stations[kind][name] = Station(name, rows[keep, time_column], fields)
    return [codes[name] for name in sorted(codes)], skipped


# --- metrics --------------------------------------------------------------------------


def nrms(a: np.ndarray, b: np.ndarray) -> float:
    """The RMS of a - b over the mean RMS of a and b, in percent: 0 to 200."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    norm_a = math.sqrt(float(np.sum(a * a)))
    norm_b = math.sqrt(float(np.sum(b * b)))
    if norm_a + norm_b == 0.0:
        return 0.0
    return 200.0 * math.sqrt(float(np.sum((a - b) ** 2))) / (norm_a + norm_b)


def sample(
    station: Station, quantity: str, times: np.ndarray, signs: dict | None = None
) -> np.ndarray | None:
    """A quantity of a station at the given times, a column per component."""
    components = VECTORS.get(quantity, (quantity,))
    if not all(name in station.fields for name in components):
        return None
    signs = signs or {}
    return np.stack(
        [
            signs.get(name, 1.0) * np.interp(times, station.time, station.fields[name])
            for name in components
        ],
        axis=1,
    )


def _window(a: Station, b: Station, end_time: float) -> tuple[float, float]:
    """The times both stations cover, without the first sample of either."""
    start = max(a.time[min(1, len(a.time) - 1)], b.time[min(1, len(b.time) - 1)])
    stop = min(a.time[-1], b.time[-1], end_time)
    return float(start), float(stop)


def _time_grid(a: Station, b: Station, start: float, stop: float) -> np.ndarray:
    steps = [float(np.median(np.diff(s.time))) for s in (a, b) if len(s.time) > 1]
    step = min(steps) if steps else (stop - start) / 1000.0
    count = int(np.clip(np.ceil((stop - start) / max(step, 1e-12)) + 1, 64, 200001))
    return np.linspace(start, stop, count)


def _shifted_nrms(
    a: Station, b: Station, quantity: str, grid: np.ndarray, shift: float, signs
) -> float:
    """nrms of a(t + shift) against b(t), over the times a covers."""
    inside = (grid + shift >= a.time[0]) & (grid + shift <= a.time[-1])
    if np.count_nonzero(inside) < 16:
        return math.inf
    values_a = sample(a, quantity, grid[inside] + shift, signs)
    values_b = sample(b, quantity, grid[inside])
    if values_a is None or values_b is None:
        return math.inf
    return nrms(values_a, values_b)


def best_shift(
    a: Station, b: Station, quantity: str, grid: np.ndarray, max_shift: float, signs
) -> float:
    """The time shift s within +-max_shift for which a(t + s) fits b(t) best.

    A scan over the shifts on a thinned grid gives the candidate, a golden-section
    search around it on the full grid the shift itself; no shift is kept unless a shift
    fits better. The fit is the metric itself, not a cross-correlation: the series of a
    station need not decay within the window, and a correlation of such series favours
    the shifts that keep the most of them in the window.
    """
    if sample(b, quantity, grid[:1]) is None or sample(a, quantity, grid[:1]) is None:
        return 0.0
    if max_shift <= 0.0 or len(grid) < 16:
        return 0.0
    step = grid[1] - grid[0]
    thinned = grid[:: max(1, len(grid) // 4000)]
    coarse = max(step, max_shift / 100.0)
    shifts = np.arange(-max_shift, max_shift + 0.5 * coarse, coarse)
    costs = [_shifted_nrms(a, b, quantity, thinned, shift, signs) for shift in shifts]
    candidate = float(shifts[int(np.argmin(costs))])

    def cost(shift):
        return _shifted_nrms(a, b, quantity, grid, shift, signs)

    low, high = candidate - coarse, candidate + coarse
    ratio = (math.sqrt(5.0) - 1.0) / 2.0
    x1, x2 = high - ratio * (high - low), low + ratio * (high - low)
    f1, f2 = cost(x1), cost(x2)
    for _ in range(40):
        if f1 <= f2:
            high, x2, f2 = x2, x1, f1
            x1 = high - ratio * (high - low)
            f1 = cost(x1)
        else:
            low, x1, f1 = x1, x2, f2
            x2 = low + ratio * (high - low)
            f2 = cost(x2)
    shift = x1 if f1 <= f2 else x2
    shift = float(np.clip(shift, -max_shift, max_shift))
    return shift if cost(shift) < cost(0.0) else 0.0


def compare_stations(
    a: Station,
    b: Station,
    kind: str,
    end_time: float,
    max_shift: float,
    signs: dict | None = None,
) -> list[tuple[str, str, float]]:
    """(quantity, metric, value) of a station pair; a is shifted against b."""
    start, stop = _window(a, b, end_time)
    if not stop > start:
        return []
    grid = _time_grid(a, b, start, stop)
    shift = best_shift(a, b, PRIMARY[kind], grid, max_shift, signs)
    results = []
    for quantity in QUANTITIES[kind]:
        values_a = sample(a, quantity, grid, signs)
        values_b = sample(b, quantity, grid)
        if values_a is None or values_b is None:
            continue
        results.append((quantity, "nrms", nrms(values_a, values_b)))
        shifted = _shifted_nrms(a, b, quantity, grid, shift, signs)
        results.append((quantity, "nrms_shifted", shifted))
    if results:
        results.append((PRIMARY[kind], "tshift", shift))
    return results


@dataclass
class Grid:
    """The points along strike and down dip at which the fault data is compared."""

    along: np.ndarray
    down: np.ndarray
    spacing: float

    @classmethod
    def over(cls, region, spacing: float) -> "Grid":
        along = np.arange(region[0], region[1] + 0.5 * spacing, spacing)
        down = np.arange(region[2], region[3] + 0.5 * spacing, spacing)
        mesh_along, mesh_down = np.meshgrid(along, down)
        return cls(mesh_along.ravel(), mesh_down.ravel(), spacing)

    @property
    def shape(self):
        return (len(np.unique(self.down)), len(np.unique(self.along)))


def _interpolate(along, down, values, grid: Grid):
    from scipy.interpolate import LinearNDInterpolator, NearestNDInterpolator
    from scipy.spatial import QhullError

    points = np.c_[along, down]
    if len(points) == 0:
        return np.full(grid.along.shape, np.nan)
    try:
        interpolator = LinearNDInterpolator(points, values)
    except (QhullError, ValueError):
        interpolator = NearestNDInterpolator(points, values)
    return np.asarray(interpolator(np.c_[grid.along, grid.down]), dtype=float)


def _coverage(along, down, grid: Grid) -> tuple[np.ndarray, np.ndarray]:
    """Whether a point of the data lies near each grid point, and which one does."""
    from scipy.spatial import cKDTree

    points = np.c_[along, down]
    tree = cKDTree(points)
    distance, index = tree.query(np.c_[grid.along, grid.down])
    if len(points) > 1:
        spacing = float(np.median(tree.query(points, k=2)[0][:, 1]))
    else:
        spacing = 0.0
    return distance <= 2.0 * max(grid.spacing, spacing), index


@dataclass
class GriddedFault:
    rupture: np.ndarray | None = None  # rupture time, nan where none or unknown
    known: np.ndarray | None = None  # whether the rupture data covers the point
    maps: dict = field(default_factory=dict)  # slip and peak slip rate
    maps_known: np.ndarray | None = None


def grid_fault(solution: Solution, grid: Grid) -> GriddedFault:
    gridded = GriddedFault()
    if solution.rupture is not None:
        along, down, time = (np.asarray(x, dtype=float) for x in solution.rupture)
        known, nearest = _coverage(along, down, grid)
        ruptured = np.isfinite(time) & (time > 0.0) & (time < NO_RUPTURE)
        values = np.full(grid.along.shape, np.nan)
        if np.count_nonzero(ruptured) >= 3:
            values = _interpolate(along[ruptured], down[ruptured], time[ruptured], grid)
        gridded.rupture = np.where(known & ruptured[nearest], values, np.nan)
        gridded.known = known
    if solution.maps is not None:
        along, down, values = solution.maps
        known, _ = _coverage(along, down, grid)
        gridded.maps_known = known
        for name in ("slip", "psr"):
            if name in values:
                gridded.maps[name] = np.where(
                    known, _interpolate(along, down, values[name], grid), np.nan
                )
    return gridded


def compare_rupture(
    a: GriddedFault, b: GriddedFault, until: float
) -> list[tuple[str, str, float]]:
    """The rupture time of a against b, over the points both cover, up to a time."""
    if a.rupture is None or b.rupture is None:
        return []
    known = a.known & b.known
    with np.errstate(invalid="ignore"):
        ruptured_a = known & np.isfinite(a.rupture) & (a.rupture <= until)
        ruptured_b = known & np.isfinite(b.rupture) & (b.rupture <= until)
    both = ruptured_a & ruptured_b
    either = ruptured_a | ruptured_b
    if not np.any(either):
        return []
    results = [
        (
            "rt",
            "area_mismatch",
            float(np.count_nonzero(either & ~both)) / np.count_nonzero(either),
        )
    ]
    if np.any(both):
        difference = (a.rupture - b.rupture)[both]
        results += [
            ("rt", "rms_ms", 1000.0 * math.sqrt(float(np.mean(difference**2)))),
            ("rt", "bias_ms", 1000.0 * float(np.mean(difference))),
            ("rt", "max_ms", 1000.0 * float(np.max(np.abs(difference)))),
        ]
    return results


def compare_maps(a: GriddedFault, b: GriddedFault) -> list[tuple[str, str, float]]:
    if a.maps_known is None or b.maps_known is None:
        return []
    results = []
    for name in ("slip", "psr"):
        if name not in a.maps or name not in b.maps:
            continue
        valid = (
            a.maps_known
            & b.maps_known
            & np.isfinite(a.maps[name])
            & np.isfinite(b.maps[name])
        )
        if not np.any(valid):
            continue
        values_a, values_b = a.maps[name][valid], b.maps[name][valid]
        results.append((name, "nrms", nrms(values_a, values_b)))
        results.append(
            (name, "rms", math.sqrt(float(np.mean((values_a - values_b) ** 2))))
        )
    return results


def compare_moment(
    a: Station, b: Station, end_time: float, max_shift: float
) -> list[tuple[str, str, float]]:
    until = min(a.time[-1], b.time[-1], end_time)
    moment_a = float(np.interp(until, a.time, a.fields["moment"]))
    moment_b = float(np.interp(until, b.time, b.fields["moment"]))
    results = []
    if moment_b != 0.0:
        results.append(("m0", "rel", abs(moment_a - moment_b) / abs(moment_b)))
    # a moment rate needs a few samples of the moment to be one
    if len(a.time) >= 5 and len(b.time) >= 5:
        results += compare_stations(a, b, "moment", end_time, max_shift)
    return results


# --- the comparisons ------------------------------------------------------------------


@dataclass
class Settings:
    end_time: float
    max_shift: float
    signs: dict


def compare_solutions(
    a: Solution,
    b: Solution,
    gridded: dict,
    settings: Settings,
    signs: dict | None = None,
) -> list[dict]:
    """Every metric of a measured against b."""
    rows = []

    def add(category, station, results):
        for quantity, metric, value in results:
            rows.append(
                {
                    "category": category,
                    "station": station,
                    "quantity": quantity,
                    "metric": metric,
                    "value": float(value),
                }
            )

    until = min(a.fault_time, b.fault_time, settings.end_time)
    add("rupture", "", compare_rupture(gridded[a.id], gridded[b.id], until))
    if not (a.scec or b.scec):
        if math.isclose(a.fault_time, b.fault_time, rel_tol=1e-6, abs_tol=1e-6):
            add("fault-map", "", compare_maps(gridded[a.id], gridded[b.id]))
    for kind in ("fault", "body"):
        for name in sorted(set(a.stations[kind]) & set(b.stations[kind])):
            add(
                kind,
                name,
                compare_stations(
                    a.stations[kind][name],
                    b.stations[kind][name],
                    kind,
                    settings.end_time,
                    settings.max_shift,
                    signs,
                ),
            )
    if a.moment is not None and b.moment is not None:
        add(
            "moment",
            "",
            compare_moment(a.moment, b.moment, settings.end_time, settings.max_shift),
        )
    return rows


def choose_signs(
    reference: Solution | None, codes: list[Solution], fixed: dict, end_time: float
) -> tuple[dict, dict]:
    """The sign of each station field of SeisSol against the uploads, and how it was
    chosen: fixed by the manifest, or the one under which the reference agrees better
    with the uploads, summed over all stations and codes."""
    signs, origin = {}, {}
    for name, value in (fixed or {}).items():
        signs[name] = float(value)
        origin[name] = "manifest"
    if reference is None or not codes:
        return signs, origin
    scores: dict[str, list[float]] = {}
    for kind, fields in (("fault", tuple(FAULT_FIELDS)), ("body", BODY_FIELDS)):
        for code in codes:
            for name in set(reference.stations[kind]) & set(code.stations[kind]):
                a, b = reference.stations[kind][name], code.stations[kind][name]
                start, stop = _window(a, b, end_time)
                if not stop > start:
                    continue
                grid = _time_grid(a, b, start, stop)
                for quantity in fields:
                    if quantity in signs and origin.get(quantity) == "manifest":
                        continue
                    values_a = sample(a, quantity, grid)
                    values_b = sample(b, quantity, grid)
                    if values_a is None or values_b is None:
                        continue
                    score = scores.setdefault(quantity, [0.0, 0.0])
                    score[0] += nrms(values_a, values_b)
                    score[1] += nrms(-values_a, values_b)
    for quantity, (plus, minus) in scores.items():
        signs[quantity] = -1.0 if minus < plus else 1.0
        origin[quantity] = "chosen"
    return signs, origin


def _summary_value(rows: list[dict], category, quantity, metric, reduce):
    values = [
        row["value"]
        for row in rows
        if row["category"] == category
        and row["quantity"] == quantity
        and row["metric"] == metric
        and math.isfinite(row["value"])
    ]
    return float(reduce(values)) if values else None


def summarize(rows: list[dict], kind: str) -> dict:
    """The few numbers of a comparison that the report and the plots show."""
    metric = "nrms_shifted" if kind.startswith("scec") else "nrms"
    summary = {
        "rt_rms_ms": _summary_value(rows, "rupture", "rt", "rms_ms", np.max),
        "rt_bias_ms": _summary_value(rows, "rupture", "rt", "bias_ms", np.max),
        "rt_area_mismatch": _summary_value(
            rows, "rupture", "rt", "area_mismatch", np.max
        ),
        "slip_nrms": _summary_value(rows, "fault-map", "slip", "nrms", np.max),
        "psr_nrms": _summary_value(rows, "fault-map", "psr", "nrms", np.max),
        "m0_rel": _summary_value(rows, "moment", "m0", "rel", np.max),
    }
    for key, category, quantity in (
        ("fault_rate", "fault", "2d-rate"),
        ("fault_psi", "fault", "psi"),
        ("fault_stress", "fault", "2d-stress"),
        ("body_vel", "body", "vel"),
    ):
        summary[f"{key}_median"] = _summary_value(
            rows, category, quantity, metric, np.median
        )
        summary[f"{key}_max"] = _summary_value(rows, category, quantity, metric, np.max)
    return summary


# --- report ---------------------------------------------------------------------------


def _format(value, digits=3) -> str:
    if value is None or (isinstance(value, float) and not math.isfinite(value)):
        return "–"
    if value == 0:
        return "0"
    return f"{value:.{digits}g}"


def _pair(summary: dict, key: str) -> str:
    median, largest = summary.get(f"{key}_median"), summary.get(f"{key}_max")
    if median is None and largest is None:
        return "–"
    return f"{_format(median)} / {_format(largest)}"


_COLUMNS = (
    ("RT RMS [ms]", lambda s: _format(s.get("rt_rms_ms"))),
    ("RT bias [ms]", lambda s: _format(s.get("rt_bias_ms"))),
    (
        "RT area mismatch [%]",
        lambda s: _format(
            None if s.get("rt_area_mismatch") is None else 100 * s["rt_area_mismatch"]
        ),
    ),
    ("slip rate nrms (med / max)", lambda s: _pair(s, "fault_rate")),
    ("psi nrms (med / max)", lambda s: _pair(s, "fault_psi")),
    ("velocity nrms (med / max)", lambda s: _pair(s, "body_vel")),
    ("final slip nrms", lambda s: _format(s.get("slip_nrms"))),
    ("PSR nrms", lambda s: _format(s.get("psr_nrms"))),
    (
        "M0 diff [%]",
        lambda s: _format(None if s.get("m0_rel") is None else 100 * s["m0_rel"]),
    ),
)


def _table(head: list[str], lines: list[list[str]]) -> list[str]:
    out = ["| " + " | ".join(head) + " |", "|" + "|".join("---" for _ in head) + "|"]
    out += ["| " + " | ".join(line) + " |" for line in lines]
    return out


def write_report(
    path: str,
    title: str,
    runs: list[Solution],
    codes: list[Solution],
    reference: str,
    comparisons: list[dict],
    settings: Settings,
    sign_origin: dict,
    region,
    plots: list[str],
    skipped: list[str],
) -> None:
    lines = [f"# {title}", ""]
    lines += [
        f"Reference: `{reference}`. Region of the rupture time comparison: "
        f"{region[0] / 1000:g} to {region[1] / 1000:g} km along strike, "
        f"{region[2] / 1000:g} to {region[3] / 1000:g} km down dip. "
        + (
            f"Compared up to t = {settings.end_time:g} s."
            if math.isfinite(settings.end_time)
            else "Compared over the time both sides cover."
        ),
        "",
        "Time series (nrms): RMS of the difference over the mean RMS of both, in "
        "percent (0 identical, 200 opposite). Between runs without a time shift, "
        "against the uploads after the shift that fits best (as in the CVWS); "
        "med / max over the stations. Rupture time: over the area that ruptured in "
        "both; the RMS is the CVWS metric, the bias (mean difference, first run minus "
        "second) shows a shift of the front as a whole and does not respond to "
        "differences at the scale of the elements, which also come from how the "
        "output samples a cell.",
        "",
        "## Runs",
        "",
    ]
    rows = []
    for run in runs:
        stations = f"{len(run.stations['fault'])} / {len(run.stations['body'])}"
        rows.append(
            [
                f"`{run.id}`",
                run.label,
                run.series,
                run.meta.get("commit", ""),
                run.meta.get("date", ""),
                _format(run.fault_time if math.isfinite(run.fault_time) else None, 6),
                stations,
                "; ".join(run.warnings) or "",
            ]
        )
    lines += _table(
        [
            "id",
            "label",
            "series",
            "commit",
            "date",
            "fault output at [s]",
            "stations on / off fault",
            "notes",
        ],
        rows,
    )
    head = ["run", "against"] + [name for name, _ in _COLUMNS]
    sections = (
        ("previous", "Jumps: each run against the one before it in its series"),
        ("reference", f"Against the reference `{reference}`"),
        ("extra", "Further comparisons (compare_with)"),
        ("scec", "Against the uploads of other codes"),
        ("scec-pair", "Spread among the uploads (for scale)"),
    )
    for kind, heading in sections:
        selected = [c for c in comparisons if c["kind"] == kind]
        if not selected:
            continue
        lines += ["", f"## {heading}", ""]
        lines += _table(
            head,
            [
                [f"`{c['run']}`", f"`{c['against']}`"]
                + [render(c["summary"]) for _, render in _COLUMNS]
                for c in selected
            ],
        )
    if codes:
        lines += ["", "## Signs of SeisSol against the uploads", ""]
        lines += [
            f"- `{name}`: {'+' if settings.signs.get(name, 1.0) > 0 else '-'}1 "
            f"({sign_origin.get(name, 'default')})"
            for name in sorted(settings.signs)
        ]
    if skipped:
        lines += ["", "## Files of the uploads that were left out", ""]
        lines += [f"- {entry}" for entry in skipped]
    if plots:
        lines += ["", "## Plots", ""]
        lines += [f"![{os.path.basename(p)}]({os.path.basename(p)})" for p in plots]
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


# --- plots ----------------------------------------------------------------------------

# the categorical slots, the ordinal ramp (light to dark) and the neutrals of the plots
SLOTS = ("#2a78d6", "#eb6834", "#1baf7a")
RAMP = (
    "#86b6ef",
    "#6da7ec",
    "#5598e7",
    "#3987e5",
    "#2a78d6",
    "#256abf",
    "#1c5cab",
    "#184f95",
    "#104281",
    "#0d366b",
)
DIVERGING = (
    "#184f95",
    "#3987e5",
    "#9ec5f4",
    "#f0efec",
    "#f3aaa9",
    "#e34948",
    "#a32f2e",
)
SURFACE, INK, INK_SECONDARY, GRID_LINE, OTHER_CODES = (
    "#fcfcfb",
    "#0b0b0b",
    "#52514e",
    "#e6e5e1",
    "#a3a29b",
)


def _style(plt) -> None:
    plt.rcParams.update(
        {
            "figure.facecolor": SURFACE,
            "axes.facecolor": SURFACE,
            "savefig.facecolor": SURFACE,
            "axes.edgecolor": GRID_LINE,
            "axes.labelcolor": INK_SECONDARY,
            "axes.titlecolor": INK,
            "axes.grid": True,
            "grid.color": GRID_LINE,
            "grid.linewidth": 0.6,
            "grid.linestyle": "-",
            "xtick.color": INK_SECONDARY,
            "ytick.color": INK_SECONDARY,
            "text.color": INK,
            "font.size": 8,
            "axes.titlesize": 9,
            "legend.frameon": False,
        }
    )


def _ramp(count: int) -> list[str]:
    if count <= 1:
        return [RAMP[-1]] * count
    return [RAMP[int(round(i))] for i in np.linspace(0, len(RAMP) - 1, count)]


# the metrics of the jump plot; a signed one is drawn on a linear scale around zero
_PANELS = (
    ("rt_rms_ms", "Rupture time: RMS difference [ms]", False),
    ("rt_bias_ms", "Rupture time: mean difference [ms]", True),
    ("fault_rate_max", "On-fault slip rate: largest nrms over the stations [%]", False),
    ("fault_psi_max", "State variable: largest nrms over the stations [%]", False),
    ("slip_nrms", "Final slip on the fault: nrms [%]", False),
    ("body_vel_max", "Off-fault velocity: largest nrms over the stations [%]", False),
    ("m0_rel", "Seismic moment: relative difference", False),
)


def plot_jumps(path: str, runs: list[Solution], comparisons: list[dict]) -> bool:
    """Per metric, a run's distance to the one before it, to the reference and to the
    uploads (their median), along the series."""
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, NullFormatter

    def values(kind, key):
        out = []
        for run in runs:
            selected = [
                c["summary"].get(key)
                for c in comparisons
                if c["run"] == run.id and c["kind"] == kind
            ]
            selected = [v for v in selected if v is not None]
            out.append(float(np.median(selected)) if selected else np.nan)
        return np.array(out)

    panels = []
    for key, title, signed in _PANELS:
        lines = {kind: values(kind, key) for kind in ("previous", "reference", "scec")}
        if any(np.any(np.isfinite(v)) for v in lines.values()):
            panels.append((title, signed, lines))
    if not panels:
        return False
    _style(plt)
    count = len(runs)
    figure, axes = plt.subplots(
        len(panels),
        1,
        sharex=True,
        figsize=(max(7.0, 0.55 * count + 3.0), 1.9 * len(panels) + 1.4),
        squeeze=False,
    )
    x = np.arange(count)
    labels = {
        "previous": "against the run before it",
        "reference": "against the reference",
        "scec": "against the uploads (median)",
    }
    for ax, (title, signed, lines) in zip(axes[:, 0], panels):
        previous = lines["previous"]
        shown = previous if signed else np.where(previous > 0, previous, np.nan)
        ax.bar(x, shown, width=0.6, color=SLOTS[0], label=labels["previous"])
        for i in [] if signed else np.flatnonzero(previous == 0):
            ax.annotate(
                "0",
                (x[i], 0),
                xycoords=("data", "axes fraction"),
                xytext=(0, 3),
                textcoords="offset points",
                ha="center",
                color=INK_SECONDARY,
            )
        for kind, color, marker in (
            ("reference", SLOTS[1], "o"),
            ("scec", SLOTS[2], "s"),
        ):
            line = (
                lines[kind]
                if signed
                else np.where(lines[kind] > 0, lines[kind], np.nan)
            )
            if np.any(np.isfinite(line)):
                ax.plot(
                    x,
                    line,
                    color=color,
                    marker=marker,
                    lw=1.5,
                    ms=4,
                    label=labels[kind],
                )
        if signed:
            ax.axhline(0.0, color=INK_SECONDARY, lw=0.8)
        else:
            finite = np.concatenate(
                [v[np.isfinite(v) & (v > 0)] for v in lines.values()]
            )
            if finite.size and finite.max() / finite.min() > 20:
                ax.set_yscale("log")
                ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
                ax.yaxis.set_minor_formatter(NullFormatter())
        ax.set_title(title, loc="left")
        ax.set_axisbelow(True)
    axes[-1, 0].set_xticks(x)
    axes[-1, 0].set_xticklabels([run.label for run in runs], rotation=40, ha="right")
    handles, names = [], []
    for ax in axes[:, 0]:
        for handle, name in zip(*ax.get_legend_handles_labels()):
            if name not in names:
                handles.append(handle)
                names.append(name)
    figure.legend(handles, names, loc="upper center", ncol=len(names))
    figure.tight_layout(rect=(0, 0, 1, 0.97))
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return True


def plot_rupture(
    path: str,
    grid: Grid,
    gridded: dict,
    solutions: list[Solution],
    reference: Solution,
    until: float,
) -> bool:
    """The rupture time of every run and upload minus that of the reference."""
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap

    base = gridded[reference.id].rupture
    if base is None:
        return False
    shape = grid.shape
    panels = []
    for solution in solutions:
        other = gridded[solution.id].rupture
        if solution.id == reference.id or other is None:
            continue
        cut = min(until, solution.fault_time)
        with np.errstate(invalid="ignore"):
            valid = (base <= cut) & (other <= cut)
        difference = np.where(valid, 1000.0 * (other - base), np.nan)
        panels.append((solution, difference.reshape(shape)))
    if not panels:
        return False
    _style(plt)
    finite = np.concatenate([np.abs(d[np.isfinite(d)]) for _, d in panels])
    limit = float(np.percentile(finite, 99)) if finite.size else 1.0
    limit = limit if limit > 0 else 1.0
    columns = min(3, len(panels))
    rows = int(math.ceil(len(panels) / columns))
    extent = (
        grid.along.min() / 1000,
        grid.along.max() / 1000,
        grid.down.max() / 1000,
        grid.down.min() / 1000,
    )
    figure, axes = plt.subplots(
        rows,
        columns,
        figsize=(3.6 * columns + 1.0, 2.2 * rows + 0.9),
        squeeze=False,
        sharex=True,
        sharey=True,
    )
    colormap = LinearSegmentedColormap.from_list("rupture-difference", DIVERGING)
    contour_times = base.reshape(shape)
    image = None
    for ax, (solution, difference) in zip(axes.ravel(), panels):
        image = ax.imshow(
            difference,
            extent=extent,
            origin="upper",
            cmap=colormap,
            vmin=-limit,
            vmax=limit,
            aspect="auto",
            interpolation="nearest",
        )
        if np.any(np.isfinite(contour_times)):
            ax.contour(
                np.unique(grid.along) / 1000,
                np.unique(grid.down) / 1000,
                contour_times,
                levels=np.arange(0.5, np.nanmax(contour_times) + 0.5, 0.5),
                colors=INK_SECONDARY,
                linewidths=0.5,
            )
        prefix = "upload " if solution.scec else ""
        ax.set_title(f"{prefix}{solution.label}", loc="left")
        ax.grid(False)
    for ax in axes.ravel()[len(panels) :]:
        ax.set_visible(False)
    for ax in axes[-1, :]:
        ax.set_xlabel("along strike [km]")
    for ax in axes[:, 0]:
        ax.set_ylabel("down dip [km]")
    figure.suptitle(
        f"Rupture time minus that of {reference.label} [ms]; "
        f"contours: {reference.label}, every 0.5 s",
        x=0.01,
        ha="left",
    )
    figure.tight_layout(rect=(0, 0, 0.9, 0.95))
    bar = figure.add_axes((0.92, 0.15, 0.015, 0.7))
    figure.colorbar(image, cax=bar)
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return True


def plot_stations(
    path: str,
    kind: str,
    quantity: str,
    unit: str,
    runs: list[Solution],
    codes: list[Solution],
    reference: Solution,
    signs: dict,
    end_time: float,
) -> bool:
    """One quantity at every station: the runs from light (early) to dark (late), the
    reference in orange, the uploads in grey."""
    import matplotlib.pyplot as plt

    names = sorted(set().union(*(run.stations[kind] for run in runs)))
    if not names:
        return False
    _style(plt)
    columns = min(3, len(names))
    rows = int(math.ceil(len(names) / columns))
    figure, axes = plt.subplots(
        rows, columns, figsize=(3.6 * columns, 2.1 * rows + 1.0), squeeze=False
    )
    colors = _ramp(len(runs))
    for ax, name in zip(axes.ravel(), names):
        for code in codes:
            station = code.stations[kind].get(name)
            if station is not None and quantity in station.fields:
                keep = station.time <= end_time
                ax.plot(
                    station.time[keep],
                    station.fields[quantity][keep],
                    color=OTHER_CODES,
                    lw=0.8,
                    label="uploads",
                )
        for run, color in zip(runs, colors):
            station = run.stations[kind].get(name)
            if station is None or quantity not in station.fields:
                continue
            keep = station.time <= end_time
            is_reference = run.id == reference.id
            values = signs.get(quantity, 1.0) * station.fields[quantity][keep]
            ax.plot(
                station.time[keep],
                values,
                color=SLOTS[1] if is_reference else color,
                lw=1.6 if is_reference else 1.0,
                label=run.label,
                zorder=3 if is_reference else 2,
            )
        ax.set_title(name, loc="left")
    for ax in axes.ravel()[len(names) :]:
        ax.set_visible(False)
    handles, labels = [], []
    for handle, label in zip(*axes[0, 0].get_legend_handles_labels()):
        if label not in labels:
            handles.append(handle)
            labels.append(label)
    figure.legend(
        handles, labels, loc="lower center", ncol=min(6, len(labels)), fontsize=7
    )
    figure.suptitle(f"{quantity} [{unit}] over time [s]", x=0.01, ha="left")
    figure.tight_layout(rect=(0, 0.06 + 0.02 * math.ceil(len(labels) / 6), 1, 0.96))
    figure.savefig(path, dpi=150)
    plt.close(figure)
    return True


def make_plots(out, grid, gridded, runs, codes, reference, comparisons, settings):
    import matplotlib

    matplotlib.use("Agg")
    plots = []
    path = os.path.join(out, "jumps.png")
    if plot_jumps(path, runs, comparisons):
        plots.append(path)
    path = os.path.join(out, "rupture-time.png")
    if plot_rupture(
        path,
        grid,
        gridded,
        runs + codes,
        reference,
        min(reference.fault_time, settings.end_time),
    ):
        plots.append(path)
    for kind, quantity, unit in (
        ("fault", "h-slip-rate", "m/s"),
        ("fault", "psi", "1"),
        ("body", "h-vel", "m/s"),
    ):
        path = os.path.join(out, f"stations-{kind}-{quantity}.png")
        if plot_stations(
            path,
            kind,
            quantity,
            unit,
            runs,
            codes,
            reference,
            settings.signs,
            settings.end_time,
        ):
            plots.append(path)
    return plots


# --- main -----------------------------------------------------------------------------


def load_manifest(path: str) -> dict:
    with open(path) as handle:
        text = handle.read()
    if path.endswith(".json"):
        return json.loads(text)
    import yaml

    return yaml.safe_load(text)


def _resolve(path: str | None, base: str) -> str | None:
    if not path:
        return None
    return path if os.path.isabs(path) else os.path.join(base, path)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__.split("\n\n")[0],
        epilog="The manifest format is described at the top of this script.",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("manifest", help="YAML or JSON file listing the runs")
    parser.add_argument("--out", default="history", help="directory of the results")
    parser.add_argument("--reference", help="id of the reference run")
    parser.add_argument("--scec", help="directory of the uploads of other codes")
    parser.add_argument("--end-time", type=float, help="compare up to this time [s]")
    parser.add_argument(
        "--region",
        type=float,
        nargs=4,
        metavar=("ALONG_MIN", "ALONG_MAX", "DOWN_MIN", "DOWN_MAX"),
        help="region of the rupture time comparison [m]",
    )
    parser.add_argument(
        "--grid-spacing",
        type=float,
        default=100.0,
        help="spacing of the points the fault data is compared at [m]",
    )
    parser.add_argument(
        "--max-shift",
        type=float,
        default=0.5,
        help="largest time shift tried between two time series [s]",
    )
    parser.add_argument(
        "--scec-pairs",
        action="store_true",
        help="compare the uploads with each other too, for scale",
    )
    parser.add_argument("--no-plots", action="store_true", help="write no plots")
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    manifest = load_manifest(args.manifest)
    base = os.path.dirname(os.path.abspath(args.manifest))
    entries = manifest.get("runs") or []
    if not entries:
        print("error: the manifest lists no runs", file=sys.stderr)
        return 2
    ids = [str(entry.get("id", "")) for entry in entries]
    if "" in ids or len(set(ids)) != len(ids):
        print("error: every run needs an id of its own", file=sys.stderr)
        return 2

    axes = Axes.from_spec(manifest.get("axes"))
    runs = [load_run(entry, axes, base) for entry in entries]
    # the default reference: the last run with data of the series "main", if any
    candidates = [run for run in runs if run.has_data] or runs
    main_series = [run for run in candidates if run.series == "main"] or candidates
    reference_id = str(
        args.reference or manifest.get("reference") or main_series[-1].id
    )
    by_id = {run.id: run for run in runs}
    if reference_id not in by_id:
        print(f"error: no run with the id {reference_id}", file=sys.stderr)
        return 2
    reference = by_id[reference_id]

    scec_dir = _resolve(args.scec or manifest.get("scec"), base)
    codes, skipped = ([], [])
    if scec_dir:
        if os.path.isdir(scec_dir):
            codes, skipped = load_scec(scec_dir)
        else:
            skipped = [f"{scec_dir}: no such directory"]
    clash = {code.id for code in codes} & set(by_id)
    for code in codes:
        if code.id in clash:
            code.id = f"scec:{code.id}"

    end_time = args.end_time or manifest.get("end_time") or math.inf
    region = args.region or manifest.get("region")
    if region is None and str(manifest.get("scenario", "")).lower() in SCENARIOS:
        region = SCENARIOS[str(manifest["scenario"]).lower()]
    if region is None and reference.rupture is not None:
        along, down, _ = reference.rupture
        region = (
            float(along.min()),
            float(along.max()),
            float(down.min()),
            float(down.max()),
        )
    region = tuple(float(x) for x in (region or (-15e3, 15e3, 0.0, 15e3)))
    grid = Grid.over(region, args.grid_spacing)

    # the signs follow from the reference, or else from the first run with stations
    sign_source = next(
        (run for run in [reference] + runs if any(run.stations.values())), None
    )
    signs, sign_origin = choose_signs(
        sign_source, codes, manifest.get("scec_signs"), end_time
    )
    settings = Settings(float(end_time), args.max_shift, signs)
    gridded = {solution.id: grid_fault(solution, grid) for solution in runs + codes}

    pairs = []
    previous: dict[str, Solution] = {}
    for run in runs:
        if not run.has_data:
            continue
        if run.series in previous:
            pairs.append((run, previous[run.series], "previous"))
        previous[run.series] = run
        if run.id != reference.id:
            pairs.append((run, reference, "reference"))
        for other in run.compare_with:
            if other in by_id and other != run.id:
                pairs.append((run, by_id[other], "extra"))
            else:
                run.warnings.append(f"compare_with: no run with the id {other}")
        pairs += [(run, code, "scec") for code in codes]
    if args.scec_pairs:
        pairs += [
            (codes[i], codes[j], "scec-pair")
            for i in range(len(codes))
            for j in range(i + 1, len(codes))
        ]

    records, comparisons = [], []
    for a, b, kind in pairs:
        rows = compare_solutions(
            a, b, gridded, settings, signs if (b.scec and not a.scec) else None
        )
        for row in rows:
            records.append({"run": a.id, "against": b.id, "kind": kind, **row})
        comparisons.append(
            {
                "run": a.id,
                "against": b.id,
                "kind": kind,
                "summary": summarize(rows, kind),
            }
        )

    os.makedirs(args.out, exist_ok=True)
    columns = [
        "run",
        "against",
        "kind",
        "category",
        "station",
        "quantity",
        "metric",
        "value",
    ]
    pd.DataFrame(records, columns=columns).to_csv(
        os.path.join(args.out, "metrics.csv"), index=False
    )
    plots = []
    if not args.no_plots:
        try:
            plots = make_plots(
                args.out, grid, gridded, runs, codes, reference, comparisons, settings
            )
        except ImportError as error:
            print(f"warning: no plots ({error})", file=sys.stderr)
    title = manifest.get("title") or (
        f"{manifest.get('scenario', 'Benchmark')}: {len(runs)} runs compared"
    )
    write_report(
        os.path.join(args.out, "report.md"),
        title,
        runs,
        codes,
        reference.id,
        comparisons,
        settings,
        sign_origin,
        region,
        plots,
        skipped,
    )
    summary = {
        "reference": reference.id,
        "region": region,
        "end_time": None if not math.isfinite(settings.end_time) else settings.end_time,
        "runs": [
            {
                "id": run.id,
                "label": run.label,
                "series": run.series,
                **run.meta,
                "fault_time": run.fault_time if math.isfinite(run.fault_time) else None,
                "warnings": run.warnings,
            }
            for run in runs
        ],
        "uploads": [{"code": code.id, "warnings": code.warnings} for code in codes],
        "signs": {name: settings.signs[name] for name in sorted(settings.signs)},
        "comparisons": comparisons,
        "skipped": skipped,
    }
    with open(os.path.join(args.out, "summary.json"), "w") as handle:
        json.dump(summary, handle, indent=1)

    jumps = [c for c in comparisons if c["kind"] == "previous"]
    if jumps:
        print("jumps (run against the one before it):")
        for c in jumps:
            s = c["summary"]
            print(
                f"  {c['against']:>24} -> {c['run']:<24} "
                f"RT RMS {_format(s['rt_rms_ms'])} ms, "
                f"slip rate nrms max {_format(s['fault_rate_max'])} %, "
                f"final slip nrms {_format(s['slip_nrms'])} %"
            )
    print(f"written to {os.path.abspath(args.out)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
