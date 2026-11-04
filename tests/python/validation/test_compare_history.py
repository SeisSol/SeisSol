# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Tests for scripts/validate/compare-history.py

The script compares runs of one rupture benchmark made with different versions of
SeisSol, and the uploads of other codes to the CVWS. What it reads comes in the formats
of several generations of SeisSol, and of the CVWS; the tests write small files of each.
"""

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pytest

SEISSOL_ROOT = Path(__file__).resolve().parents[3]
# Hyphenated filename — load via spec; dataclasses need the module registered
_spec = importlib.util.spec_from_file_location(
    "compare_history",
    SEISSOL_ROOT / "scripts" / "validate" / "compare-history.py",
)
ch = importlib.util.module_from_spec(_spec)
sys.modules["compare_history"] = ch
_spec.loader.exec_module(ch)


# ============================================================================
# writers of the formats the script reads
# ============================================================================


def write_text_receiver(path, variables, point, rows, fortran=False):
    """A receiver text file as SeisSol writes it."""
    with open(path, "w") as handle:
        handle.write('TITLE = "Temporal Signal for receiver number 00001"\n')
        handle.write("VARIABLES = " + ",".join(f'"{v}"' for v in variables) + "\n")
        for axis, value in zip(("x1", "x2", "x3"), point):
            text = f"{value:.12E}".replace("E", "D") if fortran else f"{value:.12e}"
            handle.write(f"# {axis}  {text}\n")
        for row in rows:
            handle.write("  ".join(f"{v:.15e}" for v in row) + "\n")


def write_hdf5_receivers(path, group_name, id_name, points, times, columns):
    """An HDF5 receiver file: one table, a column per receiver."""
    h5py = pytest.importorskip("h5py")
    count = len(points)
    dtype = np.dtype([(name, "<f8") for name in ["Time", *columns]])
    table = np.zeros((len(times), count), dtype=dtype)
    for receiver in range(count):
        table["Time"][:, receiver] = times
        for name, values in columns.items():
            table[name][:, receiver] = values[receiver]
    with h5py.File(path, "w") as handle:
        group = handle.create_group(group_name)
        group["Coordinates"] = np.asarray(points, dtype=float)
        group["Index"] = np.array([[0, i] for i in range(count)], dtype=np.uint64)
        group[id_name] = np.arange(count, dtype=np.uint64)
        group["SimulationIndex"] = np.zeros(count, dtype=np.uint64)
        group["group0"] = table


def write_fault_output(directory, prefix, vertices, cells, times, fields):
    """A fault output as Xdmf with one HDF5 file for the cells and one for the
    vertices, a row of each quantity per output time."""
    h5py = pytest.importorskip("h5py")
    cells = np.asarray(cells, dtype=np.uint64)
    count = len(cells)
    with h5py.File(directory / f"{prefix}-fault_cell.h5", "w") as handle:
        handle["mesh0/connect"] = cells
        for name, values in fields.items():
            handle[f"mesh0/{name}"] = np.asarray(values, dtype=float)
    with h5py.File(directory / f"{prefix}-fault_vertex.h5", "w") as handle:
        handle["mesh0/geometry"] = np.asarray(vertices, dtype=float)
    grids = []
    for step, time in enumerate(times):
        attributes = "".join(
            f'<Attribute Name="{name}" Center="Cell">'
            f'<DataItem ItemType="HyperSlab" Dimensions="{count}">'
            f'<DataItem NumberType="UInt" Precision="4" Format="XML" Dimensions="3 2">'
            f"{step} 0 1 1 1 {count}</DataItem>"
            f'<DataItem NumberType="Float" Precision="8" Format="HDF" '
            f'Dimensions="{len(times)} {count}">{prefix}-fault_cell.h5:/mesh0/{name}'
            f"</DataItem></DataItem></Attribute>"
            for name in fields
        )
        grids.append(
            f'<Grid Name="step_{step:012d}" GridType="Uniform">'
            f'<Topology TopologyType="Triangle" NumberOfElements="{count}">'
            f'<DataItem NumberType="Int" Precision="8" Format="HDF" '
            f'Dimensions="{count} 3">{prefix}-fault_cell.h5:/mesh0/connect</DataItem>'
            f"</Topology>"
            f'<Geometry name="geo" GeometryType="XYZ" NumberOfElements="{len(vertices)}">'
            f'<DataItem NumberType="Float" Precision="8" Format="HDF" '
            f'Dimensions="{len(vertices)} 3">{prefix}-fault_vertex.h5:/mesh0/geometry'
            f"</DataItem></Geometry>"
            f'<Time Value="{time}"/>{attributes}</Grid>'
        )
    xdmf = (
        '<?xml version="1.0" ?>\n<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>\n'
        '<Xdmf Version="2.0"><Domain>'
        '<Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
        + "".join(grids)
        + "</Grid></Domain></Xdmf>\n"
    )
    (directory / f"{prefix}-fault.xdmf").write_text(xdmf)


def fault_mesh(spacing=1000.0, along=(-6000.0, 6000.0), down=(0.0, 12000.0)):
    """Triangles on the vertical plane y = 0, z = -depth."""
    xs = np.arange(along[0], along[1] + 1, spacing)
    ds = np.arange(down[0], down[1] + 1, spacing)
    vertices = np.array([(x, 0.0, -d) for d in ds for x in xs])
    columns = len(xs)
    cells = []
    for row in range(len(ds) - 1):
        for column in range(columns - 1):
            v = row * columns + column
            cells += [(v, v + 1, v + columns), (v + 1, v + columns + 1, v + columns)]
    return vertices, np.array(cells)


def rupture_times(vertices, cells, delay=0.0, speed=3000.0, hypocentre=(0.0, 6000.0)):
    centroids = vertices[cells].mean(axis=1)
    distance = np.hypot(
        centroids[:, 0] - hypocentre[0], -centroids[:, 2] - hypocentre[1]
    )
    return distance / speed + delay


def pulse(time, onset, width=0.1):
    return np.exp(-(((time - onset) / width) ** 2))


def write_run(directory, delay=0.0, moment_scale=1.0):
    """A run of SeisSol: fault output, on-fault and off-fault receivers in text
    files, energy output. Everything starts later by delay."""
    directory.mkdir(parents=True)
    vertices, cells = fault_mesh()
    rt = rupture_times(vertices, cells, delay)
    times = [0.0, 2.0]
    rt_rows = np.vstack([np.zeros_like(rt), np.where(rt <= 2.0, rt, 0.0)])
    slip = np.vstack([np.zeros_like(rt), np.clip(2.0 - rt, 0.0, None)])
    write_fault_output(
        directory,
        "tpv",
        vertices,
        cells,
        times,
        {"RT": rt_rows, "ASl": slip, "PSR": slip},
    )
    time = np.arange(0.0, 2.0 + 1e-9, 0.002)
    variables = ["Time", "SRs", "SRd", "StV", "Ts0", "Td0", "Pn0", "Sls", "Sld", "ASI"]
    for number, (x, depth) in enumerate(((0.0, 3000.0), (3000.0, 6000.0)), start=1):
        onset = 0.3 + np.hypot(x, depth - 6000.0) / 3000.0 + delay
        rate = 2.0 * pulse(time, onset)
        slip_x = np.cumsum(rate) * 0.002
        rows = np.column_stack(
            [
                time,
                rate,
                0.01 * rate,
                0.56 - 0.1 * pulse(time, onset),
                40e6 - 5e6 * pulse(time, onset),
                np.zeros_like(time),
                -120e6 + np.zeros_like(time),
                slip_x,
                0.01 * slip_x,
                slip_x,
            ]
        )
        write_text_receiver(
            directory / f"tpv-faultreceiver-{number:05d}-00000.dat",
            variables,
            (x, 0.0, -depth),
            rows,
            fortran=True,
        )
    for number, point in enumerate(
        ((0.0, 2000.0, -0.1), (0.0, -2000.0, -0.1)), start=1
    ):
        onset = 0.8 + delay
        velocity = pulse(time, onset, 0.2)
        rows = np.column_stack(
            [time]
            + [np.zeros_like(time)] * 6
            + [velocity, 0.5 * velocity, 0.1 * velocity]
        )
        write_text_receiver(
            directory / f"tpv-receiver-{number:05d}-00000.dat",
            ["Time", "xx", "yy", "zz", "xy", "yz", "xz", "u", "v", "w"],
            point,
            rows,
        )
    energy_times = np.arange(0.0, 2.0 + 1e-9, 0.1)
    moment = moment_scale * 1e18 * np.clip(energy_times - delay, 0.0, None) ** 2
    with open(directory / "tpv-energy.csv", "w") as handle:
        handle.write("time,variable,simulation_index,measurement\n")
        for t, m in zip(energy_times, moment):
            handle.write(f"{t},seismic_moment,0,{m}\n")
            handle.write(f"{t},elastic_kinetic_energy,0,{0.5 * m}\n")


def write_upload(directory, code, delay=0.0, flip_normal_stress=True):
    """The CVWS files of a code: two on-fault stations, one off-fault station and the
    rupture contour data, in MPa, with the normal stress positive in compression."""
    folder = directory / code
    folder.mkdir(parents=True)
    time = np.arange(0.0, 3.0 + 1e-9, 0.005)
    for x, depth in ((0.0, 3000.0), (3000.0, 6000.0)):
        onset = 0.3 + np.hypot(x, depth - 6000.0) / 3000.0 + delay
        rate = 2.0 * pulse(time, onset)
        slip = np.cumsum(rate) * 0.005
        sign = -1.0 if flip_normal_stress else 1.0
        rows = np.column_stack(
            [
                time,
                slip,
                rate,
                40.0 - 5.0 * pulse(time, onset),
                0.01 * slip,
                0.01 * rate,
                np.zeros_like(time),
                sign * -120.0 + np.zeros_like(time),
                0.56 - 0.1 * pulse(time, onset),
            ]
        )
        name = ch.fault_station(x, depth)
        with open(folder / f"{code}_{name}.txt", "w") as handle:
            handle.write(f"# problem=TPV104\n# code={code}\n")
            handle.write(
                "t h-slip h-slip-rate h-shear-stress v-slip v-slip-rate "
                "v-shear-stress n-stress psi\n"
            )
            np.savetxt(handle, rows, fmt="%.7E")
    velocity = pulse(time, 0.8 + delay, 0.2)
    rows = np.column_stack(
        [time, np.zeros_like(time), velocity, np.zeros_like(time), 0.1 * velocity]
        + [np.zeros_like(time), 0.5 * velocity]
    )
    with open(folder / "body020st000dp000.dat", "w") as handle:
        handle.write("# location= 2.0 km off fault, 0.0 km along strike\n")
        handle.write("t h-disp h-vel v-disp v-vel n-disp n-vel\n")
        np.savetxt(handle, rows, fmt="%.7E")
    vertices, cells = fault_mesh()
    centroids = vertices[cells].mean(axis=1)
    rt = rupture_times(vertices, cells, delay)
    with open(folder / "cplot.dat", "w") as handle:
        handle.write("# Column #3 = rupture time (s)\nj k t\n")
        np.savetxt(handle, np.column_stack([centroids[:, 0], -centroids[:, 2], rt]))


# ============================================================================
# station names and geometry
# ============================================================================


class TestStationNames:
    def test_fault_station_names_follow_the_cvws(self):
        assert ch.fault_station(0.0, 7500.0) == "faultst000dp075"
        assert ch.fault_station(-9000.0, 7500.0) == "faultst-090dp075"
        assert ch.fault_station(12000.0, 12000.0) == "faultst120dp120"

    def test_body_station_names_follow_the_cvws(self):
        assert ch.body_station(9000.0, 0.0, 0.0) == "body090st000dp000"
        assert ch.body_station(-6000.0, -12000.0, 0.1) == "body-060st-120dp000"

    def test_names_in_file_names_are_found(self):
        assert ch.station_from_name("barall_faultst-090dp075.txt") == (
            "fault",
            "faultst-090dp075",
        )
        assert ch.station_from_name("x/body-060st120dp000.dat") == (
            "body",
            "body-060st120dp000",
        )
        assert ch.station_from_name("cplot.dat") is None

    def test_location_lines_describe_stations(self):
        on_fault = ["# location= on fault, 9.0 km along strike, 7.5km down-dip"]
        assert ch.station_from_location(on_fault) == ("fault", "faultst090dp075")
        off_fault = ["# location= 6.0 km off fault, -12.0 km along strike, 0 km depth"]
        assert ch.station_from_location(off_fault) == ("body", "body060st-120dp000")
        assert ch.station_from_location(["# author=nobody"]) is None


class TestAxes:
    def test_default_frame_is_that_of_the_examples(self):
        axes = ch.Axes.from_spec(None)
        strike, depth = axes.on_fault([[1000.0, 0.0, -7500.0]])
        assert strike[0] == 1000.0 and depth[0] == 7500.0
        normal, strike, depth = axes.off_fault([[0.0, -9000.0, -0.1]])
        assert (normal[0], strike[0], depth[0]) == (-9000.0, 0.0, 0.1)

    def test_axes_can_be_permuted_and_flipped(self):
        axes = ch.Axes.from_spec({"strike": "-y", "normal": "x", "down": "-z"})
        strike, depth = axes.on_fault([[5.0, 2000.0, -3000.0]])
        assert (strike[0], depth[0]) == (-2000.0, 3000.0)

    def test_an_unknown_axis_is_refused(self):
        with pytest.raises(ValueError):
            ch.parse_axis("w")


# ============================================================================
# metrics
# ============================================================================


class TestNrms:
    def test_identical_series_are_zero(self):
        a = np.sin(np.linspace(0, 3, 50))
        assert ch.nrms(a, a) == 0.0

    def test_opposite_series_are_two_hundred(self):
        a = np.sin(np.linspace(0, 3, 50))
        assert ch.nrms(a, -a) == pytest.approx(200.0)

    def test_a_zero_series_is_two_hundred(self):
        a = np.sin(np.linspace(0, 3, 50))
        assert ch.nrms(a, np.zeros_like(a)) == pytest.approx(200.0)
        assert ch.nrms(np.zeros(3), np.zeros(3)) == 0.0

    def test_is_symmetric_and_bounded(self):
        a, b = np.linspace(0, 1, 20), np.linspace(0, 2, 20)
        assert ch.nrms(a, b) == pytest.approx(ch.nrms(b, a))
        assert ch.nrms(2 * b, b) == pytest.approx(200.0 / 3.0)


class TestStationComparison:
    def station(self, onset, step=0.002):
        time = np.arange(0.0, 3.0, step)
        rate = pulse(time, onset)
        return ch.Station(
            "faultst000dp075",
            time,
            {"h-slip-rate": rate, "v-slip-rate": 0.1 * rate, "psi": 0.5 - 0.1 * rate},
        )

    def test_the_time_shift_is_found(self):
        a, b = self.station(1.0), self.station(1.123, step=0.005)
        results = {
            (q, m): v for q, m, v in ch.compare_stations(a, b, "fault", 10.0, 0.5)
        }
        assert results[("2d-rate", "tshift")] == pytest.approx(-0.123, abs=1e-3)
        assert results[("2d-rate", "nrms_shifted")] < 1.0
        assert results[("2d-rate", "nrms")] > 50.0

    def test_no_shift_without_a_range(self):
        a, b = self.station(1.0), self.station(1.1)
        results = {
            (q, m): v for q, m, v in ch.compare_stations(a, b, "fault", 10.0, 0.0)
        }
        assert results[("2d-rate", "tshift")] == 0.0
        assert results[("2d-rate", "nrms_shifted")] == pytest.approx(
            results[("2d-rate", "nrms")]
        )

    def test_the_comparison_ends_at_the_end_time(self):
        a, b = self.station(1.0), self.station(1.0)
        late = pulse(b.time, 2.5)
        b.fields["h-slip-rate"] = b.fields["h-slip-rate"] + late
        full = {(q, m): v for q, m, v in ch.compare_stations(a, b, "fault", 10.0, 0.0)}
        early = {(q, m): v for q, m, v in ch.compare_stations(a, b, "fault", 2.0, 0.0)}
        assert full[("2d-rate", "nrms")] > 50.0
        assert early[("2d-rate", "nrms")] < 1e-6

    def test_signs_apply_to_the_first_station(self):
        a, b = self.station(1.0), self.station(1.0)
        b.fields["v-slip-rate"] = -b.fields["v-slip-rate"]
        flipped = {"v-slip-rate": -1.0}
        results = {
            (q, m): v
            for q, m, v in ch.compare_stations(a, b, "fault", 10.0, 0.0, flipped)
        }
        assert results[("2d-rate", "nrms")] < 1e-9


class TestRupture:
    def gridded(self, delay, region=(-6000.0, 6000.0, 0.0, 12000.0)):
        vertices, cells = fault_mesh()
        centroids = vertices[cells].mean(axis=1)
        solution = ch.Solution(id=f"d{delay}", label="")
        solution.rupture = (
            centroids[:, 0],
            -centroids[:, 2],
            rupture_times(vertices, cells, delay),
        )
        return ch.grid_fault(solution, ch.Grid.over(region, 250.0))

    def test_a_delay_shows_as_rms_and_bias(self):
        results = {
            m: v
            for _, m, v in ch.compare_rupture(
                self.gridded(0.0), self.gridded(0.01), 10.0
            )
        }
        assert results["rms_ms"] == pytest.approx(10.0, rel=1e-6)
        assert results["bias_ms"] == pytest.approx(-10.0, rel=1e-6)
        assert results["area_mismatch"] == 0.0

    def test_only_what_ruptured_by_then_counts(self):
        early, late = self.gridded(0.0), self.gridded(0.5)
        results = {m: v for _, m, v in ch.compare_rupture(early, late, 1.0)}
        assert 0.0 < results["area_mismatch"] < 1.0
        assert results["rms_ms"] == pytest.approx(500.0, rel=1e-6)

    def test_points_that_did_not_rupture_are_left_out(self):
        vertices, cells = fault_mesh()
        centroids = vertices[cells].mean(axis=1)
        rt = rupture_times(vertices, cells)
        rt[centroids[:, 0] > 0] = 0.0
        solution = ch.Solution(id="half", label="")
        solution.rupture = (centroids[:, 0], -centroids[:, 2], rt)
        grid = ch.Grid.over((-6000.0, 6000.0, 0.0, 12000.0), 250.0)
        gridded = ch.grid_fault(solution, grid)
        assert np.all(np.isnan(gridded.rupture[grid.along > 500.0]))
        assert np.any(np.isfinite(gridded.rupture[grid.along < -500.0]))


# ============================================================================
# reading the outputs of SeisSol
# ============================================================================


class TestReadingRuns:
    def test_text_receivers_of_the_fortran_writer(self, tmp_path):
        time = np.linspace(0.0, 1.0, 11)
        rows = np.column_stack(
            [time, time, 0 * time, 40e6 + 0 * time, -120e6 + 0 * time, 2 * time]
        )
        write_text_receiver(
            tmp_path / "tpv-faultreceiver-00001-00000.dat",
            ["Time", "SRs", "SRd", "Ts0", "Pn0", "ASI"],
            (-9000.0, 0.0, -7500.0),
            rows,
            fortran=True,
        )
        stations = ch.read_fault_receivers(
            str(tmp_path), "tpv", ch.Axes.from_spec(None)
        )
        station = stations["faultst-090dp075"]
        assert station.fields["h-shear-stress"][0] == pytest.approx(40.0)
        assert station.fields["n-stress"][0] == pytest.approx(-120.0)
        assert "h-slip" not in station.fields

    def test_hdf5_receivers(self, tmp_path):
        times = np.linspace(0.0, 1.0, 11)
        write_hdf5_receivers(
            tmp_path / "tpv-receivers.h5",
            "receivers",
            "PointId",
            [(12000.0, 6000.0, -0.1), (12000.0, -6000.0, -0.1)],
            times,
            {
                **{
                    name: [np.zeros_like(times)] * 2
                    for name in ("s_xx", "s_yy", "s_zz")
                },
                "v1": [times, -times],
                "v2": [2 * times, 2 * times],
                "v3": [np.zeros_like(times)] * 2,
            },
        )
        stations = ch.read_receivers(str(tmp_path), "tpv", ch.Axes.from_spec(None))
        assert sorted(stations) == ["body-060st120dp000", "body060st120dp000"]
        station = stations["body060st120dp000"]
        assert station.fields["h-vel"] == pytest.approx(times)
        assert station.fields["n-vel"] == pytest.approx(2 * times)
        # displacement is the integral of the velocity
        assert station.fields["h-disp"][-1] == pytest.approx(0.5, rel=1e-2)

    def test_fault_output(self, tmp_path):
        pytest.importorskip("seissolxdmf")
        vertices, cells = fault_mesh()
        rt = rupture_times(vertices, cells)
        write_fault_output(
            tmp_path,
            "tpv",
            vertices,
            cells,
            [0.0, 1.5],
            {"RT": np.vstack([0 * rt, rt]), "ASl": np.vstack([0 * rt, 2 * rt])},
        )
        time, (strike, depth, values) = ch.read_fault_output(
            str(tmp_path / "tpv-fault.xdmf"), ch.Axes.from_spec(None)
        )
        assert time == 1.5
        assert depth.min() > 0.0
        assert values["rt"] == pytest.approx(rt)
        assert values["slip"] == pytest.approx(2 * rt)
        assert "psr" not in values

    def test_energy_output_keeps_fractional_times(self):
        import pandas as pd

        frame = pd.DataFrame(
            {
                "time": [0.0, 0.05, 0.1],
                "variable": ["seismic_moment"] * 3,
                "simulation_index": [0, 0, 0],
                "measurement": [0.0, 1.0, 3.0],
            }
        )
        table = ch.energy_table(frame)
        assert table["time"].tolist() == [0.0, 0.05, 0.1]
        assert table["seismic_moment"].tolist() == [0.0, 1.0, 3.0]

    def test_output_given_as_prefix_or_directory(self, tmp_path):
        (tmp_path / "output").mkdir()
        assert ch.split_output("output/tpv", None, str(tmp_path)) == (
            str(tmp_path / "output"),
            "tpv",
        )
        assert ch.split_output("output", "run", str(tmp_path)) == (
            str(tmp_path / "output"),
            "run",
        )


# ============================================================================
# reading the uploads to the CVWS
# ============================================================================


class TestReadingUploads:
    def test_table_with_comments_names_and_fortran_numbers(self, tmp_path):
        path = tmp_path / "faultst000dp075.txt"
        path.write_text(
            "# problem=TPV104\n# author=x\nt h-slip psi\n#\n"
            "0.0 0.0 0.5\n0.1D+00 1.0D-01 0.4\nbroken line\n0.2 0.2\n"
        )
        names, rows, header = ch.read_cvws_table(str(path))
        assert names == ["t", "h-slip", "psi"]
        assert rows.shape == (2, 3)
        assert rows[1, 1] == pytest.approx(0.1)
        assert header[0] == "# problem=TPV104"

    def test_codes_stations_and_rupture_data(self, tmp_path):
        write_upload(tmp_path, "alpha")
        (tmp_path / "beta_faultst000dp030.txt").write_text("t h-slip\n0 0\n1 1\n")
        (tmp_path / "notes.txt").write_text("nothing to see\n")
        codes, skipped = ch.load_scec(str(tmp_path))
        by_code = {code.id: code for code in codes}
        assert sorted(by_code) == ["alpha", "beta"]
        alpha = by_code["alpha"]
        assert sorted(alpha.stations["fault"]) == ["faultst000dp030", "faultst030dp060"]
        assert sorted(alpha.stations["body"]) == ["body020st000dp000"]
        assert alpha.rupture is not None and len(alpha.rupture[0]) > 100
        assert list(by_code["beta"].stations["fault"]) == ["faultst000dp030"]
        assert any("notes.txt" in entry for entry in skipped)


# ============================================================================
# the whole script
# ============================================================================


def test_end_to_end(tmp_path):
    pytest.importorskip("seissolxdmf")
    write_run(tmp_path / "runs" / "old", delay=0.02, moment_scale=1.01)
    write_run(tmp_path / "runs" / "new")
    write_upload(tmp_path / "scec", "alpha", delay=0.01)
    write_upload(tmp_path / "scec", "beta", delay=-0.01)
    manifest = {
        "title": "test",
        "scenario": "tpv104",
        "region": [-6000.0, 6000.0, 0.0, 12000.0],
        "scec": "scec",
        "runs": [
            {"id": "old", "output": "runs/old/tpv", "commit": "abc"},
            {"id": "new", "output": "runs/new/tpv"},
            {"id": "gone", "output": "runs/gone/tpv"},
        ],
    }
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    out = tmp_path / "out"
    code = ch.main(
        [str(tmp_path / "manifest.json"), "--out", str(out), "--grid-spacing", "250"]
        + ["--scec-pairs"]
    )
    assert code == 0
    summary = json.loads((out / "summary.json").read_text())
    # with no data, "gone" is not the reference, and takes part in no comparison
    assert summary["reference"] == "new"
    comparisons = {
        (c["run"], c["against"], c["kind"]): c["summary"]
        for c in summary["comparisons"]
    }
    assert not any(run == "gone" for run, _, _ in comparisons)

    jump = comparisons[("new", "old", "previous")]
    assert jump["rt_bias_ms"] == pytest.approx(-20.0, rel=0.02)
    assert jump["rt_rms_ms"] == pytest.approx(20.0, rel=0.02)
    old_moment = 1.01 * (2.0 - 0.02) ** 2
    assert jump["m0_rel"] == pytest.approx(abs(4.0 - old_moment) / old_moment)
    assert jump["fault_rate_max"] > 1.0

    against_upload = comparisons[("new", "alpha", "scec")]
    assert against_upload["rt_bias_ms"] == pytest.approx(-10.0, rel=0.05)
    # after the shift that fits best, the slip rates agree
    assert against_upload["fault_rate_max"] < 1.0
    # the uploads have the normal stress positive in compression
    assert summary["signs"]["n-stress"] == -1.0
    assert summary["signs"]["psi"] == 1.0
    assert ("alpha", "beta", "scec-pair") in comparisons

    report = (out / "report.md").read_text()
    assert "## Jumps" in report and "`abc`" not in report and "abc" in report
    metrics = (out / "metrics.csv").read_text().splitlines()
    assert metrics[0] == "run,against,kind,category,station,quantity,metric,value"
    if importlib.util.find_spec("matplotlib") is not None:
        assert (out / "jumps.png").is_file()
        assert (out / "rupture-time.png").is_file()


def test_an_unknown_reference_fails(tmp_path):
    manifest = {"reference": "nope", "runs": [{"id": "a", "output": "a/tpv"}]}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    assert ch.main([str(tmp_path / "manifest.json"), "--out", str(tmp_path / "o")]) == 2


def test_runs_need_ids_of_their_own(tmp_path):
    manifest = {"runs": [{"id": "a", "output": "x"}, {"id": "a", "output": "y"}]}
    (tmp_path / "manifest.json").write_text(json.dumps(manifest))
    assert ch.main([str(tmp_path / "manifest.json"), "--out", str(tmp_path / "o")]) == 2
