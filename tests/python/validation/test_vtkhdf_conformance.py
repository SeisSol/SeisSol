# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Read a VTKHDF output back with VTK itself.

The unit tests check the file against the specification as we read it; this checks it against the
reader that has to make sense of it. Point it at an output directory:

    pytest --vtkhdf=<dir>/<prefix>-volume.vtkhdf tests/python/validation

Skipped when the vtk module is not installed.
"""

import pytest
from vtkmodules.util.numpy_support import vtk_to_numpy
from vtkmodules.vtkCommonExecutionModel import vtkStreamingDemandDrivenPipeline as sddp

vtk_io = pytest.importorskip("vtkmodules.vtkIOHDF")


def read(path):
    reader = vtk_io.vtkHDFReader()
    reader.SetFileName(path)
    reader.UpdateInformation()
    return reader


def time_steps(reader):
    info = reader.GetOutputInformation(0)
    if not info.Has(sddp.TIME_STEPS()):
        return []
    return list(info.Get(sddp.TIME_STEPS()))


def step(reader, time=None):
    if time is not None:
        reader.GetOutputInformation(0).Set(sddp.UPDATE_TIME_STEP(), time)
    reader.Update()
    return reader.GetOutput()


def test_file_is_readable(vtkhdf_path):
    """The grid arrives with cells, points and at least one cell array, at every time step."""
    reader = read(vtkhdf_path)
    times = time_steps(reader)

    for time in times or [None]:
        grid = step(reader, time)
        assert grid.GetNumberOfCells() > 0
        assert grid.GetNumberOfPoints() > 0

        cell_data = grid.GetCellData()
        assert cell_data.GetNumberOfArrays() > 0
        for index in range(cell_data.GetNumberOfArrays()):
            array = vtk_to_numpy(cell_data.GetArray(index))
            # one tuple per cell, not one per cell and step
            assert array.shape[0] == grid.GetNumberOfCells()


def test_time_steps_are_distinct(vtkhdf_path):
    """A time series has to hand out different data per step, not the same block every time."""
    reader = read(vtkhdf_path)
    times = time_steps(reader)
    if len(times) < 2:
        pytest.skip("not a time series")

    assert times == sorted(times)
    assert len(set(times)) == len(times)

    seen = []
    for time in times:
        grid = step(reader, time)
        cell_data = grid.GetCellData()
        values = []
        for index in range(cell_data.GetNumberOfArrays()):
            values.extend(vtk_to_numpy(cell_data.GetArray(index)).ravel().tolist())
        seen.append(tuple(values))
    # constant arrays are expected to repeat; something has to change
    assert len(set(seen)) > 1, "every step returned the same values"
