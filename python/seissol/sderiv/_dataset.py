# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The datasets a model reads: a dataset of the core format, declared by the model and sampled
as a function of its coordinates. In a trace, a sample is an operation of the program; outside of
one, it is computed by the program SeisSol would compute it with, through the module sderiv of a
SeisSol build, so that a model gives the same values in Python as in SeisSol."""

from __future__ import annotations

import os

import numpy as np

from ._graph import TracingError, Value, where


class Dataset:
    """A dataset: `kind` "asagi" (a NetCDF file in the layout of ASAGI), "table" (a text file of
    columns) or "inline" (the knots and the values given here). A relative path is relative to
    the file of the model that declares the dataset, as a path in a core file is to that file.

    ``g.rho(x, y, z)`` samples the component rho, ``g(x, y, z)`` every component, a dict.
    """

    def __init__(
        self,
        path=None,
        components=None,
        *,
        kind=None,
        variable=None,
        interpolation="linear",
        boundary="clamp",
        axes=None,
        time=None,
        coordinates=None,
        knots=None,
        values=None,
    ):
        if components is None or isinstance(components, str):
            components = [components] if isinstance(components, str) else None
        if not components:
            raise TracingError("a dataset names its components", where())
        if kind is None:
            kind = (
                "inline"
                if knots is not None
                else "table" if coordinates is not None else "asagi"
            )
        self.kind = kind
        self.components = list(components)
        self.interpolation = interpolation
        self.boundary = boundary
        self.axes = axes
        self.time = time
        self.coordinates = coordinates
        self.variable = variable if variable is not None or kind != "asagi" else "data"
        self.path = None
        if kind == "inline":
            if knots is None or values is None:
                raise TracingError(
                    "an inline dataset gives its knots and its values", where()
                )
            self.knots = [list(map(float, axis)) for axis in knots]
            data = np.asarray(values, dtype=float)
            shape = tuple(len(axis) for axis in self.knots)
            if data.shape == shape and len(self.components) == 1:
                data = data[..., np.newaxis]
            if data.shape != shape + (len(self.components),):
                raise TracingError(
                    f"the values of the inline dataset are of the shape {data.shape}, not of its "
                    f"knots and components {shape + (len(self.components),)}",
                    where(),
                )
            self.values = data
        else:
            if path is None:
                raise TracingError(
                    f"a dataset of the kind {kind} names its file", where()
                )
            self.knots = None
            self.values = None
            path = os.fspath(path)
            if not os.path.isabs(path):
                at = where()
                base = os.path.dirname(at[0]) if at is not None else os.getcwd()
                path = os.path.join(base, path)
            self.path = os.path.normpath(path)
        self._native = {}

    def sample(self, component, *coordinates):
        """The component `component` at `coordinates`."""
        if component not in self.components:
            raise TracingError(
                f"the dataset has no component {component}; it has {', '.join(self.components)}",
                where(),
            )
        graph = next((c.graph for c in coordinates if isinstance(c, Value)), None)
        if graph is not None:
            nodes = []
            for c in coordinates:
                if isinstance(c, Value):
                    nodes.append(c.node)
                else:
                    nodes.append(graph.constant(float(c)))
            return Value(graph, graph.sample(graph.grid(self), component, nodes))
        return self._sample(component, coordinates)

    def __call__(self, *coordinates):
        return {c: self.sample(c, *coordinates) for c in self.components}

    def __getattr__(self, name):
        components = self.__dict__.get("components", ())
        if name in components:
            return lambda *coordinates: self.sample(name, *coordinates)
        raise AttributeError(name)

    def _sample(self, component, coordinates):
        """A sample outside of a trace, by the program of SeisSol: that of this dataset and its
        component, through the module sderiv."""
        try:
            import sderiv
        except ImportError as error:
            raise RuntimeError(
                "a dataset is sampled outside of a trace by the program SeisSol samples it with, "
                "which needs the module sderiv of a SeisSol build with -DSDERIV_PYTHON=ON"
            ) from error
        from ._emit import spell_dataset, spell_name

        key = (component, len(coordinates))
        program = self._native.get(key)
        if program is None:
            names = [f"c{i}" for i in range(len(coordinates))]
            text = (
                "#!sderiv 1\n"
                f"in {', '.join(names)}\n"
                f"{spell_dataset('g', self, None)}\n"
                f"out value = g.{spell_name(component)}({', '.join(names)})\n"
            )
            program = sderiv.compile(text)
            self._native[key] = program
        inputs = {f"c{i}": c for i, c in enumerate(coordinates)}
        values = program.evaluate(inputs)["value"]
        if all(np.ndim(c) == 0 for c in coordinates):
            return float(values[0])
        return values


def dataset(path=None, components=None, **attributes):
    """A dataset the model reads; see Dataset."""
    return Dataset(path, components, **attributes)
