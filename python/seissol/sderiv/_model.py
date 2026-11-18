# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""A model: a function of the inputs of a point that returns its outputs by name, which runs on
NumPy arrays and is traced into a program of sderiv."""

from __future__ import annotations

import contextlib
import functools
import inspect
import math
import os
import sys
import threading
from dataclasses import dataclass, field

import numpy as np

from . import _branches
from ._branches import Undefined
from ._convert import convert
from ._emit import core_text
from ._graph import Graph, TracingError, Value, _number, where

VERSION = "0.1.0"
GENERATOR = f"seissol-python {VERSION}"

_LOCAL = threading.local()


# ------------------------------------------------------------------ the files a model reads


def _audit(event, args):
    """Records a file opened for reading while a trace records them: a constant of the model."""
    if event != "open":
        return
    files = getattr(_LOCAL, "files", None)
    if files is None:
        return
    path, mode, flags = (tuple(args) + (None, None, None))[:3]
    if path is None or isinstance(path, int):
        return
    if isinstance(mode, str):
        reading = not any(c in mode for c in "wax+")
    else:
        reading = flags is None or (flags & (os.O_WRONLY | os.O_RDWR)) == 0
    if reading:
        files.append(os.path.abspath(os.fsdecode(path)))


sys.addaudithook(_audit)

_IGNORED = tuple(
    os.path.join(os.path.abspath(p), "")
    for p in {sys.prefix, sys.base_prefix, sys.exec_prefix, "/proc", "/sys", "/dev"}
)


def _data_files(files, source=None):
    """The files of `files` that are data of the model: not the source, not a module of Python,
    not a file of the system."""
    kept = []
    for path in files:
        if path == source or path.endswith((".py", ".pyc", ".pyd", ".so")):
            continue
        if (
            path.startswith(_IGNORED)
            or "site-packages" in path
            or not os.path.isfile(path)
        ):
            continue
        if path not in kept:
            kept.append(path)
    return kept


@contextlib.contextmanager
def recording():
    """The files Python opens for reading in this block, absolute: the data a model reads while it
    is imported or traced, which a core file records with their hashes (`meta depends`). Files
    that a C extension opens itself (netCDF4, h5py) are not seen; the model declares those with
    depends()."""
    files = []
    previous = getattr(_LOCAL, "files", None)
    _LOCAL.files = files
    try:
        yield files
    finally:
        _LOCAL.files = previous
        if previous is not None:
            previous.extend(files)


@contextlib.contextmanager
def _tracing(graph):
    previous = getattr(_LOCAL, "graph", None)
    _LOCAL.graph = graph
    try:
        yield
    finally:
        _LOCAL.graph = previous


def _active():
    return getattr(_LOCAL, "graph", None)


def depends(*paths):
    """Declares files the model reads that Python does not open itself (through netCDF4 or
    h5py, say): a core file records them with their hashes, so that SeisSol notices when they
    change. Outside of a trace, it does nothing."""
    graph = _active()
    if graph is not None:
        graph.depends.extend(os.path.abspath(os.fspath(p)) for p in paths)


# ------------------------------------------------------------------ functions of the program


def check(condition, message):
    """A condition every point meets: SeisSol stops at the first point where it is zero (a NaN
    counts as true), with the message. Outside of a trace, ValueError at such a point. In a branch
    on a value of the point, it applies to the points of the branch."""
    _check(condition, str(message), ValueError)


def _check(condition, message, error):
    """check(), with `error` outside of a trace."""
    if isinstance(condition, np.ndarray) and condition.dtype == object:
        # a vector of values: a check of each
        for element in condition.flat:
            _check(element, message, error)
        return
    number = None if isinstance(condition, Value) else _number(condition)
    if number is not None and number != 0.0:
        return
    current = _branches.path()
    if current is not None:
        condition = np.logical_or(np.logical_not(current), condition)
    if isinstance(condition, Value):
        condition.graph.checks.append((condition.node, message, where()))
        return
    if _active() is not None:
        if _number(condition) == 0.0:
            raise TracingError(f"the check fails at every point: {message}", where())
        return
    values = np.asarray(condition, dtype=float)
    failing = np.flatnonzero(values.ravel() == 0.0)
    if failing.size:
        raise error(f"{message} -- point {int(failing[0])} does not meet the check")


def percentile(values, p):
    """The percentile p (0 to 100) of `values` at every point, as numpy.percentile (its method
    linear) gives it: an operation of the program, of 1 to 1024 values."""
    values = list(values)
    graph = next((v.graph for v in values + [p] if isinstance(v, Value)), None)
    if graph is None:
        stacked = np.stack([np.asarray(v, dtype=float) for v in values])
        return np.percentile(stacked, p, axis=0)
    if not 1 <= len(values) <= 1024:
        raise TracingError(
            f"a percentile of {len(values)} values; it takes 1 to 1024", where()
        )
    nodes = [
        v.node if isinstance(v, Value) else graph.constant(float(v))
        for v in [p] + values
    ]
    return Value(graph, graph.op("percentile", *nodes))


def median(values):
    """The median of `values` at every point: percentile(values, 50)."""
    return percentile(values, 50.0)


def erf(x):
    """The error function, as C's erf."""
    if isinstance(x, Value):
        return Value(x.graph, x.graph.op("erf", x.node))
    if np.ndim(x) == 0:
        return math.erf(float(x))
    return np.vectorize(math.erf, otypes=[float])(x)


# ------------------------------------------------------------------------------ the model


@dataclass
class Trace:
    """What a trace of a model made: the graph, the inputs, the outputs (name, node, position),
    the states (name, initial value, node of the next value, position), the source of the model
    and the files it read."""

    graph: Graph
    inputs: list
    outputs: list
    states: list
    source: str = None
    files: list = field(default_factory=list)


class Model:
    """A function of the inputs of a point (its parameters, but those with a default, which are
    constants of the model) that returns its outputs: a dict by name, or a value or a tuple of
    values that `outputs` names. A parameter of the name of a state reads its value of the
    previous call, and the value the function returns under that name is its next one.

    An if, and, or, not, a conditional expression, a chained comparison, min(), max(), any(),
    all(), assert and raise on a value of the point take both ways, merged at every point as
    Python takes them there (the function is converted, see _convert.py).

    Called, it runs on numbers and NumPy arrays; trace() translates it into a program, core()
    into the core text SeisSol reads."""

    def __init__(self, function, *, outputs=None, state=None):
        self.function = function
        self.outputs = [outputs] if isinstance(outputs, str) else outputs
        self.state = dict(state or {})
        self._converted = None
        functools.update_wrapper(self, function)

    @property
    def converted(self):
        """The function the model runs: the function with its branches converted, or the
        function itself where it has none."""
        if self._converted is None:
            self._converted = convert(self.function) or self.function
        return self._converted

    def __call__(self, *args, **kwargs):
        result = self.converted(*args, **kwargs)
        if isinstance(result, dict):
            values = list(result.values())
        elif isinstance(result, (tuple, list)):
            values = list(result)
        else:
            values = [result]
        for value in values:
            if isinstance(value, Undefined):
                # an output without a value: as Python fails where it reads such a name
                raise value.error()
        return result

    @property
    def inputs(self):
        return [
            name
            for name, parameter in inspect.signature(self.function).parameters.items()
            if parameter.default is inspect.Parameter.empty and name not in self.state
        ]

    def trace(self):
        """The program of the model, as a Trace."""
        graph = Graph()
        positional = []
        keywords = {}
        inputs = []
        for name, parameter in inspect.signature(self.function).parameters.items():
            if parameter.kind in (parameter.VAR_POSITIONAL, parameter.VAR_KEYWORD):
                raise TracingError(
                    f"the model takes *{name}; its inputs are the names of its parameters",
                    _definition(self.function),
                )
            if parameter.default is not inspect.Parameter.empty:
                continue
            if name in self.state:
                value = Value(graph, graph.state(name))
            else:
                value = Value(graph, graph.input(name))
                inputs.append(name)
            if parameter.kind == parameter.POSITIONAL_ONLY:
                positional.append(value)
            else:
                keywords[name] = value
        function = self.converted
        with _tracing(graph), recording() as files:
            result = function(*positional, **keywords)

        definition = _definition(self.function)
        if isinstance(result, Undefined):
            raise TracingError(result.reason(), definition)
        named = self._named(result)
        nodes = {}
        for name, value in named:
            if isinstance(value, Undefined):
                raise TracingError(f"the output {name}: {value.reason()}", definition)
            if isinstance(value, np.ndarray) and value.ndim == 0:
                value = value.item()
            if isinstance(value, Value):
                if value.graph is not graph:
                    raise TracingError(f"the output {name} is of another trace")
                nodes[name] = value.node
                continue
            number = _number(value)
            if number is None:
                raise TracingError(
                    f"the output {name} is a {type(value).__name__}, not a number at every "
                    "point",
                    _definition(self.function),
                )
            nodes[name] = graph.constant(number)

        states = []
        for name, initial in self.state.items():
            if name not in nodes:
                raise TracingError(
                    f"the state {name} gets no next value: return it by name",
                    definition,
                )
            states.append(
                (name, float(initial), nodes[name], graph.positions.get(nodes[name]))
            )
        wanted = self.outputs if self.outputs is not None else [n for n, _ in named]
        outputs = []
        for name in wanted:
            if name not in nodes:
                raise TracingError(f"the model returns no output {name}", definition)
            outputs.append(
                (name, nodes[name], graph.positions.get(nodes[name], definition))
            )
        source = inspect.getsourcefile(self.function)
        source = os.path.abspath(source) if source is not None else None
        return Trace(
            graph,
            inputs,
            outputs,
            states,
            source,
            _data_files(files + graph.depends, source),
        )

    def _named(self, result):
        if isinstance(result, dict):
            for name in result:
                if not isinstance(name, str):
                    raise TracingError(
                        f"the outputs are named by strings, not by {name!r}",
                        _definition(self.function),
                    )
            return list(result.items())
        if self.outputs is None:
            raise TracingError(
                "the model returns a dict of its outputs by name, or names them: "
                "@sd.model(outputs=...)",
                _definition(self.function),
            )
        values = list(result) if isinstance(result, (tuple, list)) else [result]
        if len(values) != len(self.outputs):
            raise TracingError(
                f"the model returns {len(values)} values for the {len(self.outputs)} outputs "
                f"{', '.join(self.outputs)}",
                _definition(self.function),
            )
        return list(zip(self.outputs, values))

    def core(self, path=None, *, command=None, depends=()):
        """The core text of the model. Where `path` names the core file, with the provenance that
        SeisSol checks: the source, the files the model read (and `depends`, absolute paths), the
        `command` that writes the file again, and the hash of the body."""
        trace = self.trace()
        directory = os.path.dirname(os.path.abspath(path)) if path is not None else None
        files = _data_files(list(trace.files) + list(depends), trace.source)
        return core_text(
            trace, directory, generator=GENERATOR, command=command, depends=files
        )

    def write(self, path, *, command=None, depends=()):
        """Writes the core file `path`; returns its text."""
        text = self.core(path, command=command, depends=depends)
        with open(path, "w", encoding="utf-8") as file:
            file.write(text)
        return text

    def program(self):
        """The program, through the module sderiv of a SeisSol build: what SeisSol evaluates,
        with its evaluate() and step()."""
        import sderiv

        return sderiv.compile(self.core())

    def evaluate(self, inputs=None, **named):
        """The outputs at points by the program SeisSol evaluates (program()), as arrays."""
        return self.program().evaluate(inputs, **named)


def _definition(function):
    code = getattr(function, "__code__", None)
    if code is None:
        return None
    return (os.path.abspath(code.co_filename), code.co_firstlineno, 1)


def model(function=None, *, outputs=None, state=None):
    """A model of the function; see Model. As @model, or @model(outputs=..., state=...)."""
    if function is None:
        return lambda f: Model(f, outputs=outputs, state=state)
    return Model(function, outputs=outputs, state=state)


def function(f):
    """A function the model calls, converted as the model is: an if, and, or, not, min() and
    max() on a value of the point in it take both ways. The functions defined in the model are
    converted with it; one it calls from elsewhere needs @sd.function."""
    return convert(f) or f
