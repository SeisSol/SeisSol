# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The values of a trace: a graph of the operations of the core format of sderiv, and the
values of the point that build it, which NumPy computes with as with numbers.

A function of NumPy maps to the operations of the core format that give its bits: np.minimum
returns the second argument where both are equal and a NaN where either is one, so it is a
select of comparisons rather than the operation min, which is std::fmin. What the core format
cannot give to the bit (np.fmod, np.round, np.hypot, ...) is refused with the line of the
model.
"""

from __future__ import annotations

import math
import os
import struct
import sys
from dataclasses import dataclass

import numpy as np

_PACKAGE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "")
_NUMPY = os.path.join(os.path.dirname(os.path.abspath(np.__file__)), "")


class TracingError(Exception):
    """A model the tracer cannot translate: the message, with the file and the line of the
    model that does what the program cannot."""

    def __init__(self, message, at=None):
        self.file, self.line, self.column = at if at is not None else (None, None, None)
        self.message = message
        prefix = f"{self.file}:{self.line}: " if self.file is not None else ""
        super().__init__(prefix + message)


_COLUMNS = {}


def _column(frame):
    """The column of the instruction `frame` executes, from 1; 1 where Python does not say."""
    code = frame.f_code
    positions = _COLUMNS.get(code)
    if positions is None:
        positions = list(code.co_positions()) if hasattr(code, "co_positions") else []
        _COLUMNS[code] = positions
    index = frame.f_lasti // 2
    if 0 <= index < len(positions) and positions[index][2] is not None:
        return positions[index][2] + 1
    return 1


def where():
    """The file, the line and the column of the model that runs: the first frame outside of this
    package, of NumPy and of Python's own modules."""
    frame = sys._getframe(1)
    while frame is not None:
        filename = frame.f_code.co_filename
        if not (
            filename.startswith(_PACKAGE)
            or filename.startswith(_NUMPY)
            or filename.startswith("<frozen")
        ):
            return (filename, frame.f_lineno, _column(frame))
        frame = frame.f_back
    return None


@dataclass(frozen=True)
class Node:
    """An operation of the core format: `op` its name ("add", "sqrt", "select", ...), or "const",
    "input", "state", "sample"; `args` the nodes it reads."""

    op: str
    args: tuple = ()
    value: float = 0.0
    name: str = ""
    grid: int = -1


class Graph:
    """The operations a trace made, each once: a node is the same node wherever it is made again.
    Where a node was made first is its position, which the core file records (`meta line`).
    """

    def __init__(self):
        self.nodes = []
        self.index = {}
        self.positions = {}
        self.grids = []
        self.checks = []
        self.depends = []

    def add(self, node, at=None):
        key = ("const", struct.pack("<d", node.value)) if node.op == "const" else node
        found = self.index.get(key)
        if found is not None:
            return found
        self.nodes.append(node)
        self.index[key] = len(self.nodes) - 1
        if at is not None:
            self.positions[len(self.nodes) - 1] = at
        return len(self.nodes) - 1

    def constant(self, value):
        return self.add(Node("const", value=float(value)))

    def input(self, name):
        return self.add(Node("input", name=name))

    def state(self, name):
        return self.add(Node("state", name=name))

    def op(self, op, *args, at=None):
        return self.add(Node(op, tuple(args)), at=at if at is not None else where())

    def sample(self, grid, component, coordinates):
        return self.add(
            Node("sample", tuple(coordinates), name=component, grid=grid), at=where()
        )

    def cell(self, grid, component):
        """The component `component` of the data on the mesh `grid` (an index of grid()): the value
        of the cell of the point, an input channel the consumer binds."""
        return self.add(Node("cell", name=component, grid=grid))

    def grid(self, declaration):
        """The index of the dataset `declaration` (a Dataset), declared once per graph."""
        for index, known in enumerate(self.grids):
            if known is declaration:
                return index
        self.grids.append(declaration)
        return len(self.grids) - 1


def _number(value):
    """`value` as a float, where it is a number of Python or NumPy; else None."""
    if isinstance(value, (bool, int, float, np.integer, np.floating, np.bool_)):
        return float(value)
    if isinstance(value, np.ndarray) and value.ndim == 0 and value.dtype.kind in "biuf":
        return float(value)
    return None


class Value:
    """A value of the point, while a model is traced: what NumPy and the operators of Python
    compute with as with a number, making the operations of the program. A branch on it, or a
    conversion into a number of Python, fails with the line of the model."""

    __slots__ = ("graph", "node")
    __array_priority__ = 1000

    def __init__(self, graph, node):
        self.graph = graph
        self.node = node

    # ------------------------------------------------------------------ the node of an operand

    def _of(self, other):
        """The node of the operand `other` in the graph of this value; None for an array."""
        if isinstance(other, Value):
            if other.graph is not self.graph:
                raise TracingError("a value of another trace", where())
            return other.node
        value = _number(other)
        if value is None:
            return None
        return self.graph.constant(value)

    def _binary(self, other, build, reflected=False):
        node = self._of(other)
        if node is None:
            return NotImplemented
        a, b = (node, self.node) if reflected else (self.node, node)
        return Value(self.graph, build(self.graph, a, b))

    # ------------------------------------------------------------------- arithmetic

    def __add__(self, other):
        return self._binary(other, lambda g, a, b: g.op("add", a, b))

    def __radd__(self, other):
        return self._binary(other, lambda g, a, b: g.op("add", a, b), True)

    def __sub__(self, other):
        return self._binary(other, lambda g, a, b: g.op("sub", a, b))

    def __rsub__(self, other):
        return self._binary(other, lambda g, a, b: g.op("sub", a, b), True)

    def __mul__(self, other):
        return self._binary(other, lambda g, a, b: g.op("mul", a, b))

    def __rmul__(self, other):
        return self._binary(other, lambda g, a, b: g.op("mul", a, b), True)

    def __truediv__(self, other):
        return self._binary(other, lambda g, a, b: g.op("div", a, b))

    def __rtruediv__(self, other):
        return self._binary(other, lambda g, a, b: g.op("div", a, b), True)

    def __pow__(self, other):
        # `**`, as NumPy computes it for an array: the exponents 2, 0.5, -1, 1 and 0 as a square,
        # a square root, a reciprocal, the value itself and 1, others by pow
        exponent = _number(other)
        if exponent is not None:
            return Value(self.graph, _power_of(self.graph, self.node, exponent))
        return self._binary(other, lambda g, a, b: g.op("pow", a, b))

    def __rpow__(self, other):
        return self._binary(other, lambda g, a, b: g.op("pow", a, b), True)

    def __mod__(self, other):
        return self._binary(other, _remainder)

    def __rmod__(self, other):
        return self._binary(other, _remainder, True)

    def __floordiv__(self, other):
        raise TracingError(
            "`//` has no exact counterpart in the program; write np.floor(a / b) where "
            "that is meant",
            where(),
        )

    __rfloordiv__ = __floordiv__

    def __divmod__(self, other):
        raise TracingError("divmod() has no counterpart in the program", where())

    __rdivmod__ = __divmod__

    def __neg__(self):
        return Value(self.graph, self.graph.op("neg", self.node))

    def __pos__(self):
        return self

    def __abs__(self):
        return Value(self.graph, self.graph.op("abs", self.node))

    # ------------------------------------------------------------- comparisons and logic

    def __lt__(self, other):
        return self._binary(other, lambda g, a, b: g.op("lt", a, b))

    def __le__(self, other):
        return self._binary(other, lambda g, a, b: g.op("le", a, b))

    def __gt__(self, other):
        return self._binary(other, lambda g, a, b: g.op("lt", b, a))

    def __ge__(self, other):
        return self._binary(other, lambda g, a, b: g.op("le", b, a))

    def __eq__(self, other):
        return self._binary(other, lambda g, a, b: g.op("eq", a, b))

    def __ne__(self, other):
        return self._binary(other, _not_equal)

    def __and__(self, other):
        return self._binary(other, lambda g, a, b: g.op("land", a, b))

    def __rand__(self, other):
        return self._binary(other, lambda g, a, b: g.op("land", a, b), True)

    def __or__(self, other):
        return self._binary(other, lambda g, a, b: g.op("lor", a, b))

    def __ror__(self, other):
        return self._binary(other, lambda g, a, b: g.op("lor", a, b), True)

    def __xor__(self, other):
        return self._binary(other, _xor)

    def __rxor__(self, other):
        return self._binary(other, _xor, True)

    def __invert__(self):
        return Value(self.graph, _not(self.graph, self.node))

    __hash__ = object.__hash__

    # ---------------------------------------------------------- what a program cannot do

    def __bool__(self):
        raise TracingError(
            "a branch on a value of the point that the trace cannot take both ways: a while "
            "loop on it, bool() of it, or an if, and, or, not, min() or max() in a function "
            "that is not converted (the model and the functions defined in it are, another "
            "function needs @sd.function); np.where(condition, a, b), np.minimum and "
            "np.maximum, and & | ~ choose at every point",
            where(),
        )

    def _number(self, *args):
        raise TracingError(
            "a value of the point as a number of Python (math.*, float(), int(), round(), "
            "range(), an index): use the functions of NumPy on it, which the program has",
            where(),
        )

    __float__ = __int__ = __index__ = __complex__ = __round__ = __trunc__ = _number
    __floor__ = __ceil__ = _number

    def __repr__(self):
        node = self.graph.nodes[self.node]
        return f"<value of the point: {node.op} {node.name}>".replace(" >", ">")

    # ---------------------------------------------------------------- NumPy

    def __array_ufunc__(self, ufunc, method, *inputs, **kwargs):
        if method != "__call__" or kwargs:
            raise TracingError(
                f"np.{ufunc.__name__}.{method} with {sorted(kwargs)}: only the plain call of a "
                "function of NumPy maps to the program",
                where(),
            )
        build = _UFUNCS.get(ufunc.__name__)
        if build is None:
            raise TracingError(_refusal(ufunc.__name__), where())
        if any(isinstance(x, np.ndarray) and x.ndim > 0 for x in inputs):
            # an array of values: the function element by element, an array of objects
            return _elementwise(lambda *xs: _apply(ufunc, build, xs), inputs)
        return _apply(ufunc, build, inputs)

    def __array_function__(self, func, types, args, kwargs):
        build = _FUNCTIONS.get(func.__name__)
        if build is None:
            raise TracingError(_refusal(func.__name__), where())
        return build(*args, **kwargs)

    def conjugate(self):
        return self


# ---------------------------------------------------------------------- operations of NumPy


def _not(g, a):
    return g.op("eq", a, g.constant(0.0))


def _isnan(g, a):
    return _not(g, g.op("eq", a, a))


def _not_equal(g, a, b):
    return _not(g, g.op("eq", a, b))


def _xor(g, a, b):
    return _not_equal(g, _not(g, a), _not(g, b))


def _minimum(g, a, b):
    # np.minimum: a where it is below b or NaN, else b -- also where both are equal
    return g.op("select", g.op("lor", g.op("lt", a, b), _isnan(g, a)), a, b)


def _maximum(g, a, b):
    return g.op("select", g.op("lor", g.op("lt", b, a), _isnan(g, a)), a, b)


def _fmin(g, a, b):
    # np.fmin: a where it is below b or b is NaN, else b -- also where both are equal
    return g.op("select", g.op("lor", g.op("lt", a, b), _isnan(g, b)), a, b)


def _fmax(g, a, b):
    return g.op("select", g.op("lor", g.op("lt", b, a), _isnan(g, b)), a, b)


def _sign(g, a):
    # np.sign: NaN where the value is NaN, which sign gives as 0
    return g.op("select", g.op("eq", a, a), g.op("sign", a), a)


def _trunc(g, a):
    return g.op(
        "select", g.op("lt", a, g.constant(0.0)), g.op("ceil", a), g.op("floor", a)
    )


def _remainder(g, a, b):
    # np.remainder and `%`: the sign of the divisor, as mod, and a zero of that sign too
    m = g.op("mod", a, b)
    zero = g.op(
        "select", g.op("lt", b, g.constant(0.0)), g.constant(-0.0), g.constant(0.0)
    )
    return g.op("select", g.op("eq", m, g.constant(0.0)), zero, m)


def _power_of(g, a, exponent):
    if exponent == 2.0:
        return g.op("mul", a, a)
    if exponent == 0.5:
        return g.op("sqrt", a)
    if exponent == -1.0:
        return g.op("rcp", a)
    if exponent == 1.0:
        return a
    if exponent == 0.0:
        return g.constant(1.0)
    return g.op("pow", a, g.constant(exponent))


def _clip(g, a, low, high):
    return _minimum(g, _maximum(g, a, low), high)


def _unary(op):
    return lambda g, a: g.op(op, a)


def _binary(op):
    return lambda g, a, b: g.op(op, a, b)


_UFUNCS = {
    "add": _binary("add"),
    "subtract": _binary("sub"),
    "multiply": _binary("mul"),
    "divide": _binary("div"),
    "true_divide": _binary("div"),
    "negative": _unary("neg"),
    "positive": lambda g, a: a,
    "absolute": _unary("abs"),
    "fabs": _unary("abs"),
    "sqrt": _unary("sqrt"),
    "square": lambda g, a: g.op("mul", a, a),
    "reciprocal": _unary("rcp"),
    "floor": _unary("floor"),
    "ceil": _unary("ceil"),
    "trunc": _trunc,
    "sign": _sign,
    "exp": _unary("exp"),
    "log": _unary("log"),
    "log2": _unary("log2"),
    "log10": _unary("log10"),
    "sin": _unary("sin"),
    "cos": _unary("cos"),
    "tan": _unary("tan"),
    "arcsin": _unary("asin"),
    "arccos": _unary("acos"),
    "arctan": _unary("atan"),
    "sinh": _unary("sinh"),
    "cosh": _unary("cosh"),
    "tanh": _unary("tanh"),
    "arcsinh": _unary("asinh"),
    "arccosh": _unary("acosh"),
    "arctanh": _unary("atanh"),
    "power": _binary("pow"),
    "float_power": _binary("pow"),
    "arctan2": _binary("atan2"),
    "minimum": _minimum,
    "maximum": _maximum,
    "fmin": _fmin,
    "fmax": _fmax,
    "remainder": _remainder,
    "less": _binary("lt"),
    "less_equal": _binary("le"),
    "greater": lambda g, a, b: g.op("lt", b, a),
    "greater_equal": lambda g, a, b: g.op("le", b, a),
    "equal": _binary("eq"),
    "not_equal": _not_equal,
    "logical_and": _binary("land"),
    "logical_or": _binary("lor"),
    "logical_xor": _xor,
    "logical_not": _not,
    "isnan": _isnan,
    "isinf": lambda g, a: g.op("eq", g.op("abs", a), g.constant(math.inf)),
    "isfinite": lambda g, a: g.op("lt", g.op("abs", a), g.constant(math.inf)),
    "deg2rad": lambda g, a: g.op("mul", a, g.constant(math.pi / 180.0)),
    "radians": lambda g, a: g.op("mul", a, g.constant(math.pi / 180.0)),
    "rad2deg": lambda g, a: g.op("mul", a, g.constant(180.0 / math.pi)),
    "degrees": lambda g, a: g.op("mul", a, g.constant(180.0 / math.pi)),
    "clip": _clip,
    "conjugate": lambda g, a: a,
}

_ROUNDING = (
    "it rounds halfway cases to even, which the program does not; write np.floor(x + 0.5) "
    "where that is meant"
)

_REFUSED = {
    "fmod": "it has the sign of the dividend, which mod (the sign of the divisor) does not give to "
    "the bit; np.remainder and `%` have the sign of the divisor",
    "floor_divide": "write np.floor(a / b) where that is meant",
    "rint": _ROUNDING,
    "round": _ROUNDING,
    "around": _ROUNDING,
    "hypot": "write np.sqrt(x * x + y * y), which does not guard against an overflow as hypot does",
    "percentile": "use sd.percentile(values, p)",
    "median": "use sd.median(values)",
}


def _refusal(name):
    instead = _REFUSED.get(name)
    return f"np.{name} has no counterpart in the program" + (
        f"; {instead}" if instead else ""
    )


def _elementwise(function, inputs):
    """`function` of the elements of `inputs`, arrays and values, an array of objects: the values
    in arrays of their own, so that NumPy does not hand the call back to them."""
    wrapped = [np.array(x, dtype=object) if isinstance(x, Value) else x for x in inputs]
    return np.frompyfunc(function, len(inputs), 1)(*wrapped)


def _graph_of(values):
    for value in values:
        if isinstance(value, Value):
            return value.graph
    return None


def _apply(ufunc, build, inputs):
    """`ufunc` of `inputs`, numbers or values: a value where one of them is, else the number."""
    graph = _graph_of(inputs)
    if graph is None:
        return ufunc(*inputs)
    nodes = []
    for x in inputs:
        if isinstance(x, Value):
            if x.graph is not graph:
                raise TracingError("a value of another trace", where())
            nodes.append(x.node)
            continue
        number = _number(x)
        if number is None:
            raise TracingError(
                f"np.{ufunc.__name__} of {type(x).__name__}, which is no number",
                where(),
            )
        nodes.append(graph.constant(number))
    return Value(graph, build(graph, *nodes))


def _nodes_of(graph, values):
    nodes = []
    for x in values:
        if isinstance(x, Value):
            nodes.append(x.node)
            continue
        number = _number(x)
        if number is None:
            raise TracingError(f"{type(x).__name__} is no number", where())
        nodes.append(graph.constant(number))
    return nodes


def _where(condition, x=None, y=None):
    if x is None or y is None:
        raise TracingError("np.where of a condition alone has no counterpart", where())
    if any(isinstance(v, np.ndarray) and v.ndim > 0 for v in (condition, x, y)):
        return _elementwise(_where, (condition, x, y))
    graph = _graph_of((condition, x, y))
    if graph is None:
        return np.where(condition, x, y)
    c, a, b = _nodes_of(graph, (condition, x, y))
    return Value(graph, graph.op("select", c, a, b))


def _select(condlist, choicelist, default=0.0):
    # the first condition that holds, from the last on nested into the default
    result = default
    for condition, choice in reversed(list(zip(condlist, choicelist))):
        result = _where(condition, choice, result)
    return result


def _clip_function(a, a_min=None, a_max=None, **kwargs):
    if kwargs:
        raise TracingError(f"np.clip with {sorted(kwargs)}", where())
    graph = _graph_of((a, a_min, a_max))
    if graph is None:
        return np.clip(a, a_min, a_max)
    value = a
    if a_min is not None:
        value = np.maximum(value, a_min)
    if a_max is not None:
        value = np.minimum(value, a_max)
    return value


_FUNCTIONS = {
    "where": _where,
    "select": _select,
    "clip": _clip_function,
}


def _object_method(name):
    def method(self):
        return Value(self.graph, _UFUNCS[name](self.graph, self.node))

    method.__name__ = name
    return method


# NumPy calls the method of the name of a function on every element of an array of objects
for _name in (
    "sqrt",
    "exp",
    "log",
    "log2",
    "log10",
    "sin",
    "cos",
    "tan",
    "arcsin",
    "arccos",
    "arctan",
    "sinh",
    "cosh",
    "tanh",
    "arcsinh",
    "arccosh",
    "arctanh",
    "floor",
    "ceil",
    "trunc",
    "deg2rad",
    "rad2deg",
):
    setattr(Value, _name, _object_method(_name))
del _name
