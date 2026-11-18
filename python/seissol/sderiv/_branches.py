# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""What the converted branches of a model call (_convert.py): an if, the rest of a block after a
return, a break or a continue at some points, a loop the points leave, a conditional expression,
and, or, not, a chained comparison, min(), max(), any(), all(), assert, raise and the value of
the function.

A condition Python knows -- a bool, a number, a string -- takes its branch as Python does. A
condition of the point -- a value of a trace, or an array of NumPy -- runs both branches, each on
the points it holds for (the path, which a check in the branch applies to), and merges what they
leave behind where it holds (np.where): at every point, what Python computes there.
"""

from __future__ import annotations

import builtins
import contextlib
import operator
import struct
import threading
import types
from functools import reduce

import numpy as np

from ._graph import TracingError, Value, _elementwise, where

RETURN = "_sd__retval"

_STATE = threading.local()


# ------------------------------------------------------------------ names without a value


class Undefined:
    """What a name holds where it has no value: before it is assigned (unbound, `why` None), or
    after a branch on a value of the point that does not give it one value per point (`why` says
    how). Used, it fails."""

    __slots__ = ("name", "why")

    def __init__(self, name, why=None):
        self.name = name
        self.why = why

    def reason(self):
        if self.why is None:
            return f"{_display(self.name)} is used before it has a value"
        return self.why

    def error(self):
        if self.why is None:
            return UnboundLocalError(
                f"cannot access local variable '{self.name}' where it is not associated with a "
                "value"
            )
        return TracingError(self.why, where())

    def _fail(self, *args, **kwargs):
        raise self.error()

    def __getattr__(self, name):
        if name.startswith("__") and name.endswith("__"):
            raise AttributeError(name)
        self._fail()

    def __repr__(self):
        return f"<no value: {self.name}>"


for _method in (
    "__add__",
    "__radd__",
    "__sub__",
    "__rsub__",
    "__mul__",
    "__rmul__",
    "__truediv__",
    "__rtruediv__",
    "__floordiv__",
    "__rfloordiv__",
    "__mod__",
    "__rmod__",
    "__pow__",
    "__rpow__",
    "__matmul__",
    "__rmatmul__",
    "__neg__",
    "__pos__",
    "__abs__",
    "__lt__",
    "__le__",
    "__gt__",
    "__ge__",
    "__eq__",
    "__ne__",
    "__and__",
    "__rand__",
    "__or__",
    "__ror__",
    "__xor__",
    "__rxor__",
    "__invert__",
    "__bool__",
    "__float__",
    "__int__",
    "__index__",
    "__iter__",
    "__len__",
    "__getitem__",
    "__contains__",
    "__call__",
    "__array__",
    "__array_ufunc__",
    "__array_function__",
):
    setattr(Undefined, _method, Undefined._fail)
del _method
Undefined.__hash__ = object.__hash__

UNDEFINED = Undefined(RETURN)


def unbound(name):
    """The value of `name` before it is assigned."""
    return Undefined(name)


def _display(name):
    if name == RETURN:
        return "the value the function returns"
    return f"`{name}`" if name.isidentifier() else name


# ----------------------------------------------------------------------- conditions and paths


def _static(value):
    """Whether Python knows `value` as a condition: neither a value of the point nor an array of
    NumPy, which holds a value per point. A vector of values has no truth, as an array of
    NumPy has none."""
    if isinstance(value, Value):
        return False
    if isinstance(value, np.ndarray) and value.ndim > 0:
        if value.dtype == object:
            raise TracingError(
                "the truth of an array of values of the point is ambiguous, as NumPy has it: "
                "compare its elements",
                where(),
            )
        return False
    return True


_TRUTHS = frozenset({"lt", "le", "eq", "land", "lor"})


def _truthlike(value):
    """Whether `value` is a truth value of the points already: 1 where it holds, else 0."""
    if isinstance(value, Value):
        return value.graph.nodes[value.node].op in _TRUTHS
    return isinstance(value, np.ndarray) and value.dtype == bool


def _truth(condition):
    """`condition` as a truth value: 1 where it holds (where it is not 0, a NaN holds), else 0."""
    if _truthlike(condition):
        return condition
    return np.where(condition, True, False)


def _falsity(condition):
    return np.logical_not(condition)


def path():
    """The points the code runs for where a branch on a value of the point runs it: a truth
    value of the points, or None for all of them."""
    return getattr(_STATE, "path", None)


@contextlib.contextmanager
def _on(condition):
    previous = path()
    _STATE.path = condition if previous is None else np.logical_and(previous, condition)
    try:
        yield
    finally:
        _STATE.path = previous


def _any(flags):
    """Whether one of `flags` holds: True or False where Python knows, else a truth value."""
    pending = []
    for flag in flags:
        if _static(flag):
            if flag:
                return True
        else:
            pending.append(flag)
    if not pending:
        return False
    return reduce(np.logical_or, pending)


# --------------------------------------------------------------------------- the merges


def _same(a, b):
    """Whether `a` and `b` are the same value at every point."""
    if a is b:
        return True
    if isinstance(a, Value) and isinstance(b, Value):
        return a.graph is b.graph and a.node == b.node
    if type(a) is type(b):
        if type(a) is float:
            return struct.pack("<d", a) == struct.pack("<d", b)
        if type(a) in (bool, int, str, bytes):
            return a == b
    return False


def _number_like(value):
    if isinstance(value, (bool, int, float, np.number, np.bool_, Value)):
        return True
    return isinstance(value, np.ndarray) and value.dtype.kind in "biufO"


def _choose(condition, then, otherwise):
    """`then` where `condition` holds, else `otherwise`: np.where, element by element where an
    array holds values of the point."""
    if any(
        isinstance(v, np.ndarray) and v.dtype == object
        for v in (condition, then, otherwise)
    ):
        return _elementwise(_choose, (condition, then, otherwise))
    return np.where(condition, then, otherwise)


def _is(value, truth):
    return isinstance(value, (bool, np.bool_)) and bool(value) is truth


def merge(name, condition, then, otherwise):
    """What `name` holds after a branch on `condition`, a value of the point: `then` where it
    holds, else `otherwise`."""
    if _same(then, otherwise):
        return then
    if isinstance(then, Undefined) or isinstance(otherwise, Undefined):
        for side in (then, otherwise):
            if isinstance(side, Undefined) and side.why is not None:
                return side
        if isinstance(then, Undefined) and isinstance(otherwise, Undefined):
            return then
        if name == RETURN:
            # the value of the function where it returned; on the other side it has not
            return otherwise if isinstance(then, Undefined) else then
        undefined = then if isinstance(then, Undefined) else otherwise
        return Undefined(
            undefined.name,
            f"{_display(undefined.name)} has a value at some points only -- a branch on a value "
            "of the point assigns it on one side and not on the other -- and it is used after "
            "the branch; give it a value before the branch",
        )
    display = _display(name)
    if isinstance(then, dict) and isinstance(otherwise, dict):
        keys = list(then) + [key for key in otherwise if key not in then]
        return {
            key: merge(
                f"{display}[{key!r}]",
                condition,
                then.get(key, Undefined(f"{display}[{key!r}]")),
                otherwise.get(key, Undefined(f"{display}[{key!r}]")),
            )
            for key in keys
        }
    if (
        isinstance(then, (tuple, list))
        and type(then) is type(otherwise)
        and len(then) == len(otherwise)
    ):
        merged = [
            merge(f"{display}[{i}]", condition, a, b)
            for i, (a, b) in enumerate(zip(then, otherwise))
        ]
        if isinstance(then, list):
            return merged
        return type(then)(*merged) if hasattr(then, "_fields") else tuple(merged)
    if not (_number_like(then) and _number_like(otherwise)):
        kinds = (type(then).__name__, type(otherwise).__name__)
        what = (
            f"different objects of the type {kinds[0]}"
            if kinds[0] == kinds[1]
            else f"a {kinds[0]} and a {kinds[1]}"
        )
        if isinstance(then, (tuple, list)) and isinstance(otherwise, (tuple, list)):
            what = f"sequences of {len(then)} and {len(otherwise)} elements"
        return Undefined(
            name,
            f"{display} is {what} on the two sides of a branch on a value of the point, which "
            "have no value per point to choose between, and it is used after the branch",
        )
    if isinstance(condition, Value):
        if _is(then, True) and _is(otherwise, False):
            return _truth(condition)
        if _is(then, False) and _is(otherwise, True):
            return _falsity(condition)
    return _choose(condition, then, otherwise)


# ------------------------------------------------- what a branch changes in place


class _Saved:
    """What a container held when it was saved, to see whether a branch changed it and to put
    it back: the items of a dict or a list (saved in turn), the elements of an array of NumPy,
    the attributes of an object."""

    __slots__ = ("value", "kind", "held")

    def __init__(self, value, kind, held=None):
        self.value = value
        self.kind = kind
        self.held = held


def _has_attributes(value):
    return hasattr(value, "__dict__") and not isinstance(
        value,
        (type, types.ModuleType, types.FunctionType, types.BuiltinFunctionType),
    )


def _saved(value, changed):
    """Whether a branch saves `value`, a name it binds (or `changed`, changes in place)."""
    if isinstance(value, (dict, list)):
        return True
    return changed and (isinstance(value, np.ndarray) or _has_attributes(value))


def _save(value, top=False, seen=None):
    seen = set() if seen is None else seen
    if isinstance(value, (dict, list)) and id(value) not in seen:
        seen.add(id(value))
        items = value.items() if isinstance(value, dict) else enumerate(value)
        return _Saved(
            value,
            "dict" if isinstance(value, dict) else "list",
            [(key, _save(item, False, seen)) for key, item in items],
        )
    if isinstance(value, np.ndarray):
        return _Saved(value, "array", value.copy())
    if top and _has_attributes(value):
        return _Saved(value, "object", _save(vars(value), False, seen))
    return _Saved(value, "leaf")


def _restore(saved):
    value = saved.value
    if saved.kind == "dict":
        items = [(key, _restore(item)) for key, item in saved.held]
        value.clear()
        value.update(items)
    elif saved.kind == "list":
        value[:] = [_restore(item) for _, item in saved.held]
    elif saved.kind == "array":
        if value.flags.writeable and value.shape == saved.held.shape:
            value[...] = saved.held
    elif saved.kind == "object":
        _restore(saved.held)
    return value


def _equal(a, b):
    """Whether two savings of a value hold the same."""
    if a.value is not b.value or a.kind != b.kind:
        return False
    if a.kind == "leaf":
        return True
    if a.kind == "array":
        x, y = a.held, b.held
        if x.shape != y.shape or x.dtype != y.dtype:
            return False
        if x.dtype == object:
            return all(p is q for p, q in zip(x.flat, y.flat))
        return x.tobytes() == y.tobytes()
    if a.kind == "object":
        return _equal(a.held, b.held)
    return len(a.held) == len(b.held) and all(
        k == l and _equal(p, q) for (k, p), (l, q) in zip(a.held, b.held)
    )


def _plain(saved):
    """What a saving holds, as dicts, lists and arrays of their own: what the branches merge."""
    if saved.kind == "dict":
        return {key: _plain(item) for key, item in saved.held}
    if saved.kind == "list":
        return [_plain(item) for _, item in saved.held]
    if saved.kind == "array":
        return saved.held
    if saved.kind == "object":
        return _plain(saved.held)
    return saved.value


def _assign(target, value, top=False):
    """Writes `value`, what the branches merged, into `target` in place; returns what holds it."""
    if isinstance(target, dict) and isinstance(value, dict):
        current = dict(target)
        target.clear()
        target.update(
            {
                key: _assign(current[key], item) if key in current else item
                for key, item in value.items()
            }
        )
        return target
    if isinstance(target, list) and isinstance(value, list):
        current = list(target)
        target[:] = [
            _assign(current[i], item) if i < len(current) else item
            for i, item in enumerate(value)
        ]
        return target
    if (
        isinstance(target, np.ndarray)
        and isinstance(value, np.ndarray)
        and target.flags.writeable
        and target.shape == value.shape
        and (target.dtype == object or value.dtype != object)
    ):
        target[...] = value
        return target
    if top and _has_attributes(target) and isinstance(value, dict):
        _assign(vars(target), value)
        return target
    return value


# ------------------------------------------------------------------------- the statements


def _nothing():
    pass


def _both(condition, then, otherwise, get, set_, changed):
    """Both branches, each on its points, and what they leave behind merged."""
    before = get()
    saved = {
        name: _save(value, top=True)
        for name, value in before.items()
        if _saved(value, name in changed)
    }
    sides = []
    for branch, on in ((then, _truth), (otherwise, _falsity)):
        if sides:
            for item in saved.values():
                _restore(item)
            set_(before)
        with _on(on(condition)):
            branch()
        after = get()
        held = {
            name: _save(after[name], top=True)
            for name in saved
            if after[name] is before[name]
        }
        sides.append((after, held))
    (after_then, held_then), (after_otherwise, held_otherwise) = sides
    merged = {}
    for name, value in before.items():
        if name in held_then and name in held_otherwise:
            # changed in place, if at all, on both sides: what it holds is merged in place
            a, b = held_then[name], held_otherwise[name]
            if not (_equal(a, saved[name]) and _equal(b, saved[name])):
                _assign(value, merge(name, condition, _plain(a), _plain(b)), top=True)
            merged[name] = value
            continue
        a = _plain(held_then[name]) if name in held_then else after_then[name]
        b = (
            _plain(held_otherwise[name])
            if name in held_otherwise
            else after_otherwise[name]
        )
        merged[name] = merge(name, condition, a, b)
    set_(merged)


def if_stmt(condition, then, otherwise, get, set_, changed=()):
    """An if: the branch Python takes, or both, merged where `condition` holds. get() gives the
    names the branches bind or change, set_() binds them; `changed` names those they change in
    place."""
    if _static(condition):
        if condition:
            then()
        else:
            otherwise()
        return
    _both(condition, then, otherwise, get, set_, changed)


def unless(flags, body, get, set_, changed=()):
    """The rest of a block after a return, a break or a continue at some points: where none of
    `flags` holds."""
    condition = _any(flags)
    if _static(condition):
        if not condition:
            body()
        return
    _both(condition, _nothing, body, get, set_, changed)


def done(*flags):
    """Whether every point has left a loop: one of `flags` holds, where Python knows it."""
    condition = _any(flags)
    return _static(condition) and bool(condition)


def going(flags, test):
    """The condition of a while loop that points may leave: false where every point has left,
    else the test, which Python has to know -- a loop that goes on at some points only is no
    program."""
    if done(*flags):
        return False
    value = test()
    if _static(value):
        return bool(value)
    raise TracingError(
        "the condition of a while loop is a value of the point (after a break or a return at "
        "some points, a name it reads may be one): the loop would go on at some points only, "
        "which the program cannot do -- a for loop over a range, with a break, can",
        where(),
    )


def loop(iterable, flags):
    """The items of `iterable` until done(*flags()): a loop every point has left takes no more."""
    iterator = iter(iterable)
    while not done(*flags()):
        try:
            item = next(iterator)
        except StopIteration:
            return
        yield item


# ------------------------------------------------------------------------ the expressions


def ifexp(condition, then, otherwise):
    """`then() if condition else otherwise()`."""
    if _static(condition):
        return then() if condition else otherwise()
    with _on(_truth(condition)):
        a = then()
    with _on(_falsity(condition)):
        b = otherwise()
    return merge("the conditional expression", condition, a, b)


def and_(value, *rest):
    """`value and rest[0]() and ...`: the first operand that is false, else the last."""
    if not rest:
        return value
    if _static(value):
        return and_(rest[0](), *rest[1:]) if value else value
    with _on(_truth(value)):
        later = and_(rest[0](), *rest[1:])
    if _truthlike(value) and _truthlike(later):
        return np.logical_and(value, later)
    return merge("`and`", value, later, value)


def or_(value, *rest):
    """`value or rest[0]() or ...`: the first operand that is true, else the last."""
    if not rest:
        return value
    if _static(value):
        return value if value else or_(rest[0](), *rest[1:])
    with _on(_falsity(value)):
        later = or_(rest[0](), *rest[1:])
    if _truthlike(value) and _truthlike(later):
        return np.logical_or(value, later)
    return merge("`or`", value, value, later)


def not_(value):
    """`not value`."""
    if _static(value):
        return not value
    return _falsity(value)


_COMPARISONS = {
    "Eq": operator.eq,
    "NotEq": operator.ne,
    "Lt": operator.lt,
    "LtE": operator.le,
    "Gt": operator.gt,
    "GtE": operator.ge,
    "Is": operator.is_,
    "IsNot": operator.is_not,
    "In": lambda a, b: a in b,
    "NotIn": lambda a, b: a not in b,
}


def compare(left, *rest):
    """A chained comparison: `left`, then the operators and the thunks of the operands in
    `rest` -- `and` of the comparisons in turn, each operand evaluated once."""
    right = rest[1]()
    value = _COMPARISONS[rest[0]](left, right)
    if len(rest) == 2:
        return value
    if _static(value):
        return compare(right, *rest[2:]) if value else value
    with _on(_truth(value)):
        later = compare(right, *rest[2:])
    if _truthlike(value) and _truthlike(later):
        return np.logical_and(value, later)
    return merge("the chained comparison", value, later, value)


def _single(name, value):
    if isinstance(value, Value):
        raise TracingError(
            f"{name}() of one value of the point, which goes over the points: the program "
            f"computes every point alike; {name}(a, b, ...) or {name}([a, b, ...]) of several "
            "values works at every point",
            where(),
        )


def _extremum(function, builtin, better, args, kwargs):
    if function is not builtin or not args or set(kwargs) - {"key", "default"}:
        return function(*args, **kwargs)
    if len(args) == 1:
        _single(builtin.__name__, args[0])
        values = list(args[0])
    elif "default" in kwargs:
        return function(*args, **kwargs)
    else:
        values = list(args)
    if not values:
        return function(values, **kwargs)
    key = kwargs.get("key")
    keys = values if key is None else [key(value) for value in values]
    # as Python: the first value of the least (greatest) key
    best, best_key = values[0], keys[0]
    for value, k in zip(values[1:], keys[1:]):
        condition = better(k, best_key)
        if _static(condition):
            if condition:
                best, best_key = value, k
            continue
        name = f"{builtin.__name__}()"
        best = merge(name, condition, value, best)
        best_key = best if key is None else merge(name, condition, k, best_key)
    return best


def call_min(function, *args, **kwargs):
    """min(), where it is the builtin: at every point, as Python computes it."""
    return _extremum(function, builtins.min, operator.lt, args, kwargs)


def call_max(function, *args, **kwargs):
    """max(), where it is the builtin: at every point, as Python computes it."""
    return _extremum(function, builtins.max, operator.gt, args, kwargs)


def _truths(function, builtin, args, kwargs, stop):
    if function is not builtin or len(args) != 1 or kwargs:
        return function(*args, **kwargs)
    _single(builtin.__name__, args[0])
    pending = []
    for value in args[0]:
        if _static(value):
            if bool(value) is stop:
                return stop
        else:
            pending.append(_truth(value))
    if not pending:
        return not stop
    return reduce(np.logical_or if stop else np.logical_and, pending)


def call_any(function, *args, **kwargs):
    """any(), where it is the builtin: at every point, as Python computes it."""
    return _truths(function, builtins.any, args, kwargs, True)


def call_all(function, *args, **kwargs):
    """all(), where it is the builtin: at every point, as Python computes it."""
    return _truths(function, builtins.all, args, kwargs, False)


# --------------------------------------------------- assert, raise and the value returned


def assert_(condition, message, source):
    """assert: on a value of the point, a check at the points of the path (natively, an
    AssertionError at the first point that fails it)."""
    from ._model import _check

    if _static(condition):
        if condition:
            return
        if path() is None:
            if message is None:
                raise AssertionError
            raise AssertionError(message())
    text = str(message()) if message is not None else f"the assertion {source} fails"
    _check(condition, text, AssertionError)


def raising(exception):
    """Whether to raise `exception`: where the path is every point, as Python does, and natively
    where it has points. In a trace on a path, a check that fails at its points instead.
    """
    current = path()
    if current is None:
        return True
    if isinstance(exception, type) and issubclass(exception, BaseException):
        instance = exception()
    elif isinstance(exception, BaseException):
        instance = exception
    else:
        return True
    if isinstance(current, Value):
        from ._model import _check

        text = str(instance)
        _check(
            False,
            f"{type(instance).__name__}: {text}" if text else type(instance).__name__,
            type(instance),
        )
        return False
    points = np.flatnonzero(np.asarray(current, dtype=bool).ravel())
    if points.size == 0:
        return False
    if instance is exception and hasattr(instance, "add_note"):
        instance.add_note(f"(at point {int(points[0])} of the arrays)")
    return True


def result(returned, value, name):
    """The value of the function `name`: `value` where it returned. A point where it did not
    fails a check: Python would return None there."""
    if _static(returned):
        return value if returned else None
    from ._model import _check

    _check(
        returned,
        f"{name}() returns no value at this point: it ends there without a return",
        ValueError,
    )
    return value
