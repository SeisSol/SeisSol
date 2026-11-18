# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Branches on values of the point in a model, which the conversion makes take both ways and
merge at every point as Python takes them there. A model is checked against its function as
written, run point by point on scalars of NumPy -- Python's own semantics: on arrays, and with
the module sderiv of a SeisSol build (-DSDERIV_PYTHON=ON) in the program SeisSol evaluates.
"""

import math
import re
import sys

import numpy as np
import pytest
from seissol import sderiv as sd

try:
    import sderiv as compiled

    if not hasattr(compiled, "compile"):
        compiled = None
except ImportError:
    compiled = None

needs_seissol = pytest.mark.skipif(
    compiled is None, reason="needs the module sderiv of a SeisSol build"
)

VALUES = np.array(
    [-0.0, 0.0, -1.5, 2.5, np.inf, -np.inf, np.nan, 1e-300, -3.25, 7.0, -2000.0, 4.0]
)


def points(count):
    grids = np.meshgrid(*([VALUES] * count), indexing="ij")
    return [grid.ravel() for grid in grids]


def same_bits(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return a.shape == b.shape and bool(
        np.all((a.view(np.uint64) == b.view(np.uint64)) | (np.isnan(a) & np.isnan(b)))
    )


def agrees(model):
    """Checks that `model` computes on arrays, and in the program, what its function computes
    point by point; returns its core text."""
    names = model.inputs
    arrays = points(len(names))
    expected = {}
    with np.errstate(all="ignore"):
        for i in range(arrays[0].size):
            values = model.function(*[np.float64(a[i]) for a in arrays])
            for name, value in values.items():
                expected.setdefault(name, []).append(float(value))
        native = model(*arrays)
    for name, values in expected.items():
        got = np.broadcast_to(np.asarray(native[name], dtype=float), arrays[0].shape)
        assert same_bits(got, values), f"{name} on arrays"
    if compiled is not None:
        traced = model.evaluate(**dict(zip(names, arrays)))
        for name, values in expected.items():
            assert same_bits(traced[name], values), f"{name} in the program"
    return model.core()


def body_of(text):
    return [line for line in text.splitlines()[1:] if not line.startswith("meta ")]


LAYERS = [(-1000.0, 2500.0), (-5000.0, 2700.0)]


def test_if_elif_else():
    @sd.model
    def layered(z):
        depth = -z
        if depth < 1000.0:
            vp, vs, rho = 4000.0, 2300.0, 2600.0
        elif depth < 5000.0:
            vp, vs, rho = 5500.0, 3200.0, 2700.0
        else:
            vp, vs, rho = 6500.0, 3700.0, 2900.0
        mu = rho * vs**2
        return {"rho": rho, "mu": mu, "lambda": rho * vp**2 - 2.0 * mu}

    text = agrees(layered)
    assert "def t4 = select(t1, 2600.0, t3)" in text
    assert "check" not in text


def test_returns_breaks_and_continues_in_branches():
    @sd.model
    def early(z):
        if z > 0.0:
            return {"rho": 0.0}
        rho = 2500.0 - z
        return {"rho": rho}

    # the rest of the function where it did not return: one select
    assert body_of(agrees(early)) == [
        "in z",
        "def t0 = lt(0.0, z)",
        "def t1 = 2500.0 - z",
        "def t2 = select(t0, 0.0, t1)",
        "out rho = t2",
    ]

    @sd.model
    def first_layer(z):
        rho = 3000.0
        for top, value in LAYERS:
            if z > top:
                rho = value
                break
        return {"rho": rho}

    @sd.model
    def returned_in_a_loop(z):
        for top, value in LAYERS:
            if z > top:
                return {"rho": value}
        return {"rho": 3000.0}

    @sd.model
    def skipping(x):
        total = 0.0
        for scale in (1.0, 2.0, 3.0):
            if x * scale > 4.0:
                continue
            total = total + x * scale
        else:
            total = -total
        return {"total": total}

    for model in (first_layer, returned_in_a_loop, skipping):
        assert "check" not in agrees(model)


def test_expressions():
    @sd.model
    def expressions(x, y):
        return {
            "ifexp": x if x > y else y,
            "andor": (x > 0.0 and y > 0.0) or x < -5.0,
            "values": (x and y) + 2.0 * (x or y),
            "not": not x,
            "chained": -1.0 < x <= 2.0 < y,
            "min": min(x, y, 1.0),
            "max": max([y, x], key=lambda v: -v),
            "anyall": any([x > 1.0, y > 1.0]) + 2.0 * all([x > 1.0, y > 1.0]),
        }

    text = agrees(expressions)
    # truth values combine as truth values
    assert "land(" in text and "lor(" in text


def test_nested_branches_functions_and_containers():
    @sd.model
    def nested(x, y):
        def clamp(v, low, high):
            if v < low:
                return low
            if v > high:
                return high
            return v

        u = 0.0
        if x > 0.0:
            if y:
                u = x + y
            else:
                u = x - y
            v = 2.0 * u
        else:
            v = -1.0
        return {"u": u, "v": v, "w": clamp(x, -1.0, 2.0) + clamp(y, 0.0, 1.0)}

    @sd.model
    def changed(x):
        out = {"rho": 1.0, "vs": [1.0, 2.0]}
        if x > 0.0:
            out["rho"] = x
            out["vs"].append(3.0)
        else:
            out["vs"].append(-x)
        out["vs"][0] += 1.0
        return {"rho": out["rho"], "vs": out["vs"][2] + out["vs"][0]}

    @sd.model
    def counted(x):
        k = 0
        total = 0.0
        while k < 3:
            if x > k:
                total = total + 1.0
            k += 1
        return {"total": total}

    for model in (nested, changed, counted):
        agrees(model)


def test_statements_and_loops_in_branches():
    @sd.model
    def statements(x, y):
        import contextlib

        if x > 0.0:
            try:
                u = x / 2.0
            except ZeroDivisionError:
                u = 0.0
            with contextlib.nullcontext():
                v = u + 1.0
        else:
            v = 2.0
        match 3:
            case 3:
                w = v if y > 0.0 else -v
            case _:
                w = 0.0
        t = 2.0 * x
        del t
        if y > 0.0:
            t = 3.0
        else:
            t = 4.0
        pick = lambda a, b: a if a > b else b  # noqa: E731
        return {"w": w + t + pick(x, y)}

    @sd.model
    def loops(x):
        found = -1.0
        for i, top in enumerate((5.0, 2.0, -1.0)):
            if x > top:
                found = float(i)
                break
        else:
            found = 99.0
        total = 0.0
        for i in range(3):
            for j in range(3):
                if x > i + j:
                    continue
                if x < -i:
                    break
                total = total + 1.0
        return {"found": found, "total": total}

    @sd.model
    def nested_returns(x):
        import contextlib

        for i in range(3):
            for j in range(3):
                if x > 10.0 * i + j:
                    with contextlib.nullcontext():
                        return {"u": float(10 * i + j)}
        return {"u": -1.0}

    for model in (statements, loops, nested_returns):
        agrees(model)

    @sd.model
    def counting(x):
        k = 0
        while k < 5:
            k += 1
            if x * k > 6.0:
                break
        return {"k": k}

    # a while loop some points leave goes on as the others: its condition is Python's
    assert counting(2.0) == {"k": 4}
    with pytest.raises(
        sd.TracingError, match="the condition of a while loop is a value"
    ):
        counting.trace()


def test_a_condition_python_knows_takes_one_way():
    @sd.model
    def chosen(z, layered=True):
        if layered:
            return {"rho": np.where(z < -1000.0, 2700.0, 2600.0)}
        return {"rho": np.fmod(z, 2.0)}  # refused, were it traced

    assert "select(" in chosen.core()
    for count in range(4):

        @sd.model
        def stepped(x, count=count):
            total = 0.0
            for i in range(10):
                if i == count:
                    break
                total = total + x
            return {"total": total}

        total = stepped(np.array([1.0, 2.0]))["total"]
        assert np.broadcast_to(total, (2,)).tolist() == [count, 2.0 * count]
        assert body_of(stepped.core())[-1].startswith("out total = ")


def test_checks_and_assertions_apply_on_their_branch():
    @sd.model
    def checked(z):
        if z > 0.0:
            sd.check(z < 10.0, "below 10 above the surface")
            assert z != 5.0, "not at 5"
        else:
            assert z > -100.0
        return {"u": z}

    np.testing.assert_array_equal(checked(np.array([-20.0, 3.0]))["u"], [-20.0, 3.0])
    with pytest.raises(ValueError, match="below 10 above the surface -- point 1"):
        checked(np.array([-20.0, 12.0]))
    with pytest.raises(AssertionError, match="not at 5 -- point 0"):
        checked(np.array([5.0]))
    with pytest.raises(AssertionError):
        checked(-200.0)
    text = checked.core()
    assert [line for line in text.splitlines() if line.startswith("check ")][
        -1
    ].endswith('"the assertion z > -100.0 fails"')
    if compiled is not None:
        program = checked.program()
        program.evaluate(z=[-20.0, 3.0])
        with pytest.raises(compiled.EvaluationError, match="not at 5"):
            program.evaluate(z=[1.0, 5.0])


def test_a_raise_in_a_branch_and_a_point_without_a_return():
    @sd.model
    def guarded(z):
        if z > 0.0:
            raise ValueError("no material above the surface")
        return {"rho": 2700.0 - z}

    assert guarded(-1.0) == {"rho": 2701.0}
    with pytest.raises(ValueError, match="no material above the surface") as caught:
        guarded(np.array([-1.0, 2.0]))
    if sys.version_info >= (3, 11):
        assert "(at point 1 of the arrays)" in caught.value.__notes__
    assert '"ValueError: no material above the surface"' in guarded.core()

    @sd.model
    def partial(z):
        if z < 0.0:
            return {"rho": 2700.0}

    assert partial(1.0) is None
    with pytest.raises(
        ValueError, match="returns no value at this point: .* -- point 1"
    ):
        partial(np.array([-1.0, 1.0]))
    assert "returns no value at this point" in partial.core()
    if compiled is not None:
        with pytest.raises(compiled.EvaluationError, match="no material above"):
            guarded.evaluate(z=[-1.0, 2.0])
        with pytest.raises(compiled.EvaluationError, match="returns no value"):
            partial.evaluate(z=[-1.0, 1.0])


def test_a_name_with_a_value_at_some_points_only():
    @sd.model
    def onesided(z):
        if z < 0.0:
            rho = 2700.0
        return {"rho": rho}

    assert onesided(-1.0) == {"rho": 2700.0}
    with pytest.raises(UnboundLocalError, match="'rho'"):
        onesided(1.0)
    with pytest.raises(sd.TracingError, match="`rho` has a value at some points only"):
        onesided.trace()
    with pytest.raises(sd.TracingError, match="`rho` has a value at some points only"):
        onesided(np.array([-1.0, 1.0]))

    @sd.model
    def different(z):
        if z < 0.0:
            label = "below"  # noqa: F841
        else:
            label = "above"  # noqa: F841
        return {"rho": 2700.0}  # unused, a label of each side is no matter

    assert body_of(different.core()) == ["in z", "out rho = 2700.0"]


def test_the_lines_of_the_model():
    @sd.model
    def deep(x):
        if x > 0.0:
            y = math.sqrt(x)  # the line of the error
        else:
            y = 0.0
        return {"y": y}

    first = deep.function.__code__.co_firstlineno
    with pytest.raises(sd.TracingError) as caught:
        deep.trace()
    assert caught.value.file == __file__
    assert caught.value.line == first + 3

    @sd.model
    def merged(x):
        if x > 0.0:
            y = 2.0 * x
        else:
            y = -x
        return {"y": y}

    first = merged.function.__code__.co_firstlineno
    # the select is where the if is (its column where Python has columns, from 3.11)
    column = 9 if sys.version_info >= (3, 11) else 1
    assert f'meta line = "y", "{__file__}", {first + 2}, {column}' in merged.core()


def test_what_the_conversion_does_not_follow():
    @sd.model
    def walrus(x):
        if x > 0.0 and (y := 2.0 * x) > 1.0:
            return {"u": y}
        return {"u": x}

    with pytest.warns(UserWarning, match=re.escape("assignment expression (:=)")):
        with pytest.raises(sd.TracingError, match="@sd.function"):
            walrus.trace()

    @sd.model
    def looping(z):
        while z > 0.0:  # the line of the error
            z = z - 1.0
        return {"z": z}

    with pytest.raises(sd.TracingError, match="a while loop on it") as caught:
        looping.trace()
    assert caught.value.line == looping.function.__code__.co_firstlineno + 2


def clamped(v):
    if v < 0.0:
        return 0.0
    return v


@sd.function
def converted_clamped(v):
    if v < 0.0:
        return 0.0
    return v


def test_functions_of_their_own():
    @sd.model
    def calling(x):
        return {"u": converted_clamped(x)}

    agrees(calling)
    with pytest.raises(sd.TracingError, match="@sd.function"):
        sd.model(lambda x: {"u": clamped(x)}).trace()
    # a lambda as the model
    lambda_model = sd.model(lambda z: {"rho": 2600.0 if z > -1000.0 else 2700.0})
    assert "select(" in lambda_model.core()


def test_states_and_vectors_in_branches():
    @sd.model(state={"peak": 0.0})
    def peak(v, peak):
        if abs(v) > peak:
            peak = abs(v)
        return {"peak": peak}

    @sd.model
    def turned(x, y):
        v = np.array([x, y]) if x > 0.0 else np.array([y, x])
        return {"a": v[0], "b": v[1]}

    agrees(turned)
    assert (
        len([line for line in peak.core().splitlines() if line.startswith("next ")])
        == 1
    )
    if compiled is not None:
        program = peak.program()
        _, states = program.step({"v": [1.0, -3.0]})
        _, states = program.step({"v": [-2.0, 1.0]}, states)
        np.testing.assert_array_equal(states["peak"], [2.0, 3.0])
