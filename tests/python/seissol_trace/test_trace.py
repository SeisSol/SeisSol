# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Models in Python traced into sderiv (seissol.sderiv, python -m seissol.trace). The tests that
evaluate a program need the module sderiv of a SeisSol build (-DSDERIV_PYTHON=ON) and are skipped
without it."""

import hashlib
import math
import os
import re
import subprocess
import sys

import numpy as np
import pytest
from conftest import PYTHON
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


@sd.model
def material(x, y, z):
    depth = -z
    rho = np.where(depth > 2000.0, 2700.0, 2600.0)
    vs = 3000.0 + 0.5 * depth
    return {"rho": rho, "mu": rho * vs**2}


def lines_of(text, prefix):
    return [line for line in text.splitlines() if line.startswith(prefix)]


def test_a_model_traces_into_core_text():
    text = material.core()
    assert text.startswith("#!sderiv 1\n")
    assert lines_of(text, "in ") == ["in x, y, z"]
    assert lines_of(text, "out ") == ["out rho = t2", "out mu = t6"]
    # `**2` as NumPy computes it for an array: a square
    assert "def t5 = t4 * t4" in text
    # where every definition and output comes from, in this file
    assert f'meta line = "rho", "{__file__}", ' in text
    # called, it runs as written
    values = material(x=0.0, y=0.0, z=np.array([-3000.0, -100.0]))
    np.testing.assert_array_equal(values["rho"], [2700.0, 2600.0])


def test_a_branch_on_a_value_of_the_point_names_the_line():
    @sd.model
    def branching(z):
        if z > 0.0:  # the line of the error
            return {"u": z}
        return {"u": -z}

    with pytest.raises(sd.TracingError) as caught:
        branching.trace()
    error = caught.value
    assert error.file == __file__
    assert error.line == branching.function.__code__.co_firstlineno + 2
    assert "np.where" in str(error)


@pytest.mark.parametrize(
    "body, fragment",
    [
        (lambda x: {"u": math.sqrt(x)}, "number of Python"),
        (lambda x: {"u": np.round(x)}, "halfway cases to even"),
        (lambda x: {"u": np.fmod(x, 2.0)}, "sign of the dividend"),
        (lambda x: {"u": x // 2.0}, "np.floor(a / b)"),
        (lambda x: {"u": max(x, 1.0)}, "np.maximum"),
        (lambda x: {"u": np.hypot(x, 1.0)}, "hypot"),
    ],
)
def test_what_the_program_cannot_do_is_refused(body, fragment):
    with pytest.raises(sd.TracingError, match=re.escape(fragment)):
        sd.model(body).trace()


def test_outputs_named_by_the_model():
    single = sd.model(lambda x: 2.0 * x, outputs="u")
    assert lines_of(single.core(), "out ") == ["out u = t0"]
    pair = sd.model(lambda x: (x, -x), outputs=["a", "b"])
    assert len(lines_of(pair.core(), "out ")) == 2
    with pytest.raises(sd.TracingError, match="dict of its outputs"):
        sd.model(lambda x: x).trace()
    with pytest.raises(sd.TracingError, match="no output v"):
        sd.model(lambda x: {"u": x}, outputs=["v"]).trace()


def test_names_and_constants():
    # names the core format reserves, in backquotes; a parameter with a default is a constant
    @sd.model
    def named(type, next, scale=3.0):
        return {"out": scale * type + next, "g": 1.0}

    text = named.core()
    assert lines_of(text, "in ") == ["in `type`, `next`"]
    assert "def t0 = 3.0 * `type`" in text
    assert "out `out` = t1" in text
    assert "out g = 1.0" in text


def test_states_carry_over_calls():
    @sd.model(state={"peak": 0.0})
    def peak(v, peak):
        return {"peak": np.maximum(peak, np.abs(v))}

    text = peak.core()
    assert lines_of(text, "in ") == ["in v"]
    assert lines_of(text, "state ") == ["state peak = 0.0"]
    assert len(lines_of(text, "next peak = ")) == 1
    assert len(lines_of(text, "out peak = ")) == 1
    with pytest.raises(sd.TracingError, match="no next value"):
        sd.model(lambda v, s: {"u": v + s}, state={"s": 0.0}).trace()


def test_checks():
    @sd.model
    def checked(z):
        sd.check(z <= 0.0, "below the surface")
        return {"u": -z}

    assert any(
        '"below the surface"' in line for line in lines_of(checked.core(), "check ")
    )
    np.testing.assert_array_equal(checked(np.array([-1.0, 0.0]))["u"], [1.0, -0.0])
    with pytest.raises(ValueError, match="point 1 does not meet the check"):
        checked(np.array([-1.0, 2.0]))


def test_percentile_and_median():
    @sd.model
    def statistics(a, b, c):
        return {"p": sd.percentile([a, b, c], 25.0), "m": sd.median([a, b, c])}

    text = statistics.core()
    assert "percentile(25.0, a, b, c)" in text
    assert "percentile(50.0, a, b, c)" in text
    values = statistics(
        np.array([1.0, 3.0]), np.array([2.0, 1.0]), np.array([3.0, 2.0])
    )
    np.testing.assert_array_equal(values["m"], [2.0, 2.0])


def write_model(directory, text, name="model.py"):
    path = directory / name
    path.write_text(text)
    return path


MODEL = """
import numpy as np
from seissol import sderiv as sd

layers = np.loadtxt("layers.csv", delimiter=",")

@sd.model
def material(z):
    rho = np.where(z < layers[0], layers[1], layers[2])
    return {"rho": rho}
"""


def run_trace(directory, *args):
    environment = dict(os.environ, PYTHONPATH=os.pathsep.join([str(PYTHON)] + sys.path))
    return subprocess.run(
        [sys.executable, "-m", "seissol.trace", *args],
        cwd=directory,
        env=environment,
        capture_output=True,
        text=True,
    )


def test_the_command_writes_the_core_file_with_its_provenance(tmp_path):
    (tmp_path / "layers.csv").write_text("-2000.0, 2700.0, 2600.0\n")
    source = write_model(tmp_path, MODEL)
    done = run_trace(tmp_path, "model.py")
    assert done.returncode == 0, done.stderr
    text = (tmp_path / "model.sderiv").read_text()
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    assert f'meta source = "model.py", "sha256:{digest}"' in text
    csv = hashlib.sha256((tmp_path / "layers.csv").read_bytes()).hexdigest()
    assert f'meta depends = "layers.csv", "sha256:{csv}"' in text
    assert 'meta command = "python -m seissol.trace model.py -o model.sderiv"' in text
    # the hash of the body: every line but the header and the meta lines
    body = "".join(
        line + "\n" for line in text.splitlines()[1:] if not line.startswith("meta ")
    )
    assert f"sha256:{hashlib.sha256(body.encode()).hexdigest()}" in text
    # the data the model read became constants of the program
    assert "lt(z, -2000.0)" in text


def test_the_command_chooses_among_models(tmp_path):
    write_model(
        tmp_path,
        "from seissol import sderiv as sd\n"
        "a = sd.model(lambda x: {'u': x})\n"
        "b = sd.model(lambda x: {'u': -x})\n",
    )
    assert run_trace(tmp_path, "model.py").returncode != 0
    done = run_trace(tmp_path, "model.py", "--model", "b", "-o", "b.sderiv")
    assert done.returncode == 0, done.stderr
    assert "def t0 = -x" in (tmp_path / "b.sderiv").read_text()
    write_model(
        tmp_path,
        "import numpy as np\nfrom seissol import sderiv as sd\n"
        "@sd.model\ndef m(x):\n    if x > 0:\n        return {'u': x}\n    return {'u': x}\n",
        "branch.py",
    )
    done = run_trace(tmp_path, "branch.py")
    assert done.returncode == 1
    assert "branch.py:5:" in done.stderr


# ----------------------------------------------- with the program SeisSol evaluates


def ladder():
    values = np.array(
        [-0.0, 0.0, -1.5, 2.5, np.inf, -np.inf, np.nan, 1e-300, -3.25, 7.0]
    )
    x, y, z = np.meshgrid(values, values, values, indexing="ij")
    return x.ravel(), y.ravel(), z.ravel()


def same_bits(a, b):
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    return bool(
        np.all((a.view(np.uint64) == b.view(np.uint64)) | (np.isnan(a) & np.isnan(b)))
    )


@needs_seissol
def test_the_program_gives_the_bits_of_numpy():
    @sd.model
    def exact(x, y, z):
        return {
            "minmax": np.minimum(x, y)
            + np.maximum(y, x)
            - np.fmin(x, z) * np.fmax(z, y),
            "rest": x % 1.5 + np.sign(y) + np.trunc(z) + np.clip(x, -1.0, 1.0),
            "powers": np.sqrt(np.abs(x)) ** 0.5 + y**-1 + z**2 + x**1 + y**0,
            "logic": np.where((x < y) & ~(z >= 0.0) | (x != x), x, -y),
            "select": np.select([x > 1.0, y > 1.0], [y, z], 0.5),
            "tests": np.isnan(x) + 2.0 * np.isinf(y) + 4.0 * np.isfinite(z),
            "angles": np.deg2rad(x) + np.rad2deg(y) + np.floor(z) - np.ceil(x),
        }

    x, y, z = ladder()
    with np.errstate(all="ignore"):
        native = exact(x, y, z)
    traced = exact.evaluate(x=x, y=y, z=z)
    for name, values in native.items():
        assert same_bits(traced[name], values), name


@needs_seissol
def test_functions_of_the_math_library_within_a_few_ulp():
    @sd.model
    def bounded(x, y):
        return {
            "u": np.exp(x) + np.log(y) + np.sin(x) * np.cos(y) + np.arctan2(x, y),
            "v": np.tanh(x) + np.arcsinh(y) + np.power(y, 1.5) + sd.erf(x),
        }

    rng = np.random.default_rng(7)
    x = rng.uniform(-3.0, 3.0, 2000)
    y = rng.uniform(0.1, 5.0, 2000)
    native = bounded(x, y)
    traced = bounded.evaluate(x=x, y=y)
    for name, values in native.items():
        np.testing.assert_allclose(traced[name], values, rtol=1e-14, atol=1e-14)


@needs_seissol
def test_vectors_of_values():
    @sd.model
    def traction(s_xx, s_yy, s_xy):
        stress = np.array([[s_xx, s_xy], [s_xy, s_yy]])
        normal = np.array([0.6, 0.8])
        t = stress @ normal
        return {"normal": t @ normal, "shear": np.sqrt(t @ t - (t @ normal) ** 2)}

    values = traction.evaluate(s_xx=1.0, s_yy=2.0, s_xy=0.5)
    t = np.array([[1.0, 0.5], [0.5, 2.0]]) @ np.array([0.6, 0.8])
    assert values["normal"][0] == pytest.approx(t @ np.array([0.6, 0.8]))


@needs_seissol
def test_states_and_checks_in_the_program():
    @sd.model(state={"peak": 0.0})
    def peak(v, peak):
        sd.check(np.isfinite(v), "a finite velocity")
        return {"peak": np.maximum(peak, np.abs(v))}

    program = peak.program()
    outputs, states = program.step({"v": [1.0, -3.0]})
    outputs, states = program.step({"v": [-2.0, 1.0]}, states)
    np.testing.assert_array_equal(states["peak"], [2.0, 3.0])
    with pytest.raises(compiled.EvaluationError, match="a finite velocity"):
        program.step({"v": [np.inf]})


@needs_seissol
def test_a_grid_samples_as_in_seissol():
    layers = sd.grid(
        components="rho",
        knots=[[-3000.0, -1000.0, 0.0]],
        values=[2700.0, 2600.0, 2500.0],
        interpolation="linear",
    )

    @sd.model
    def layered(z):
        return {"rho": layers.rho(z)}

    text = layered.core()
    assert lines_of(text, "grid g0 ")
    z = np.array([-4000.0, -2000.0, -500.0, 100.0])
    traced = layered.evaluate(z=z)["rho"]
    np.testing.assert_allclose(traced, [2700.0, 2650.0, 2550.0, 2500.0])
    # outside of a trace, the same program samples it
    assert same_bits(layered(z)["rho"], traced)


@needs_seissol
def test_the_command_checks_the_core_file(tmp_path):
    (tmp_path / "layers.csv").write_text("-2000.0, 2700.0, 2600.0\n")
    write_model(tmp_path, MODEL)
    done = run_trace(tmp_path, "model.py", "--check")
    assert done.returncode == 0, done.stderr
    assert "program " in done.stdout and "z -> rho" in done.stdout


@needs_seissol
def test_seissol_refuses_a_core_file_whose_model_changed(tmp_path):
    (tmp_path / "layers.csv").write_text("-2000.0, 2700.0, 2600.0\n")
    source = write_model(tmp_path, MODEL)
    assert run_trace(tmp_path, "model.py").returncode == 0
    core = tmp_path / "model.sderiv"
    text = core.read_text()
    assert compiled.compile(text, file=str(core)).outputs == ["rho"]
    # the model edited since: an error that gives the command which writes the file again
    source.write_text(MODEL + "# edited\n")
    with pytest.raises(
        compiled.Error, match="python -m seissol.trace model.py -o model.sderiv"
    ):
        compiled.compile(text, file=str(core))
