# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The Python module of sderiv: load, print, validate and evaluate programs."""

import math
import os
import tempfile
import unittest

import numpy as np
import sderiv


class Programs(unittest.TestCase):
    def test_evaluate_at_points(self):
        program = sderiv.compile("out def u = x * y + 1\nout def v = sqrt(abs(x))\n")
        self.assertEqual(program.inputs, ["x", "y"])
        self.assertEqual(program.outputs, ["u", "v"])
        self.assertEqual(program.compute_type, "f64")
        self.assertEqual(len(program.fingerprint), 16)
        # arrays, a number at every point, and inputs by name
        values = program.evaluate({"x": np.array([1.0, -4.0, 9.0])}, y=2.0)
        np.testing.assert_array_equal(values["u"], [3.0, -7.0, 19.0])
        np.testing.assert_array_equal(values["v"], [1.0, 2.0, 3.0])
        # one point, of numbers; inputs the program does not read are left
        self.assertEqual(program.evaluate(x=2.0, y=3.0, z=5.0)["u"].tolist(), [7.0])

    def test_inputs_given_wrong(self):
        program = sderiv.compile("out def u = x + y\n")
        with self.assertRaises(KeyError):
            program.evaluate(x=1.0)
        with self.assertRaises(ValueError):
            program.evaluate(x=[1.0, 2.0], y=[1.0, 2.0, 3.0])
        with self.assertRaises(ValueError):
            program.evaluate(x=np.zeros((2, 2)), y=1.0)

    def test_normal_form(self):
        program = sderiv.compile("out def rho = if z < -2000 then 2700 else 2600\n")
        core = program.core()
        self.assertTrue(core.startswith("#!sderiv 1\n"))
        self.assertEqual(str(program), core)
        # the core format reads back into the same program
        again = sderiv.compile(core)
        self.assertEqual(again.fingerprint, program.fingerprint)
        self.assertIsNone(again.validate())
        self.assertIn(program.fingerprint, repr(program))

    def test_file(self):
        example = os.environ.get("SDERIV_EXAMPLE")
        if example is None:
            self.skipTest("SDERIV_EXAMPLE names no example")
        program = sderiv.load(example)
        values = program.evaluate(z=np.array([-3000.0, 0.0]))
        np.testing.assert_array_equal(values["rho"], [2700.0, 2600.0])

    def test_errors_name_the_line_and_column(self):
        with self.assertRaises(sderiv.Error) as caught:
            sderiv.compile("def a = 1\nout def u = a +\n", file="model.sderiv")
        error = caught.exception
        self.assertIsInstance(error, ValueError)
        self.assertEqual(error.file, "model.sderiv")
        self.assertEqual(error.line, 3)
        self.assertIn("model.sderiv:3:", str(error))
        with self.assertRaises(sderiv.Error) as caught:
            sderiv.compile("#!sderiv 1\nin x\nout u = y\n")
        self.assertEqual(caught.exception.line, 3)
        with self.assertRaises(sderiv.Error):
            sderiv.load("/nonexistent/model.sderiv")

    def test_state(self):
        program = sderiv.compile("state s = 0\nnext s = s + x\nout u = s\n")
        self.assertEqual(program.states, {"s": 0.0})
        with self.assertRaises(ValueError):
            program.evaluate(x=1.0)
        # the state of the previous call, from the value it starts from
        outputs, states = program.step({"x": [1.0, 2.0]})
        np.testing.assert_array_equal(outputs["u"], [0.0, 0.0])
        np.testing.assert_array_equal(states["s"], [1.0, 2.0])
        outputs, states = program.step({"x": [1.0, 2.0]}, states)
        np.testing.assert_array_equal(outputs["u"], [1.0, 2.0])
        np.testing.assert_array_equal(states["s"], [2.0, 4.0])
        # the states give the points where the inputs are numbers
        outputs, states = program.step({"x": 1.0}, {"s": [5.0, 6.0, 7.0]})
        np.testing.assert_array_equal(states["s"], [6.0, 7.0, 8.0])

    def test_check(self):
        program = sderiv.compile('check x < 1, "x below 1"\nout def u = x\n')
        self.assertEqual(program.checks, ["x below 1"])
        np.testing.assert_array_equal(program.evaluate(x=[0.0, 0.5])["u"], [0.0, 0.5])
        with self.assertRaises(sderiv.EvaluationError) as caught:
            program.evaluate(x=[0.0, 2.0])
        self.assertIn("x below 1", str(caught.exception))
        self.assertIn("point 1", str(caught.exception))

    def test_g_is_an_input(self):
        program = sderiv.compile("out def p = 2600 * g * -z\n")
        self.assertEqual(sorted(program.inputs), ["g", "z"])
        self.assertEqual(
            program.evaluate(g=9.81e-3, z=-2.0)["p"].tolist(), [2600 * 9.81e-3 * 2.0]
        )

    def test_meta(self):
        program = sderiv.compile('meta generator = "test", 1\nout u = x\n')
        self.assertEqual(program.meta, [("generator", ["test", 1])])

    def test_compute_type(self):
        program = sderiv.compile("#!sderiv 1\nin x\nout u = real(f32(x) * f32(3))\n")
        value = program.evaluate(x=0.1)["u"][0]
        self.assertEqual(value, float(np.float32(0.1) * np.float32(3)))
        self.assertNotEqual(value, 0.1 * 3)

    def test_lua_and_easi(self):
        with tempfile.TemporaryDirectory() as directory:
            lua = os.path.join(directory, "model.lua")
            with open(lua, "w") as file:
                file.write(
                    "local M = {}\n"
                    'M.output_parameters = { "rho" }\n'
                    "function M.evaluate(fields, x, z)\n"
                    "  if z < -2000 then return 2700 + x end\n"
                    "  return 2600\n"
                    "end\n"
                    "return M\n"
                )
            program = sderiv.load(lua)
            values = program.evaluate(x=[1.0, 1.0], z=[-3000.0, 0.0])
            np.testing.assert_array_equal(values["rho"], [2701.0, 2600.0])
            easi = os.path.join(directory, "material.yaml")
            with open(easi, "w") as file:
                file.write(
                    "!FunctionMap\n"
                    "map:\n"
                    "  rho: return 2600 - 0.5 * z;\n"
                    "  mu: return 10000000000.0;\n"
                )
            with self.assertRaises(sderiv.Error):
                sderiv.load(easi)
            program = sderiv.load_easi(easi, ["rho"])
            self.assertEqual(program.outputs, ["rho"])
            values = program.evaluate(x=0.0, y=0.0, z=[-10.0, 0.0], group=1)
            np.testing.assert_array_equal(values["rho"], [2605.0, 2600.0])
            self.assertEqual(sorted(sderiv.load_easi(easi).outputs), ["mu", "rho"])

    def test_grid(self):
        # a grid the program gives itself; the module starts MPI for it, and ends it with Python
        program = sderiv.compile(
            'grid layers (kind = "inline", components = ("rho"), knots = ((-3000, -1000)),\n'
            '             values = ((2700), (2600)), interpolation = "previous")\n'
            "out def rho = layers.rho(z)\n"
        )
        self.assertEqual(len(program.grids), 1)
        values = program.evaluate(z=[-4000.0, -2000.0, 0.0])["rho"]
        np.testing.assert_array_equal(values, [2700.0, 2700.0, 2600.0])

    def test_nan_and_infinity(self):
        program = sderiv.compile("out def u = x / y\n")
        values = program.evaluate(x=[1.0, 0.0], y=[0.0, 0.0])["u"]
        self.assertTrue(math.isinf(values[0]))
        self.assertTrue(math.isnan(values[1]))


if __name__ == "__main__":
    unittest.main()
