# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""Models in Python, translated into sderiv for SeisSol.

A model is a function of the inputs of a point -- its parameters, x, y, z, t or the parameters
of a material, as the consumer gives them -- that returns its outputs by name::

    import numpy as np
    from seissol import sderiv as sd

    @sd.model
    def material(x, y, z):
        depth = -z
        rho = np.where(depth > 2000.0, 2700.0, 2600.0)
        vs = 3000.0 + 0.5 * depth
        return {"rho": rho, "mu": rho * vs**2}

It runs as written on NumPy arrays, and a trace translates it into a program: every operation on
a value of the point becomes an operation of the program, with the bits NumPy gives (np.minimum,
np.remainder and `**` included). A branch on a value of the point (if, and, or, min(), max()) or
its use as a number of Python (math.*, float()) fails with the line of the model; np.where,
np.select, np.minimum and & | ~ express them. `python -m seissol.trace model.py` writes the core
file model.sderiv next to the model, which SeisSol reads as python:model.py.
"""

from ._graph import TracingError
from ._grid import Grid, grid
from ._model import (
    GENERATOR,
    VERSION,
    Model,
    Trace,
    check,
    depends,
    erf,
    median,
    model,
    percentile,
    recording,
)

__version__ = VERSION

__all__ = [
    "GENERATOR",
    "Grid",
    "Model",
    "Trace",
    "TracingError",
    "check",
    "depends",
    "erf",
    "grid",
    "median",
    "model",
    "percentile",
    "recording",
]
