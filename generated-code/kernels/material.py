# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
# SPDX-FileContributor: Carsten Uphoff

import numpy as np
from yateto import Tensor
from yateto.input import parseJSONMatrixFile


def addKernels(
    generator, aderdg, matricesDir, PlasticityMethod, materialorder, includeTensors
):
    if materialorder is None:
        materialorder = 1

    # Load matrices
    db = parseJSONMatrixFile(
        f"{matricesDir}/hom-{PlasticityMethod}-h{materialorder}.json",
        clones=dict(),
        alignStride=aderdg.alignStride,
    )

    numberOfNodes = db.homproject.shape()[1]

    starM = Tensor(
        "starM", shape=tuple([numberOfNodes] + list(aderdg.starMatrix(0).shape())[1:])
    )
    godLocalM = Tensor(
        "godLocalM", shape=tuple([numberOfNodes] + list(aderdg.QgodLocal.shape())[1:])
    )
    godNeighborM = Tensor(
        "godNeighborM",
        shape=tuple([numberOfNodes] + list(aderdg.QgodNeighbor.shape())[1:]),
    )

    generator.add(
        "homStar", aderdg.starMatrix(0)["Mij"] <= db.homproject["MK"] * starM["Kij"]
    )
    generator.add(
        "homFluxsolver",
        [
            aderdg.QgodLocal["Mij"] <= db.homproject["MK"] * godLocalM["Kij"],
            aderdg.QgodNeighbor["Mij"] <= db.homproject["MK"] * godNeighborM["Kij"],
        ],
    )

    for x in db.__dict__:
        if isinstance(db.__dict__[x], dict):
            for y in db.__dict__[x]:
                includeTensors.add(db.__dict__[x][y])
        else:
            includeTensors.add(db.__dict__[x])
