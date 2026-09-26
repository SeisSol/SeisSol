# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The point set the material is sampled at inside a cell.

Separate from the one the plastic strain lives on, because the two answer
different questions. The plastic strain wants points to evaluate a yield
criterion at; the material wants points whose spacing carries the operator.
A set whose face traces are the two-dimensional nodal set gives the flux its
samples for free but only integrates a product exactly while the material is
linear; the conical-product set integrates far beyond that and has no point on
a face at all.
"""

from kernels.multsim import OptionalDimTensor
from yateto import Tensor
from yateto.input import parseJSONMatrixFile

#: The sets a build can choose between, as the matrix files name them.
SETS = ("nb", "ip")


def _db(matricesDir, pointSet, order, alignStride):
    return parseJSONMatrixFile(
        f"{matricesDir}/plasticity-{pointSet}-matrices-{order}.json",
        clones=dict(),
        alignStride=alignStride,
    )


def tensors(matricesDir, aderdg, pointSet):
    """The three tensors a nodal material needs: where its samples sit, how a
    modal field is read at those points, and how a field given there comes
    back to the modal basis.

    Renamed off the plasticity matrices so that nothing downstream has to know
    which set a build chose, and so that the two can differ.
    """
    db = _db(matricesDir, pointSet, aderdg.order, aderdg.alignStride)

    def renamed(name, source):
        # der Tensor traegt seine Werte im Sparsity-Muster, also von dort holen
        return Tensor(name, source.shape(), spp=dict(source.values()),
                      alignStride=aderdg.alignStride(name))

    return {
        name: renamed(name, source)
        for name, source in (("materialNodes", db.vNodes),
                             ("materialEval", db.v),
                             ("materialProject", db.vInv))
    }


def pointCount(matricesDir, aderdg, pointSet):
    return _db(matricesDir, pointSet, aderdg.order, aderdg.alignStride).vNodes.shape()[0]


def includeTensors(matricesDir, aderdg, pointSet, include):
    """The sample points are read by the host, which builds the query that
    fills them, so they have to reach the generated code even where no kernel
    names them."""
    for tensor in tensors(matricesDir, aderdg, pointSet).values():
        include.add(tensor)


def addKernels(generator, aderdg, matricesDir, pointSet):
    """The one kernel a nodal material needs before anything reads it: the
    samples, as the modal field the operator carries.

    Whoever reads the material -- the output, a diagnostic, a later operator --
    reads it modally and so does not have to know which point set a build
    sampled at.
    """
    mats = tensors(matricesDir, aderdg, pointSet)
    npoints = pointCount(matricesDir, aderdg, pointSet)

    samples = OptionalDimTensor(
        "materialSamples",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (npoints,),
        alignStride=True,
    )
    modal = OptionalDimTensor(
        "modalVar",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.num3DBasisFunctions(),),
        alignStride=True,
    )
    generator.add(
        "projectMaterialToModal",
        modal["k"] <= mats["materialProject"]["kn"] * samples["n"],
    )
