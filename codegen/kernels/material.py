# SPDX-FileCopyrightText: 2026 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

"""The point set the material is sampled at inside a cell.

Separate from the one the plastic strain lives on, because the two answer
different questions. The plastic strain wants points to evaluate a yield
criterion at; the material wants points to form the operator at. Both sets a
build can sample the material at are quadrature rules -- the symmetric rule of
degree 2N with positive weights and interior points (sq), or the conical
product of Gauss-Jacobi rules (ip) -- so the operator is formed at the samples
themselves and projected back with the weights of the rule: the Galerkin
projection of the product of material and field, exact where the material is
constant. Neither set has a point on a face.
"""

import numpy as np
from kernels.multsim import OptionalDimTensor
from yateto import Tensor
from yateto.input import parseJSONMatrixFile

#: The sets a build can choose between, as the matrix files name them; the
#: first one is the default.
SETS = ("sq", "ip")


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
        # the tensor carries its values in the sparsity pattern, so take them there
        return Tensor(
            name,
            source.shape(),
            spp=dict(source.values()),
            alignStride=aderdg.alignStride(name),
        )

    return {
        name: renamed(name, source)
        for name, source in (
            ("materialNodes", db.vNodes),
            ("materialEval", db.v),
            ("materialProject", db.vInv),
        )
    }


def _interpolation(name, aderdg, samples, points, same):
    """The values the samples give at the points of another set.

    The samples give a modal field, and the field has a value at every point;
    at the samples themselves that is the samples again.
    """

    def dense(source):
        matrix = np.zeros(source.shape())
        for index, value in source.values().items():
            matrix[index] = float(value)
        return matrix

    if same:
        interpolation = np.eye(samples.vNodes.shape()[0])
    else:
        interpolation = dense(points.v) @ dense(samples.vInv)
    interpolation[np.abs(interpolation) < 1e-14] = 0.0
    return Tensor(
        name,
        interpolation.shape,
        spp={
            index: repr(float(value))
            for index, value in np.ndenumerate(interpolation)
            if value != 0.0
        },
        alignStride=aderdg.alignStride(name),
    )


def quadratureInterpolation(matricesDir, aderdg, pointSet):
    """The material at the points of the volume quadrature, the conical-product
    set a modal field is integrated over.

    A quantity integrated over the cell with a material that varies inside it,
    such as an energy, reads the material there: the samples themselves where
    they are that set, and what they interpolate there where they are not.
    """
    quadrature = _db(matricesDir, "ip", aderdg.order, aderdg.alignStride)
    samples = _db(matricesDir, pointSet, aderdg.order, aderdg.alignStride)
    return _interpolation(
        "materialToQuadrature", aderdg, samples, quadrature, pointSet == "ip"
    )


def pointCount(matricesDir, aderdg, pointSet):
    return _db(matricesDir, pointSet, aderdg.order, aderdg.alignStride).vNodes.shape()[
        0
    ]


def includeTensors(matricesDir, aderdg, pointSet, include):
    """The sample points are read by the host, which builds the query that
    fills them, and so is their interpolation to the volume quadrature, which
    the energies are integrated with; they have to reach the generated code even
    where no kernel names them."""
    for tensor in tensors(matricesDir, aderdg, pointSet).values():
        include.add(tensor)
    include.add(quadratureInterpolation(matricesDir, aderdg, pointSet))


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


#: The three reparametrisations a shared face can have, as permutations of the
#: barycentric coordinates of the reference triangle. All three are odd: two
#: tetrahedra see the face between them with opposite orientation, so the map
#: from one parametrisation to the other always reverses it, and which of the
#: three it is depends on which vertices pair up.
FACE_REFLECTIONS = ((0, 2, 1), (1, 0, 2), (2, 1, 0))


def faceOrientationPermutations(matricesDir, aderdg):
    """How the points of a face are renumbered between the two cells sharing it.

    A value given at the points of a face reaches the other cell's numbering by
    being reordered, nothing more: the rule the flux is integrated with is
    symmetric under the triangle's reflections, so each one maps it onto
    itself. The generated fP carries the same map in the modal basis with a
    mass factor, which is what the weak form needs and what a value does not.
    """
    db = parseJSONMatrixFile(
        f"{matricesDir}/face-sq-matrices-{aderdg.order}.json",
        clones=dict(),
    )
    shape = db.faceQuadPoints.shape()
    points = np.zeros(shape)
    for idx, value in db.faceQuadPoints.values().items():
        points[idx] = float(value)

    bary = np.stack(
        [1.0 - points[:, 0] - points[:, 1], points[:, 0], points[:, 1]], axis=1
    )

    permutations = []
    for sigma in FACE_REFLECTIONS:
        moved = np.stack([bary[:, sigma[1]], bary[:, sigma[2]]], axis=1)
        order = []
        for target in moved:
            distance = np.linalg.norm(points - target, axis=1)
            nearest = int(np.argmin(distance))
            if distance[nearest] > 1e-10:
                raise RuntimeError(
                    "the points of the face rule are not symmetric under the "
                    "reflections of the triangle, so a face cannot be renumbered"
                )
            order.append(nearest)
        permutations.append(tuple(order))
    return tuple(permutations)
