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
from yateto import Tensor, simpleParameterSpace
from yateto.input import parseJSONMatrixFile
from yateto.util import tensor_collection_from_constant_expression

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


def pointCount(matricesDir, aderdg, pointSet):
    return _db(matricesDir, pointSet, aderdg.order, aderdg.alignStride).vNodes.shape()[
        0
    ]


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


def addFaceKernels(generator, aderdg, matricesDir, pointSet):
    """The material where a face needs it.

    A flux is built from the material of both cells at the points of the face,
    and those points are the two-dimensional nodal set. Which points a cell
    sampled its material at does not have to be among them: the samples give a
    modal field, and the modal field has a value everywhere. Folding the two
    steps into one matrix per side is what makes that free -- the samples go
    straight to the face, without the modal coefficients ever being written.
    """
    mats = tensors(matricesDir, aderdg, pointSet)
    npoints = pointCount(matricesDir, aderdg, pointSet)

    folded = tensor_collection_from_constant_expression(
        base_name="materialToFace",
        expressions=lambda side: aderdg.db.V3mTo2nFace[side][aderdg.t("kl")]
        * mats["materialProject"]["ln"],
        group_indices=simpleParameterSpace(4),
        target_indices="kn",
    )
    aderdg.db.update(folded)

    samples = OptionalDimTensor(
        "materialSamples",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (npoints,),
        alignStride=True,
    )
    faceValues = OptionalDimTensor(
        "materialAtFace",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.num2DBasisFunctions(),),
        alignStride=True,
    )
    generator.addFamily(
        "projectMaterialToFace",
        simpleParameterSpace(4),
        lambda side: faceValues["k"]
        <= aderdg.db.materialToFace[side]["kn"] * samples["n"],
    )


def addFaultKernels(generator, aderdg, matricesDir, pointSet, faultDb):
    """The material where a fault needs it.

    A fault face reads its Riemann problem at the quadrature points of the
    dynamic rupture rule, which is a different set from the nodal one the flux
    uses and has its own matrix per side and reparametrisation. The route is the
    same as for a face: the samples give a modal field and the field is read
    wherever it is wanted, so the two matrices fold into one and the modal
    coefficients are never written.
    """
    mats = tensors(matricesDir, aderdg, pointSet)
    npoints = pointCount(matricesDir, aderdg, pointSet)

    folded = tensor_collection_from_constant_expression(
        base_name="materialToFault",
        expressions=lambda side, relation: faultDb.V3mTo2n[side, relation][
            aderdg.t("kl")
        ]
        * mats["materialProject"]["ln"],
        group_indices=simpleParameterSpace(4, 4),
        target_indices="kn",
    )
    aderdg.db.update(folded)

    samples = OptionalDimTensor(
        "materialSamples",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (npoints,),
        alignStride=True,
    )
    faultValues = OptionalDimTensor(
        "materialAtFault",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.t(faultDb.V3mTo2n[0, 0].shape())[0],),
        alignStride=True,
    )
    generator.addFamily(
        "projectMaterialToFault",
        simpleParameterSpace(4, 4),
        lambda side, relation: faultValues["k"]
        <= aderdg.db.materialToFault[side, relation]["kn"] * samples["n"],
    )


#: The three reparametrisations a shared face can have, as permutations of the
#: barycentric coordinates of the reference triangle. All three are odd: two
#: tetrahedra see the face between them with opposite orientation, so the map
#: from one parametrisation to the other always reverses it, and which of the
#: three it is depends on which vertices pair up.
FACE_REFLECTIONS = ((0, 2, 1), (1, 0, 2), (2, 1, 0))


def faceOrientationPermutations(matricesDir, aderdg):
    """How the nodes of a face are renumbered between the two cells sharing it.

    A value given at the nodes of a face reaches the other cell's numbering by
    being reordered, nothing more: the nodal set is symmetric under the
    triangle's reflections, so each one maps it onto itself. The generated fP
    carries the same map in the modal basis with a mass factor -- M2 times this
    permutation -- which is what the weak form needs and what a value does not.
    """
    import numpy as np
    from yateto.input import parseJSONMatrixFile

    db = parseJSONMatrixFile(
        f"{matricesDir}/nodal/nodalBoundary_matrices_{aderdg.order}.json",
        clones=dict(),
        alignStride=aderdg.alignStride,
    )
    shape = db.nodes2D.shape()
    points = np.zeros(shape)
    for idx, value in db.nodes2D.values().items():
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
                    "the two-dimensional nodal set is not symmetric under the "
                    "reflections of the triangle, so a face cannot be renumbered"
                )
            order.append(nearest)
        permutations.append(tuple(order))
    return tuple(permutations)


def addNeighborFaceKernels(generator, aderdg, matricesDir, pointSet):
    """The neighbour's material at the nodes of the shared face.

    The flux of a face is built from the material on both sides of it, at the
    same points and in the same order. For the cell itself that is
    materialToFace; for its neighbour the face is parametrised the other way
    round, so its own evaluation is followed by the renumbering into this
    cell's ordering -- the same fold the field takes, with the projection from
    the samples in front of it.
    """
    mats = tensors(matricesDir, aderdg, pointSet)
    npoints = pointCount(matricesDir, aderdg, pointSet)
    permutations = faceOrientationPermutations(matricesDir, aderdg)
    nodes = len(permutations[0])

    renumber = [
        Tensor(
            f"materialRenumber({h})",
            (nodes, nodes),
            spp={(row, permutation[row]): "1.0" for row in range(nodes)},
        )
        for h, permutation in enumerate(permutations)
    ]

    folded = tensor_collection_from_constant_expression(
        base_name="materialNeighborToFace",
        expressions=lambda h, j: renumber[h]["km"]
        * aderdg.db.V3mTo2nFace[j][aderdg.t("ml")]
        * mats["materialProject"]["ln"],
        group_indices=simpleParameterSpace(3, 4),
        target_indices="kn",
    )
    aderdg.db.update(folded)

    samples = OptionalDimTensor(
        "materialSamples",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (npoints,),
        alignStride=True,
    )
    faceValues = OptionalDimTensor(
        "materialAtFace",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (nodes,),
        alignStride=True,
    )
    generator.addFamily(
        "projectMaterialToNeighborFace",
        simpleParameterSpace(3, 4),
        lambda h, j: faceValues["k"]
        <= aderdg.db.materialNeighborToFace[h, j]["kn"] * samples["n"],
    )


def addNeighborFaceMatrices(aderdg, matricesDir):
    """Reading a neighbour's field at the nodes of the shared face.

    Its own face evaluation, then the renumbering into this cell's ordering --
    folded into one matrix per reparametrisation and neighbour side, so that the
    kernel does one product instead of two.
    """
    permutations = faceOrientationPermutations(matricesDir, aderdg)
    nodes = len(permutations[0])

    renumber = [
        Tensor(
            f"faceRenumber({h})",
            (nodes, nodes),
            spp={(row, permutation[row]): "1.0" for row in range(nodes)},
        )
        for h, permutation in enumerate(permutations)
    ]

    folded = tensor_collection_from_constant_expression(
        base_name="neighborToFace",
        expressions=lambda h, j: renumber[h]["nm"]
        * aderdg.db.V3mTo2nFace[j][aderdg.t("ml")],
        group_indices=simpleParameterSpace(3, 4),
        target_indices="nl",
    )
    aderdg.db.update(folded)
    return nodes
