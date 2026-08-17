# SPDX-FileCopyrightText: 2024 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

from kernels.common import generate_kernel_name_prefix
from yateto import Tensor, simpleParameterSpace
from yateto.input import parseJSONMatrixFile


def volumePointCount(degree):
    return ((degree + 1) * (degree + 2) * (degree + 3)) // 6


def facePointCount(degree):
    return ((degree + 1) * (degree + 2)) // 2


def addKernels(generator, aderdg, PlasticityMethod, matricesDir, targets=["cpu"]):
    # The projection matrices are built at run time (cf. seissol::numerical::projection), since
    # they depend on the output refinement. They are therefore declared as dense tensors of the
    # right shape only -- deriving a sparsity pattern from the unrefined case would be wrong.
    #
    # The nodal-to-modal transforms (vInv, MV2nTo2m) are folded into the run-time matrices, so
    # the nodal kernels only differ from the modal ones in the number of columns.
    plasticityDB = parseJSONMatrixFile(
        f"{matricesDir}/plasticity-{PlasticityMethod}-matrices-{aderdg.order}.json",
        clones=dict(),
        alignStride=aderdg.alignStride,
    )

    order = aderdg.order
    volumeBasisCount = aderdg.numberOf3DBasisFunctions()
    faceBasisCount = aderdg.numberOf2DBasisFunctions()
    volumeNodeCount = plasticityDB.v.shape()[0]
    # the face displacement is stored at the nodes2D points; these are unisolvent, hence the
    # node count coincides with the number of 2D basis functions
    faceNodeCount = faceBasisCount

    maxOrder = 8
    rangeLimit = maxOrder + 1

    def collection(name, points, columns):
        return [
            Tensor(f"{name}({order},{i})", (points(i), columns), alignStride=True)
            for i in range(rangeLimit)
        ]

    # volume basis -> volume points
    collvv = collection("collvv", volumePointCount, volumeBasisCount)
    # face basis -> face points
    collff = collection("collff", facePointCount, faceBasisCount)
    # volume basis -> face points (the face is selected via the run-time matrix)
    collvf = collection("collvf", facePointCount, volumeBasisCount)
    # volume nodes -> volume points
    collnv = collection("collnv", volumePointCount, volumeNodeCount)
    # face nodes -> face points
    collnf = collection("collnf", facePointCount, faceNodeCount)

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)

        simcount = aderdg.multipleSimulations

        # the following is due to a shortcut in Yateto,
        # where 1-column matrices are interpreted as rank-1 vectors

        qb = Tensor("qb", (simcount, volumeBasisCount))
        qn = Tensor("qn", (simcount, volumeNodeCount))
        pb = Tensor("pb", (simcount, faceBasisCount))
        pn = Tensor("pn", (simcount, faceNodeCount))
        xv = [Tensor(f"xv({i})", (volumePointCount(i),)) for i in range(rangeLimit)]
        xf = [Tensor(f"xf({i})", (facePointCount(i),)) for i in range(rangeLimit)]

        simselect = Tensor("simselect", (simcount,))

        generator.addFamily(
            f"{name_prefix}projectBasisToVtkVolume",
            simpleParameterSpace(rangeLimit),
            lambda i: xv[i]["p"] <= collvv[i]["pb"] * simselect["s"] * qb["sb"],
            target=target,
        )
        generator.addFamily(
            f"{name_prefix}projectNodalToVtkVolume",
            simpleParameterSpace(rangeLimit),
            lambda i: xv[i]["p"] <= collnv[i]["pn"] * simselect["s"] * qn["sn"],
            target=target,
        )
        generator.addFamily(
            f"{name_prefix}projectBasisToVtkFace",
            simpleParameterSpace(rangeLimit),
            lambda i: xf[i]["p"] <= collff[i]["pb"] * simselect["s"] * pb["sb"],
            target=target,
        )
        generator.addFamily(
            f"{name_prefix}projectNodalToVtkFace",
            simpleParameterSpace(rangeLimit),
            lambda i: xf[i]["p"] <= collnf[i]["pn"] * simselect["s"] * pn["sn"],
            target=target,
        )
        generator.addFamily(
            f"{name_prefix}projectBasisToVtkFaceFromVolume",
            simpleParameterSpace(rangeLimit),
            lambda i: xf[i]["p"] <= collvf[i]["pb"] * simselect["s"] * qb["sb"],
            target=target,
        )


def includeTensors(matricesDir, includeTensors):
    vtkbase = parseJSONMatrixFile(f"{matricesDir}/vtkbase.json")
    for x in vtkbase.__dict__:
        if isinstance(vtkbase.__dict__[x], dict):
            for y in vtkbase.__dict__[x]:
                includeTensors.add(vtkbase.__dict__[x][y])
        else:
            includeTensors.add(vtkbase.__dict__[x])
