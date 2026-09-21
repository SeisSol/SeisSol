# SPDX-FileCopyrightText: 2017 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from kernels.common import generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from yateto import Scalar, Tensor, simpleParameterSpace
from yateto.util import tensor_collection_from_constant_expression


def addKernels(
    generator,
    aderdg,
    include_tensors,
    matricesDir,
    dynamicRuptureMethod,
    targets,
):
    dirichlet_offset = OptionalDimTensor(
        "dirichletOffset",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.numberOfQuantities(),),
        alignStride=True,
    )

    dirichlet_map = Tensor(
        "dirichletMap",
        (aderdg.numberOfQuantities(), aderdg.numberOfQuantities()),
        alignStride=False,
    )

    # The boundary condition is given in global coordinates; the face-aligned
    # form that the flux solver absorbs is derived from it once per face.
    dirichlet_offset_global = OptionalDimTensor(
        "dirichletOffsetGlobal",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.numberOfQuantities(),),
        alignStride=True,
    )

    dirichlet_map_global = Tensor(
        "dirichletMapGlobal",
        (aderdg.numberOfQuantities(), aderdg.numberOfQuantities()),
        alignStride=False,
    )

    # The boundary condition acts on the quantities that enter the Riemann
    # problem, which is the leading block of the rotation for materials that
    # carry more quantities than that.
    nq = aderdg.numberOfQuantities()

    generator.add(
        "rotateBoundaryCondition",
        [
            dirichlet_map["ab"]
            <= aderdg.Tinv["ac"].subslice("a", 0, nq).subslice("c", 0, nq)
            * dirichlet_map_global["cd"]
            * aderdg.T["db"].subslice("d", 0, nq).subslice("b", 0, nq),
            dirichlet_offset["a"]
            <= aderdg.Tinv["am"].subslice("a", 0, nq).subslice("m", 0, nq)
            * dirichlet_offset_global["m"],
        ],
    )

    rho = Tensor("rho", ())

    mainstresscnt = 3 if aderdg.velocityOffset() > 1 else 1

    averageNormalDisplacement = OptionalDimTensor(
        "averageNormalDisplacement",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.numberOf2DBasisFunctions(),),
        alignStride=True,
    )

    g2m = Scalar("g2m")  # -2 * g
    dt = Scalar("dt")

    main_stress_select = np.zeros(aderdg.numberOfQuantities())
    main_stress_select[0:mainstresscnt] = 1.0
    main_stress_select = Tensor(
        "mainStressSelect", main_stress_select.shape, main_stress_select
    )

    # The free-surface-gravity map is diag(-1, ..., -1, 1, ..., 1) in the
    # face-aligned basis and hence constant over the face. Folding it into the
    # local flux solver turns the boundary into an ordinary local flux; only the
    # displacement-driven offset is left over, and that one is rank one.
    fsg_map = np.eye(aderdg.numberOfQuantities())
    for i in range(mainstresscnt):
        fsg_map[i, i] = -1.0
    fsg_map = Tensor("fsgMap", fsg_map.shape, fsg_map)

    # The Dirichlet map is constant over the face as well, so it folds the same
    # way; the offset then lifts a constant nodal function, which is a fixed
    # vector per face.
    face_node_sum = Tensor(
        "faceNodeSum",
        (aderdg.numberOf2DBasisFunctions(),),
        np.ones(aderdg.numberOf2DBasisFunctions()),
    )
    dirichlet_lift = tensor_collection_from_constant_expression(
        base_name="dirichletLift",
        expressions=lambda i: aderdg.db.project2nFaceTo3m[i]["kn"] * face_node_sum["n"],
        group_indices=simpleParameterSpace(4),
        target_indices="k",
    )
    aderdg.db.update(dirichlet_lift)

    fold_dirichlet = (
        aderdg.AplusT["mp"]
        <= aderdg.AplusT["mp"]
        + aderdg.Tinv["bm"].subslice("b", 0, nq).subslice("m", 0, nq)
        * dirichlet_map["ab"]
        * aderdg.AminusT["ap"]
    )
    generator.add("foldDirichlet", fold_dirichlet)

    fold_free_surface_gravity = (
        aderdg.AplusT["mp"]
        <= aderdg.AplusT["mp"]
        + aderdg.Tinv["om"].subslice("o", 0, nq).subslice("m", 0, nq)
        * fsg_map["oq"]
        * aderdg.AminusT["qp"]
    )
    generator.add("foldFreeSurfaceGravity", fold_free_surface_gravity)

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        dirichlet_flux = (
            lambda i: aderdg.extendedQTensor()["kp"]
            <= aderdg.extendedQTensor()["kp"]
            + dt
            * aderdg.db.dirichletLift[i]["k"]
            * dirichlet_offset["o"]
            * aderdg.AminusT["op"]
        )

        fsg_flux = (
            lambda i: aderdg.extendedQTensor()["kp"]
            <= aderdg.extendedQTensor()["kp"]
            + g2m
            * rho[""]
            * aderdg.db.project2nFaceTo3m[i]["kn"]
            * averageNormalDisplacement["n"]
            * main_stress_select["o"]
            * aderdg.AminusT["op"]
        )

        generator.addFamily(
            f"{name_prefix}dirichletFlux",
            simpleParameterSpace(4),
            dirichlet_flux,
            target=target,
        )

        generator.addFamily(
            f"{name_prefix}fsgFlux",
            simpleParameterSpace(4),
            fsg_flux,
            target=target,
        )

    # To be used as Tinv in flux solver - this way we can save two rotations
    # for the Dirichlet boundary, as ghost cell dofs are already rotated
    identity_rotation = np.double(aderdg.transformation_spp())
    quantities = aderdg.numberOfQuantities()
    identity_rotation[0:quantities, 0:quantities] = np.eye(quantities)
    identity_rotation = Tensor(
        "identityT",
        aderdg.transformation_spp().shape,
        identity_rotation,
    )
    include_tensors.add(identity_rotation)

    aderdg.INodalUpdate = OptionalDimTensor(
        "INodalUpdate",
        aderdg.INodal.optName(),
        aderdg.INodal.optSize(),
        aderdg.INodal.optPos(),
        (aderdg.numberOf2DBasisFunctions(), aderdg.numberOfQuantities()),
        alignStride=True,
    )

    factor = Scalar("factor")
    updateINodal = (
        aderdg.INodal["kp"] <= aderdg.INodal["kp"] + factor * aderdg.INodalUpdate["kp"]
    )
    generator.add("updateINodal", updateINodal)
