# SPDX-FileCopyrightText: 2017 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

import numpy as np
from kernels.common import cold_kernel_attrs, generate_kernel_name_prefix
from kernels.multsim import OptionalDimTensor
from kernels.quantities import layout, total_extent
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
        (aderdg.numQuantities(),),
        alignStride=True,
    )

    dirichlet_map = Tensor(
        "dirichletMap",
        (aderdg.numQuantities(), aderdg.numQuantities()),
        alignStride=False,
    )

    # The boundary condition is given in global coordinates; the face-aligned
    # form that the flux solver absorbs is derived from it once per face.
    dirichlet_offset_global = OptionalDimTensor(
        "dirichletOffsetGlobal",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.numQuantities(),),
        alignStride=True,
    )

    dirichlet_map_global = Tensor(
        "dirichletMapGlobal",
        (aderdg.numQuantities(), aderdg.numQuantities()),
        alignStride=False,
    )

    # The boundary condition acts on the quantities that enter the Riemann
    # problem, which is the leading block of the rotation for materials that
    # carry more quantities than that.
    nq = aderdg.numQuantities()

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
        attrs=cold_kernel_attrs(),
    )

    rho = Tensor("rho", ())

    mainstresscnt = 3 if aderdg.velocityOffset() > 1 else 1

    averageNormalDisplacement = OptionalDimTensor(
        "averageNormalDisplacement",
        aderdg.Q.optName(),
        aderdg.Q.optSize(),
        aderdg.Q.optPos(),
        (aderdg.num2DBasisFunctions(),),
        alignStride=True,
    )

    g2m = Scalar("g2m")  # -2 * g
    dt = Scalar("dt")

    main_stress_select = np.zeros(aderdg.numQuantities())
    main_stress_select[0:mainstresscnt] = 1.0
    main_stress_select = Tensor(
        "mainStressSelect", main_stress_select.shape, main_stress_select
    )

    # The free-surface-gravity map is diag(-1, ..., -1, 1, ..., 1) in the
    # face-aligned basis and hence constant over the face. Folding it into the
    # local flux solver turns the boundary into an ordinary local flux; only the
    # displacement-driven offset is left over, and that one is rank one.
    fsg_map = np.eye(aderdg.numQuantities())
    for i in range(mainstresscnt):
        fsg_map[i, i] = -1.0
    fsg_map = Tensor("fsgMap", fsg_map.shape, fsg_map)

    # The Dirichlet map is constant over the face as well, so it folds the same
    # way; the offset then lifts a constant nodal function, which is a fixed
    # vector per face.
    face_node_sum = Tensor(
        "faceNodeSum",
        (aderdg.num2DBasisFunctions(),),
        np.ones(aderdg.num2DBasisFunctions()),
    )
    dirichlet_lift = tensor_collection_from_constant_expression(
        base_name="dirichletLift",
        expressions=lambda i: aderdg.db.project2nFaceTo3m[i]["kn"] * face_node_sum["n"],
        group_indices=simpleParameterSpace(4),
        target_indices="k",
    )
    aderdg.db.update(dirichlet_lift)

    # The fold lands in the rows of the local flux solver, and only the
    # quantities of the Riemann problem have rows there: a material that keeps
    # its relaxation in Q (the fused anelastic solver) holds the rest
    # structurally zero, so the map can only act on that block.
    nr = total_extent(layout(aderdg.primaryGroups()))
    fold_dirichlet = aderdg.AplusT["mp"].subslice("m", 0, nr) <= (
        aderdg.AplusT["mp"].subslice("m", 0, nr)
        + aderdg.Tinv["bm"].subslice("b", 0, nr).subslice("m", 0, nr)
        * dirichlet_map["ab"].subslice("a", 0, nr).subslice("b", 0, nr)
        * aderdg.AminusT["ap"].subslice("a", 0, nr)
    )
    generator.add("foldDirichlet", fold_dirichlet, attrs=cold_kernel_attrs())

    fold_free_surface_gravity = (
        aderdg.AplusT["mp"]
        <= aderdg.AplusT["mp"]
        + aderdg.Tinv["om"].subslice("o", 0, nq).subslice("m", 0, nq)
        * fsg_map["oq"]
        * aderdg.AminusT["qp"]
    )
    generator.add(
        "foldFreeSurfaceGravity",
        fold_free_surface_gravity,
        attrs=cold_kernel_attrs(),
    )

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
    # for the Dirichlet boundary, as ghost cell dofs are already rotated.
    # It stands in for Tinv, so it has to be shaped like Tinv: the two need not
    # agree, and where they do not, the flux solver would read this with the
    # wrong stride.
    inv_spp = aderdg.transformation_inv_spp()
    identity_rotation = np.double(inv_spp)
    quantities = min(aderdg.numQuantities(), inv_spp.shape[0])
    identity_rotation[0:quantities, 0:quantities] = np.eye(quantities)
    identity_rotation = Tensor(
        "identityT",
        inv_spp.shape,
        identity_rotation,
    )
    include_tensors.add(identity_rotation)

    aderdg.INodalUpdate = OptionalDimTensor(
        "INodalUpdate",
        aderdg.INodal.optName(),
        aderdg.INodal.optSize(),
        aderdg.INodal.optPos(),
        (aderdg.num2DBasisFunctions(), aderdg.numQuantities()),
        alignStride=True,
    )

    factor = Scalar("factor")
    updateINodal = (
        aderdg.INodal["kp"] <= aderdg.INodal["kp"] + factor * aderdg.INodalUpdate["kp"]
    )
    generator.add("updateINodal", updateINodal)
