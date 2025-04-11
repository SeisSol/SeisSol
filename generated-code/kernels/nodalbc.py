import numpy as np
import kernels.equations.viscoelastic2
from kernels.multsim import OptionalDimTensor
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.util import tensor_collection_from_constant_expression
from kernels.common import generate_kernel_name_prefix

def addKernels(generator, aderdg, include_tensors, matricesDir, dynamicRuptureMethod, targets):
    easi_ident_map = np.stack([np.eye(aderdg.numberOfQuantities())] * aderdg.numberOf2DBasisFunctions(), axis=2)
    assert(easi_ident_map.shape ==
           (aderdg.numberOfQuantities(), aderdg.numberOfQuantities(), aderdg.numberOf2DBasisFunctions()))
    if aderdg.multipleSimulations > 1:
        easi_ident_map_spp = np.stack([easi_ident_map] * aderdg.multipleSimulations, axis=0)
    else:
        easi_ident_map_spp = easi_ident_map

    easi_ident_map = OptionalDimTensor('easiIdentMap', 
            aderdg.Q.optName(), 
            aderdg.multipleSimulations, 
            0, 
            easi_ident_map.shape, 
            spp=easi_ident_map_spp,
            alignStride=False)
    easi_boundary_constant = Tensor('easiBoundaryConstant',
                                    (aderdg.numberOfQuantities(), aderdg.numberOf2DBasisFunctions()),
                                    alignStride=False)
    easi_boundary_map = Tensor('easiBoundaryMap',
                               (aderdg.numberOfQuantities(), aderdg.numberOfQuantities(), aderdg.numberOf2DBasisFunctions(),),
                               alignStride=False)
    create_easi_boundary_ghost_cells = (
            aderdg.INodal['la'] <= easi_boundary_map['abl'] * aderdg.INodal['lb'] + easi_ident_map['abl'] * easi_boundary_constant['bl']
    )
    generator.add('createEasiBoundaryGhostCells', create_easi_boundary_ghost_cells)

    projectToNodalBoundary = (
        lambda j: aderdg.INodal["kp"]
        <= aderdg.db.V3mTo2nFace[j][aderdg.t("km")] * aderdg.I["mp"]
    )

    generator.addFamily('projectToNodalBoundary',
                        simpleParameterSpace(4),
                        projectToNodalBoundary)

    for target in targets:
        name_prefix = generate_kernel_name_prefix(target)
        projectToNodalBoundaryRotated = (
            lambda j: aderdg.INodal["kp"]
            <= aderdg.db.V3mTo2nFace[j][aderdg.t("kl")]
            * aderdg.I["lm"]
            * aderdg.Tinv["pm"]
        )
        
        generator.addFamily(f'{name_prefix}projectToNodalBoundaryRotated',
                          simpleParameterSpace(4),
                          projectToNodalBoundaryRotated,
                          target=target)

        projectDerivativeToNodalBoundaryRotated = (
            lambda i, j: aderdg.INodal["kp"]
            <= aderdg.db.V3mTo2nFace[j][aderdg.t("kl")]
            * aderdg.dQs[i]["lm"]
            * aderdg.Tinv["pm"]
        )

        generator.addFamily(f'{name_prefix}projectDerivativeToNodalBoundaryRotated',
                          simpleParameterSpace(aderdg.order, 4),
                          projectDerivativeToNodalBoundaryRotated,
                          target=target)

    # To be used as Tinv in flux solver - this way we can save two rotations for
    # Dirichlet boundary, as ghost cell dofs are already rotated

    identity_rotation = np.double(aderdg.transformation_spp())
    quantities = aderdg.numberOfQuantities()
    identity_rotation[0:quantities, 0:quantities] = np.eye(quantities)
    identity_rotation = Tensor('identityT',
                               aderdg.transformation_spp().shape,
                               identity_rotation,
                               )
    include_tensors.add(identity_rotation)

    aderdg.INodalUpdate = OptionalDimTensor('INodalUpdate',
                                          aderdg.INodal.optName(),
                                          aderdg.INodal.optSize(),
                                          aderdg.INodal.optPos(),
                                          (aderdg.numberOf2DBasisFunctions(), aderdg.numberOfQuantities()),
                                          alignStride=True)

    factor = Scalar('factor')
    updateINodal = aderdg.INodal['kp'] <= aderdg.INodal['kp'] + factor * aderdg.INodalUpdate['kp']
    generator.add('updateINodal', updateINodal)
