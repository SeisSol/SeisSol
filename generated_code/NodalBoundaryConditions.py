import numpy as np
import viscoelastic2
from multSim import OptionalDimTensor
from yateto import Tensor, Scalar, simpleParameterSpace
from yateto.memory import CSCMemoryLayout
from yateto.util import tensor_collection_from_constant_expression

def addKernels(generator, aderdg, include_tensors, matricesDir, dynamicRuptureMethod):
    easi_ident_map = np.stack([np.eye(aderdg.numberOfQuantities())] * aderdg.numberOf2DBasisFunctions(), axis=2)
    assert(easi_ident_map.shape ==
           (aderdg.numberOfQuantities(), aderdg.numberOfQuantities(), aderdg.numberOf2DBasisFunctions()))
    easi_ident_map = Tensor('easiIdentMap',
                            easi_ident_map.shape,
                            easi_ident_map,
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

    projectToNodalBoundary = lambda j: aderdg.INodal['kp'] <= aderdg.db.V3mTo2nFace[j]['km'] * aderdg.I['mp']

    generator.addFamily('projectToNodalBoundary',
                        simpleParameterSpace(4),
                        projectToNodalBoundary)

    projectToNodalBoundaryRotated = lambda j: aderdg.INodal['kp'] <= aderdg.db.V3mTo2nFace[j]['kl'] \
                                              * aderdg.I['lm'] \
                                              * aderdg.Tinv['pm']

    generator.addFamily('projectToNodalBoundaryRotated',
                        simpleParameterSpace(4),
                        projectToNodalBoundaryRotated)

    # To be used as Tinv in flux solver - this way we can save two rotations for
    # Dirichlet boundary, as ghost cell dofs are already rotated
    identity_rotation = np.double(aderdg.transformation_spp())
    identity_rotation[0:9, 0:9] = np.eye(9)
    identity_rotation = Tensor('identityT',
                               aderdg.transformation_spp().shape,
                               identity_rotation,
                               )
    include_tensors.add(identity_rotation)

    rDivM_mult_V2nTo2m = tensor_collection_from_constant_expression(
        base_name='rDivMMultV2nTo2m',
        expressions=lambda i: aderdg.db.rDivM[i]['jk'] * aderdg.db.V2nTo2m['kl'],
        group_indices=range(4),
        target_indices='jl')

    aderdg.db.update(rDivM_mult_V2nTo2m)

    if not isinstance(aderdg, viscoelastic2.Viscoelastic2ADERDG):
        localFluxNodal = lambda i: aderdg.Q['kp'] <= aderdg.Q['kp'] + aderdg.db.rDivMMultV2nTo2m[i]['kn'] * aderdg.INodal['no'] * aderdg.AminusT['op']
        localFluxNodalPrefetch = lambda i: aderdg.I if i == 0 else (aderdg.Q if i == 1 else None)
        generator.addFamily('localFluxNodal', simpleParameterSpace(4), localFluxNodal, localFluxNodalPrefetch)
    else:
        # Nodal bc not supported for visc2, but we  need to generate rDivM_mult_V2nTo2m.
        # include_tensors doesnt allow tensor families directly.
        include_tensors.update([
            rDivM_mult_V2nTo2m["rDivMMultV2nTo2m"][i] for i in range(4)
        ])

    selectZDisplacement = np.zeros((aderdg.numberOfQuantities(), 1))
    selectZDisplacement[8, 0] = 1
    selectZDisplacement = Tensor('selectZDisplacement',
                                 selectZDisplacement.shape,
                                 selectZDisplacement,
                                 CSCMemoryLayout)

    selectZDisplacementFromDisplacements = np.zeros((3, 1))
    selectZDisplacementFromDisplacements[2, 0] = 1
    selectZDisplacementFromDisplacements = Tensor('selectZDisplacementFromDisplacements',
                                                  selectZDisplacementFromDisplacements.shape,
                                                  selectZDisplacementFromDisplacements,
                                                  CSCMemoryLayout)

    aderdg.INodalDisplacement = OptionalDimTensor('INodalDisplacement',
                                                aderdg.Q.optName(),
                                                aderdg.Q.optSize(),
                                                aderdg.Q.optPos(),
                                                (aderdg.numberOf2DBasisFunctions(), 1),
                                                alignStride=True)

    displacement = OptionalDimTensor('displacement',
                                     aderdg.Q.optName(),
                                     aderdg.Q.optSize(),
                                     aderdg.Q.optPos(),
                                     (aderdg.numberOf3DBasisFunctions(), 3),
                                     alignStride=True)

    dt = Scalar('dt')
    displacementAvgNodal = lambda side: aderdg.INodalDisplacement['ip'] <= aderdg.db.V3mTo2nFace[side]['ij'] * aderdg.I['jk'] * selectZDisplacement['kp'] \
                                        + dt * aderdg.db.V3mTo2nFace[side]['ij'] * displacement['jk'] * selectZDisplacementFromDisplacements['kp']

    generator.addFamily('displacementAvgNodal',
                        simpleParameterSpace(4),
                        displacementAvgNodal)

    aderdg.INodalUpdate = OptionalDimTensor('INodalUpdate',
                                          aderdg.INodal.optName(),
                                          aderdg.INodal.optSize(),
                                          aderdg.INodal.optPos(),
                                          (aderdg.numberOf2DBasisFunctions(), aderdg.numberOfQuantities()),
                                          alignStride=True)

    factor = Scalar('factor')
    updateINodal = aderdg.INodal['kp'] <= aderdg.INodal['kp'] + factor * aderdg.INodalUpdate['kp']
    generator.add('updateINodal', updateINodal)