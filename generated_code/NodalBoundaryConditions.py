import numpy as np
import viscoelastic2
from multSim import OptionalDimTensor
from yateto import Tensor, Scalar, simpleParameterSpace
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

    project2nFaceTo3m = tensor_collection_from_constant_expression(
        base_name='project2nFaceTo3m',
        expressions=lambda i: aderdg.db.rDivM[i]['jk'] * aderdg.db.V2nTo2m['kl'],
        group_indices=range(4),
        target_indices='jl')

    aderdg.db.update(project2nFaceTo3m)

    selectZDisplacementFromQuantities = np.zeros(aderdg.numberOfQuantities())
    selectZDisplacementFromQuantities[8] = 1
    selectZDisplacementFromQuantities= Tensor('selectZDisplacementFromQuantities',
                                              selectZDisplacementFromQuantities.shape,
                                              selectZDisplacementFromQuantities,
                                              )

    selectZDisplacementFromDisplacements = np.zeros(3)
    selectZDisplacementFromDisplacements[2] = 1
    selectZDisplacementFromDisplacements = Tensor('selectZDisplacementFromDisplacements',
                                                  selectZDisplacementFromDisplacements.shape,
                                                  selectZDisplacementFromDisplacements,
                                                  )

    aderdg.INodalDisplacement = OptionalDimTensor('INodalDisplacement',
                                                aderdg.Q.optName(),
                                                aderdg.Q.optSize(),
                                                aderdg.Q.optPos(),
                                                (aderdg.numberOf2DBasisFunctions(),),
                                                alignStride=True)

    displacement = OptionalDimTensor('displacement',
                                     aderdg.Q.optName(),
                                     aderdg.Q.optSize(),
                                     aderdg.Q.optPos(),
                                     (aderdg.numberOf3DBasisFunctions(), 3),
                                     alignStride=True)

    dt = Scalar('dt')
    displacementAvgNodal = lambda side: aderdg.INodalDisplacement['i'] <= \
                                        aderdg.db.V3mTo2nFace[side]['ij'] * aderdg.I['jk'] * selectZDisplacementFromQuantities['k'] \
                                        + dt * aderdg.db.V3mTo2nFace[side]['ij'] * displacement['jk'] * selectZDisplacementFromDisplacements['k']

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