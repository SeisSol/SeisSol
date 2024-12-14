
def plasticity(generator, aderdg):
    cohesionTimesCosAngularFriction = forge.scalar()
    sinAngularFriction = forge.scalar()
    initialLoading = forge.scalar()

    qplastic = forge.tensor((aderdg.basisFunctions(), 7))

    qstress = qplastic[:,:6]
    qeta = qplastic[:,6]

    # forge.layout_like(qplastic, 'qplastic')

    # forge.layout(qplastic, 'qplastic')

    nodalStressTensors = vandermondeMatrix @ qstress
    etaNodal = vandermondeMatrix @ qeta

    # m2n
    # plComputeMean
    # plSubtractMean
    # plComputeSecondInvariant
    # tau = sqrt(secondInvariant)
    # taulim[ip] = std::max((real) 0.0, plasticityData->cohesionTimesCosAngularFriction - meanStress[ip] * plasticityData->sinAngularFriction);
    # bool adjust = false;
    # scalar: if (tau[ip] > taulim[ip]) { yieldFactor[ip] = (taulim[ip] / tau[ip] - 1.0) * oneMinusIntegratingFactor; } else {yieldFactor[ip] = 0.0;}
    # adjust = reduction: any (tau[ip] > taulim[ip])

    meanStress = forge.sum(nodalStressTensors[0:3,:]) / 3
    localStresses = nodalStressTensors[0:3,:] - meanStress

    secondInvariant = 0.5 * forge.sum(localStresses**2) + forge.sum(nodalStressTensors[3:6,:]**2)

    tau = forge.sqrt(secondInvariant)
    taulim = forge.max(0, cohesionTimesCosAngularFriction - meanStress * sinAngularFriction)

    localAdjust = tau > taulim
    yieldFactor = forge.where(localAdjust, (taulim / tau - 1) * oneMinusIntegratingFactor, 0)
    adjust = forge.any(localAdjust)

    # 0.5 + forge.sqrt(forge.sum(TODO))

    with forge.block_conditional(adjust):
        pass

    factor = mufactor / (T_v * oneMinusIntegratingFactor)
    nodeDuDtPstrain = factor * (localPrevDofs['PS'] - localDofs['PS'])

    # adj
    # plastic strain
    # plConvertToNodalNoLoading
    # plConvertEtaModal2Nodal
    # compute eta
    # plConvertEtaNodal2Modal
