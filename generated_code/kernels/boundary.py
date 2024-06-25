

def freeSurfaceGravity():
    elementBoundaryDofs['CP'] <= -2 * rho[''] * g * elementDisplacement['P'] - elementBoundaryDofs['CP']

def extractRotationMatrices():
    displacementToFaceNormal['ij'] <= Tinv.slice(0, 6, 9).slice(1, 6, 9)['ij']
    displacementToGlobalData['ij'] <= T.slice(0, 6, 9).slice(1, 6, 9)['ij']

def computeInvAcousticImpedance():
    forge.scalarblock(forge.assign(forge.tensor(invImpedances['']), 1 / forge.sqrt(forge.tensor(lambdas['']), forge.tensor(rhos['']))))

def initializeTaylorSeriesForGravitationalBoundary():
    prevCoefficients['P'] = rotatedFaceDisplacement['P']
    integratedDisplacementNodal['P'] = deltaTInt * rotatedFaceDisplacement['P']

def updateRotatedFaceDisplacement(elastic):
    # TODO: factorEvaluated, factorInt become family

    uInside = dofsFaceNodal.slice(1, 6)['IP']
    vInside = dofsFaceNodal.slice(1, 7)['IP']
    wInside = dofsFaceNodal.slice(1, 8)['IP']
    pressureInside = dofsFaceNodal.slice(1, 0)['IP']
    prevCoefficients
    if elastic:
        curCoeff = uInside - invImpedance[''] * (rho[''] * g * prevCoefficients['P'] + pressureInside)
    else:
        curCoeff = uInside
    
    prevCoefficients['P'] <= curCoeff

    rotatedFaceDisplacement.slice(0, 0)['IP'] <= rotatedFaceDisplacement.slice(0, 0)['IP'] + factorEvaluated * curCoeff
    rotatedFaceDisplacement.slice(0, 1)['IP'] <= rotatedFaceDisplacement.slice(0, 1)['IP'] + factorEvaluated * vInside
    rotatedFaceDisplacement.slice(0, 2)['IP'] <= rotatedFaceDisplacement.slice(0, 2)['IP'] + factorEvaluated * wInside

    integratedDisplacementNodal <= integratedDisplacementNodal + factorInt * curCoeff['P']
