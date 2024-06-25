def adjustDeviatoricTensors():
  localStresses3 = nodalStressTensors.slice(1, 0, 3)['PS']

  # 1. Compute the mean stress for each node
  meanStress['P'] <= localStresses3['PS'] / 3.0

  # 2. Compute deviatoric stress tensor
  localStresses3['PS'] <= localStresses3['PS'] - meanStress['P']

  localStresses6 = nodalStressTensors.slice(1, 3, 6)['PZ']

  # 3. Compute the second invariant for each node
  tau['P'] <= 0.5 * localStresses3['PS'] * localStresses3['PS'] + localStresses6['PZ'] * localStresses6['PZ']

  forge.scalarblock(forge.assign(tau, forge.sqrt(tau)))

  # 4. Compute the plasticity criteria
  cohesionTimesCosAngularFriction = cohesionTimesCosAngularFriction['']
  sinAngularFriction = sinAngularFriction['']
  forge.scalarblock(forge.assign(taulim, forge.max(0, cohesionTimesCosAngularFriction - meanStress * sinAngularFriction)))

  # 5. Compute the yield factor
  forge.scalarblock(forge.assign(factor, forge.ternary(tau['P'] > taulim, ((taulim / tau) - 1.0) * oneMinusIntegratingFactor, 0)))

  isAdjusted[''] <= tau['P'] > taulim
  isAdjustableVector[''] <= isAdjusted

  # 6. Adjust deviatoric stress tensor if a node within a node exceeds the elasticity region
  forge.conditionalblock(
    isAdjusted[''],
    [
      elementTensors.slice(1, 0, 3)['PS'] <= localStresses3['PS'] * factor,
      elementTensors.slice(1, 3, 6)['PS'] <= localStresses6['PS'] * factor
    ]
  )

def computePstrains():
  factor = mufactor / (T_v * oneMinusIntegratingFactor)
  nodeDuDtPstrain = factor * (localPrevDofs['PS'] - localDofs['PS'])
  forge.conditionalblock(
    isAdjustableVector[''],
    [
      localPstream <= localPstrain + timeStepWidth * nodeDuDtPstrain,
      localDuDtPstrain['PS'] <= nodeDuDtPstrain,
    ]
  )

def pstrainToQEtaModal():
  forge.conditionalblock(
    isAdjustableVector[''],
    [
      localQEtaModal['P'] <= localPstream.slice(0, 6)['PI']
    ]
  )

def qEtaModalToPstrain():
  forge.conditionalblock(
    isAdjustableVector[''],
    [
      localPstream.slice(0, 6)['PI'] <= localQEtaModal['P']
    ]
  )

def updateQEtaNodal():
  factor = localQStressNodal['PS'] * localQStressNodal['PS']
  forge.conditionalblock(
    isAdjustableVector[''],
    [
      forge.scalarblock(forge.assign(localQEtaNodal, forge.max(0, localQEtaNodal) + timeStepWidth * forge.sqrt(0.5 * factor)))
    ]
  )
