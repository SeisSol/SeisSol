import tensorforge.functions as forge
from yateto.type import Tensor, Scalar

def smoothStep(currentTime, t0):
    tau = currentTime - t0
    return forge.ternary(currentTime <= 0, 0, forge.ternary(currentTime < t0, forge.exp(tau**2 / (currentTime * (currentTime - 2 * t0))), 1))

def smoothStepIncrement(currentTime, dt, t0):
    return smoothStep(currentTime, t0) - smoothStep(currentTime - dt, t0)

class FrictionLawBase:
    def __init__(self, aderdg, numberOfPoints):
        self.numberOfPoints = numberOfPoints
        self.plusQ = self.drvarcoll('plusQ', (aderdg.order, aderdg.numberOfQuantities(),))
        self.minusQ = self.drvarcoll('minusQ', (aderdg.order, aderdg.numberOfQuantities(),))
        self.plusImposedState = self.drvarcoll('plusImposedState', (aderdg.numberOfQuantities(),))
        self.minusImposedState = self.drvarcoll('minusImposedState', (aderdg.numberOfQuantities(),))
        self.aderdg = aderdg

        self.fullUpdateTime = self.drscalar('fullUpdateTime')
        self.t0 = self.drscalar('t0')
        self.dt = []
        self.timeWeights = []
        for i in range(self.aderdg.order):
            self.dt += [self.drscalar(f'dt({i})')]
            self.timeWeights += [self.drscalar(f'timeWeights({i})')]
        self.sumDt = self.drscalar('sumDt')
        
        self.etaP = self.drscalar('etaP')
        self.etaS = self.drscalar('etaS')
        self.zp = self.drscalar('zp')
        self.zs = self.drscalar('zs')
        self.zpNeig = self.drscalar('zpNeig')
        self.zsNeig = self.drscalar('zsNeig')
        self.invEtaS = self.drscalar('invEtaS')
        self.invZp = self.drscalar('invZp')
        self.invZs = self.drscalar('invZs')
        self.invZpNeig = self.drscalar('invZpNeig')
        self.invZsNeig = self.drscalar('invZsNeig')
        
        self.initialStressInFaultCS = []
        self.nucleationStressInFaultCS = []
        for i in range(6):
            self.initialStressInFaultCS += [self.drvar(f'initialStressInFaultCS({i})')]
            self.nucleationStressInFaultCS += [self.drvar(f'nucleationStressInFaultCS({i})')]
        self.initialPressure = self.drvar('initialPressure')
        self.nucleationPressure = self.drvar('nucleationPressure')
        self.mu = self.drvar('mu')
        self.stateVariable = self.drvar('stateVariable')
        self.dynStressTime = self.drvar('dynStressTime')
        self.dynStressTimePending = self.drvar('dynStressTimePending')
        self.peakSlipRate = self.drvar('peakSlipRate')
        self.slipRateMagnitude = self.drvar('slipRateMagnitude')
        self.accumulatedSlipMagnitude = self.drvar('accumulatedSlipMagnitude')
        self.accumulatedSlip = self.drvar('accumulatedSlip')
        self.frictionalEnergy = self.drvar('frictionalEnergy')
        self.slip1 = self.drvar('slip1')
        self.slip2 = self.drvar('slip2')
        self.traction1 = self.drvar('traction1')
        self.traction2 = self.drvar('traction2')
        self.slipRate1 = self.drvar('slipRate1')
        self.slipRate2 = self.drvar('slipRate2')
        self.strength = self.drvar('strength')

        self.faultStressesNormalStress = self.drvarcoll('faultStressesNormalStress', (aderdg.order,))
        self.faultStressesTraction1 = self.drvarcoll('faultStressesTraction1', (aderdg.order,))
        self.faultStressesTraction2 = self.drvarcoll('faultStressesTraction2', (aderdg.order,))

        self.tractionResultsTraction1 = self.drvarcoll('tractionResultsTraction1', (aderdg.order,))
        self.tractionResultsTraction2 = self.drvarcoll('tractionResultsTraction2', (aderdg.order,))

        self.ruptureTime = self.drvar('ruptureTime')
        self.ruptureTimePending = self.drvar('ruptureTimePending')

        self.terminatorSlipRateThreshold = self.drscalar('terminatorSlipRateThreshold')

    def drvar(self, name):
        return forge.tensor(Tensor(name, tuple([self.numberOfPoints]))['P']) # TODO: align stride
    
    def drvarcoll(self, name, subshape):
        tensor = Tensor(name, tuple(list(subshape) + [self.numberOfPoints]))
        spareIndices = 'ZYXVUW'
        def handleDim(curr, indices, dimid):
            if dimid == len(subshape):
                return forge.tensor(curr[indices + 'P'])
            else:
                return [handleDim(curr.slice(dimid, i), indices + spareIndices[dimid], dimid+1) for i in range(subshape[dimid])]
        return handleDim(tensor, '', 0)
    
    def drscalar(self, name):
        return forge.scalar(Scalar(name))
    
    def generate(self, routineName, generator):
        routine = []
        self.precomputeStressFromQInterpolated(routine)
        self.preHook(routine)
        for i in range(self.aderdg.order):
            self.adjustInitialStress(routine, i)
            self.updateFrictionAndSlip(routine, i)
        self.postHook(routine)
        self.saveRuptureFrontOutput(routine)
        self.saveDynamicStressOutput(routine)
        self.savePeakSlipRateOutput(routine)
        self.postcomputeImposedStateFromNewStress(routine)
        self.computeFrictionEnergy(routine)

        generator.add(routineName, forge.scalarblock(routine), target='gpu')
    
    def precomputeStressFromQInterpolated(self, routine):
        U,V,W,N,T1,T2 = 6,7,8,0,3,5
        for i in range(self.aderdg.order):
            routine += [
                forge.assign(self.faultStressesNormalStress[i], self.etaP * (self.minusQ[i][U] - self.plusQ[i][U] + self.plusQ[i][N] * self.invZp + self.minusQ[i][N] * self.invZpNeig)),
                forge.assign(self.faultStressesTraction1[i], self.etaS * (self.minusQ[i][V] - self.plusQ[i][V] + self.plusQ[i][T1] * self.invZs + self.minusQ[i][T1] * self.invZsNeig)),
                forge.assign(self.faultStressesTraction2[i], self.etaS * (self.minusQ[i][W] - self.plusQ[i][W] + self.plusQ[i][T2] * self.invZs + self.minusQ[i][T2] * self.invZsNeig))
            ]

    def adjustInitialStress(self, routine, i):
        gNuc = smoothStepIncrement(self.fullUpdateTime, self.dt[i], self.t0)
        routine += [forge.conditional(self.fullUpdateTime <= self.t0, [
            forge.assign(self.initialStressInFaultCS[0], self.initialStressInFaultCS[0] + gNuc * self.nucleationStressInFaultCS[0]),
            forge.assign(self.initialPressure, self.initialPressure + gNuc * self.nucleationPressure)
        ])]

    def saveRuptureFrontOutput(self, routine):
        ruptureFrontThreshold = 0.001
        routine += [
            forge.conditional(self.ruptureTimePending and abs(self.slipRateMagnitude) >= ruptureFrontThreshold, [
                forge.assign(self.ruptureTime, self.fullUpdateTime),
                forge.assign(self.ruptureTimePending, False)
            ])
        ]

    def savePeakSlipRateOutput(self, routine):
        routine += [
            forge.assign(self.peakSlipRate, forge.max(self.peakSlipRate, self.slipRateMagnitude))
        ]

    def postcomputeImposedStateFromNewStress(self, routine):
        U,V,W,N,T1,T2 = 6,7,8,0,3,5
        # TODO: init to zero
        minusImposedState = [0] * self.aderdg.numberOfQuantities()
        plusImposedState = [0] * self.aderdg.numberOfQuantities()
        for i in range(self.aderdg.order):
            minusImposedState[N] = minusImposedState[N] + self.timeWeights[i] * self.faultStressesNormalStress[i]
            minusImposedState[T1] = minusImposedState[T1] + self.timeWeights[i] * self.tractionResultsTraction1[i]
            minusImposedState[T2] = minusImposedState[T2] + self.timeWeights[i] * self.tractionResultsTraction2[i]
            minusImposedState[U] = minusImposedState[U] + self.timeWeights[i] * (self.minusQ[i][U] - self.invZpNeig * (self.faultStressesNormalStress[i] - self.minusQ[i][N]))
            minusImposedState[V] = minusImposedState[V] + self.timeWeights[i] * (self.minusQ[i][V] - self.invZsNeig * (self.tractionResultsTraction1[i] - self.minusQ[i][T1]))
            minusImposedState[W] = minusImposedState[W] + self.timeWeights[i] * (self.minusQ[i][W] - self.invZsNeig * (self.tractionResultsTraction2[i] - self.minusQ[i][T2]))

            plusImposedState[N] = plusImposedState[N] + self.timeWeights[i] * self.faultStressesNormalStress[i]
            plusImposedState[T1] = plusImposedState[T1] + self.timeWeights[i] * self.tractionResultsTraction1[i]
            plusImposedState[T2] = plusImposedState[T2] + self.timeWeights[i] * self.tractionResultsTraction2[i]
            plusImposedState[U] = plusImposedState[U] + self.timeWeights[i] * (self.plusQ[i][U] + self.invZp * (self.faultStressesNormalStress[i] - self.plusQ[i][N]))
            plusImposedState[V] = plusImposedState[V] + self.timeWeights[i] * (self.plusQ[i][V] + self.invZs * (self.tractionResultsTraction1[i] - self.plusQ[i][T1]))
            plusImposedState[W] = plusImposedState[W] + self.timeWeights[i] * (self.plusQ[i][W] + self.invZs * (self.tractionResultsTraction2[i] - self.plusQ[i][T2]))
        routine += [
            forge.assign(self.minusImposedState[N], minusImposedState[N]),
            forge.assign(self.minusImposedState[T1], minusImposedState[T1]),
            forge.assign(self.minusImposedState[T2], minusImposedState[T2]),
            forge.assign(self.minusImposedState[U], minusImposedState[U]),
            forge.assign(self.minusImposedState[V], minusImposedState[V]),
            forge.assign(self.minusImposedState[W], minusImposedState[W]),
            forge.assign(self.plusImposedState[N], plusImposedState[N]),
            forge.assign(self.plusImposedState[T1], plusImposedState[T1]),
            forge.assign(self.plusImposedState[T2], plusImposedState[T2]),
            forge.assign(self.plusImposedState[U], plusImposedState[U]),
            forge.assign(self.plusImposedState[V], plusImposedState[V]),
            forge.assign(self.plusImposedState[W], plusImposedState[W]),
        ]

    def computeFrictionEnergy(self, routine):
        """
        U,V,W,N,T1,T2 = 6,7,8,0,3,5
        accumulatedSlip = self.accumulatedSlip
        frictionalEnergy = self.frictionalEnergy
        bPlus = self.etaS * self.invZs
        bMinus = self.etaS * self.invZsNeig
        energyDataSlip = 0
        for i in range(self.aderdg.order):
            interpolatedSlipRate1 = self.minusQ[i][U] - self.plusQ[i][U]
            interpolatedSlipRate2 = self.minusQ[i][V] - self.plusQ[i][V]
            interpolatedSlipRate3 = self.minusQ[i][W] - self.plusQ[i][W]

            interpolatedSlipRateMagnitude = forge.sqrt(interpolatedSlipRate1**2 + interpolatedSlipRate2**2 + interpolatedSlipRate3**2)
            accumulatedSlip = self.accumulatedSlip + self.timeWeights[i] * interpolatedSlipRateMagnitude
            energyDataSlip[0] = energyDataSlip[0] + self.timeWeights[i] * interpolatedSlipRate1
            energyDataSlip[1] = energyDataSlip[1] + self.timeWeights[i] * interpolatedSlipRate2
            energyDataSlip[2] = energyDataSlip[2] + self.timeWeights[i] * interpolatedSlipRate3

            interpolatedTraction12 = bPlus * self.minusQ[i][T1] + bMinus * self.plusQ[i][T1]
            interpolatedTraction13 = bPlus * self.minusQ[i][T2] + bMinus * self.plusQ[i][T2]

            weight = -self.timeWeights[i] * self.spaceWeights * doubledSurfaceArea
            frictionalEnergy = frictionalEnergy + weight * (interpolatedTraction12 * interpolatedSlipRate2 + interpolatedTraction13 * interpolatedSlipRate3)
        routine += [
            self.assign(self.accumulatedSlip, accumulatedSlip),
            self.assign(self.frictionalEnergy, frictionalEnergy)
        ]
        """

    def updateTimeSinceSlipRateBelowThreshold(self, routine):
        forge.assign(self.timeSinceSlipRateBelowThreshold, forge.ternary(
            self.ruptureTimePending,
            1e1000,
            forge.ternary(
                self.slipRateMagnitude < self.terminatorSlipRateThreshold,
                self.sumDt,
                0
            )
        ))

    def preHook(self, routine):
        raise NotImplementedError()
    
    def postHook(self, routine):
        raise NotImplementedError()

    def updateFrictionAndSlip(self, routine, i):
        raise NotImplementedError()
    
    def saveDynamicStressOutput(self, routine):
        raise NotImplementedError()

class LinearSlipWeakeningBase(FrictionLawBase):
    def __init__(self, aderdg, numberOfPoints):
        super().__init__(aderdg, numberOfPoints)
        self.muS = self.drvar('muS')
        self.muD = self.drvar('muD')
        self.dC = self.drvar('dC')
        self.cohesion = self.drvar('cohesion')
        self.forcedRuptureTime = self.drvar('forcedRuptureTime')

    def preHook(self, routine):
        pass
    
    def postHook(self, routine):
        pass

    def updateFrictionAndSlip(self, routine, i):
        self.calcStrengthHook(routine, i)
        self.calcSlipRateAndTraction(routine, i)
        self.calcStateVariableHook(routine, i)
        self.frictionFunctionHook(routine)
    
    def saveDynamicStressOutput(self, routine):
        routine += [
            forge.conditional(self.dynStressTimePending and abs(self.accumulatedSlipMagnitude) >= self.dC, [
                forge.assign(self.dynStressTime, self.fullUpdateTime),
                forge.assign(self.dynStressTimePending, False)
            ])
        ]

    def calcStrengthHook(self, routine, i):
        totalNormalStress = self.initialStressInFaultCS[0] + self.faultStressesNormalStress[i]
        strength = -self.cohesion - self.mu * forge.min(totalNormalStress, 0)
        # only FL=16 for now
        routine += [
            forge.assign(self.strength, strength)
        ]

    def calcSlipRateAndTraction(self, routine, i):
        totalStress1 = self.initialStressInFaultCS[3] + self.faultStressesTraction1[i]
        totalStress2 = self.initialStressInFaultCS[5] + self.faultStressesTraction2[i]
        absoluteShearStress = forge.sqrt(totalStress1**2 + totalStress2**2)

        slipRateMagnitude = forge.max(0, (absoluteShearStress - self.strength) * self.invEtaS)
        divisor = self.strength + self.etaS * self.slipRateMagnitude
        slipRate1 = slipRateMagnitude * totalStress1 / divisor
        slipRate2 = slipRateMagnitude * totalStress2 / divisor

        traction1 = self.faultStressesTraction1[i] - self.etaS * slipRate1
        traction2 = self.faultStressesTraction2[i] - self.etaS * slipRate2

        slip1 = self.slip1 + slipRate1 * self.dt[i]
        slip2 = self.slip2 + slipRate2 * self.dt[i]

        routine += [
            forge.assign(self.slipRateMagnitude, slipRateMagnitude),
            forge.assign(self.slipRate1, slipRate1),
            forge.assign(self.slipRate2, slipRate2),
            forge.assign(self.traction1, traction1),
            forge.assign(self.traction2, traction2),
            forge.assign(self.slip1, slip1),
            forge.assign(self.slip2, slip2)
        ]

    def calcStateVariableHook(self, routine, i):
        # resampledSlipRate = self.resampleMatrix @ self.slipRateMagnitude
        # self.accumulatedSlipMagnitude.pretensor <= self.accumulatedSlipMagnitude.pretensor + self.resampleMatrix['PS'] * self.slipRateMagnitude.pretensor.tensor['S']
        resampledSlipRate = self.slipRateMagnitude
        accumulatedSlipMagnitude = self.accumulatedSlipMagnitude + resampledSlipRate * self.dt[i]
        localStateVariable = forge.min(abs(accumulatedSlipMagnitude) / self.dC, 1)
        tn = self.fullUpdateTime + self.dt[i]
        clampfactor = (tn - self.forcedRuptureTime) / self.t0
        f2 = forge.ternary(self.t0 == 0, 
            forge.ternary(tn >= self.forcedRuptureTime, 1, 0),
            forge.min(forge.max(clampfactor, 0), 1)
        )
        stateVariable = forge.max(localStateVariable, f2)

        routine += [
            forge.assign(self.accumulatedSlipMagnitude, accumulatedSlipMagnitude),
            forge.assign(self.stateVariable, stateVariable)
        ]
    
    def frictionFunctionHook(self, routine):
        routine += [
            forge.assign(self.mu, self.muS - (self.muS - self.muD) * self.stateVariable)
        ]


class LinearSlipWeakening(LinearSlipWeakeningBase):
    def strengthHook(self, routine):
        pass

    def stateVariableHook(self, routine):
        pass

    def resampleSlipRate(self, routine):
        pass

class BiMaterialFault(LinearSlipWeakeningBase):
    def strengthHook(self, routine):
        pass

    def stateVariableHook(self, routine):
        pass

    def resampleSlipRate(self, routine):
        pass

class TPApprox(LinearSlipWeakeningBase):
    def strengthHook(self, routine):
        pass

    def stateVariableHook(self, routine):
        pass

    def resampleSlipRate(self, routine):
        pass

"""
class RateAndState(FrictionLawBase):
    def __init__(self, aderdg, numberOfPoints, tpmethod):
        super().__init__(aderdg, numberOfPoints)
        self.a = self.drvar('a')
        self.sl0 = self.drvar('sl0')
        self.newtonTolerance = self.drscalar('newtonTolerance')
        self.maxNumberSlipRateUpdates = self.drscalar('maxNumberSlipRateUpdates')
        self.stateVariableBuffer = forge.temp()
        self.tpmethod = tpmethod

        self.initialVariablesNormalStress = self.drvar('initialVariablesNormalStress')

    def preHook(self, routine):
        routine += [
            forge.assign(self.stateVariableBuffer, self.stateVariable)
        ]
    
    def postHook(self, routine):
        pass

    def calcInitialVariables(self, routine, i):
        routine += [
            
        ]
    
    def updateStateVariableIterative(self, routine, i):
        slipRateTest = forge.temp()
        iterator = forge.temp()
        routine += [
            forge.assign(slipRateTest, self.slipRateMagnitude),
            forge.assign(iterator, 0)
        ]

        almostZero = 1e-30 # TODO: adjust

        muF = self.updateMu(routine)
        dMuF = self.updateMuDerivative(routine)
        g = -self.invEtaS * (abs(self.normalStress) * muF - self.absoluteShearStress) - slipRateTest
        dG = -self.invEtaS * (abs(self.normalStress) * dMuF) - 1.0
        slipRateTestNew = forge.max(almostZero, slipRateTest - (g / dG))
        iteration = forge.loop(
            abs(g) >= self.newtonTolerance and iterator < self.maxNumberSlipRateUpdates,
            [
                forge.assign(slipRateTest, slipRateTestNew),
                forge.assign(iterator, iterator + 1)
            ]
        )

        routine += [iteration]

    def updateMu(self, routine):
        raise NotImplementedError()
    
    def updateMuDerivative(self, routine):
        raise NotImplementedError()

    def updateNormalStress(self, routine, i):
        routine += [
            forge.assign(self.initialVariablesNormalStress, forge.min(0, self.faultStressesNormalStress + self.initialStressInFaultCS[0] - self.tpmethod.calcFluidPressure()))
        ]

    def calcSlipRateAndTraction(self, routine, i):
        pass

    def updateFrictionAndSlip(self, routine, i):
        self.calcInitialVariables(routine, i)
        self.updateStateVariableIterative(routine, i)
        self.tpmethod.calcFluidPressure()
        self.updateNormalStress(routine, i)
        self.calcSlipRateAndTraction(routine, i)

class FastVelocityWeakeningLaw(RateAndState):
    def updateMu(self, routine):
        localA = forge.cast(self.a, FloatingPointType.DOUBLE)
        x = 0.5 / self.rsSr0 * forge.exp(localStateVariable / localA) * localSlipRateMagnitude
        return localA * forge.asinh(x)
    
    def updateMuDerivative(self, routine):
        localA = forge.cast(self.a, FloatingPointType.DOUBLE)
        c = 0.5 / self.rsSr0 * forge.exp(localStateVariable / localA)
        return localA * c / forge.sqrt((localSlipRateMagnitude * c)**2 + 1.0)

class SlowVelocityWeakeningLaw(RateAndState):
    def __init__(self, aderdg, numberOfPoints, tpmethod):
        super().__init__(aderdg, numberOfPoints, tpmethod)
        # self.a = self.drvar('a')
    
    def updateMu(self, routine):
        localA = forge.cast(self.a, FloatingPointType.DOUBLE)
        localSl0 = forge.cast(self.sl0, FloatingPointType.DOUBLE)
        log1 = forge.log(self.rsSr0 * localStateVariable / localSl0)
        x = 0.5 * (localSlipRateMagnitude / self.rsSr0) * forge.exp((details.rsF0 + details.rsB * log1) / localA)
        return localA * forge.asinh(x)
    
    def updateMuDerivative(self, routine):
        pass

class SlipLaw(SlowVelocityWeakeningLaw):
    def updateStateVariable(self, routine, i):
        sl0 = forge.cast(self.sl0, FloatingPointType.DOUBLE)
        localSlipRate = forge.cast(self.localSlipRate, FloatingPointType.DOUBLE)
        exp1 = forge.exp(-localSlipRate * (timeIncrement / sl0))

        stateVarReference = forge.cast(self.stateVariable, FloatingPointType.DOUBLE)
        self.stateVariableBuffer = forge.cast(sl0 / localSlipRate * forge.pow(localSlipRate * stateVarReference / sl0, exp1), self.fptype)

class AgingLaw(SlowVelocityWeakeningLaw):
    def updateStateVariable(self, routine, i):
        sl0 = forge.cast(self.sl0, FloatingPointType.DOUBLE)
        localSlipRate = forge.cast(self.localSlipRate, FloatingPointType.DOUBLE)
        exp1 = forge.exp(-localSlipRate * (timeIncrement / sl0))

        stateVarReference = forge.cast(self.stateVariable, FloatingPointType.DOUBLE)
        self.stateVariableBuffer = forge.cast(stateVarReference * exp1 + localSl0 / localSlipRate * (1.0 - exp1), self.fptype)

class TPMethod:
    def calcFluidPressure(self):
        raise NotImplementedError()

class NoTP(TPMethod):
    def calcFluidPressure(self):
        return 0

class TP(TPMethod):
    def calcFluidPressure(self, slipRateMagnitude, routine):
        temperatureUpdate = forge.TempVar()
        pressureUpdate = forge.TempVar()

        routine += [
            forge.assign(temperatureUpdate, 0),
            forge.assign(pressureUpdate, 0)
        ]

        tauV = self.faultStrength * slipRateMagnitude
        lambdaPrime = self.undrainedTPResponse * self.thermalDiffusivity / (self.hydraulicDiffusivity - self.thermalDiffusivity)

        for i in range(TODO):
            squaredNormalizedTPGrid = pow2(tpGridPoints[i] / self.halfWidthShearZone)
            expTheta = forge.exp(-self.thermalDiffusivity * self.dt * squaredNormalizedTPGrid)
            expSigma = forge.exp(-self.hydraulicDiffusivity * self.dt * squaredNormalizedTPGrid)

            thetaDiffusion = thetaBuffer[i] * expTheta
            sigmaDiffusion = sigmaBuffer[i] * expSigma

            omega = tauV * heatSource[i]
            thetaGeneration = omega / (self.heatCapacity * squaredNormalizedTPGrid * self.thermalDiffusivity) * (1 - expTheta)
            sigmaGeneration = omega * (self.undrainedTPResponse + lambdaPrime) / (self.heatCapacity * squaredNormalizedTPGrid * self.hydraulicDiffusivity) * (1 - expSigma)

            thetaBuffer[i] = thetaBuffer[i] + thetaGeneration
            sigmaBuffer[i] = sigmaBuffer[i] + sigmaDiffusion

            scaledInverseFourierCoefficient = tpInverseFourierCoefficients[i] / self.halfWidthShearZone
            temperatureUpdate += scaledInverseFourierCoefficient * thetaBuffer[i]
            pressureUpdate += scaledInverseFourierCoefficient * sigmaBuffer[i]
        
        pressureUpdate -= lambdaPrime * temperatureUpdate

        routine += [
            forge.assign(self.temperature, temperatureUpdate + self.initialTemperature),
            forge.assign(self.pressure, .pressureUpdate + self.initialPressure)
        ]
        return 0

class ImposedSlipLaw(FrictionLawBase):
    def preHook(self, routine):
        pass
    
    def postHook(self, routine):
        pass

    def updateFrictionAndSlip(self, routine, i):
        stfEvaluated = self.evaluateSTF(currentTime, timeIncrement)

        traction1 = self.faultStressesTraction1 - self.etaS * self.imposedSlipDirection1 * stfEvaluated
        traction2 = self.faultStressesTraction2 - self.etaS * self.imposedSlipDirection2 * stfEvaluated
        slipRate1 = self.imposedSlipDirection1 * stfEvaluated
        slipRate2 = self.imposedSlipDirection2 * stfEvaluated
        slipRateMagnitude = forge.sqrt(slipRate1**2, slipRate2**2)

        slip1 = self.slip1 + slipRate1 * timeIncrement
        slip2 = self.slip2 + slipRate2 * timeIncrement
        accumulatedSlipMagnitude = self.accumulatedSlipMagnitude + slipRateMagnitude * timeIncrement

        routine += [
            forge.assign(self.traction1, traction1),
            forge.assign(self.traction2, traction2),
            forge.assign(self.slipRate1, slipRate1),
            forge.assign(self.slipRate2, slipRate2),
            forge.assign(self.slip1, slip1),
            forge.assign(self.slip2, slip2),
            forge.assign(self.slipRateMagnitude, slipRateMagnitude),
            forge.assign(self.accumulatedSlipMagnitude, accumulatedSlipMagnitude)
            # TODO traction results
        ]
    
    def saveDynamicStressOutput(self, routine):
        pass

    def evaluateSTF(self, currentTime, timeIncrement):
        raise NotImplementedError()

class GaussianSTF(ImposedSlipLaw):
    def __init__(self, aderdg, numberOfPoints):
        super().__init__(aderdg, numberOfPoints)
        self.onsetTime = self.drvar('onsetTime')
        self.riseTime = self.drvar('riseTime')

    def evaluateSTF(self, currentTime, timeIncrement):
        return smoothStepIncrement(currentTime - self.onsetTime, timeIncrement, self.riseTime) / timeIncrement

class YoffeSTF(ImposedSlipLaw):
    def __init__(self, aderdg, numberOfPoints):
        super().__init__(aderdg, numberOfPoints)
        self.onsetTime = self.drvar('onsetTime')
        self.tauR = self.drvar('tauR')
        self.tauS = self.drvar('tauS')

    def evaluateSTF(self, currentTime, timeIncrement):
        time = currentTime - self.onsetTime
        tauR = self.tauR
        tauS = self.tauS

        k = 2.0 / (np.pi * tauR * tauS * tauS)

        # c1 to c6 are analytical functions used for building the regularized Yoffe function
        c1 = (0.5 * time + 0.25 * tauR) * forge.sqrt(time * (tauR - time)) + (time * tauR - tauR * tauR) * forge.asin(forge.sqrt(time / tauR)) - 0.75 * tauR * tauR * forge.atan(forge.sqrt((tauR - time) / time))

        c2 = 0.375 * np.pi * tauR * tauR

        c3 = (tauS - time - 0.5 * tauR) * forge.sqrt((time - tauS) * (tauR - time + tauS)) + tauR * (2 * tauR - 2 * time + 2 * tauS) * forge.asin(forge.sqrt((time - tauS) / tauR)) + 1.5 * tauR * tauR * forge.atan(forge.sqrt((tauR - time + tauS) / (time - tauS)))
        
        c4 = (-tauS + 0.5 * time + 0.25 * tauR) * forge.sqrt((time - 2.0 * tauS) * (tauR - time + 2.0 * tauS)) - tauR * (tauR - time + 2.0 * tauS) * forge.asin(forge.sqrt((time - 2.0 * tauS) / tauR)) - 0.75 * tauR * tauR * forge.atan(forge.sqrt((tauR - time + 2.0 * tauS) / (time - 2.0 * tauS)))
        
        c5 = 0.5 * np.pi * tauR * (time - tauR)

        c6 = 0.5 * np.pi * tauR * (2.0 * tauS - time + tauR)

        return forge.ternary(tauR > 2.0 * tauS,
            forge.ternary(time <= 0, 0,
            forge.ternary(time <= tauS, k * (c1 + c2),
            forge.ternary(time <= 2.0 * tauS, k * (c1 - c2 + c3),
            forge.ternary(time < tauR, k * (c1 + c3 + c4),
            forge.ternary(time < tauR + tauS, k * (c3 + c4 + c5),
            forge.ternary(time < tauR + 2.0 * tauS, k * (c4 + c6),
            0
            )
            )
            )
            )
            )
            ),
            forge.ternary(time <= 0, 0,
            forge.ternary(time <= tauS, k * (c1 + c2),
            forge.ternary(time <= 2.0 * tauS, k * (c5 - c2 + c3), # only change to above
            forge.ternary(time < tauR, k * (c1 + c3 + c4),
            forge.ternary(time < tauR + tauS, k * (c3 + c4 + c5),
            forge.ternary(time < tauR + 2.0 * tauS, k * (c4 + c6),
            0
            )
            )
            )
            )
            )
            )
        )
"""
