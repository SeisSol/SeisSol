/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2016, SeisSol Group
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 * 
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 * 
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 **/
 
#ifndef INITIALIZER_DR_H_
#define INITIALIZER_DR_H_

#include <Initializer/typedefs.hpp>
#include <Initializer/tree/LTSTree.hpp>
#include <generated_code/tensor.h>
#include <DynamicRupture/Misc.h>

namespace seissol {
  namespace initializer {
    struct DynamicRupture;
    struct LTSLinearSlipWeakening;
    struct LTSLinearSlipWeakeningForcedRuptureTime;
    struct LTSLinearSlipWeakeningBimaterial;
    struct LTSRateAndState;
    struct LTSRateAndStateFastVelocityWeakening;
    struct LTSRateAndStateThermalPressurization;
    struct LTSImposedSlipRates;
    struct LTSImposedSlipRatesYoffe;
    struct LTSImposedSlipRatesGaussian;
  } // namespace initializer
} // namespace seissol

#ifndef ACL_DEVICE
#define MEMKIND_NEIGHBOUR_INTEGRATION seissol::memory::Standard
#define MEMKIND_Q_INTERPOLATED seissol::memory::Standard
#define MEMKIND_IMPOSED_STATE seissol::memory::Standard
#define MEMKIND_STANDARD seissol::memory::Standard
#else
#define MEMKIND_NEIGHBOUR_INTEGRATION seissol::memory::DeviceUnifiedMemory
#define MEMKIND_IMPOSED_STATE seissol::memory::DeviceUnifiedMemory
#define MEMKIND_STANDARD seissol::memory::DeviceUnifiedMemory
#endif

struct seissol::initializer::DynamicRupture {
public:
  virtual ~DynamicRupture() = default;
  Variable<real*>                                                   timeDerivativePlus;
  Variable<real*>                                                   timeDerivativeMinus;
  Variable<real[tensor::QInterpolated::size()]>                     imposedStatePlus;
  Variable<real[tensor::QInterpolated::size()]>                     imposedStateMinus;
  Variable<DRGodunovData>                                           godunovData;
  Variable<real[tensor::fluxSolver::size()]>                        fluxSolverPlus;
  Variable<real[tensor::fluxSolver::size()]>                        fluxSolverMinus;
  Variable<DRFaceInformation>                                       faceInformation;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsPlus;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsMinus;
  Variable<DREnergyOutput>                                          drEnergyOutput;

  Variable<seissol::dr::ImpedancesAndEta>                           impAndEta;
  Variable<seissol::dr::ImpedanceMatrices>                          impedanceMatrices;
  //size padded for vectorization
  //CS = coordinate system
  Variable<real[dr::misc::numPaddedPoints][6]> initialStressInFaultCS;
  Variable<real[dr::misc::numPaddedPoints][6]> nucleationStressInFaultCS;
  // will be always zero, if not using poroelasticity
  Variable<real[dr::misc::numPaddedPoints]> initialPressure;
  Variable<real[dr::misc::numPaddedPoints]> nucleationPressure;
  Variable<real[dr::misc::numPaddedPoints]> mu;
  Variable<real[dr::misc::numPaddedPoints]> accumulatedSlipMagnitude;
  Variable<real[dr::misc::numPaddedPoints]> slip1; // slip at given fault node along local direction 1
  Variable<real[dr::misc::numPaddedPoints]> slip2; // slip at given fault node along local direction 2
  Variable<real[dr::misc::numPaddedPoints]> slipRateMagnitude;
  Variable<real[dr::misc::numPaddedPoints]> slipRate1; // slip rate at given fault node along local direction 1
  Variable<real[dr::misc::numPaddedPoints]> slipRate2; // slip rate at given fault node along local direction 2
  Variable<real[dr::misc::numPaddedPoints]> ruptureTime;
  Variable<real[dr::misc::numPaddedPoints]> dynStressTime;
  Variable<bool[dr::misc::numPaddedPoints]> ruptureTimePending;
  Variable<bool[dr::misc::numPaddedPoints]> dynStressTimePending;
  Variable<real[dr::misc::numPaddedPoints]> peakSlipRate;
  Variable<real[dr::misc::numPaddedPoints]> traction1;
  Variable<real[dr::misc::numPaddedPoints]> traction2;
  Variable<real[CONVERGENCE_ORDER][tensor::QInterpolated::size()]> qInterpolatedPlus;
  Variable<real[CONVERGENCE_ORDER][tensor::QInterpolated::size()]> qInterpolatedMinus;

#ifdef ACL_DEVICE
  ScratchpadMemory                        idofsPlusOnDevice;
  ScratchpadMemory                        idofsMinusOnDevice;
#endif
  
  virtual void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      timeDerivativePlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(     timeDerivativeMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(        imposedStatePlus,             mask,     PAGESIZE_HEAP,      MEMKIND_IMPOSED_STATE );
    tree.addVar(       imposedStateMinus,             mask,     PAGESIZE_HEAP,      MEMKIND_IMPOSED_STATE );
    tree.addVar(             godunovData,             mask,                 1,      MEMKIND_NEIGHBOUR_INTEGRATION );
    tree.addVar(          fluxSolverPlus,             mask,                 1,      MEMKIND_NEIGHBOUR_INTEGRATION );
    tree.addVar(         fluxSolverMinus,             mask,                 1,      MEMKIND_NEIGHBOUR_INTEGRATION );
    tree.addVar(         faceInformation,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(          waveSpeedsPlus,             mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(         waveSpeedsMinus,             mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(          drEnergyOutput,             mask,         ALIGNMENT,      MEMKIND_STANDARD );
    tree.addVar(      impAndEta,                      mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      impedanceMatrices,              mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      initialStressInFaultCS,         mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      nucleationStressInFaultCS,      mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      initialPressure,                mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      nucleationPressure,             mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      ruptureTime,                    mask,                 1,      MEMKIND_STANDARD );

    tree.addVar(ruptureTimePending, mask, 1, MEMKIND_STANDARD);
    tree.addVar(dynStressTime, mask, 1, MEMKIND_STANDARD);
    tree.addVar(dynStressTimePending, mask, 1, MEMKIND_STANDARD);
    tree.addVar(mu, mask, 1, MEMKIND_STANDARD);
    tree.addVar(accumulatedSlipMagnitude, mask, 1, MEMKIND_STANDARD);
    tree.addVar(slip1, mask, 1, MEMKIND_STANDARD);
    tree.addVar(slip2, mask, 1, MEMKIND_STANDARD);
    tree.addVar(slipRateMagnitude, mask, 1, MEMKIND_STANDARD);
    tree.addVar(slipRate1, mask, 1, MEMKIND_STANDARD);
    tree.addVar(slipRate2, mask, 1, MEMKIND_STANDARD);
    tree.addVar(peakSlipRate, mask, 1, MEMKIND_STANDARD);
    tree.addVar(traction1, mask, 1, MEMKIND_STANDARD);
    tree.addVar(traction2, mask, 1, MEMKIND_STANDARD);
    tree.addVar(qInterpolatedPlus, mask, ALIGNMENT, MEMKIND_STANDARD);
    tree.addVar(qInterpolatedMinus, mask, ALIGNMENT, MEMKIND_STANDARD);

#ifdef ACL_DEVICE
    tree.addScratchpadMemory(idofsPlusOnDevice,  1, seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(idofsMinusOnDevice, 1,  seissol::memory::DeviceGlobalMemory);
#endif
  }
};

struct seissol::initializer::LTSLinearSlipWeakening : public seissol::initializer::DynamicRupture {
    Variable<real[dr::misc::numPaddedPoints]> dC;
    Variable<real[dr::misc::numPaddedPoints]> muS;
    Variable<real[dr::misc::numPaddedPoints]> muD;
    Variable<real[dr::misc::numPaddedPoints]> cohesion;
    Variable<real[dr::misc::numPaddedPoints]> forcedRuptureTime;


    virtual void addTo(initializer::LTSTree& tree) {
        seissol::initializer::DynamicRupture::addTo(tree);
        LayerMask mask = LayerMask(Ghost);
        tree.addVar(dC, mask, 1, MEMKIND_STANDARD);
        tree.addVar(muS, mask, 1, MEMKIND_STANDARD);
        tree.addVar(muD, mask, 1, MEMKIND_STANDARD);
        tree.addVar(cohesion, mask,1, MEMKIND_STANDARD);
        tree.addVar(forcedRuptureTime, mask, 1, MEMKIND_STANDARD);
    }
};

struct seissol::initializer::LTSLinearSlipWeakeningBimaterial : public seissol::initializer::LTSLinearSlipWeakening {
  Variable<real[dr::misc::numPaddedPoints]> regularisedStrength;

  virtual void addTo(initializer::LTSTree& tree) {
    seissol::initializer::LTSLinearSlipWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(regularisedStrength, mask, 1, MEMKIND_STANDARD);
  }
};

struct seissol::initializer::LTSRateAndState : public seissol::initializer::DynamicRupture {
  Variable<real[dr::misc::numPaddedPoints]> rsA;
  Variable<real[dr::misc::numPaddedPoints]> rsSl0;
  Variable<real[dr::misc::numPaddedPoints]> stateVariable;

  virtual void addTo(initializer::LTSTree& tree) {
    seissol::initializer::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsA, mask, 1, MEMKIND_STANDARD);
    tree.addVar(rsSl0, mask, 1, MEMKIND_STANDARD);
    tree.addVar(stateVariable, mask, 1, MEMKIND_STANDARD);
  }
};


struct seissol::initializer::LTSRateAndStateFastVelocityWeakening : public seissol::initializer::LTSRateAndState {
  Variable<real[dr::misc::numPaddedPoints]> rsSrW;

  virtual void addTo(initializer::LTSTree& tree) {
    seissol::initializer::LTSRateAndState::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsSrW, mask, 1, MEMKIND_STANDARD);
  }
};

struct seissol::initializer::LTSRateAndStateThermalPressurization : public seissol::initializer::LTSRateAndStateFastVelocityWeakening {

  Variable<real[dr::misc::numPaddedPoints]> temperature;
  Variable<real[dr::misc::numPaddedPoints]> pressure;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> theta;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigma;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> thetaTmpBuffer;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigmaTmpBuffer;
  Variable<real[dr::misc::numPaddedPoints]> faultStrength;
  Variable<real[dr::misc::numPaddedPoints]>halfWidthShearZone;
  Variable<real[dr::misc::numPaddedPoints]> hydraulicDiffusivity;

  virtual void addTo(initializer::LTSTree& tree) {
    seissol::initializer::LTSRateAndStateFastVelocityWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(temperature, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(pressure, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(theta, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(sigma, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(thetaTmpBuffer, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(sigmaTmpBuffer, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(faultStrength, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(halfWidthShearZone, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(hydraulicDiffusivity, mask, ALIGNMENT, seissol::memory::Standard);
  }
};


struct seissol::initializer::LTSImposedSlipRates : public seissol::initializer::DynamicRupture {
  Variable<real[dr::misc::numPaddedPoints]> imposedSlipDirection1;
  Variable<real[dr::misc::numPaddedPoints]> imposedSlipDirection2;
  Variable<real[dr::misc::numPaddedPoints]> onsetTime;

  virtual void addTo(initializer::LTSTree& tree) {
    seissol::initializer::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(imposedSlipDirection1, mask, 1, MEMKIND_STANDARD);
    tree.addVar(imposedSlipDirection2, mask, 1, MEMKIND_STANDARD);
    tree.addVar(slip2, mask, 1, MEMKIND_STANDARD);
    tree.addVar(onsetTime, mask, 1, MEMKIND_STANDARD);
  }
};


struct seissol::initializer::LTSImposedSlipRatesYoffe : public seissol::initializer::LTSImposedSlipRates {
  Variable<real[dr::misc::numPaddedPoints]> tauS;
  Variable<real[dr::misc::numPaddedPoints]> tauR;

  virtual void addTo(initializer::LTSTree& tree) {
    seissol::initializer::LTSImposedSlipRates::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(tauS, mask, 1, MEMKIND_STANDARD);
    tree.addVar(tauR, mask, 1, MEMKIND_STANDARD);
  }
};


struct seissol::initializer::LTSImposedSlipRatesGaussian : public seissol::initializer::LTSImposedSlipRates {
  Variable<real[dr::misc::numPaddedPoints]> riseTime;

  virtual void addTo(initializer::LTSTree& tree) {
    seissol::initializer::LTSImposedSlipRates::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(riseTime, mask, 1, MEMKIND_STANDARD);
  }
};

#endif // INITIALIZER_DR_H_
