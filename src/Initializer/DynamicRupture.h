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

#include "Initializer/tree/VariableContainer.hpp"

namespace seissol {
  namespace initializers {

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

template<typename Config>
struct DynamicRupture : LTSVariableContainer {
public:
  using RealT = typename Config::RealT;

  virtual ~DynamicRupture() = default;
  Variable<RealT*>                                                   timeDerivativePlus;
  Variable<RealT*>                                                   timeDerivativeMinus;
  Variable<RealT[tensor::QInterpolated::size()]>                     imposedStatePlus;
  Variable<RealT[tensor::QInterpolated::size()]>                     imposedStateMinus;
  Variable<DRGodunovData>                                           godunovData;
  Variable<RealT[tensor::fluxSolver::size()]>                        fluxSolverPlus;
  Variable<RealT[tensor::fluxSolver::size()]>                        fluxSolverMinus;
  Variable<DRFaceInformation>                                       faceInformation;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsPlus;
  Variable<model::IsotropicWaveSpeeds>                              waveSpeedsMinus;
  Variable<DREnergyOutput>                                          drEnergyOutput;

  Variable<seissol::dr::ImpedancesAndEta>                           impAndEta;
  //size padded for vectorization
  //CS = coordinate system
  Variable<RealT[dr::misc::numPaddedPoints][6]> initialStressInFaultCS;
  Variable<RealT[dr::misc::numPaddedPoints][6]> nucleationStressInFaultCS;
  Variable<RealT[dr::misc::numPaddedPoints]> mu;
  Variable<RealT[dr::misc::numPaddedPoints]> accumulatedSlipMagnitude;
  Variable<RealT[dr::misc::numPaddedPoints]> slip1; // slip at given fault node along local direction 1
  Variable<RealT[dr::misc::numPaddedPoints]> slip2; // slip at given fault node along local direction 2
  Variable<RealT[dr::misc::numPaddedPoints]> slipRateMagnitude;
  Variable<RealT[dr::misc::numPaddedPoints]> slipRate1; // slip rate at given fault node along local direction 1
  Variable<RealT[dr::misc::numPaddedPoints]> slipRate2; // slip rate at given fault node along local direction 2
  Variable<RealT[dr::misc::numPaddedPoints]> ruptureTime;
  Variable<RealT[dr::misc::numPaddedPoints]> dynStressTime;
  Variable<bool[dr::misc::numPaddedPoints]> ruptureTimePending;
  Variable<bool[dr::misc::numPaddedPoints]> dynStressTimePending;
  Variable<RealT[dr::misc::numPaddedPoints]> peakSlipRate;
  Variable<RealT[dr::misc::numPaddedPoints]> traction1;
  Variable<RealT[dr::misc::numPaddedPoints]> traction2;
  Variable<RealT[ConvergenceOrder][tensor::QInterpolated::size()]> qInterpolatedPlus;
  Variable<RealT[ConvergenceOrder][tensor::QInterpolated::size()]> qInterpolatedMinus;

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
    tree.addVar(          drEnergyOutput,             mask,         Alignment,      MEMKIND_STANDARD );
    tree.addVar(      impAndEta,                      mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      initialStressInFaultCS,         mask,                 1,      MEMKIND_STANDARD );
    tree.addVar(      nucleationStressInFaultCS,      mask,                 1,      MEMKIND_STANDARD );
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
    tree.addVar(qInterpolatedPlus, mask, Alignment, MEMKIND_STANDARD);
    tree.addVar(qInterpolatedMinus, mask, Alignment, MEMKIND_STANDARD);

#ifdef ACL_DEVICE
    tree.addScratchpadMemory(idofsPlusOnDevice,  1, seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(idofsMinusOnDevice, 1,  seissol::memory::DeviceGlobalMemory);
#endif
  }
};

template<typename Config>
struct LTSLinearSlipWeakening : public DynamicRupture<Config> {
    using RealT = typename Config::RealT;
    Variable<RealT[dr::misc::numPaddedPoints]> dC;
    Variable<RealT[dr::misc::numPaddedPoints]> muS;
    Variable<RealT[dr::misc::numPaddedPoints]> muD;
    Variable<RealT[dr::misc::numPaddedPoints]> cohesion;
    Variable<RealT[dr::misc::numPaddedPoints]> forcedRuptureTime;


    virtual void addTo(initializers::LTSTree& tree) {
        DynamicRupture<Config>::addTo(tree);
        LayerMask mask = LayerMask(Ghost);
        tree.addVar(dC, mask, 1, MEMKIND_STANDARD);
        tree.addVar(muS, mask, 1, MEMKIND_STANDARD);
        tree.addVar(muD, mask, 1, MEMKIND_STANDARD);
        tree.addVar(cohesion, mask,1, MEMKIND_STANDARD);
        tree.addVar(forcedRuptureTime, mask, 1, MEMKIND_STANDARD);
    }
};

template<typename Config>
struct LTSLinearSlipWeakeningBimaterial : public LTSLinearSlipWeakening<Config> {
  using RealT = typename Config::RealT;
  Variable<RealT[dr::misc::numPaddedPoints]> regularisedStrength;

  virtual void addTo(initializers::LTSTree& tree) {
    LTSLinearSlipWeakening<Config>::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(regularisedStrength, mask, 1, MEMKIND_STANDARD);
  }
};

template<typename Config>
struct LTSRateAndState : public DynamicRupture<Config> {
  using RealT = typename Config::RealT;
  Variable<RealT[dr::misc::numPaddedPoints]> rsA;
  Variable<RealT[dr::misc::numPaddedPoints]> rsSl0;
  Variable<RealT[dr::misc::numPaddedPoints]> stateVariable;

  virtual void addTo(initializers::LTSTree& tree) {
    DynamicRupture<Config>::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsA, mask, 1, MEMKIND_STANDARD);
    tree.addVar(rsSl0, mask, 1, MEMKIND_STANDARD);
    tree.addVar(stateVariable, mask, 1, MEMKIND_STANDARD);
  }
};

template<typename Config>
struct LTSRateAndStateFastVelocityWeakening : public LTSRateAndState<Config> {
  using RealT = typename Config::RealT;
  Variable<RealT[dr::misc::numPaddedPoints]> rsSrW;

  virtual void addTo(initializers::LTSTree& tree) {
    LTSRateAndState<Config>::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsSrW, mask, 1, MEMKIND_STANDARD);
  }
};

template<typename Config>
struct LTSRateAndStateThermalPressurization : public LTSRateAndStateFastVelocityWeakening<Config> {
using RealT = typename Config::RealT;
  Variable<RealT[dr::misc::numPaddedPoints]> temperature;
  Variable<RealT[dr::misc::numPaddedPoints]> pressure;
  Variable<RealT[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> theta;
  Variable<RealT[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigma;
  Variable<RealT[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> thetaTmpBuffer;
  Variable<RealT[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigmaTmpBuffer;
  Variable<RealT[dr::misc::numPaddedPoints]> faultStrength;
  Variable<RealT[dr::misc::numPaddedPoints]>halfWidthShearZone;
  Variable<RealT[dr::misc::numPaddedPoints]> hydraulicDiffusivity;

  virtual void addTo(initializers::LTSTree& tree) {
    LTSRateAndStateFastVelocityWeakening<Config>::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(temperature, mask, Alignment, seissol::memory::Standard);
    tree.addVar(pressure, mask, Alignment, seissol::memory::Standard);
    tree.addVar(theta, mask, Alignment, seissol::memory::Standard);
    tree.addVar(sigma, mask, Alignment, seissol::memory::Standard);
    tree.addVar(thetaTmpBuffer, mask, Alignment, seissol::memory::Standard);
    tree.addVar(sigmaTmpBuffer, mask, Alignment, seissol::memory::Standard);
    tree.addVar(faultStrength, mask, Alignment, seissol::memory::Standard);
    tree.addVar(halfWidthShearZone, mask, Alignment, seissol::memory::Standard);
    tree.addVar(hydraulicDiffusivity, mask, Alignment, seissol::memory::Standard);
  }
};

template<typename Config>
struct LTSImposedSlipRates : public DynamicRupture<Config> {
  using RealT = typename Config::RealT;
  Variable<RealT[dr::misc::numPaddedPoints]> imposedSlipDirection1;
  Variable<RealT[dr::misc::numPaddedPoints]> imposedSlipDirection2;
  Variable<RealT[dr::misc::numPaddedPoints]> onsetTime;

  virtual void addTo(initializers::LTSTree& tree) {
    DynamicRupture<Config>::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(imposedSlipDirection1, mask, 1, seissol::memory::Standard);
    tree.addVar(imposedSlipDirection2, mask, 1, seissol::memory::Standard);
    tree.addVar(onsetTime, mask, 1, seissol::memory::Standard);
  }
};

template<typename Config>
struct LTSImposedSlipRatesYoffe : public LTSImposedSlipRates<Config> {
  using RealT = typename Config::RealT;
  Variable<RealT[dr::misc::numPaddedPoints]> tauS;
  Variable<RealT[dr::misc::numPaddedPoints]> tauR;

  virtual void addTo(initializers::LTSTree& tree) {
    LTSImposedSlipRates<Config>::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(tauS, mask, 1, seissol::memory::Standard);
    tree.addVar(tauR, mask, 1, seissol::memory::Standard);
  }
};

template<typename Config>
struct LTSImposedSlipRatesGaussian : public LTSImposedSlipRates<Config> {
  using RealT = typename Config::RealT;
  Variable<RealT[dr::misc::numPaddedPoints]> riseTime;

  virtual void addTo(initializers::LTSTree& tree) {
    LTSImposedSlipRates<Config>::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(riseTime, mask, 1, seissol::memory::Standard);
  }
};

  }}

#endif // INITIALIZER_DR_H_
