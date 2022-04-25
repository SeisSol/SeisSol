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
  namespace initializers {
    struct DynamicRupture;
    struct LTS_LinearSlipWeakening;
    struct LTS_LinearSlipWeakeningForcedRuptureTime;
    struct LTS_LinearSlipWeakeningBimaterial;
    struct LTS_RateAndState;
    struct LTS_RateAndStateFastVelocityWeakening;
    struct LTS_RateAndStateThermalPressurization;
    struct LTS_ImposedSlipRates;
    struct LTS_ImposedSlipRatesYoffe;
    struct LTS_ImposedSlipRatesGaussian;
  } // namespace initializers
} // namespace seissol

#ifndef ACL_DEVICE
#	define MEMKIND_NEIGHBOUR_INTEGRATION seissol::memory::Standard
#	define MEMKIND_Q_INTERPOLATED seissol::memory::Standard
#	define MEMKIND_IMPOSED_STATE seissol::memory::Standard
#else
#	define MEMKIND_NEIGHBOUR_INTEGRATION seissol::memory::DeviceUnifiedMemory
#	define MEMKIND_Q_INTERPOLATED seissol::memory::PinnedMemory
#	define MEMKIND_IMPOSED_STATE seissol::memory::DeviceGlobalMemory
#endif


struct seissol::initializers::DynamicRupture {
public:
  virtual ~DynamicRupture() {}
  bool isFaultParameterizedByTraction;        //true if Traction T_n , T_s, T_d is stored in iniBulk_XX, iniShearXY, iniShearXZ
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

  Variable<seissol::dr::ImpedancesAndEta>                           impAndEta;
  //size padded for vectorization
  //CS = coordinate system
  Variable<real[dr::misc::numPaddedPoints][6]> initialStressInFaultCS;
  Variable<real[dr::misc::numPaddedPoints][6]> nucleationStressInFaultCS;
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
  Variable<real> averagedSlip;
  Variable<real[CONVERGENCE_ORDER][tensor::QInterpolated::size()]> qInterpolatedPlus;
  Variable<real[CONVERGENCE_ORDER][tensor::QInterpolated::size()]> qInterpolatedMinus;

#ifdef ACL_DEVICE
  ScratchpadMemory                        idofsPlusOnDevice;
  ScratchpadMemory                        idofsMinusOnDevice;
  ScratchpadMemory                        QInterpolatedPlusOnDevice;
  ScratchpadMemory                        QInterpolatedMinusOnDevice;

  ScratchpadMemory                        QInterpolatedPlusOnHost;
  ScratchpadMemory                        QInterpolatedMinusOnHost;
  ScratchpadMemory                        imposedStatePlusOnHost;
  ScratchpadMemory                        imposedStateMinusOnHost;
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
    tree.addVar(          waveSpeedsPlus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(         waveSpeedsMinus,             mask,                 1,      seissol::memory::Standard );
    tree.addVar(      impAndEta,                      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      initialStressInFaultCS,         mask,                 1,      seissol::memory::Standard );
    tree.addVar(      nucleationStressInFaultCS,      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      ruptureTime,                    mask,                 1,      seissol::memory::Standard );

    tree.addVar(ruptureTimePending, mask, 1, seissol::memory::Standard);
    tree.addVar(dynStressTime, mask, 1, seissol::memory::Standard);
    tree.addVar(dynStressTimePending, mask, 1, seissol::memory::Standard);
    tree.addVar(mu, mask, 1, seissol::memory::Standard);
    tree.addVar(accumulatedSlipMagnitude, mask, 1, seissol::memory::Standard);
    tree.addVar(slip1, mask, 1, seissol::memory::Standard);
    tree.addVar(slip2, mask, 1, seissol::memory::Standard);
    tree.addVar(slipRateMagnitude, mask, 1, seissol::memory::Standard);
    tree.addVar(slipRate1, mask, 1, seissol::memory::Standard);
    tree.addVar(slipRate2, mask, 1, seissol::memory::Standard);
    tree.addVar(peakSlipRate, mask, 1, seissol::memory::Standard);
    tree.addVar(traction1, mask, 1, seissol::memory::Standard);
    tree.addVar(traction2, mask, 1, seissol::memory::Standard);
    tree.addVar(averagedSlip, mask, 1, seissol::memory::Standard);
    tree.addVar(qInterpolatedPlus, mask, ALIGNMENT, seissol::memory::Standard);
    tree.addVar(qInterpolatedMinus, mask, ALIGNMENT, seissol::memory::Standard);

#ifdef ACL_DEVICE
    tree.addScratchpadMemory(  idofsPlusOnDevice,              1,      seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(  idofsMinusOnDevice,             1,      seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(  QInterpolatedPlusOnDevice,      1,      seissol::memory::DeviceGlobalMemory);
    tree.addScratchpadMemory(  QInterpolatedMinusOnDevice,     1,      seissol::memory::DeviceGlobalMemory);

    tree.addScratchpadMemory(  QInterpolatedPlusOnHost,         1,      seissol::memory::PinnedMemory);
    tree.addScratchpadMemory(  QInterpolatedMinusOnHost,        1,      seissol::memory::PinnedMemory);
    tree.addScratchpadMemory(  imposedStatePlusOnHost,          1,      seissol::memory::PinnedMemory);
    tree.addScratchpadMemory(  imposedStateMinusOnHost,         1,      seissol::memory::PinnedMemory);
#endif
  }
};

struct seissol::initializers::LTS_LinearSlipWeakening : public seissol::initializers::DynamicRupture {
    Variable<real[dr::misc::numPaddedPoints]> dC;
    Variable<real[dr::misc::numPaddedPoints]> muS;
    Variable<real[dr::misc::numPaddedPoints]> muD;
    Variable<real[dr::misc::numPaddedPoints]> cohesion;


    virtual void addTo(initializers::LTSTree& tree) {
        seissol::initializers::DynamicRupture::addTo(tree);
        LayerMask mask = LayerMask(Ghost);
        tree.addVar(dC, mask, 1, seissol::memory::Standard);
        tree.addVar(muS, mask, 1, seissol::memory::Standard);
        tree.addVar(muD, mask, 1, seissol::memory::Standard);
        tree.addVar(cohesion, mask,1, seissol::memory::Standard);
    }
};

struct seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime : public seissol::initializers::LTS_LinearSlipWeakening {
  Variable<real[dr::misc::numPaddedPoints]> forcedRuptureTime;
  Variable<real> tn;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_LinearSlipWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(forcedRuptureTime, mask, 1, seissol::memory::Standard);
    tree.addVar(tn, mask, 1, seissol::memory::Standard);
  }
};

struct seissol::initializers::LTS_LinearSlipWeakeningBimaterial : public seissol::initializers::LTS_LinearSlipWeakening {
  Variable<real[dr::misc::numPaddedPoints]> regularisedStrength;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_LinearSlipWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(regularisedStrength, mask, 1, seissol::memory::Standard);
  }
};

struct seissol::initializers::LTS_RateAndState : public seissol::initializers::DynamicRupture {
  Variable<real[dr::misc::numPaddedPoints]> rsA;
  Variable<real[dr::misc::numPaddedPoints]> rsSl0;
  Variable<real[dr::misc::numPaddedPoints]> stateVariable;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsA, mask, 1, seissol::memory::Standard);
    tree.addVar(rsSl0, mask, 1, seissol::memory::Standard);
    tree.addVar(stateVariable, mask, 1, seissol::memory::Standard);
  }
};


struct seissol::initializers::LTS_RateAndStateFastVelocityWeakening : public seissol::initializers::LTS_RateAndState {
  Variable<real[dr::misc::numPaddedPoints]> rsSrW;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_RateAndState::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsSrW, mask, 1, seissol::memory::Standard);
  }
};

struct seissol::initializers::LTS_RateAndStateThermalPressurization : public seissol::initializers::LTS_RateAndStateFastVelocityWeakening {

  Variable<real[dr::misc::numPaddedPoints]> temperature;  //this is TP[1] in fortran
  Variable<real[dr::misc::numPaddedPoints]> pressure;     //this is TP[2] in fortran
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> theta;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigma;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> thetaTmpBuffer;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigmaTmpBuffer;
  Variable<real[dr::misc::numPaddedPoints]> faultStrength;
  Variable<real[dr::misc::numPaddedPoints]>halfWidthShearZone;
  Variable<real[dr::misc::numPaddedPoints]> hydraulicDiffusivity;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_RateAndStateFastVelocityWeakening::addTo(tree);
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


struct seissol::initializers::LTS_ImposedSlipRates : public seissol::initializers::DynamicRupture {
  Variable<real[dr::misc::numPaddedPoints]> strikeSlip;
  Variable<real[dr::misc::numPaddedPoints]> dipSlip;
  Variable<real[dr::misc::numPaddedPoints]> onsetTime;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(strikeSlip, mask, 1, seissol::memory::Standard);
    tree.addVar(dipSlip, mask, 1, seissol::memory::Standard);
    tree.addVar(onsetTime, mask, 1, seissol::memory::Standard);
  }
};


struct seissol::initializers::LTS_ImposedSlipRatesYoffe : public seissol::initializers::LTS_ImposedSlipRates {
  Variable<real[dr::misc::numPaddedPoints]> tauS;
  Variable<real[dr::misc::numPaddedPoints]> tauR;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(tauS, mask, 1, seissol::memory::Standard);
    tree.addVar(tauR, mask, 1, seissol::memory::Standard);
  }
};


struct seissol::initializers::LTS_ImposedSlipRatesGaussian : public seissol::initializers::LTS_ImposedSlipRates {
  Variable<real[dr::misc::numPaddedPoints]> riseTime;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(riseTime, mask, 1, seissol::memory::Standard);
  }
};


#endif
