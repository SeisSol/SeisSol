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

namespace seissol {
  namespace initializers {
    struct DynamicRupture;
    struct LTS_LinearSlipWeakening;
    struct LTS_LinearSlipWeakeningForcedRuptureTime;
    struct LTS_LinearSlipWeakeningBimaterial;
    struct LTS_RateAndState;
    struct LTS_RateAndStateFastVelocityWeakening;
    struct LTS_RateAndStateThermalPressurisation;
    struct LTS_ImposedSlipRates;
  }
}

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
  //assert(init::QInterpolated::Start[0] == 0); ?
protected:
  static constexpr int numOfPointsPadded = init::QInterpolated::Stop[0];
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
  Variable<real[numOfPointsPadded][6]> initialStressInFaultCS;
  Variable<real[numOfPointsPadded][6]> nucleationStressInFaultCS;
  Variable<real[ numOfPointsPadded ]> mu;
  Variable<real[ numOfPointsPadded ]> slip;
  Variable<real[ numOfPointsPadded ]> slipStrike; // = Slip1
  Variable<real[ numOfPointsPadded ]> slipDip;    // = Slip2
  Variable<real[ numOfPointsPadded ]> slipRateMagnitude;
  Variable<real[ numOfPointsPadded ]> slipRateStrike;  // slip rate in Y-dirction (strike) Fortran: slipRate1
  Variable<real[ numOfPointsPadded ]> slipRateDip; // slip rate in Z-direction (dip) Fortran: slipRate2
  Variable<real[ numOfPointsPadded ]> ruptureTime;
  Variable<real[ numOfPointsPadded ]> dynStressTime;
  Variable<bool[ numOfPointsPadded ]> ruptureFront;
  Variable<real[ numOfPointsPadded ]> peakSlipRate;
  Variable<real[ numOfPointsPadded ]> tractionXY;
  Variable<real[ numOfPointsPadded ]> tractionXZ;
  Variable<bool[ numOfPointsPadded ]> ds;
  Variable<real> averagedSlip;

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
    tree.addVar(dynStressTime, mask, 1, seissol::memory::Standard );
    tree.addVar(      ruptureFront,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      mu,                             mask,                 1,      seissol::memory::Standard );
    tree.addVar(      slip,                           mask,                 1,      seissol::memory::Standard );
    tree.addVar(slipStrike, mask, 1, seissol::memory::Standard );
    tree.addVar(slipDip, mask, 1, seissol::memory::Standard );
    tree.addVar(slipRateMagnitude, mask, 1, seissol::memory::Standard );
    tree.addVar(slipRateStrike, mask, 1, seissol::memory::Standard );
    tree.addVar(slipRateDip, mask, 1, seissol::memory::Standard );
    tree.addVar(peakSlipRate, mask, 1, seissol::memory::Standard );
    tree.addVar(tractionXY, mask, 1, seissol::memory::Standard );
    tree.addVar(tractionXZ, mask, 1, seissol::memory::Standard );
    tree.addVar(ds, mask, 1, seissol::memory::Standard );
    tree.addVar(averagedSlip, mask, 1, seissol::memory::Standard );

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
    Variable<real[ numOfPointsPadded ]> d_c;
    Variable<real[ numOfPointsPadded ]> mu_s;
    Variable<real[ numOfPointsPadded ]> mu_d;
    Variable<real[ numOfPointsPadded ]> cohesion;


    virtual void addTo(initializers::LTSTree& tree) {
        seissol::initializers::DynamicRupture::addTo(tree);
        LayerMask mask = LayerMask(Ghost);
        tree.addVar(d_c,                              mask,                 1,      seissol::memory::Standard );
        tree.addVar(mu_s, mask, 1, seissol::memory::Standard );
        tree.addVar(mu_d, mask, 1, seissol::memory::Standard );
        tree.addVar(cohesion, mask,1, seissol::memory::Standard );
    }
};

struct seissol::initializers::LTS_LinearSlipWeakeningForcedRuptureTime : public seissol::initializers::LTS_LinearSlipWeakening {
  Variable<real[ numOfPointsPadded ]>                   forcedRuptureTime;
  Variable<real>                                        tn;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_LinearSlipWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      forcedRuptureTime,                mask,                 1,      seissol::memory::Standard );
    tree.addVar(      tn,                               mask,                 1,      seissol::memory::Standard );
  }
};

struct seissol::initializers::LTS_LinearSlipWeakeningBimaterial : public seissol::initializers::LTS_LinearSlipWeakening {
  Variable<real[ numOfPointsPadded ]>                   regularisedStrength;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_LinearSlipWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      regularisedStrength,                mask,                 1,      seissol::memory::Standard );
  }
};

struct seissol::initializers::LTS_RateAndState : public seissol::initializers::DynamicRupture {
  Variable<real[ numOfPointsPadded ]> rs_a;
  Variable<real[ numOfPointsPadded ]> rs_sl0;
  Variable<real[ numOfPointsPadded ]> stateVariable;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rs_a, mask, 1, seissol::memory::Standard );
    tree.addVar(rs_sl0, mask, 1, seissol::memory::Standard );
    tree.addVar(stateVariable, mask, 1, seissol::memory::Standard );
  }
};


struct seissol::initializers::LTS_RateAndStateFastVelocityWeakening : public seissol::initializers::LTS_RateAndState {
  Variable<real[ numOfPointsPadded ]> rs_srW;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_RateAndState::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rs_srW, mask, 1, seissol::memory::Standard );
  }
};

struct seissol::initializers::LTS_RateAndStateThermalPressurisation : public seissol::initializers::LTS_RateAndStateFastVelocityWeakening {

  Variable<real[numOfPointsPadded]>                               temperature;  //this is TP[1] in fortran
  Variable<real[numOfPointsPadded]>                               pressure;     //this is TP[2] in fortran
  Variable<real[numOfPointsPadded][seissol::dr::numberOfTPGridPoints]>      TP_theta;
  Variable<real[numOfPointsPadded][seissol::dr::numberOfTPGridPoints]>      TP_sigma;
  Variable<real[numOfPointsPadded]>                               TP_halfWidthShearZone;
  Variable<real[numOfPointsPadded]>                               alphaHy;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::LTS_RateAndStateFastVelocityWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      temperature,                mask,                 1,      seissol::memory::Standard );
    tree.addVar(      pressure,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      TP_theta,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      TP_sigma,                   mask,                 1,      seissol::memory::Standard );
    tree.addVar(      TP_halfWidthShearZone,      mask,                 1,      seissol::memory::Standard );
    tree.addVar(      alphaHy,                    mask,                 1,      seissol::memory::Standard );
  }
};


struct seissol::initializers::LTS_ImposedSlipRates : public seissol::initializers::DynamicRupture {
  Variable<real[numOfPointsPadded][6]>                            nucleationStressInFaultCS;

  virtual void addTo(initializers::LTSTree& tree) {
    seissol::initializers::DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(      nucleationStressInFaultCS,  mask,                 1,      seissol::memory::Standard );
  }
};


#endif
