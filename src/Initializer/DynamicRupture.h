/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 *http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

#include "DynamicRupture/Misc.h"
#include "IO/Instance/Checkpoint/CheckpointManager.hpp"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Layer.hpp"
#include "Initializer/typedefs.hpp"
#include "Parallel/Helper.hpp"
#include "generated_code/tensor.h"
#include "tree/Layer.hpp"

namespace seissol::initializer {

inline auto allocationModeDR() {
#ifndef ACL_DEVICE
  return AllocationMode::HostOnly;
#else
  return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
#endif
}

struct DynamicRupture {
  public:
  virtual ~DynamicRupture() = default;
  Variable<real*> timeDerivativePlus;
  Variable<real*> timeDerivativeMinus;
  Variable<real*> timeDerivativePlusDevice;
  Variable<real*> timeDerivativeMinusDevice;
  Variable<real[tensor::QInterpolated::size()]> imposedStatePlus;
  Variable<real[tensor::QInterpolated::size()]> imposedStateMinus;
  Variable<DRGodunovData> godunovData;
  Variable<real[tensor::fluxSolver::size()]> fluxSolverPlus;
  Variable<real[tensor::fluxSolver::size()]> fluxSolverMinus;
  Variable<DRFaceInformation> faceInformation;
  Variable<model::IsotropicWaveSpeeds> waveSpeedsPlus;
  Variable<model::IsotropicWaveSpeeds> waveSpeedsMinus;
  Variable<DREnergyOutput> drEnergyOutput;

  Variable<seissol::dr::ImpedancesAndEta> impAndEta;
  Variable<seissol::dr::ImpedanceMatrices> impedanceMatrices;
  // size padded for vectorization
  // CS = coordinate system
  Variable<real[dr::misc::numPaddedPoints][6]> initialStressInFaultCS;
  Variable<real[dr::misc::numPaddedPoints][6]> nucleationStressInFaultCS;
  Variable<real[dr::misc::numPaddedPoints][6]> nucleationStressInFaultCS2;
  // will be always zero, if not using poroelasticity
  Variable<real[dr::misc::numPaddedPoints]> initialPressure;
  Variable<real[dr::misc::numPaddedPoints]> nucleationPressure;
  Variable<real[dr::misc::numPaddedPoints]> mu;
  Variable<real[dr::misc::numPaddedPoints]> accumulatedSlipMagnitude;
  Variable<real[dr::misc::numPaddedPoints]>
      slip1; // slip at given fault node along local direction 1
  Variable<real[dr::misc::numPaddedPoints]>
      slip2; // slip at given fault node along local direction 2
  Variable<real[dr::misc::numPaddedPoints]> slipRateMagnitude;
  Variable<real[dr::misc::numPaddedPoints]>
      slipRate1; // slip rate at given fault node along local direction 1
  Variable<real[dr::misc::numPaddedPoints]>
      slipRate2; // slip rate at given fault node along local direction 2
  Variable<real[dr::misc::numPaddedPoints]> ruptureTime;
  Variable<real[dr::misc::numPaddedPoints]> dynStressTime;
  Variable<bool[dr::misc::numPaddedPoints]> ruptureTimePending;
  Variable<bool[dr::misc::numPaddedPoints]> dynStressTimePending;
  Variable<real[dr::misc::numPaddedPoints]> peakSlipRate;
  Variable<real[dr::misc::numPaddedPoints]> traction1;
  Variable<real[dr::misc::numPaddedPoints]> traction2;
  Variable<real[ConvergenceOrder][tensor::QInterpolated::size()]> qInterpolatedPlus;
  Variable<real[ConvergenceOrder][tensor::QInterpolated::size()]> qInterpolatedMinus;

#ifdef ACL_DEVICE
  ScratchpadMemory idofsPlusOnDevice;
  ScratchpadMemory idofsMinusOnDevice;
#endif

  virtual void addTo(LTSTree& tree) {
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(timeDerivativePlus, mask, 1, AllocationMode::HostOnly, true);
    tree.addVar(timeDerivativeMinus, mask, 1, AllocationMode::HostOnly, true);
    tree.addVar(timeDerivativePlusDevice, mask, 1, AllocationMode::HostOnly, true);
    tree.addVar(timeDerivativeMinusDevice, mask, 1, AllocationMode::HostOnly, true);
    tree.addVar(imposedStatePlus, mask, PagesizeHeap, allocationModeDR());
    tree.addVar(imposedStateMinus, mask, PagesizeHeap, allocationModeDR());
    tree.addVar(godunovData, mask, 1, allocationModeDR());
    tree.addVar(fluxSolverPlus, mask, 1, allocationModeDR());
    tree.addVar(fluxSolverMinus, mask, 1, allocationModeDR());
    tree.addVar(faceInformation, mask, 1, AllocationMode::HostOnly, true);
    tree.addVar(waveSpeedsPlus, mask, 1, allocationModeDR(), true);
    tree.addVar(waveSpeedsMinus, mask, 1, allocationModeDR(), true);
    tree.addVar(drEnergyOutput, mask, Alignment, allocationModeDR());
    tree.addVar(impAndEta, mask, 1, allocationModeDR(), true);
    tree.addVar(impedanceMatrices, mask, 1, allocationModeDR(), true);
    tree.addVar(initialStressInFaultCS, mask, 1, allocationModeDR());
    tree.addVar(nucleationStressInFaultCS, mask, 1, allocationModeDR(), true);
    tree.addVar(nucleationStressInFaultCS2, mask, 1, allocationModeDR(), true);
    tree.addVar(initialPressure, mask, 1, allocationModeDR());
    tree.addVar(nucleationPressure, mask, 1, allocationModeDR(), true);
    tree.addVar(ruptureTime, mask, 1, allocationModeDR());

    tree.addVar(ruptureTimePending, mask, 1, allocationModeDR());
    tree.addVar(dynStressTime, mask, 1, allocationModeDR());
    tree.addVar(dynStressTimePending, mask, 1, allocationModeDR());
    tree.addVar(mu, mask, 1, allocationModeDR());
    tree.addVar(accumulatedSlipMagnitude, mask, 1, allocationModeDR());
    tree.addVar(slip1, mask, 1, allocationModeDR());
    tree.addVar(slip2, mask, 1, allocationModeDR());
    tree.addVar(slipRateMagnitude, mask, 1, allocationModeDR());
    tree.addVar(slipRate1, mask, 1, allocationModeDR());
    tree.addVar(slipRate2, mask, 1, allocationModeDR());
    tree.addVar(peakSlipRate, mask, 1, allocationModeDR());
    tree.addVar(traction1, mask, 1, allocationModeDR());
    tree.addVar(traction2, mask, 1, allocationModeDR());
    tree.addVar(qInterpolatedPlus, mask, Alignment, allocationModeDR());
    tree.addVar(qInterpolatedMinus, mask, Alignment, allocationModeDR());

#ifdef ACL_DEVICE
    tree.addScratchpadMemory(idofsPlusOnDevice, 1, AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(idofsMinusOnDevice, 1, AllocationMode::DeviceOnly);
#endif
  }

  virtual void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                           LTSTree* tree) {
    manager.registerData("initialStressInFaultCS", tree, initialStressInFaultCS);
    manager.registerData("initialPressure", tree, initialPressure);
    manager.registerData("mu", tree, mu);
    manager.registerData("slipRate1", tree, slipRate1);
    manager.registerData("slipRate2", tree, slipRate2);
    manager.registerData("accumulatedSlipMagnitude", tree, accumulatedSlipMagnitude);
    manager.registerData("slip1", tree, slip1);
    manager.registerData("slip2", tree, slip2);
    manager.registerData("peakSlipRate", tree, peakSlipRate);
    manager.registerData("ruptureTime", tree, ruptureTime);
    manager.registerData("ruptureTimePending", tree, ruptureTimePending);
    manager.registerData("dynStressTime", tree, dynStressTime);
    manager.registerData("dynStressTimePending", tree, dynStressTimePending);
    manager.registerData("drEnergyOutput", tree, drEnergyOutput);
  }
};

struct LTSLinearSlipWeakening : public DynamicRupture {
  Variable<real[dr::misc::numPaddedPoints]> dC;
  Variable<real[dr::misc::numPaddedPoints]> muS;
  Variable<real[dr::misc::numPaddedPoints]> muD;
  Variable<real[dr::misc::numPaddedPoints]> cohesion;
  Variable<real[dr::misc::numPaddedPoints]> forcedRuptureTime;

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(dC, mask, 1, allocationModeDR(), true);
    tree.addVar(muS, mask, 1, allocationModeDR(), true);
    tree.addVar(muD, mask, 1, allocationModeDR(), true);
    tree.addVar(cohesion, mask, 1, allocationModeDR(), true);
    tree.addVar(forcedRuptureTime, mask, 1, allocationModeDR(), true);
  }
};

struct LTSLinearSlipWeakeningBimaterial : public LTSLinearSlipWeakening {
  Variable<real[dr::misc::numPaddedPoints]> regularisedStrength;

  void addTo(LTSTree& tree) override {
    LTSLinearSlipWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(regularisedStrength, mask, 1, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) override {
    seissol::initializer::LTSLinearSlipWeakening::registerCheckpointVariables(manager, tree);
    manager.registerData("regularisedStrength", tree, regularisedStrength);
  }
};

struct LTSRateAndState : public DynamicRupture {
  Variable<real[dr::misc::numPaddedPoints]> rsA;
  Variable<real[dr::misc::numPaddedPoints]> rsSl0;
  Variable<real[dr::misc::numPaddedPoints]> stateVariable;

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsA, mask, 1, allocationModeDR(), true);
    tree.addVar(rsSl0, mask, 1, allocationModeDR(), true);
    tree.addVar(stateVariable, mask, 1, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) override {
    seissol::initializer::DynamicRupture::registerCheckpointVariables(manager, tree);
    manager.registerData("stateVariable", tree, stateVariable);
  }
};

struct LTSRateAndStateFastVelocityWeakening : public LTSRateAndState {
  Variable<real[dr::misc::numPaddedPoints]> rsSrW;

  void addTo(LTSTree& tree) override {
    LTSRateAndState::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsSrW, mask, 1, allocationModeDR(), true);
  }
};

struct LTSRateAndStateThermalPressurization : public LTSRateAndStateFastVelocityWeakening {

  Variable<real[dr::misc::numPaddedPoints]> temperature;
  Variable<real[dr::misc::numPaddedPoints]> pressure;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> theta;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigma;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> thetaTmpBuffer;
  Variable<real[dr::misc::numPaddedPoints][seissol::dr::misc::numberOfTPGridPoints]> sigmaTmpBuffer;
  Variable<real[dr::misc::numPaddedPoints]> faultStrength;
  Variable<real[dr::misc::numPaddedPoints]> halfWidthShearZone;
  Variable<real[dr::misc::numPaddedPoints]> hydraulicDiffusivity;

  void addTo(LTSTree& tree) override {
    LTSRateAndStateFastVelocityWeakening::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(temperature, mask, Alignment, allocationModeDR());
    tree.addVar(pressure, mask, Alignment, allocationModeDR());
    tree.addVar(theta, mask, Alignment, allocationModeDR(), true);
    tree.addVar(sigma, mask, Alignment, allocationModeDR(), true);
    tree.addVar(thetaTmpBuffer, mask, Alignment, allocationModeDR(), true);
    tree.addVar(sigmaTmpBuffer, mask, Alignment, allocationModeDR(), true);
    tree.addVar(faultStrength, mask, Alignment, allocationModeDR(), true);
    tree.addVar(halfWidthShearZone, mask, Alignment, allocationModeDR(), true);
    tree.addVar(hydraulicDiffusivity, mask, Alignment, allocationModeDR(), true);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) override {
    seissol::initializer::LTSRateAndStateFastVelocityWeakening::registerCheckpointVariables(manager,
                                                                                            tree);
    manager.registerData("temperature", tree, temperature);
    manager.registerData("pressure", tree, pressure);
  }
};

struct LTSImposedSlipRates : public DynamicRupture {
  Variable<real[dr::misc::numPaddedPoints]> imposedSlipDirection1;
  Variable<real[dr::misc::numPaddedPoints]> imposedSlipDirection2;
  Variable<real[dr::misc::numPaddedPoints]> onsetTime;

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(imposedSlipDirection1, mask, 1, allocationModeDR(), true);
    tree.addVar(imposedSlipDirection2, mask, 1, allocationModeDR(), true);
    tree.addVar(onsetTime, mask, 1, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesYoffe : public LTSImposedSlipRates {
  Variable<real[dr::misc::numPaddedPoints]> tauS;
  Variable<real[dr::misc::numPaddedPoints]> tauR;

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(tauS, mask, 1, allocationModeDR(), true);
    tree.addVar(tauR, mask, 1, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesGaussian : public LTSImposedSlipRates {
  Variable<real[dr::misc::numPaddedPoints]> riseTime;

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(riseTime, mask, 1, allocationModeDR(), true);
  }
};

} // namespace seissol::initializer

#endif // INITIALIZER_DR_H_
