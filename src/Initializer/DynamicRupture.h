// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_DYNAMICRUPTURE_H_
#define SEISSOL_SRC_INITIALIZER_DYNAMICRUPTURE_H_

#include "DynamicRupture/Misc.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Layer.h"
#include "Initializer/Typedefs.h"
#include "Parallel/Helper.h"
#include "Tree/Layer.h"
#include "generated_code/tensor.h"

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
  Variable<real[dr::misc::NumPaddedPoints][6]> initialStressInFaultCS;
  Variable<real[dr::misc::NumPaddedPoints][6]> nucleationStressInFaultCS;
  // will be always zero, if not using poroelasticity
  Variable<real[dr::misc::NumPaddedPoints]> initialPressure;
  Variable<real[dr::misc::NumPaddedPoints]> nucleationPressure;
  Variable<real[dr::misc::NumPaddedPoints]> mu;
  Variable<real[dr::misc::NumPaddedPoints]> accumulatedSlipMagnitude;
  Variable<real[dr::misc::NumPaddedPoints]>
      slip1; // slip at given fault node along local direction 1
  Variable<real[dr::misc::NumPaddedPoints]>
      slip2; // slip at given fault node along local direction 2
  Variable<real[dr::misc::NumPaddedPoints]> slipRateMagnitude;
  Variable<real[dr::misc::NumPaddedPoints]>
      slipRate1; // slip rate at given fault node along local direction 1
  Variable<real[dr::misc::NumPaddedPoints]>
      slipRate2; // slip rate at given fault node along local direction 2
  Variable<real[dr::misc::NumPaddedPoints]> ruptureTime;
  Variable<real[dr::misc::NumPaddedPoints]> dynStressTime;
  Variable<bool[dr::misc::NumPaddedPoints]> ruptureTimePending;
  Variable<bool[dr::misc::NumPaddedPoints]> dynStressTimePending;
  Variable<real[dr::misc::NumPaddedPoints]> peakSlipRate;
  Variable<real[dr::misc::NumPaddedPoints]> traction1;
  Variable<real[dr::misc::NumPaddedPoints]> traction2;
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
    tree.addVar(impedanceMatrices, mask, Alignment, allocationModeDR(), true);
    tree.addVar(initialStressInFaultCS, mask, 1, allocationModeDR());
    tree.addVar(nucleationStressInFaultCS, mask, 1, allocationModeDR(), true);
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
  Variable<real[dr::misc::NumPaddedPoints]> dC;
  Variable<real[dr::misc::NumPaddedPoints]> muS;
  Variable<real[dr::misc::NumPaddedPoints]> muD;
  Variable<real[dr::misc::NumPaddedPoints]> cohesion;
  Variable<real[dr::misc::NumPaddedPoints]> forcedRuptureTime;

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
  Variable<real[dr::misc::NumPaddedPoints]> regularisedStrength;

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
  Variable<real[dr::misc::NumPaddedPoints]> rsA;
  Variable<real[dr::misc::NumPaddedPoints]> rsSl0;
  Variable<real[dr::misc::NumPaddedPoints]> stateVariable;

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
  Variable<real[dr::misc::NumPaddedPoints]> rsSrW;

  void addTo(LTSTree& tree) override {
    LTSRateAndState::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(rsSrW, mask, 1, allocationModeDR(), true);
  }
};

struct LTSRateAndStateThermalPressurization : public LTSRateAndStateFastVelocityWeakening {

  Variable<real[dr::misc::NumPaddedPoints]> temperature;
  Variable<real[dr::misc::NumPaddedPoints]> pressure;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> theta;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> sigma;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> thetaTmpBuffer;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> sigmaTmpBuffer;
  Variable<real[dr::misc::NumPaddedPoints]> faultStrength;
  Variable<real[dr::misc::NumPaddedPoints]> halfWidthShearZone;
  Variable<real[dr::misc::NumPaddedPoints]> hydraulicDiffusivity;

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
  Variable<real[dr::misc::NumPaddedPoints]> imposedSlipDirection1;
  Variable<real[dr::misc::NumPaddedPoints]> imposedSlipDirection2;
  Variable<real[dr::misc::NumPaddedPoints]> onsetTime;

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(imposedSlipDirection1, mask, 1, allocationModeDR(), true);
    tree.addVar(imposedSlipDirection2, mask, 1, allocationModeDR(), true);
    tree.addVar(onsetTime, mask, 1, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesYoffe : public LTSImposedSlipRates {
  Variable<real[dr::misc::NumPaddedPoints]> tauS;
  Variable<real[dr::misc::NumPaddedPoints]> tauR;

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(tauS, mask, 1, allocationModeDR(), true);
    tree.addVar(tauR, mask, 1, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesGaussian : public LTSImposedSlipRates {
  Variable<real[dr::misc::NumPaddedPoints]> riseTime;

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    LayerMask mask = LayerMask(Ghost);
    tree.addVar(riseTime, mask, 1, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesDelta : public LTSImposedSlipRates {};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_DYNAMICRUPTURE_H_
