// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_DYNAMICRUPTURE_H_
#define SEISSOL_SRC_INITIALIZER_DYNAMICRUPTURE_H_

#include "DynamicRupture/Misc.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"
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
    const auto mask = LayerMask(Ghost);
    tree.addVar(timeDerivativePlus, mask, Alignment, AllocationMode::HostOnly, true);
    tree.addVar(timeDerivativeMinus, mask, Alignment, AllocationMode::HostOnly, true);
    tree.addVar(timeDerivativePlusDevice, mask, Alignment, AllocationMode::HostOnly, true);
    tree.addVar(timeDerivativeMinusDevice, mask, Alignment, AllocationMode::HostOnly, true);
    tree.addVar(imposedStatePlus, mask, PagesizeHeap, allocationModeDR());
    tree.addVar(imposedStateMinus, mask, PagesizeHeap, allocationModeDR());
    tree.addVar(godunovData, mask, Alignment, allocationModeDR());
    tree.addVar(fluxSolverPlus, mask, Alignment, allocationModeDR());
    tree.addVar(fluxSolverMinus, mask, Alignment, allocationModeDR());
    tree.addVar(faceInformation, mask, Alignment, AllocationMode::HostOnly, true);
    tree.addVar(waveSpeedsPlus, mask, Alignment, allocationModeDR(), true);
    tree.addVar(waveSpeedsMinus, mask, Alignment, allocationModeDR(), true);
    tree.addVar(drEnergyOutput, mask, Alignment, allocationModeDR());
    tree.addVar(impAndEta, mask, Alignment, allocationModeDR(), true);
    tree.addVar(impedanceMatrices, mask, Alignment, allocationModeDR(), true);
    tree.addVar(initialStressInFaultCS, mask, Alignment, allocationModeDR());
    tree.addVar(nucleationStressInFaultCS, mask, Alignment, allocationModeDR(), true);
    tree.addVar(initialPressure, mask, Alignment, allocationModeDR());
    tree.addVar(nucleationPressure, mask, Alignment, allocationModeDR(), true);
    tree.addVar(ruptureTime, mask, Alignment, allocationModeDR());

    tree.addVar(ruptureTimePending, mask, Alignment, allocationModeDR());
    tree.addVar(dynStressTime, mask, Alignment, allocationModeDR());
    tree.addVar(dynStressTimePending, mask, Alignment, allocationModeDR());
    tree.addVar(mu, mask, Alignment, allocationModeDR());
    tree.addVar(accumulatedSlipMagnitude, mask, Alignment, allocationModeDR());
    tree.addVar(slip1, mask, Alignment, allocationModeDR());
    tree.addVar(slip2, mask, Alignment, allocationModeDR());
    tree.addVar(slipRateMagnitude, mask, Alignment, allocationModeDR());
    tree.addVar(slipRate1, mask, Alignment, allocationModeDR());
    tree.addVar(slipRate2, mask, Alignment, allocationModeDR());
    tree.addVar(peakSlipRate, mask, Alignment, allocationModeDR());
    tree.addVar(traction1, mask, Alignment, allocationModeDR());
    tree.addVar(traction2, mask, Alignment, allocationModeDR());
    tree.addVar(qInterpolatedPlus, mask, Alignment, allocationModeDR());
    tree.addVar(qInterpolatedMinus, mask, Alignment, allocationModeDR());

#ifdef ACL_DEVICE
    tree.addScratchpadMemory(idofsPlusOnDevice, Alignment, AllocationMode::DeviceOnly);
    tree.addScratchpadMemory(idofsMinusOnDevice, Alignment, AllocationMode::DeviceOnly);
#endif
  }

  virtual void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                           LTSTree* tree) const {
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
    const auto mask = LayerMask(Ghost);
    tree.addVar(dC, mask, Alignment, allocationModeDR(), true);
    tree.addVar(muS, mask, Alignment, allocationModeDR(), true);
    tree.addVar(muD, mask, Alignment, allocationModeDR(), true);
    tree.addVar(cohesion, mask, Alignment, allocationModeDR(), true);
    tree.addVar(forcedRuptureTime, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSLinearSlipWeakeningBimaterial : public LTSLinearSlipWeakening {
  Variable<real[dr::misc::NumPaddedPoints]> regularizedStrength;

  void addTo(LTSTree& tree) override {
    LTSLinearSlipWeakening::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.addVar(regularizedStrength, mask, Alignment, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const override {
    seissol::initializer::LTSLinearSlipWeakening::registerCheckpointVariables(manager, tree);
    manager.registerData("regularizedStrength", tree, regularizedStrength);
  }
};

struct LTSRateAndState : public DynamicRupture {
  Variable<real[dr::misc::NumPaddedPoints]> rsA;
  Variable<real[dr::misc::NumPaddedPoints]> rsSl0;
  Variable<real[dr::misc::NumPaddedPoints]> stateVariable;

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.addVar(rsA, mask, Alignment, allocationModeDR(), true);
    tree.addVar(rsSl0, mask, Alignment, allocationModeDR(), true);
    tree.addVar(stateVariable, mask, Alignment, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const override {
    seissol::initializer::DynamicRupture::registerCheckpointVariables(manager, tree);
    manager.registerData("stateVariable", tree, stateVariable);
  }
};

struct LTSRateAndStateFastVelocityWeakening : public LTSRateAndState {
  Variable<real[dr::misc::NumPaddedPoints]> rsSrW;

  void addTo(LTSTree& tree) override {
    LTSRateAndState::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.addVar(rsSrW, mask, Alignment, allocationModeDR(), true);
  }
};

struct ThermalPressurization {
  Variable<real[dr::misc::NumPaddedPoints]> temperature;
  Variable<real[dr::misc::NumPaddedPoints]> pressure;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> theta;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> sigma;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> thetaTmpBuffer;
  Variable<real[dr::misc::NumPaddedPoints][seissol::dr::misc::NumTpGridPoints]> sigmaTmpBuffer;
  Variable<real[dr::misc::NumPaddedPoints]> faultStrength;
  Variable<real[dr::misc::NumPaddedPoints]> halfWidthShearZone;
  Variable<real[dr::misc::NumPaddedPoints]> hydraulicDiffusivity;

  void addTo(LTSTree& tree) {
    const auto mask = LayerMask(Ghost);
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
                                   LTSTree* tree) const {
    manager.registerData("temperature", tree, temperature);
    manager.registerData("pressure", tree, pressure);
  }
};

struct LTSRateAndStateThermalPressurization : public LTSRateAndState, public ThermalPressurization {
  void addTo(LTSTree& tree) override {
    LTSRateAndState::addTo(tree);
    ThermalPressurization::addTo(tree);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const override {
    LTSRateAndState::registerCheckpointVariables(manager, tree);
    ThermalPressurization::registerCheckpointVariables(manager, tree);
  }
};

struct LTSRateAndStateThermalPressurizationFastVelocityWeakening
    : public LTSRateAndStateFastVelocityWeakening,
      public ThermalPressurization {
  void addTo(LTSTree& tree) override {
    LTSRateAndStateFastVelocityWeakening::addTo(tree);
    ThermalPressurization::addTo(tree);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const override {
    LTSRateAndStateFastVelocityWeakening::registerCheckpointVariables(manager, tree);
    ThermalPressurization::registerCheckpointVariables(manager, tree);
  }
};

struct LTSImposedSlipRates : public DynamicRupture {
  Variable<real[dr::misc::NumPaddedPoints]> imposedSlipDirection1;
  Variable<real[dr::misc::NumPaddedPoints]> imposedSlipDirection2;
  Variable<real[dr::misc::NumPaddedPoints]> onsetTime;

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.addVar(imposedSlipDirection1, mask, Alignment, allocationModeDR(), true);
    tree.addVar(imposedSlipDirection2, mask, Alignment, allocationModeDR(), true);
    tree.addVar(onsetTime, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesYoffe : public LTSImposedSlipRates {
  Variable<real[dr::misc::NumPaddedPoints]> tauS;
  Variable<real[dr::misc::NumPaddedPoints]> tauR;

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.addVar(tauS, mask, Alignment, allocationModeDR(), true);
    tree.addVar(tauR, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesGaussian : public LTSImposedSlipRates {
  Variable<real[dr::misc::NumPaddedPoints]> riseTime;

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.addVar(riseTime, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesDelta : public LTSImposedSlipRates {};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_DYNAMICRUPTURE_H_
