// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_DYNAMICRUPTURE_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_DYNAMICRUPTURE_H_

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"
#include "generated_code/tensor.h"
#include <Initializer/Parameters/DRParameters.h>
#include <Kernels/Common.h>

namespace seissol::initializer {

inline auto allocationModeDR() {
  if constexpr (!isDeviceOn()) {
    return AllocationMode::HostOnly;
  } else {
    return useUSM() ? AllocationMode::HostDeviceUnified : AllocationMode::HostDeviceSplit;
  }
}

// NOTE: for the sake of GPU performance, make sure that NumPaddedPoints is always last.

struct DynamicRupture {
  public:
  DynamicRupture() = default;
  explicit DynamicRupture(const parameters::DRParameters* parameters)
      : nucleationStressInFaultCS(parameters->nucleationCount),
        nucleationPressure(parameters->nucleationCount) {}

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
  Variable<real[6][dr::misc::NumPaddedPoints]> initialStressInFaultCS;
  std::vector<Variable<real[6][dr::misc::NumPaddedPoints]>> nucleationStressInFaultCS;
  // will be always zero, if not using poroelasticity
  Variable<real[dr::misc::NumPaddedPoints]> initialPressure;
  std::vector<Variable<real[dr::misc::NumPaddedPoints]>> nucleationPressure;
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

  Scratchpad<real> idofsPlusOnDevice;
  Scratchpad<real> idofsMinusOnDevice;

  virtual void addTo(LTSTree& tree) {
    const auto mask = LayerMask(Ghost);
    tree.add(timeDerivativePlus, mask, Alignment, AllocationMode::HostOnly, true);
    tree.add(timeDerivativeMinus, mask, Alignment, AllocationMode::HostOnly, true);
    tree.add(timeDerivativePlusDevice, mask, Alignment, AllocationMode::HostOnly, true);
    tree.add(timeDerivativeMinusDevice, mask, Alignment, AllocationMode::HostOnly, true);
    tree.add(imposedStatePlus, mask, PagesizeHeap, allocationModeDR());
    tree.add(imposedStateMinus, mask, PagesizeHeap, allocationModeDR());
    tree.add(godunovData, mask, Alignment, allocationModeDR());
    tree.add(fluxSolverPlus, mask, Alignment, allocationModeDR());
    tree.add(fluxSolverMinus, mask, Alignment, allocationModeDR());
    tree.add(faceInformation, mask, Alignment, AllocationMode::HostOnly, true);
    tree.add(waveSpeedsPlus, mask, Alignment, allocationModeDR(), true);
    tree.add(waveSpeedsMinus, mask, Alignment, allocationModeDR(), true);
    tree.add(drEnergyOutput, mask, Alignment, allocationModeDR());
    tree.add(impAndEta, mask, Alignment, allocationModeDR(), true);
    tree.add(impedanceMatrices, mask, Alignment, allocationModeDR(), true);
    tree.add(initialStressInFaultCS, mask, Alignment, allocationModeDR());
    tree.add(initialPressure, mask, Alignment, allocationModeDR());
    tree.add(ruptureTime, mask, Alignment, allocationModeDR());

    for (auto& nucleation : nucleationStressInFaultCS) {
      tree.add(nucleation, mask, Alignment, allocationModeDR(), true);
    }
    for (auto& nucleation : nucleationPressure) {
      tree.add(nucleation, mask, Alignment, allocationModeDR(), true);
    }

    tree.add(ruptureTimePending, mask, Alignment, allocationModeDR());
    tree.add(dynStressTime, mask, Alignment, allocationModeDR());
    tree.add(dynStressTimePending, mask, Alignment, allocationModeDR());
    tree.add(mu, mask, Alignment, allocationModeDR());
    tree.add(accumulatedSlipMagnitude, mask, Alignment, allocationModeDR());
    tree.add(slip1, mask, Alignment, allocationModeDR());
    tree.add(slip2, mask, Alignment, allocationModeDR());
    tree.add(slipRateMagnitude, mask, Alignment, allocationModeDR());
    tree.add(slipRate1, mask, Alignment, allocationModeDR());
    tree.add(slipRate2, mask, Alignment, allocationModeDR());
    tree.add(peakSlipRate, mask, Alignment, allocationModeDR());
    tree.add(traction1, mask, Alignment, allocationModeDR());
    tree.add(traction2, mask, Alignment, allocationModeDR());
    tree.add(qInterpolatedPlus, mask, Alignment, allocationModeDR());
    tree.add(qInterpolatedMinus, mask, Alignment, allocationModeDR());

    if constexpr (isDeviceOn()) {
      tree.add(idofsPlusOnDevice, LayerMask(), Alignment, AllocationMode::DeviceOnly);
      tree.add(idofsMinusOnDevice, LayerMask(), Alignment, AllocationMode::DeviceOnly);
    }
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

  explicit LTSLinearSlipWeakening(const parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.add(dC, mask, Alignment, allocationModeDR(), true);
    tree.add(muS, mask, Alignment, allocationModeDR(), true);
    tree.add(muD, mask, Alignment, allocationModeDR(), true);
    tree.add(cohesion, mask, Alignment, allocationModeDR(), true);
    tree.add(forcedRuptureTime, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSLinearSlipWeakeningBimaterial : public LTSLinearSlipWeakening {
  Variable<real[dr::misc::NumPaddedPoints]> regularizedStrength;

  explicit LTSLinearSlipWeakeningBimaterial(const parameters::DRParameters* parameters)
      : LTSLinearSlipWeakening(parameters) {}

  void addTo(LTSTree& tree) override {
    LTSLinearSlipWeakening::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.add(regularizedStrength, mask, Alignment, allocationModeDR());
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
  Variable<real[dr::misc::NumPaddedPoints]> rsF0;
  Variable<real[dr::misc::NumPaddedPoints]> rsMuW;
  Variable<real[dr::misc::NumPaddedPoints]> rsB;

  explicit LTSRateAndState(const parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.add(rsA, mask, Alignment, allocationModeDR(), true);
    tree.add(rsSl0, mask, Alignment, allocationModeDR(), true);
    tree.add(stateVariable, mask, Alignment, allocationModeDR());
    tree.add(rsF0, mask, Alignment, allocationModeDR(), true);
    tree.add(rsMuW, mask, Alignment, allocationModeDR(), true);
    tree.add(rsB, mask, Alignment, allocationModeDR(), true);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const override {
    seissol::initializer::DynamicRupture::registerCheckpointVariables(manager, tree);
    manager.registerData("stateVariable", tree, stateVariable);
  }
};

struct LTSRateAndStateFastVelocityWeakening : public LTSRateAndState {
  Variable<real[dr::misc::NumPaddedPoints]> rsSrW;

  explicit LTSRateAndStateFastVelocityWeakening(const parameters::DRParameters* parameters)
      : LTSRateAndState(parameters) {}

  void addTo(LTSTree& tree) override {
    LTSRateAndState::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.add(rsSrW, mask, Alignment, allocationModeDR(), true);
  }
};

struct ThermalPressurization {
  Variable<real[dr::misc::NumPaddedPoints]> temperature;
  Variable<real[dr::misc::NumPaddedPoints]> pressure;
  Variable<real[seissol::dr::misc::NumTpGridPoints][dr::misc::NumPaddedPoints]> theta;
  Variable<real[seissol::dr::misc::NumTpGridPoints][dr::misc::NumPaddedPoints]> sigma;
  Variable<real[dr::misc::NumPaddedPoints]> halfWidthShearZone;
  Variable<real[dr::misc::NumPaddedPoints]> hydraulicDiffusivity;

  void addTo(LTSTree& tree) {
    const auto mask = LayerMask(Ghost);
    tree.add(temperature, mask, Alignment, allocationModeDR());
    tree.add(pressure, mask, Alignment, allocationModeDR());
    tree.add(theta, mask, Alignment, allocationModeDR());
    tree.add(sigma, mask, Alignment, allocationModeDR());
    tree.add(halfWidthShearZone, mask, Alignment, allocationModeDR(), true);
    tree.add(hydraulicDiffusivity, mask, Alignment, allocationModeDR(), true);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   LTSTree* tree) const {
    manager.registerData("temperature", tree, temperature);
    manager.registerData("pressure", tree, pressure);
    manager.registerData("theta", tree, theta);
    manager.registerData("sigma", tree, sigma);
  }
};

struct LTSRateAndStateThermalPressurization : public LTSRateAndState, public ThermalPressurization {
  explicit LTSRateAndStateThermalPressurization(const parameters::DRParameters* parameters)
      : LTSRateAndState(parameters) {}

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
  explicit LTSRateAndStateThermalPressurizationFastVelocityWeakening(
      const parameters::DRParameters* parameters)
      : LTSRateAndStateFastVelocityWeakening(parameters) {}

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

  explicit LTSImposedSlipRates(const parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(LTSTree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.add(imposedSlipDirection1, mask, Alignment, allocationModeDR(), true);
    tree.add(imposedSlipDirection2, mask, Alignment, allocationModeDR(), true);
    tree.add(onsetTime, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesYoffe : public LTSImposedSlipRates {
  Variable<real[dr::misc::NumPaddedPoints]> tauS;
  Variable<real[dr::misc::NumPaddedPoints]> tauR;

  explicit LTSImposedSlipRatesYoffe(const parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.add(tauS, mask, Alignment, allocationModeDR(), true);
    tree.add(tauR, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesGaussian : public LTSImposedSlipRates {
  Variable<real[dr::misc::NumPaddedPoints]> riseTime;

  explicit LTSImposedSlipRatesGaussian(const parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}

  void addTo(LTSTree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    const auto mask = LayerMask(Ghost);
    tree.add(riseTime, mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesDelta : public LTSImposedSlipRates {
  explicit LTSImposedSlipRatesDelta(const parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_DYNAMICRUPTURE_H_
