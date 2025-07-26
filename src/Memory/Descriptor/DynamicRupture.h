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
#include <Memory/Tree/Backmap.h>

namespace seissol {

inline auto allocationModeDR() {
  using namespace seissol::initializer;
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
  std::size_t nucleationCount{1};
  explicit DynamicRupture(const initializer::parameters::DRParameters* parameters)
      : nucleationCount(parameters->nucleationCount) {}

  virtual ~DynamicRupture() = default;
  struct TimeDerivativePlus : public initializer::Variable<real*> {};
  struct TimeDerivativeMinus : public initializer::Variable<real*> {};
  struct TimeDerivativePlusDevice : public initializer::Variable<real*> {};
  struct TimeDerivativeMinusDevice : public initializer::Variable<real*> {};
  struct ImposedStatePlus : public initializer::Variable<real[tensor::QInterpolated::size()]> {};
  struct ImposedStateMinus : public initializer::Variable<real[tensor::QInterpolated::size()]> {};
  struct GodunovData : public initializer::Variable<DRGodunovData> {};
  struct FluxSolverPlus : public initializer::Variable<real[tensor::fluxSolver::size()]> {};
  struct FluxSolverMinus : public initializer::Variable<real[tensor::fluxSolver::size()]> {};
  struct FaceInformation : public initializer::Variable<DRFaceInformation> {};
  struct WaveSpeedsPlus : public initializer::Variable<model::IsotropicWaveSpeeds> {};
  struct WaveSpeedsMinus : public initializer::Variable<model::IsotropicWaveSpeeds> {};
  struct DREnergyOutputVar : public initializer::Variable<DREnergyOutput> {};

  struct ImpAndEta : public initializer::Variable<seissol::dr::ImpedancesAndEta> {};
  struct ImpedanceMatrices : public initializer::Variable<seissol::dr::ImpedanceMatrices> {};
  // size padded for vectorization
  // CS = coordinate system
  struct InitialStressInFaultCS : public initializer::Variable<real[6][dr::misc::NumPaddedPoints]> {
  };
  struct NucleationStressInFaultCS
      : public initializer::Variable<real[6][dr::misc::NumPaddedPoints]> {};
  // will be always zero, if not using poroelasticity
  struct InitialPressure : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct NucleationPressure : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct Mu : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct AccumulatedSlipMagnitude : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {
  };
  struct Slip1 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {
  }; // slip at given fault node along local direction 1
  struct Slip2 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {
  }; // slip at given fault node along local direction 2
  struct SlipRateMagnitude : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct SlipRate1 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {
  }; // slip rate at given fault node along local direction 1
  struct SlipRate2 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {
  }; // slip rate at given fault node along local direction 2
  struct RuptureTime : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct DynStressTime : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct RuptureTimePending : public initializer::Variable<bool[dr::misc::NumPaddedPoints]> {};
  struct DynStressTimePending : public initializer::Variable<bool[dr::misc::NumPaddedPoints]> {};
  struct PeakSlipRate : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct Traction1 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct Traction2 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct QInterpolatedPlus
      : public initializer::Variable<real[ConvergenceOrder][tensor::QInterpolated::size()]> {};
  struct QInterpolatedMinus
      : public initializer::Variable<real[ConvergenceOrder][tensor::QInterpolated::size()]> {};

  struct IdofsPlusOnDevice : public initializer::Scratchpad<real> {};
  struct IdofsMinusOnDevice : public initializer::Scratchpad<real> {};

  struct DynrupVarmap : public initializer::GenericVarmap {};

  using Tree = initializer::LTSTree<DynrupVarmap>;
  using Layer = initializer::Layer<DynrupVarmap>;
  using Ref = initializer::Layer<DynrupVarmap>::CellRef;
  using Backmap = initializer::StorageBackmap<1>;

  virtual void addTo(Tree& tree) {
    using namespace seissol::initializer;
    const auto mask = LayerMask(Ghost);
    tree.add<TimeDerivativePlus>(mask, Alignment, AllocationMode::HostOnly, true);
    tree.add<TimeDerivativeMinus>(mask, Alignment, AllocationMode::HostOnly, true);
    tree.add<TimeDerivativePlusDevice>(mask, Alignment, AllocationMode::HostOnly, true);
    tree.add<TimeDerivativeMinusDevice>(mask, Alignment, AllocationMode::HostOnly, true);
    tree.add<ImposedStatePlus>(mask, PagesizeHeap, allocationModeDR());
    tree.add<ImposedStateMinus>(mask, PagesizeHeap, allocationModeDR());
    tree.add<GodunovData>(mask, Alignment, allocationModeDR());
    tree.add<FluxSolverPlus>(mask, Alignment, allocationModeDR());
    tree.add<FluxSolverMinus>(mask, Alignment, allocationModeDR());
    tree.add<FaceInformation>(mask, Alignment, AllocationMode::HostOnly, true);
    tree.add<WaveSpeedsPlus>(mask, Alignment, allocationModeDR(), true);
    tree.add<WaveSpeedsMinus>(mask, Alignment, allocationModeDR(), true);
    tree.add<DREnergyOutputVar>(mask, Alignment, allocationModeDR());
    tree.add<ImpAndEta>(mask, Alignment, allocationModeDR(), true);
    tree.add<ImpedanceMatrices>(mask, Alignment, allocationModeDR(), true);
    tree.add<InitialStressInFaultCS>(mask, Alignment, allocationModeDR());
    tree.add<InitialPressure>(mask, Alignment, allocationModeDR());
    tree.add<RuptureTime>(mask, Alignment, allocationModeDR());

    tree.add<NucleationStressInFaultCS>(mask, Alignment, allocationModeDR(), true, nucleationCount);
    tree.add<NucleationPressure>(mask, Alignment, allocationModeDR(), true, nucleationCount);

    tree.add<RuptureTimePending>(mask, Alignment, allocationModeDR());
    tree.add<DynStressTime>(mask, Alignment, allocationModeDR());
    tree.add<DynStressTimePending>(mask, Alignment, allocationModeDR());
    tree.add<Mu>(mask, Alignment, allocationModeDR());
    tree.add<AccumulatedSlipMagnitude>(mask, Alignment, allocationModeDR());
    tree.add<Slip1>(mask, Alignment, allocationModeDR());
    tree.add<Slip2>(mask, Alignment, allocationModeDR());
    tree.add<SlipRateMagnitude>(mask, Alignment, allocationModeDR());
    tree.add<SlipRate1>(mask, Alignment, allocationModeDR());
    tree.add<SlipRate2>(mask, Alignment, allocationModeDR());
    tree.add<PeakSlipRate>(mask, Alignment, allocationModeDR());
    tree.add<Traction1>(mask, Alignment, allocationModeDR());
    tree.add<Traction2>(mask, Alignment, allocationModeDR());
    tree.add<QInterpolatedPlus>(mask, Alignment, allocationModeDR());
    tree.add<QInterpolatedMinus>(mask, Alignment, allocationModeDR());

    if constexpr (isDeviceOn()) {
      tree.add<IdofsPlusOnDevice>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
      tree.add<IdofsMinusOnDevice>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
    }
  }

  virtual void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                           Tree* tree) const {
    manager.registerData<InitialStressInFaultCS>("initialStressInFaultCS", tree);
    manager.registerData<InitialPressure>("initialPressure", tree);
    manager.registerData<Mu>("mu", tree);
    manager.registerData<SlipRate1>("slipRate1", tree);
    manager.registerData<SlipRate2>("slipRate2", tree);
    manager.registerData<AccumulatedSlipMagnitude>("accumulatedSlipMagnitude", tree);
    manager.registerData<Slip1>("slip1", tree);
    manager.registerData<Slip2>("slip2", tree);
    manager.registerData<PeakSlipRate>("peakSlipRate", tree);
    manager.registerData<RuptureTime>("ruptureTime", tree);
    manager.registerData<RuptureTimePending>("ruptureTimePending", tree);
    manager.registerData<DynStressTime>("dynStressTime", tree);
    manager.registerData<DynStressTimePending>("dynStressTimePending", tree);
    manager.registerData<DREnergyOutputVar>("drEnergyOutput", tree);
  }
};

struct LTSLinearSlipWeakening : public DynamicRupture {
  struct DC : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct MuS : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct MuD : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct Cohesion : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct ForcedRuptureTime : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  explicit LTSLinearSlipWeakening(const initializer::parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(Tree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<DC>(mask, Alignment, allocationModeDR(), true);
    tree.add<MuS>(mask, Alignment, allocationModeDR(), true);
    tree.add<MuD>(mask, Alignment, allocationModeDR(), true);
    tree.add<Cohesion>(mask, Alignment, allocationModeDR(), true);
    tree.add<ForcedRuptureTime>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSLinearSlipWeakeningBimaterial : public LTSLinearSlipWeakening {
  struct RegularizedStrength : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  explicit LTSLinearSlipWeakeningBimaterial(const initializer::parameters::DRParameters* parameters)
      : LTSLinearSlipWeakening(parameters) {}

  void addTo(Tree& tree) override {
    LTSLinearSlipWeakening::addTo(tree);
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<RegularizedStrength>(mask, Alignment, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Tree* tree) const override {
    LTSLinearSlipWeakening::registerCheckpointVariables(manager, tree);
    manager.registerData<RegularizedStrength>("regularizedStrength", tree);
  }
};

struct LTSRateAndState : public DynamicRupture {
  struct RsA : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct RsSl0 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct StateVariable : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  explicit LTSRateAndState(const initializer::parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(Tree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<RsA>(mask, Alignment, allocationModeDR(), true);
    tree.add<RsSl0>(mask, Alignment, allocationModeDR(), true);
    tree.add<StateVariable>(mask, Alignment, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Tree* tree) const override {
    DynamicRupture::registerCheckpointVariables(manager, tree);
    manager.registerData<StateVariable>("stateVariable", tree);
  }
};

struct LTSRateAndStateFastVelocityWeakening : public LTSRateAndState {
  struct RsSrW : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  explicit LTSRateAndStateFastVelocityWeakening(
      const initializer::parameters::DRParameters* parameters)
      : LTSRateAndState(parameters) {}

  void addTo(Tree& tree) override {
    LTSRateAndState::addTo(tree);
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<RsSrW>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSThermalPressurization {
  struct Temperature : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct Pressure : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct Theta : public initializer::Variable<
                     real[seissol::dr::misc::NumTpGridPoints][dr::misc::NumPaddedPoints]> {};
  struct Sigma : public initializer::Variable<
                     real[seissol::dr::misc::NumTpGridPoints][dr::misc::NumPaddedPoints]> {};
  struct HalfWidthShearZone : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct HydraulicDiffusivity : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  void addTo(DynamicRupture::Tree& tree) {
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<Temperature>(mask, Alignment, allocationModeDR());
    tree.add<Pressure>(mask, Alignment, allocationModeDR());
    tree.add<Theta>(mask, Alignment, allocationModeDR());
    tree.add<Sigma>(mask, Alignment, allocationModeDR());
    tree.add<HalfWidthShearZone>(mask, Alignment, allocationModeDR(), true);
    tree.add<HydraulicDiffusivity>(mask, Alignment, allocationModeDR(), true);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   DynamicRupture::Tree* tree) const {
    manager.registerData<Temperature>("temperature", tree);
    manager.registerData<Pressure>("pressure", tree);
    manager.registerData<Theta>("theta", tree);
    manager.registerData<Sigma>("sigma", tree);
  }
};

struct LTSRateAndStateThermalPressurization : public LTSRateAndState,
                                              public LTSThermalPressurization {
  explicit LTSRateAndStateThermalPressurization(
      const initializer::parameters::DRParameters* parameters)
      : LTSRateAndState(parameters) {}

  void addTo(Tree& tree) override {
    LTSRateAndState::addTo(tree);
    LTSThermalPressurization::addTo(tree);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Tree* tree) const override {
    LTSRateAndState::registerCheckpointVariables(manager, tree);
    LTSThermalPressurization::registerCheckpointVariables(manager, tree);
  }
};

struct LTSRateAndStateThermalPressurizationFastVelocityWeakening
    : public LTSRateAndStateFastVelocityWeakening,
      public LTSThermalPressurization {
  explicit LTSRateAndStateThermalPressurizationFastVelocityWeakening(
      const initializer::parameters::DRParameters* parameters)
      : LTSRateAndStateFastVelocityWeakening(parameters) {}

  void addTo(Tree& tree) override {
    LTSRateAndStateFastVelocityWeakening::addTo(tree);
    LTSThermalPressurization::addTo(tree);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Tree* tree) const override {
    LTSRateAndStateFastVelocityWeakening::registerCheckpointVariables(manager, tree);
    LTSThermalPressurization::registerCheckpointVariables(manager, tree);
  }
};

struct LTSImposedSlipRates : public DynamicRupture {
  struct ImposedSlipDirection1 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct ImposedSlipDirection2 : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct OnsetTime : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  explicit LTSImposedSlipRates(const initializer::parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(Tree& tree) override {
    DynamicRupture::addTo(tree);
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<ImposedSlipDirection1>(mask, Alignment, allocationModeDR(), true);
    tree.add<ImposedSlipDirection2>(mask, Alignment, allocationModeDR(), true);
    tree.add<OnsetTime>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesYoffe : public LTSImposedSlipRates {
  struct TauS : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};
  struct TauR : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  explicit LTSImposedSlipRatesYoffe(const initializer::parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}

  void addTo(Tree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<TauS>(mask, Alignment, allocationModeDR(), true);
    tree.add<TauR>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesGaussian : public LTSImposedSlipRates {
  struct RiseTime : public initializer::Variable<real[dr::misc::NumPaddedPoints]> {};

  explicit LTSImposedSlipRatesGaussian(const initializer::parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}

  void addTo(Tree& tree) override {
    LTSImposedSlipRates::addTo(tree);
    const auto mask = initializer::LayerMask(Ghost);
    tree.add<RiseTime>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesDelta : public LTSImposedSlipRates {
  explicit LTSImposedSlipRatesDelta(const initializer::parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_DYNAMICRUPTURE_H_
