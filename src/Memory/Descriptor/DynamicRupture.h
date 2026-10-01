// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_MEMORY_DESCRIPTOR_DYNAMICRUPTURE_H_
#define SEISSOL_SRC_MEMORY_DESCRIPTOR_DYNAMICRUPTURE_H_

#include "Common/Real.h"
#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Typedefs.h"
#include "GeneratedCode/tensor.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Memory/Tree/Backmap.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"

namespace seissol {

inline auto allocationModeDR() {
  using namespace seissol::initializer;
  if constexpr (!isDeviceOn()) {
    return AllocationMode::HostOnly;
  } else {
    const auto modeMaybeCompress = useDeviceL2Compress() ? AllocationMode::HostDeviceCompress
                                                         : AllocationMode::HostDeviceSplit;
    return useUSM() ? AllocationMode::HostDeviceUnified : modeMaybeCompress;
  }
}

// NOTE: for the sake of GPU performance, make sure that NumPaddedPoints is always last.

struct DynamicRupture {
  public:
  DynamicRupture() = default;
  /// the stress sources of a face: every configured nucleation, then the initial state
  std::size_t stressSourceCount{1};
  explicit DynamicRupture(const initializer::parameters::DRParameters* parameters)
      : stressSourceCount(dr::stressSourceCount(*parameters)) {}

  virtual ~DynamicRupture() = default;

  // The data of a fault face are held in the reals and the layout of the configuration of its
  // layer.
  template <typename Cfg>
  using RealPtr = Real<Cfg>*;
  template <typename Cfg>
  using PointReals = Real<Cfg>[dr::misc::NumPaddedPoints<Cfg>];
  template <typename Cfg>
  using PointFlags = bool[dr::misc::NumPaddedPoints<Cfg>];
  template <typename Cfg>
  using PointStresses = Real<Cfg>[6][dr::misc::NumPaddedPoints<Cfg>];
  template <typename Cfg>
  using QInterpolatedArray = Real<Cfg>[tensor::QInterpolated<Cfg>::size()];
  template <typename Cfg>
  using QInterpolatedTimeArray =
      Real<Cfg>[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()];
  template <typename Cfg>
  using FluxSolverArray = Real<Cfg>[tensor::fluxSolver<Cfg>::size()];

  struct TimeDofsPlus : public initializer::VariantVariable<RealPtr> {};
  struct TimeDofsMinus : public initializer::VariantVariable<RealPtr> {};
  struct TimeDerivativePlus : public initializer::VariantVariable<RealPtr> {};
  struct TimeDerivativeMinus : public initializer::VariantVariable<RealPtr> {};
  struct TimeDerivativePlusDevice : public initializer::VariantVariable<RealPtr> {};
  struct TimeDerivativeMinusDevice : public initializer::VariantVariable<RealPtr> {};
  struct ImposedStatePlus : public initializer::VariantVariable<QInterpolatedArray> {};
  struct ImposedStateMinus : public initializer::VariantVariable<QInterpolatedArray> {};
  struct GodunovData : public initializer::VariantVariable<DRGodunovData> {};
  struct FluxSolverPlus : public initializer::VariantVariable<FluxSolverArray> {};
  struct FluxSolverMinus : public initializer::VariantVariable<FluxSolverArray> {};
  struct FaceInformation : public initializer::Variable<DRFaceInformation> {};
  struct WaveSpeedsPlus : public initializer::Variable<model::IsotropicWaveSpeeds> {};
  struct WaveSpeedsMinus : public initializer::Variable<model::IsotropicWaveSpeeds> {};
  struct DREnergyOutputVar : public initializer::VariantVariable<DREnergyOutput> {};

  struct ImpAndEta : public initializer::VariantVariable<seissol::dr::ImpedancesAndEta> {};
  struct ImpedanceMatrices : public initializer::VariantVariable<seissol::dr::ImpedanceMatrices> {};
  // size padded for vectorization
  // CS = coordinate system
  /// the stress of every source of this face, the initial state first; see dr::stressSourceCount
  struct StressSourceInFaultCS : public initializer::VariantVariable<PointStresses> {};
  // will be always zero, if not using poroelasticity
  struct StressSourcePressure : public initializer::VariantVariable<PointReals> {};
  /// the onset of every source of this face, per point; see dr::stressSourceCount
  struct StressSourceOnset : public initializer::VariantVariable<PointReals> {};
  /// the rise time of every source of this face, per point; see dr::stressSourceCount
  struct StressSourceRiseTime : public initializer::VariantVariable<PointReals> {};
  struct Mu : public initializer::VariantVariable<PointReals> {};
  struct AccumulatedSlipMagnitude : public initializer::VariantVariable<PointReals> {};
  // slip at given fault node along local direction 1
  struct Slip1 : public initializer::VariantVariable<PointReals> {};
  // slip at given fault node along local direction 2
  struct Slip2 : public initializer::VariantVariable<PointReals> {};
  struct SlipRateMagnitude : public initializer::VariantVariable<PointReals> {};
  // slip rate at given fault node along local direction 1
  struct SlipRate1 : public initializer::VariantVariable<PointReals> {};
  // slip rate at given fault node along local direction 2
  struct SlipRate2 : public initializer::VariantVariable<PointReals> {};
  struct RuptureTime : public initializer::VariantVariable<PointReals> {};
  struct DynStressTime : public initializer::VariantVariable<PointReals> {};
  struct RuptureTimePending : public initializer::VariantVariable<PointFlags> {};
  struct DynStressTimePending : public initializer::VariantVariable<PointFlags> {};
  struct PeakSlipRate : public initializer::VariantVariable<PointReals> {};
  struct Traction1 : public initializer::VariantVariable<PointReals> {};
  struct Traction2 : public initializer::VariantVariable<PointReals> {};
  struct QInterpolatedPlus : public initializer::VariantVariable<QInterpolatedTimeArray> {};
  struct QInterpolatedMinus : public initializer::VariantVariable<QInterpolatedTimeArray> {};

  struct IdofsPlusOnDevice : public initializer::VariantScratchpad<Real> {};
  struct IdofsMinusOnDevice : public initializer::VariantScratchpad<Real> {};

  struct DynrupVarmap : public initializer::GenericVarmap {};

  using Storage = initializer::Storage<DynrupVarmap>;
  using Layer = initializer::Layer<DynrupVarmap>;
  template <typename Cfg>
  using Ref = initializer::Layer<DynrupVarmap>::CellRef<Cfg>;
  using Backmap = initializer::StorageBackmap<1>;

  virtual void addTo(Storage& storage) {
    using namespace seissol::initializer;
    const auto mask = LayerMask(Ghost);
    storage.add<TimeDofsPlus>(mask, Alignment, AllocationMode::HostOnly, true);
    storage.add<TimeDofsMinus>(mask, Alignment, AllocationMode::HostOnly, true);
    storage.add<TimeDerivativePlus>(mask, Alignment, AllocationMode::HostOnly, true);
    storage.add<TimeDerivativeMinus>(mask, Alignment, AllocationMode::HostOnly, true);
    storage.add<TimeDerivativePlusDevice>(mask, Alignment, AllocationMode::HostOnly, true);
    storage.add<TimeDerivativeMinusDevice>(mask, Alignment, AllocationMode::HostOnly, true);
    storage.add<ImposedStatePlus>(mask, PagesizeHeap, allocationModeDR());
    storage.add<ImposedStateMinus>(mask, PagesizeHeap, allocationModeDR());
    storage.add<GodunovData>(mask, Alignment, allocationModeDR());
    storage.add<FluxSolverPlus>(mask, Alignment, allocationModeDR());
    storage.add<FluxSolverMinus>(mask, Alignment, allocationModeDR());
    storage.add<FaceInformation>(mask, Alignment, AllocationMode::HostOnly, true);
    storage.add<WaveSpeedsPlus>(mask, Alignment, allocationModeDR(), true);
    storage.add<WaveSpeedsMinus>(mask, Alignment, allocationModeDR(), true);
    storage.add<DREnergyOutputVar>(mask, Alignment, allocationModeDR());
    storage.add<ImpAndEta>(mask, Alignment, allocationModeDR(), true);
    storage.add<ImpedanceMatrices>(mask, Alignment, allocationModeDR(), true);
    storage.add<RuptureTime>(mask, Alignment, allocationModeDR());

    // NOTE: the number of stress sources per face is passed here.
    storage.add<StressSourceInFaultCS>(
        mask, Alignment, allocationModeDR(), true, stressSourceCount);
    storage.add<StressSourcePressure>(mask, Alignment, allocationModeDR(), true, stressSourceCount);
    storage.add<StressSourceOnset>(mask, Alignment, allocationModeDR(), true, stressSourceCount);
    storage.add<StressSourceRiseTime>(mask, Alignment, allocationModeDR(), true, stressSourceCount);

    storage.add<RuptureTimePending>(mask, Alignment, allocationModeDR());
    storage.add<DynStressTime>(mask, Alignment, allocationModeDR());
    storage.add<DynStressTimePending>(mask, Alignment, allocationModeDR());
    storage.add<Mu>(mask, Alignment, allocationModeDR());
    storage.add<AccumulatedSlipMagnitude>(mask, Alignment, allocationModeDR());
    storage.add<Slip1>(mask, Alignment, allocationModeDR());
    storage.add<Slip2>(mask, Alignment, allocationModeDR());
    storage.add<SlipRateMagnitude>(mask, Alignment, allocationModeDR());
    storage.add<SlipRate1>(mask, Alignment, allocationModeDR());
    storage.add<SlipRate2>(mask, Alignment, allocationModeDR());
    storage.add<PeakSlipRate>(mask, Alignment, allocationModeDR());
    storage.add<Traction1>(mask, Alignment, allocationModeDR());
    storage.add<Traction2>(mask, Alignment, allocationModeDR());
    storage.add<QInterpolatedPlus>(mask, Alignment, allocationModeDR());
    storage.add<QInterpolatedMinus>(mask, Alignment, allocationModeDR());

    if constexpr (isDeviceOn()) {
      storage.add<IdofsPlusOnDevice>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
      storage.add<IdofsMinusOnDevice>(LayerMask(), Alignment, AllocationMode::DeviceOnly);
    }
  }

  virtual void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                           Storage& storage) const {
    manager.registerData<Mu>("mu", storage);
    manager.registerData<SlipRate1>("slipRate1", storage);
    manager.registerData<SlipRate2>("slipRate2", storage);
    manager.registerData<AccumulatedSlipMagnitude>("accumulatedSlipMagnitude", storage);
    manager.registerData<Slip1>("slip1", storage);
    manager.registerData<Slip2>("slip2", storage);
    manager.registerData<PeakSlipRate>("peakSlipRate", storage);
    manager.registerData<RuptureTime>("ruptureTime", storage);
    manager.registerData<RuptureTimePending>("ruptureTimePending", storage);
    manager.registerData<DynStressTime>("dynStressTime", storage);
    manager.registerData<DynStressTimePending>("dynStressTimePending", storage);
    manager.registerData<DREnergyOutputVar>("drEnergyOutput", storage);
  }
};

struct LTSLinearSlipWeakening : public DynamicRupture {
  struct DC : public initializer::VariantVariable<PointReals> {};
  struct MuS : public initializer::VariantVariable<PointReals> {};
  struct MuD : public initializer::VariantVariable<PointReals> {};
  struct Cohesion : public initializer::VariantVariable<PointReals> {};
  struct ForcedRuptureTime : public initializer::VariantVariable<PointReals> {};

  explicit LTSLinearSlipWeakening(const initializer::parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(Storage& storage) override {
    DynamicRupture::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<DC>(mask, Alignment, allocationModeDR(), true);
    storage.add<MuS>(mask, Alignment, allocationModeDR(), true);
    storage.add<MuD>(mask, Alignment, allocationModeDR(), true);
    storage.add<Cohesion>(mask, Alignment, allocationModeDR(), true);
    storage.add<ForcedRuptureTime>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSLinearSlipWeakeningBimaterial : public LTSLinearSlipWeakening {
  struct RegularizedStrength : public initializer::VariantVariable<PointReals> {};

  explicit LTSLinearSlipWeakeningBimaterial(const initializer::parameters::DRParameters* parameters)
      : LTSLinearSlipWeakening(parameters) {}

  void addTo(Storage& storage) override {
    LTSLinearSlipWeakening::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<RegularizedStrength>(mask, Alignment, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Storage& storage) const override {
    LTSLinearSlipWeakening::registerCheckpointVariables(manager, storage);
    manager.registerData<RegularizedStrength>("regularizedStrength", storage);
  }
};

struct LTSRateAndState : public DynamicRupture {
  struct RsA : public initializer::VariantVariable<PointReals> {};
  struct RsSl0 : public initializer::VariantVariable<PointReals> {};
  struct StateVariable : public initializer::VariantVariable<PointReals> {};
  struct RsF0 : public initializer::VariantVariable<PointReals> {};
  struct RsMuW : public initializer::VariantVariable<PointReals> {};
  struct RsB : public initializer::VariantVariable<PointReals> {};
  struct ConvergenceInner : public initializer::VariantVariable<PointFlags> {};
  struct ConvergenceOuter : public initializer::VariantVariable<PointFlags> {};

  explicit LTSRateAndState(const initializer::parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(Storage& storage) override {
    DynamicRupture::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<RsA>(mask, Alignment, allocationModeDR(), true);
    storage.add<RsSl0>(mask, Alignment, allocationModeDR(), true);
    storage.add<StateVariable>(mask, Alignment, allocationModeDR());
    storage.add<RsF0>(mask, Alignment, allocationModeDR(), true);
    storage.add<RsMuW>(mask, Alignment, allocationModeDR(), true);
    storage.add<RsB>(mask, Alignment, allocationModeDR(), true);
    storage.add<ConvergenceInner>(mask, Alignment, allocationModeDR());
    storage.add<ConvergenceOuter>(mask, Alignment, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Storage& storage) const override {
    DynamicRupture::registerCheckpointVariables(manager, storage);
    manager.registerData<StateVariable>("stateVariable", storage);
    manager.registerData<ConvergenceInner>("convergenceInner", storage);
    manager.registerData<ConvergenceOuter>("convergenceOuter", storage);
  }
};

struct LTSRateAndStateFastVelocityWeakening : public LTSRateAndState {
  struct RsSrW : public initializer::VariantVariable<PointReals> {};

  explicit LTSRateAndStateFastVelocityWeakening(
      const initializer::parameters::DRParameters* parameters)
      : LTSRateAndState(parameters) {}

  void addTo(Storage& storage) override {
    LTSRateAndState::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<RsSrW>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSThermalPressurization {
  template <typename Cfg>
  using PointReals = DynamicRupture::PointReals<Cfg>;
  template <typename Cfg>
  using GridPointReals = Real<Cfg>[dr::misc::NumTpGridPoints][dr::misc::NumPaddedPoints<Cfg>];

  struct Temperature : public initializer::VariantVariable<PointReals> {};
  struct Pressure : public initializer::VariantVariable<PointReals> {};
  struct Theta : public initializer::VariantVariable<GridPointReals> {};
  struct Sigma : public initializer::VariantVariable<GridPointReals> {};
  struct HalfWidthShearZone : public initializer::VariantVariable<PointReals> {};
  struct HydraulicDiffusivity : public initializer::VariantVariable<PointReals> {};

  void addTo(DynamicRupture::Storage& storage) {
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<Temperature>(mask, Alignment, allocationModeDR());
    storage.add<Pressure>(mask, Alignment, allocationModeDR());
    storage.add<Theta>(mask, Alignment, allocationModeDR());
    storage.add<Sigma>(mask, Alignment, allocationModeDR());
    storage.add<HalfWidthShearZone>(mask, Alignment, allocationModeDR(), true);
    storage.add<HydraulicDiffusivity>(mask, Alignment, allocationModeDR(), true);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   DynamicRupture::Storage& storage) const {
    manager.registerData<Temperature>("temperature", storage);
    manager.registerData<Pressure>("pressure", storage);
    manager.registerData<Theta>("theta", storage);
    manager.registerData<Sigma>("sigma", storage);
  }
};

struct LTSRateAndStateThermalPressurization : public LTSRateAndState,
                                              public LTSThermalPressurization {
  explicit LTSRateAndStateThermalPressurization(
      const initializer::parameters::DRParameters* parameters)
      : LTSRateAndState(parameters) {}

  void addTo(Storage& storage) override {
    LTSRateAndState::addTo(storage);
    LTSThermalPressurization::addTo(storage);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Storage& storage) const override {
    LTSRateAndState::registerCheckpointVariables(manager, storage);
    LTSThermalPressurization::registerCheckpointVariables(manager, storage);
  }
};

struct LTSRateAndStateThermalPressurizationFastVelocityWeakening
    : public LTSRateAndStateFastVelocityWeakening,
      public LTSThermalPressurization {
  explicit LTSRateAndStateThermalPressurizationFastVelocityWeakening(
      const initializer::parameters::DRParameters* parameters)
      : LTSRateAndStateFastVelocityWeakening(parameters) {}

  void addTo(Storage& storage) override {
    LTSRateAndStateFastVelocityWeakening::addTo(storage);
    LTSThermalPressurization::addTo(storage);
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Storage& storage) const override {
    LTSRateAndStateFastVelocityWeakening::registerCheckpointVariables(manager, storage);
    LTSThermalPressurization::registerCheckpointVariables(manager, storage);
  }
};

struct LTSImposedSlipRates : public DynamicRupture {
  struct ImposedSlipDirection1 : public initializer::VariantVariable<PointReals> {};
  struct ImposedSlipDirection2 : public initializer::VariantVariable<PointReals> {};
  struct OnsetTime : public initializer::VariantVariable<PointReals> {};

  explicit LTSImposedSlipRates(const initializer::parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(Storage& storage) override {
    DynamicRupture::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<ImposedSlipDirection1>(mask, Alignment, allocationModeDR(), true);
    storage.add<ImposedSlipDirection2>(mask, Alignment, allocationModeDR(), true);
    storage.add<OnsetTime>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesYoffe : public LTSImposedSlipRates {
  struct TauS : public initializer::VariantVariable<PointReals> {};
  struct TauR : public initializer::VariantVariable<PointReals> {};

  explicit LTSImposedSlipRatesYoffe(const initializer::parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}

  void addTo(Storage& storage) override {
    LTSImposedSlipRates::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<TauS>(mask, Alignment, allocationModeDR(), true);
    storage.add<TauR>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesGaussian : public LTSImposedSlipRates {
  struct RiseTime : public initializer::VariantVariable<PointReals> {};

  explicit LTSImposedSlipRatesGaussian(const initializer::parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}

  void addTo(Storage& storage) override {
    LTSImposedSlipRates::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<RiseTime>(mask, Alignment, allocationModeDR(), true);
  }
};

struct LTSImposedSlipRatesDelta : public LTSImposedSlipRates {
  explicit LTSImposedSlipRatesDelta(const initializer::parameters::DRParameters* parameters)
      : LTSImposedSlipRates(parameters) {}
};

} // namespace seissol

#endif // SEISSOL_SRC_MEMORY_DESCRIPTOR_DYNAMICRUPTURE_H_
