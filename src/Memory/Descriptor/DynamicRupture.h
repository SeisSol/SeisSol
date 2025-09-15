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
#include "GeneratedCode/tensor.h"
#include "IO/Instance/Checkpoint/CheckpointManager.h"
#include "Initializer/Typedefs.h"
#include "Memory/Tree/LTSTree.h"
#include "Memory/Tree/Layer.h"
#include "Parallel/Helper.h"
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

// NOTE: for the sake of GPU performance, make sure that NumPaddedPoints<Cfg> is always last.

template <typename Cfg>
using DrDataArray = Real<Cfg>[dr::misc::NumPaddedPoints<Cfg>];

template <typename Cfg>
using RealPtr = Real<Cfg>*;

struct DynamicRupture {
  public:
  DynamicRupture() = default;
  std::size_t nucleationCount{1};
  explicit DynamicRupture(const initializer::parameters::DRParameters* parameters)
      : nucleationCount(parameters->nucleationCount) {}

  virtual ~DynamicRupture() = default;
  struct TimeDerivativePlus : public initializer::VariantVariable<RealPtr> {};
  struct TimeDerivativeMinus : public initializer::VariantVariable<RealPtr> {};
  struct TimeDerivativePlusDevice : public initializer::VariantVariable<RealPtr> {};
  struct TimeDerivativeMinusDevice : public initializer::VariantVariable<RealPtr> {};
  struct ImposedStatePlus : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[tensor::QInterpolated<Cfg>::size()];
  };
  struct ImposedStateMinus : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[tensor::QInterpolated<Cfg>::size()];
  };
  struct GodunovData : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = DRGodunovData<Cfg>;
  };
  struct FluxSolverPlus : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[tensor::fluxSolver<Cfg>::size()];
  };
  struct FluxSolverMinus : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[tensor::fluxSolver<Cfg>::size()];
  };
  struct FaceInformation : public initializer::Variable<DRFaceInformation> {};
  struct WaveSpeedsPlus : public initializer::Variable<model::IsotropicWaveSpeeds> {};
  struct WaveSpeedsMinus : public initializer::Variable<model::IsotropicWaveSpeeds> {};
  struct DREnergyOutputVar : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = DREnergyOutput<Cfg>;
  };

  struct ImpAndEta : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = seissol::dr::ImpedancesAndEta<Real<Cfg>>;
  };
  struct ImpedanceMatrices : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = seissol::dr::ImpedanceMatrices<Cfg>;
  };
  // size padded for vectorization
  // CS = coordinate system
  struct InitialStressInFaultCS : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[6][dr::misc::NumPaddedPoints<Cfg>];
  };
  struct NucleationStressInFaultCS : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[6][dr::misc::NumPaddedPoints<Cfg>];
  };
  // will be always zero, if not using poroelasticity
  struct InitialPressure : public initializer::VariantVariable<DrDataArray> {};
  struct NucleationPressure : public initializer::VariantVariable<DrDataArray> {};
  struct Mu : public initializer::VariantVariable<DrDataArray> {};
  struct AccumulatedSlipMagnitude : public initializer::VariantVariable<DrDataArray> {};
  // slip at given fault node along local direction 1
  struct Slip1 : public initializer::VariantVariable<DrDataArray> {};
  // slip at given fault node along local direction 2
  struct Slip2 : public initializer::VariantVariable<DrDataArray> {};
  struct SlipRateMagnitude : public initializer::VariantVariable<DrDataArray> {};
  // slip rate at given fault node along local direction 1
  struct SlipRate1 : public initializer::VariantVariable<DrDataArray> {};
  // slip rate at given fault node along local direction 2
  struct SlipRate2 : public initializer::VariantVariable<DrDataArray> {};
  struct RuptureTime : public initializer::VariantVariable<DrDataArray> {};
  struct DynStressTime : public initializer::VariantVariable<DrDataArray> {};
  struct RuptureTimePending : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = bool[dr::misc::NumPaddedPoints<Cfg>];
  };
  struct DynStressTimePending : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = bool[dr::misc::NumPaddedPoints<Cfg>];
  };
  struct PeakSlipRate : public initializer::VariantVariable<DrDataArray> {};
  struct Traction1 : public initializer::VariantVariable<DrDataArray> {};
  struct Traction2 : public initializer::VariantVariable<DrDataArray> {};
  struct QInterpolatedPlus : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()];
  };
  struct QInterpolatedMinus : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>[Cfg::ConvergenceOrder][tensor::QInterpolated<Cfg>::size()];
  };

  struct IdofsPlusOnDevice : public initializer::Scratchpad<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>;
  };
  struct IdofsMinusOnDevice : public initializer::Scratchpad<void> {
    template <typename Cfg>
    using VariantType = Real<Cfg>;
  };

  struct DynrupVarmap : public initializer::GenericVarmap {};

  using Storage = initializer::Storage<DynrupVarmap>;
  using Layer = initializer::Layer<DynrupVarmap>;
  template <typename Config>
  using Ref = initializer::Layer<DynrupVarmap>::CellRef<Config>;
  using Backmap = initializer::StorageBackmap<1>;

  virtual void addTo(Storage& storage) {
    using namespace seissol::initializer;
    const auto mask = LayerMask(Ghost);
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
    storage.add<InitialStressInFaultCS>(mask, Alignment, allocationModeDR());
    storage.add<InitialPressure>(mask, Alignment, allocationModeDR());
    storage.add<RuptureTime>(mask, Alignment, allocationModeDR());

    // NOTE: nucleation count (multi-nucleation support) is passed here.
    storage.add<NucleationStressInFaultCS>(
        mask, Alignment, allocationModeDR(), true, nucleationCount);
    storage.add<NucleationPressure>(mask, Alignment, allocationModeDR(), true, nucleationCount);

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
    manager.registerData<InitialStressInFaultCS>("initialStressInFaultCS", storage);
    manager.registerData<InitialPressure>("initialPressure", storage);
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
  struct DC : public initializer::VariantVariable<DrDataArray> {};
  struct MuS : public initializer::VariantVariable<DrDataArray> {};
  struct MuD : public initializer::VariantVariable<DrDataArray> {};
  struct Cohesion : public initializer::VariantVariable<DrDataArray> {};
  struct ForcedRuptureTime : public initializer::VariantVariable<DrDataArray> {};

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
  struct RegularizedStrength : public initializer::VariantVariable<DrDataArray> {};

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
  struct RsA : public initializer::VariantVariable<DrDataArray> {};
  struct RsSl0 : public initializer::VariantVariable<DrDataArray> {};
  struct StateVariable : public initializer::VariantVariable<DrDataArray> {};

  explicit LTSRateAndState(const initializer::parameters::DRParameters* parameters)
      : DynamicRupture(parameters) {}

  void addTo(Storage& storage) override {
    DynamicRupture::addTo(storage);
    const auto mask = initializer::LayerMask(Ghost);
    storage.add<RsA>(mask, Alignment, allocationModeDR(), true);
    storage.add<RsSl0>(mask, Alignment, allocationModeDR(), true);
    storage.add<StateVariable>(mask, Alignment, allocationModeDR());
  }

  void registerCheckpointVariables(io::instance::checkpoint::CheckpointManager& manager,
                                   Storage& storage) const override {
    DynamicRupture::registerCheckpointVariables(manager, storage);
    manager.registerData<StateVariable>("stateVariable", storage);
  }
};

struct LTSRateAndStateFastVelocityWeakening : public LTSRateAndState {
  struct RsSrW : public initializer::VariantVariable<DrDataArray> {};

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
  struct Temperature : public initializer::VariantVariable<DrDataArray> {};
  struct Pressure : public initializer::VariantVariable<DrDataArray> {};
  struct Theta : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType =
        Real<Cfg>[seissol::dr::misc::NumTpGridPoints][dr::misc::NumPaddedPoints<Cfg>];
  };
  struct Sigma : public initializer::Variable<void> {
    template <typename Cfg>
    using VariantType =
        Real<Cfg>[seissol::dr::misc::NumTpGridPoints][dr::misc::NumPaddedPoints<Cfg>];
  };
  struct HalfWidthShearZone : public initializer::VariantVariable<DrDataArray> {};
  struct HydraulicDiffusivity : public initializer::VariantVariable<DrDataArray> {};

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
  struct ImposedSlipDirection1 : public initializer::VariantVariable<DrDataArray> {};
  struct ImposedSlipDirection2 : public initializer::VariantVariable<DrDataArray> {};
  struct OnsetTime : public initializer::VariantVariable<DrDataArray> {};

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
  struct TauS : public initializer::VariantVariable<DrDataArray> {};
  struct TauR : public initializer::VariantVariable<DrDataArray> {};

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
  struct RiseTime : public initializer::VariantVariable<DrDataArray> {};

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
