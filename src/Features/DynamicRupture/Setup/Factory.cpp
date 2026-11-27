// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Factory.h"

#include "Core/Data/DataReader.h"
#include "Features/DynamicRupture/DRParameters.h"
#include "Features/DynamicRupture/Misc.h"
#include "Features/DynamicRupture/Output/Output.h"
#include "Features/DynamicRupture/Output/ReceiverBasedOutput.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionLaws.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"
#include "Initializer/Initializers.h"

#include <memory>
#include <utility>
#include <utils/logger.h>

namespace friction_law_cpu = seissol::dr::friction_law::cpu;

// for now, fake the friction laws to be on CPU instead, if we have no GPU
#ifdef ACL_DEVICE
namespace friction_law_gpu = seissol::dr::friction_law::gpu;
#else
namespace friction_law_gpu = seissol::dr::friction_law::cpu;
#endif

namespace seissol::dr::factory {

namespace {

/// The friction law `LawT` specialized by `TailT` (its thermal pressurization, its specialization
/// or its source time function), both for the configuration `Cfg`.
template <template <typename, typename> typename LawT, template <typename> typename TailT>
struct Bound {
  template <typename Cfg>
  using Type = LawT<Cfg, TailT<Cfg>>;
};

/**
 * The parts of a friction law: the storage `StorageT` of its faces, its initializer
 * `InitializerT`, its friction solvers `LawT` on the host and `LawDeviceT` on the device (both
 * templates of the configuration), and its output `OutputT`.
 */
template <typename StorageT,
          typename InitializerT,
          template <typename> typename LawT,
          template <typename> typename LawDeviceT,
          typename OutputT>
class LawFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;

  FrictionLawSolver solver() override {
    return {.storage = std::make_unique<StorageT>(drParameters_.get()),
            .frictionLaw = friction_law::makeFrictionSolverFactory<LawT>(*drParameters_),
            .frictionLawDevice =
                friction_law::makeFrictionSolverFactory<LawDeviceT>(*drParameters_)};
  }

  std::unique_ptr<initializer::BaseDRInitializer> initializer() override {
    return makeInitializer<InitializerT>();
  }

  std::unique_ptr<output::ReceiverOutput> output() override { return std::make_unique<OutputT>(); }
};

using NoFaultFactory = LawFactory<seissol::DynamicRupture,
                                  initializer::NoFaultInitializer,
                                  friction_law_cpu::NoFault,
                                  friction_law_gpu::NoFault,
                                  output::NoFault>;

using LinearSlipWeakeningFactory = LawFactory<
    seissol::LTSLinearSlipWeakening,
    initializer::LinearSlipWeakeningInitializer,
    Bound<friction_law_cpu::LinearSlipWeakeningLaw, friction_law_cpu::NoSpecialization>::Type,
    Bound<friction_law_gpu::LinearSlipWeakeningLaw, friction_law_gpu::NoSpecialization>::Type,
    output::LinearSlipWeakening>;

using LinearSlipWeakeningBimaterialFactory = LawFactory<
    seissol::LTSLinearSlipWeakeningBimaterial,
    initializer::LinearSlipWeakeningBimaterialInitializer,
    Bound<friction_law_cpu::LinearSlipWeakeningLaw, friction_law_cpu::BiMaterialFault>::Type,
    Bound<friction_law_gpu::LinearSlipWeakeningLaw, friction_law_gpu::BiMaterialFault>::Type,
    output::LinearSlipWeakeningBimaterial>;

using LinearSlipWeakeningTPApproxFactory =
    LawFactory<seissol::LTSLinearSlipWeakening,
               initializer::LinearSlipWeakeningInitializer,
               Bound<friction_law_cpu::LinearSlipWeakeningLaw, friction_law_cpu::TPApprox>::Type,
               Bound<friction_law_gpu::LinearSlipWeakeningLaw, friction_law_gpu::TPApprox>::Type,
               output::LinearSlipWeakening>;

using RateAndStateAgingFactory =
    LawFactory<seissol::LTSRateAndState,
               initializer::RateAndStateInitializer,
               Bound<friction_law_cpu::AgingLaw, friction_law_cpu::NoTP>::Type,
               Bound<friction_law_gpu::AgingLaw, friction_law_gpu::NoTP>::Type,
               output::RateAndState>;

using RateAndStateAgingTPFactory =
    LawFactory<seissol::LTSRateAndStateThermalPressurization,
               initializer::RateAndStateThermalPressurizationInitializer,
               Bound<friction_law_cpu::AgingLaw, friction_law_cpu::ThermalPressurization>::Type,
               Bound<friction_law_gpu::AgingLaw, friction_law_gpu::ThermalPressurization>::Type,
               output::RateAndStateThermalPressurization>;

using RateAndStateSlipFactory =
    LawFactory<seissol::LTSRateAndState,
               initializer::RateAndStateInitializer,
               Bound<friction_law_cpu::SlipLaw, friction_law_cpu::NoTP>::Type,
               Bound<friction_law_gpu::SlipLaw, friction_law_gpu::NoTP>::Type,
               output::RateAndState>;

using RateAndStateSlipTPFactory =
    LawFactory<seissol::LTSRateAndStateThermalPressurization,
               initializer::RateAndStateThermalPressurizationInitializer,
               Bound<friction_law_cpu::SlipLaw, friction_law_cpu::ThermalPressurization>::Type,
               Bound<friction_law_gpu::SlipLaw, friction_law_gpu::ThermalPressurization>::Type,
               output::RateAndStateThermalPressurization>;

using RateAndStateFastVelocityWeakeningFactory =
    LawFactory<seissol::LTSRateAndStateFastVelocityWeakening,
               initializer::RateAndStateFastVelocityInitializer,
               Bound<friction_law_cpu::FastVelocityWeakeningLaw, friction_law_cpu::NoTP>::Type,
               Bound<friction_law_gpu::FastVelocityWeakeningLaw, friction_law_gpu::NoTP>::Type,
               output::RateAndState>;

using RateAndStateFastVelocityWeakeningTPFactory =
    LawFactory<seissol::LTSRateAndStateThermalPressurizationFastVelocityWeakening,
               initializer::RateAndStateFastVelocityThermalPressurizationInitializer,
               Bound<friction_law_cpu::FastVelocityWeakeningLaw,
                     friction_law_cpu::ThermalPressurization>::Type,
               Bound<friction_law_gpu::FastVelocityWeakeningLaw,
                     friction_law_gpu::ThermalPressurization>::Type,
               output::RateAndStateThermalPressurization>;

using RateAndStateSevereVelocityWeakeningFactory =
    LawFactory<seissol::LTSRateAndState,
               initializer::RateAndStateSevereVelocityInitializer,
               Bound<friction_law_cpu::SevereVelocityWeakeningLaw, friction_law_cpu::NoTP>::Type,
               Bound<friction_law_gpu::SevereVelocityWeakeningLaw, friction_law_gpu::NoTP>::Type,
               output::RateAndState>;

using RateAndStateSevereVelocityWeakeningTPFactory =
    LawFactory<seissol::LTSRateAndStateThermalPressurization,
               initializer::RateAndStateSevereVelocityThermalPressurizationInitializer,
               Bound<friction_law_cpu::SevereVelocityWeakeningLaw,
                     friction_law_cpu::ThermalPressurization>::Type,
               Bound<friction_law_gpu::SevereVelocityWeakeningLaw,
                     friction_law_gpu::ThermalPressurization>::Type,
               output::RateAndStateThermalPressurization>;

using ImposedSlipRatesYoffeFactory =
    LawFactory<seissol::LTSImposedSlipRatesYoffe,
               initializer::ImposedSlipRatesYoffeInitializer,
               Bound<friction_law_cpu::ImposedSlipRates, friction_law_cpu::YoffeSTF>::Type,
               Bound<friction_law_gpu::ImposedSlipRates, friction_law_gpu::YoffeSTF>::Type,
               output::ImposedSlipRates>;

using ImposedSlipRatesGaussianFactory =
    LawFactory<seissol::LTSImposedSlipRatesGaussian,
               initializer::ImposedSlipRatesGaussianInitializer,
               Bound<friction_law_cpu::ImposedSlipRates, friction_law_cpu::GaussianSTF>::Type,
               Bound<friction_law_gpu::ImposedSlipRates, friction_law_gpu::GaussianSTF>::Type,
               output::ImposedSlipRates>;

using ImposedSlipRatesDeltaFactory =
    LawFactory<seissol::LTSImposedSlipRatesDelta,
               initializer::ImposedSlipRatesDeltaInitializer,
               Bound<friction_law_cpu::ImposedSlipRates, friction_law_cpu::DeltaSTF>::Type,
               Bound<friction_law_gpu::ImposedSlipRates, friction_law_gpu::DeltaSTF>::Type,
               output::ImposedSlipRates>;

/// The factory `WithoutTP`, or `WithTP` if the parameters turn on thermal pressurization.
template <typename WithoutTP, typename WithTP, typename... Args>
std::unique_ptr<AbstractFactory> makeTPFactory(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
    const Args&... args) {
  if (drParameters->isThermalPressureOn) {
    return std::make_unique<WithTP>(drParameters, args...);
  }
  return std::make_unique<WithoutTP>(drParameters, args...);
}

/// The factory of the friction law of `drParameters`, made with `args` after the parameters.
template <typename... Args>
std::unique_ptr<AbstractFactory>
    makeFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
                const Args&... args) {
  switch (drParameters->frictionLawType) {
  case seissol::dr::misc::FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesDelta:
    return std::make_unique<ImposedSlipRatesDeltaFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakeningTPApprox:
    return std::make_unique<LinearSlipWeakeningTPApproxFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::RateAndStateAgingLaw:
    return makeTPFactory<RateAndStateAgingFactory, RateAndStateAgingTPFactory>(drParameters,
                                                                               args...);
  case seissol::dr::misc::FrictionLawType::RateAndStateSlipLaw:
    return makeTPFactory<RateAndStateSlipFactory, RateAndStateSlipTPFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::RateAndStateSevereVelocityWeakening:
    return makeTPFactory<RateAndStateSevereVelocityWeakeningFactory,
                         RateAndStateSevereVelocityWeakeningTPFactory>(drParameters, args...);
  case seissol::dr::misc::FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return {nullptr};
  case seissol::dr::misc::FrictionLawType::RateAndStateFastVelocityWeakening:
    return makeTPFactory<RateAndStateFastVelocityWeakeningFactory,
                         RateAndStateFastVelocityWeakeningTPFactory>(drParameters, args...);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

} // namespace

DynamicRuptureTuple AbstractFactory::produce() {
  auto parts = solver();
  return {.storage = std::move(parts.storage),
          .initializer = initializer(),
          .frictionLaw = std::move(parts.frictionLaw),
          .frictionLawDevice = std::move(parts.frictionLawDevice),
          .output = output()};
}

std::unique_ptr<AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
               const initializer::FaultMesh& faultMesh,
               const data::ModelLoader& models) {
  return makeFactory(drParameters, faultMesh, models);
}

std::unique_ptr<AbstractFactory> getFactory(
    const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters) {
  return makeFactory(drParameters);
}

} // namespace seissol::dr::factory
