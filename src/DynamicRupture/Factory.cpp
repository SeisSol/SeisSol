// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Factory.h"

#include "FrictionLaws/FrictionLaws.h"
#include "Initializer/Initializers.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Output/Output.h"
#include <memory>
#include <utils/logger.h>

namespace friction_law_cpu = seissol::dr::friction_law::cpu;

// for now, fake the friction laws to be on CPU instead, if we have no GPU
#ifdef ACL_DEVICE
namespace friction_law_gpu = seissol::dr::friction_law::gpu;
#else
namespace friction_law_gpu = seissol::dr::friction_law::cpu;
#endif

namespace {
using namespace seissol::dr;
using namespace seissol::dr::factory;

class NoFaultFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::initializer::DynamicRupture>(),
            std::make_unique<seissol::dr::initializer::NoFaultInitializer>(drParameters,
                                                                           seissolInstance),
            std::make_unique<friction_law_cpu::NoFault>(drParameters.get()),
            std::make_unique<friction_law_gpu::NoFault>(drParameters.get()),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::NoFault>(), seissolInstance)};
  }
};

class LinearSlipWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::initializer::LTSLinearSlipWeakening>(),
            std::make_unique<seissol::dr::initializer::LinearSlipWeakeningInitializer>(
                drParameters, seissolInstance),
            std::make_unique<
                friction_law_cpu::LinearSlipWeakeningLaw<friction_law_cpu::NoSpecialization>>(
                drParameters.get()),
            std::make_unique<
                friction_law_gpu::LinearSlipWeakeningLaw<friction_law_gpu::NoSpecialization>>(
                drParameters.get()),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::LinearSlipWeakening>(), seissolInstance)};
  }
};

class RateAndStateAgingFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters->isThermalPressureOn) {
      return {
          std::make_unique<seissol::initializer::LTSRateAndStateThermalPressurization>(),
          std::make_unique<seissol::dr::initializer::RateAndStateThermalPressurizationInitializer>(
              drParameters, seissolInstance),
          std::make_unique<friction_law_cpu::AgingLaw<friction_law_cpu::ThermalPressurization>>(
              drParameters.get()),
          std::make_unique<friction_law_gpu::AgingLaw<friction_law_gpu::ThermalPressurization>>(
              drParameters.get()),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance)};
    } else {
      return {
          std::make_unique<seissol::initializer::LTSRateAndState>(),
          std::make_unique<seissol::dr::initializer::RateAndStateInitializer>(drParameters,
                                                                              seissolInstance),
          std::make_unique<friction_law_cpu::AgingLaw<friction_law_cpu::NoTP>>(drParameters.get()),
          std::make_unique<friction_law_gpu::AgingLaw<friction_law_gpu::NoTP>>(drParameters.get()),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance)};
    }
  }
};

class RateAndStateSlipFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters->isThermalPressureOn) {
      return {
          std::make_unique<seissol::initializer::LTSRateAndStateThermalPressurization>(),
          std::make_unique<seissol::dr::initializer::RateAndStateThermalPressurizationInitializer>(
              drParameters, seissolInstance),
          std::make_unique<friction_law_cpu::SlipLaw<friction_law_cpu::ThermalPressurization>>(
              drParameters.get()),
          std::make_unique<friction_law_gpu::SlipLaw<friction_law_gpu::ThermalPressurization>>(
              drParameters.get()),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance)};
    } else {
      return {
          std::make_unique<seissol::initializer::LTSRateAndState>(),
          std::make_unique<seissol::dr::initializer::RateAndStateInitializer>(drParameters,
                                                                              seissolInstance),
          std::make_unique<friction_law_cpu::SlipLaw<friction_law_cpu::NoTP>>(drParameters.get()),
          std::make_unique<friction_law_gpu::SlipLaw<friction_law_gpu::NoTP>>(drParameters.get()),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance)};
    }
  }
};

class LinearSlipWeakeningBimaterialFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    using Specialization = friction_law_cpu::BiMaterialFault;
    using FrictionLawType = friction_law_cpu::LinearSlipWeakeningLaw<Specialization>;
    using SpecializationGpu = friction_law_gpu::BiMaterialFault;
    using FrictionLawTypeGpu = friction_law_gpu::LinearSlipWeakeningLaw<SpecializationGpu>;

    return {std::make_unique<seissol::initializer::LTSLinearSlipWeakeningBimaterial>(),
            std::make_unique<seissol::dr::initializer::LinearSlipWeakeningBimaterialInitializer>(
                drParameters, seissolInstance),
            std::make_unique<FrictionLawType>(drParameters.get()),
            std::make_unique<FrictionLawTypeGpu>(drParameters.get()),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::LinearSlipWeakeningBimaterial>(),
                seissolInstance)};
  }
};

class LinearSlipWeakeningTPApproxFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    using Specialization = friction_law_cpu::TPApprox;
    using FrictionLawType = friction_law_cpu::LinearSlipWeakeningLaw<Specialization>;
    using SpecializationGpu = friction_law_gpu::TPApprox;
    using FrictionLawTypeGpu = friction_law_gpu::LinearSlipWeakeningLaw<SpecializationGpu>;

    return {std::make_unique<seissol::initializer::LTSLinearSlipWeakening>(),
            std::make_unique<seissol::dr::initializer::LinearSlipWeakeningInitializer>(
                drParameters, seissolInstance),
            std::make_unique<FrictionLawType>(drParameters.get()),
            std::make_unique<FrictionLawTypeGpu>(drParameters.get()),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::LinearSlipWeakening>(), seissolInstance)};
  }
};

class ImposedSlipRatesYoffeFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::initializer::LTSImposedSlipRatesYoffe>(),
            std::make_unique<seissol::dr::initializer::ImposedSlipRatesYoffeInitializer>(
                drParameters, seissolInstance),
            std::make_unique<friction_law_cpu::ImposedSlipRates<friction_law_cpu::YoffeSTF>>(
                drParameters.get()),
            std::make_unique<friction_law_gpu::ImposedSlipRates<friction_law_gpu::YoffeSTF>>(
                drParameters.get()),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::ImposedSlipRates>(), seissolInstance)};
  }
};

class ImposedSlipRatesGaussianFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::initializer::LTSImposedSlipRatesGaussian>(),
            std::make_unique<seissol::dr::initializer::ImposedSlipRatesGaussianInitializer>(
                drParameters, seissolInstance),
            std::make_unique<friction_law_cpu::ImposedSlipRates<friction_law_cpu::GaussianSTF>>(
                drParameters.get()),
            std::make_unique<friction_law_gpu::ImposedSlipRates<friction_law_gpu::GaussianSTF>>(
                drParameters.get()),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::ImposedSlipRates>(), seissolInstance)};
  }
};

class ImposedSlipRatesDeltaFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::initializer::LTSImposedSlipRatesDelta>(),
            std::make_unique<seissol::dr::initializer::ImposedSlipRatesDeltaInitializer>(
                drParameters, seissolInstance),
            std::make_unique<friction_law_cpu::ImposedSlipRates<friction_law_cpu::DeltaSTF>>(
                drParameters.get()),
            std::make_unique<friction_law_gpu::ImposedSlipRates<friction_law_gpu::DeltaSTF>>(
                drParameters.get()),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::ImposedSlipRates>(), seissolInstance)};
  }
};

class RateAndStateFastVelocityWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters->isThermalPressureOn) {
      return {
          std::make_unique<
              seissol::initializer::LTSRateAndStateThermalPressurizationFastVelocityWeakening>(),
          std::make_unique<
              seissol::dr::initializer::RateAndStateFastVelocityThermalPressurizationInitializer>(
              drParameters, seissolInstance),
          std::make_unique<
              friction_law_cpu::FastVelocityWeakeningLaw<friction_law_cpu::ThermalPressurization>>(
              drParameters.get()),
          std::make_unique<
              friction_law_gpu::FastVelocityWeakeningLaw<friction_law_gpu::ThermalPressurization>>(
              drParameters.get()),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance)};
    } else {
      return {std::make_unique<seissol::initializer::LTSRateAndStateFastVelocityWeakening>(),
              std::make_unique<seissol::dr::initializer::RateAndStateFastVelocityInitializer>(
                  drParameters, seissolInstance),
              std::make_unique<friction_law_cpu::FastVelocityWeakeningLaw<friction_law_cpu::NoTP>>(
                  drParameters.get()),
              std::make_unique<friction_law_gpu::FastVelocityWeakeningLaw<friction_law_gpu::NoTP>>(
                  drParameters.get()),
              std::make_unique<seissol::dr::output::OutputManager>(
                  std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance)};
    }
  }
};

class RateAndStateSevereVelocityWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters->isThermalPressureOn) {
      return {
          std::make_unique<seissol::initializer::LTSRateAndStateThermalPressurization>(),
          std::make_unique<seissol::dr::initializer::RateAndStateThermalPressurizationInitializer>(
              drParameters, seissolInstance),
          std::make_unique<friction_law_cpu::SevereVelocityWeakeningLaw<
              friction_law_cpu::ThermalPressurization>>(drParameters.get()),
          std::make_unique<friction_law_gpu::SevereVelocityWeakeningLaw<
              friction_law_gpu::ThermalPressurization>>(drParameters.get()),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance)};
    } else {
      return {
          std::make_unique<seissol::initializer::LTSRateAndState>(),
          std::make_unique<seissol::dr::initializer::RateAndStateInitializer>(drParameters,
                                                                              seissolInstance),
          std::make_unique<friction_law_cpu::SevereVelocityWeakeningLaw<friction_law_cpu::NoTP>>(
              drParameters.get()),
          std::make_unique<friction_law_gpu::SevereVelocityWeakeningLaw<friction_law_gpu::NoTP>>(
              drParameters.get()),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance)};
    }
  }
};
} // namespace

namespace seissol::dr::factory {
std::unique_ptr<AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
               seissol::SeisSol& seissolInstance) {
  switch (drParameters->frictionLawType) {
  case seissol::initializer::parameters::FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesDelta:
    return std::make_unique<ImposedSlipRatesDeltaFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningTPApprox:
    return std::make_unique<LinearSlipWeakeningTPApproxFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingLaw:
    return std::make_unique<RateAndStateAgingFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateSlipLaw:
    return std::make_unique<RateAndStateSlipFactory>(drParameters, seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateSevereVelocityWeakening:
    return std::make_unique<RateAndStateSevereVelocityWeakeningFactory>(drParameters,
                                                                        seissolInstance);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return {nullptr};
  case seissol::initializer::parameters::FrictionLawType::RateAndStateFastVelocityWeakening:
    return std::make_unique<RateAndStateFastVelocityWeakeningFactory>(drParameters,
                                                                      seissolInstance);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}
} // namespace seissol::dr::factory
