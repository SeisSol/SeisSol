// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Factory.h"

#include "DynamicRupture/Misc.h"
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

namespace seissol::dr::factory {

namespace {

class NoFaultFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::DynamicRupture>(drParameters_.get()),
            std::make_unique<seissol::dr::initializer::NoFaultInitializer>(drParameters_,
                                                                           seissolInstance_),
            std::make_unique<friction_law_cpu::NoFault>(FrictionLawParameters(*drParameters_)),
            std::make_unique<friction_law_gpu::NoFault>(FrictionLawParameters(*drParameters_)),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::NoFault>(), seissolInstance_)};
  }
};

class LinearSlipWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::LTSLinearSlipWeakening>(drParameters_.get()),
            std::make_unique<seissol::dr::initializer::LinearSlipWeakeningInitializer>(
                drParameters_, seissolInstance_),
            std::make_unique<
                friction_law_cpu::LinearSlipWeakeningLaw<friction_law_cpu::NoSpecialization>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<
                friction_law_gpu::LinearSlipWeakeningLaw<friction_law_gpu::NoSpecialization>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::LinearSlipWeakening>(), seissolInstance_)};
  }
};

class RateAndStateAgingFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters_->isThermalPressureOn) {
      return {
          std::make_unique<seissol::LTSRateAndStateThermalPressurization>(drParameters_.get()),
          std::make_unique<seissol::dr::initializer::RateAndStateThermalPressurizationInitializer>(
              drParameters_, seissolInstance_),
          std::make_unique<friction_law_cpu::AgingLaw<friction_law_cpu::ThermalPressurization>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<friction_law_gpu::AgingLaw<friction_law_gpu::ThermalPressurization>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance_)};
    } else {
      return {std::make_unique<seissol::LTSRateAndState>(drParameters_.get()),
              std::make_unique<seissol::dr::initializer::RateAndStateInitializer>(drParameters_,
                                                                                  seissolInstance_),
              std::make_unique<friction_law_cpu::AgingLaw<friction_law_cpu::NoTP>>(
                  FrictionLawParameters(*drParameters_)),
              std::make_unique<friction_law_gpu::AgingLaw<friction_law_gpu::NoTP>>(
                  FrictionLawParameters(*drParameters_)),
              std::make_unique<seissol::dr::output::OutputManager>(
                  std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance_)};
    }
  }
};

class RateAndStateSlipFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters_->isThermalPressureOn) {
      return {
          std::make_unique<seissol::LTSRateAndStateThermalPressurization>(drParameters_.get()),
          std::make_unique<seissol::dr::initializer::RateAndStateThermalPressurizationInitializer>(
              drParameters_, seissolInstance_),
          std::make_unique<friction_law_cpu::SlipLaw<friction_law_cpu::ThermalPressurization>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<friction_law_gpu::SlipLaw<friction_law_gpu::ThermalPressurization>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance_)};
    } else {
      return {std::make_unique<seissol::LTSRateAndState>(drParameters_.get()),
              std::make_unique<seissol::dr::initializer::RateAndStateInitializer>(drParameters_,
                                                                                  seissolInstance_),
              std::make_unique<friction_law_cpu::SlipLaw<friction_law_cpu::NoTP>>(
                  FrictionLawParameters(*drParameters_)),
              std::make_unique<friction_law_gpu::SlipLaw<friction_law_gpu::NoTP>>(
                  FrictionLawParameters(*drParameters_)),
              std::make_unique<seissol::dr::output::OutputManager>(
                  std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance_)};
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

    return {std::make_unique<seissol::LTSLinearSlipWeakeningBimaterial>(drParameters_.get()),
            std::make_unique<seissol::dr::initializer::LinearSlipWeakeningBimaterialInitializer>(
                drParameters_, seissolInstance_),
            std::make_unique<FrictionLawType>(FrictionLawParameters(*drParameters_)),
            std::make_unique<FrictionLawTypeGpu>(FrictionLawParameters(*drParameters_)),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::LinearSlipWeakeningBimaterial>(),
                seissolInstance_)};
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

    return {std::make_unique<seissol::LTSLinearSlipWeakening>(drParameters_.get()),
            std::make_unique<seissol::dr::initializer::LinearSlipWeakeningInitializer>(
                drParameters_, seissolInstance_),
            std::make_unique<FrictionLawType>(FrictionLawParameters(*drParameters_)),
            std::make_unique<FrictionLawTypeGpu>(FrictionLawParameters(*drParameters_)),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::LinearSlipWeakening>(), seissolInstance_)};
  }
};

class ImposedSlipRatesYoffeFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::LTSImposedSlipRatesYoffe>(drParameters_.get()),
            std::make_unique<seissol::dr::initializer::ImposedSlipRatesYoffeInitializer>(
                drParameters_, seissolInstance_),
            std::make_unique<friction_law_cpu::ImposedSlipRates<friction_law_cpu::YoffeSTF>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<friction_law_gpu::ImposedSlipRates<friction_law_gpu::YoffeSTF>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::ImposedSlipRates>(), seissolInstance_)};
  }
};

class ImposedSlipRatesGaussianFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::LTSImposedSlipRatesGaussian>(drParameters_.get()),
            std::make_unique<seissol::dr::initializer::ImposedSlipRatesGaussianInitializer>(
                drParameters_, seissolInstance_),
            std::make_unique<friction_law_cpu::ImposedSlipRates<friction_law_cpu::GaussianSTF>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<friction_law_gpu::ImposedSlipRates<friction_law_gpu::GaussianSTF>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::ImposedSlipRates>(), seissolInstance_)};
  }
};

class ImposedSlipRatesDeltaFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    return {std::make_unique<seissol::LTSImposedSlipRatesDelta>(drParameters_.get()),
            std::make_unique<seissol::dr::initializer::ImposedSlipRatesDeltaInitializer>(
                drParameters_, seissolInstance_),
            std::make_unique<friction_law_cpu::ImposedSlipRates<friction_law_cpu::DeltaSTF>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<friction_law_gpu::ImposedSlipRates<friction_law_gpu::DeltaSTF>>(
                FrictionLawParameters(*drParameters_)),
            std::make_unique<seissol::dr::output::OutputManager>(
                std::make_unique<seissol::dr::output::ImposedSlipRates>(), seissolInstance_)};
  }
};

class RateAndStateFastVelocityWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters_->isThermalPressureOn) {
      return {
          std::make_unique<seissol::LTSRateAndStateThermalPressurizationFastVelocityWeakening>(
              drParameters_.get()),
          std::make_unique<
              seissol::dr::initializer::RateAndStateFastVelocityThermalPressurizationInitializer>(
              drParameters_, seissolInstance_),
          std::make_unique<
              friction_law_cpu::FastVelocityWeakeningLaw<friction_law_cpu::ThermalPressurization>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<
              friction_law_gpu::FastVelocityWeakeningLaw<friction_law_gpu::ThermalPressurization>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance_)};
    } else {
      return {std::make_unique<seissol::LTSRateAndStateFastVelocityWeakening>(drParameters_.get()),
              std::make_unique<seissol::dr::initializer::RateAndStateFastVelocityInitializer>(
                  drParameters_, seissolInstance_),
              std::make_unique<friction_law_cpu::FastVelocityWeakeningLaw<friction_law_cpu::NoTP>>(
                  FrictionLawParameters(*drParameters_)),
              std::make_unique<friction_law_gpu::FastVelocityWeakeningLaw<friction_law_gpu::NoTP>>(
                  FrictionLawParameters(*drParameters_)),
              std::make_unique<seissol::dr::output::OutputManager>(
                  std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance_)};
    }
  }
};

class RateAndStateSevereVelocityWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override {
    if (drParameters_->isThermalPressureOn) {
      return {
          std::make_unique<seissol::LTSRateAndStateThermalPressurization>(drParameters_.get()),
          std::make_unique<seissol::dr::initializer::RateAndStateThermalPressurizationInitializer>(
              drParameters_, seissolInstance_),
          std::make_unique<friction_law_cpu::SevereVelocityWeakeningLaw<
              friction_law_cpu::ThermalPressurization>>(FrictionLawParameters(*drParameters_)),
          std::make_unique<friction_law_gpu::SevereVelocityWeakeningLaw<
              friction_law_gpu::ThermalPressurization>>(FrictionLawParameters(*drParameters_)),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndStateThermalPressurization>(),
              seissolInstance_)};
    } else {
      return {
          std::make_unique<seissol::LTSRateAndState>(drParameters_.get()),
          std::make_unique<seissol::dr::initializer::RateAndStateInitializer>(drParameters_,
                                                                              seissolInstance_),
          std::make_unique<friction_law_cpu::SevereVelocityWeakeningLaw<friction_law_cpu::NoTP>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<friction_law_gpu::SevereVelocityWeakeningLaw<friction_law_gpu::NoTP>>(
              FrictionLawParameters(*drParameters_)),
          std::make_unique<seissol::dr::output::OutputManager>(
              std::make_unique<seissol::dr::output::RateAndState>(), seissolInstance_)};
    }
  }
};
} // namespace

std::unique_ptr<AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
               seissol::SeisSol& seissolInstance) {
  switch (drParameters->frictionLawType) {
  case seissol::dr::misc::FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::ImposedSlipRatesDelta:
    return std::make_unique<ImposedSlipRatesDeltaFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::LinearSlipWeakeningTPApprox:
    return std::make_unique<LinearSlipWeakeningTPApproxFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::RateAndStateAgingLaw:
    return std::make_unique<RateAndStateAgingFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::RateAndStateSlipLaw:
    return std::make_unique<RateAndStateSlipFactory>(drParameters, seissolInstance);
  case seissol::dr::misc::FrictionLawType::RateAndStateSevereVelocityWeakening:
    return std::make_unique<RateAndStateSevereVelocityWeakeningFactory>(drParameters,
                                                                        seissolInstance);
  case seissol::dr::misc::FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return {nullptr};
  case seissol::dr::misc::FrictionLawType::RateAndStateFastVelocityWeakening:
    return std::make_unique<RateAndStateFastVelocityWeakeningFactory>(drParameters,
                                                                      seissolInstance);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}
} // namespace seissol::dr::factory
