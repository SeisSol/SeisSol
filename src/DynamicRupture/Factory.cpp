#include "Factory.h"

#include "FrictionLaws/FrictionLaws.h"
#include "FrictionLaws/ThermalPressurization/ThermalPressurization.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/Initializers.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Output/Output.hpp"
#include <memory>
#include <utils/logger.h>

#ifdef ACL_DEVICE
namespace friction_law_impl = seissol::dr::friction_law::gpu;
#else
namespace friction_law_impl = seissol::dr::friction_law;
#endif

namespace seissol::dr::factory {
std::unique_ptr<AbstractFactory>
    getFactory(std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters,
               seissol::SeisSol& seissolInstance,
               unsigned int numFused) {
  switch (drParameters->frictionLawType) {
  case seissol::initializer::parameters::FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory>(drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory>(drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory>(
        drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory>(drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory>(
        drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::LinearSlipWeakeningTPApprox:
    return std::make_unique<LinearSlipWeakeningTPApproxFactory>(
        drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingLaw:
    return std::make_unique<RateAndStateAgingFactory>(drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateSlipLaw:
    return std::make_unique<RateAndStateSlipFactory>(drParameters, seissolInstance, numFused);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateVelocityWeakening:
    logError() << "friction law 7 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return std::unique_ptr<AbstractFactory>(nullptr);
  case seissol::initializer::parameters::FrictionLawType::RateAndStateFastVelocityWeakening:
    return std::make_unique<RateAndStateFastVelocityWeakeningFactory>(
        drParameters, seissolInstance, numFused);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

DynamicRuptureTuple NoFaultFactory::produce() {
  return {std::make_unique<seissol::initializer::DynamicRupture>(),
          std::make_unique<initializer::NoFaultInitializer>(drParameters, seissolInstance),
          std::make_unique<friction_law::NoFault>(drParameters.get()),
          std::make_unique<output::OutputManager>(
              std::make_unique<output::NoFault>(), seissolInstance, numFused)};
}

DynamicRuptureTuple LinearSlipWeakeningFactory::produce() {
  return {
      std::make_unique<seissol::initializer::LTSLinearSlipWeakening>(),
      std::make_unique<initializer::LinearSlipWeakeningInitializer>(drParameters, seissolInstance),
      std::make_unique<
          friction_law_impl::LinearSlipWeakeningLaw<friction_law_impl::NoSpecialization>>(
          drParameters.get()),
      std::make_unique<output::OutputManager>(
          std::make_unique<output::LinearSlipWeakening>(), seissolInstance, numFused)};
}

DynamicRuptureTuple RateAndStateAgingFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSRateAndState>(),
            std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
            std::make_unique<friction_law::AgingLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(),
                seissolInstance,
                numFused)};
  } else {
    return {
        std::make_unique<seissol::initializer::LTSRateAndState>(),
        std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
        std::make_unique<friction_law_impl::AgingLaw<friction_law_impl::NoTP>>(drParameters.get()),
        std::make_unique<output::OutputManager>(
            std::make_unique<output::RateAndState>(), seissolInstance, numFused)};
  }
}

DynamicRuptureTuple RateAndStateSlipFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSRateAndState>(),
            std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
            std::make_unique<friction_law::SlipLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(),
                seissolInstance,
                numFused)};
  } else {
    return {
        std::make_unique<seissol::initializer::LTSRateAndState>(),
        std::make_unique<initializer::RateAndStateInitializer>(drParameters, seissolInstance),
        std::make_unique<friction_law_impl::SlipLaw<friction_law_impl::NoTP>>(drParameters.get()),
        std::make_unique<output::OutputManager>(
            std::make_unique<output::RateAndState>(), seissolInstance, numFused)};
  }
}

DynamicRuptureTuple LinearSlipWeakeningBimaterialFactory::produce() {
  using Specialization = friction_law_impl::BiMaterialFault;
  using FrictionLawType = friction_law_impl::LinearSlipWeakeningLaw<Specialization>;

  return {
      std::make_unique<seissol::initializer::LTSLinearSlipWeakeningBimaterial>(),
      std::make_unique<initializer::LinearSlipWeakeningBimaterialInitializer>(drParameters,
                                                                              seissolInstance),
      std::make_unique<FrictionLawType>(drParameters.get()),
      std::make_unique<output::OutputManager>(
          std::make_unique<output::LinearSlipWeakeningBimaterial>(), seissolInstance, numFused)};
}

DynamicRuptureTuple LinearSlipWeakeningTPApproxFactory::produce() {
  using Specialization = friction_law_impl::TPApprox;
  using FrictionLawType = friction_law_impl::LinearSlipWeakeningLaw<Specialization>;

  return {
      std::make_unique<seissol::initializer::LTSLinearSlipWeakening>(),
      std::make_unique<initializer::LinearSlipWeakeningInitializer>(drParameters, seissolInstance),
      std::make_unique<FrictionLawType>(drParameters.get()),
      std::make_unique<output::OutputManager>(
          std::make_unique<output::LinearSlipWeakening>(), seissolInstance, numFused)};
}

DynamicRuptureTuple ImposedSlipRatesYoffeFactory::produce() {
  return {
      std::make_unique<seissol::initializer::LTSImposedSlipRatesYoffe>(),
      std::make_unique<initializer::ImposedSlipRatesYoffeInitializer>(drParameters,
                                                                      seissolInstance),
      std::make_unique<friction_law::ImposedSlipRates<friction_law::YoffeSTF>>(drParameters.get()),
      std::make_unique<output::OutputManager>(
          std::make_unique<output::ImposedSlipRates>(), seissolInstance, numFused)};
}

DynamicRuptureTuple ImposedSlipRatesGaussianFactory::produce() {
  return {std::make_unique<seissol::initializer::LTSImposedSlipRatesGaussian>(),
          std::make_unique<initializer::ImposedSlipRatesGaussianInitializer>(drParameters,
                                                                             seissolInstance),
          std::make_unique<friction_law::ImposedSlipRates<friction_law::GaussianSTF>>(
              drParameters.get()),
          std::make_unique<output::OutputManager>(
              std::make_unique<output::ImposedSlipRates>(), seissolInstance, numFused)};
}

DynamicRuptureTuple RateAndStateFastVelocityWeakeningFactory::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializer::LTSRateAndStateThermalPressurization>(),
            std::make_unique<initializer::RateAndStateThermalPressurizationInitializer>(
                drParameters, seissolInstance),
            std::make_unique<
                friction_law::FastVelocityWeakeningLaw<friction_law::ThermalPressurization>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>(),
                seissolInstance,
                numFused)};
  } else {
    return {std::make_unique<seissol::initializer::LTSRateAndStateFastVelocityWeakening>(),
            std::make_unique<initializer::RateAndStateFastVelocityInitializer>(drParameters,
                                                                               seissolInstance),
            std::make_unique<friction_law_impl::FastVelocityWeakeningLaw<friction_law_impl::NoTP>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndState>(), seissolInstance, numFused)};
  }
}
} // namespace seissol::dr::factory
