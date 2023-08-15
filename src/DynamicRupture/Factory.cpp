#include "Factory.h"

#include "FrictionLaws/FrictionLaws.h"
#include "FrictionLaws/ThermalPressurization/ThermalPressurization.h"

#ifdef ACL_DEVICE
namespace friction_law_impl = seissol::dr::friction_law::gpu;
#else
namespace friction_law_impl = seissol::dr::friction_law;
#endif

namespace seissol::dr::factory {
template <typename Config>
class NoFaultFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
class LinearSlipWeakeningFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
class RateAndStateAgingFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
class RateAndStateSlipFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
class LinearSlipWeakeningBimaterialFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
class ImposedSlipRatesYoffeFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
class ImposedSlipRatesGaussianFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
class RateAndStateFastVelocityWeakeningFactory : public AbstractFactory<Config> {
  public:
  using AbstractFactory<Config>::AbstractFactory;
  DynamicRuptureTuple<Config> produce() override;
};

template <typename Config>
std::unique_ptr<AbstractFactory<Config>>
    getFactory(std::shared_ptr<dr::DRParameters> drParameters) {
  switch (drParameters->frictionLawType) {
  case FrictionLawType::NoFault:
    return std::make_unique<NoFaultFactory<Config>>(drParameters);
  case FrictionLawType::ImposedSlipRatesYoffe:
    return std::make_unique<ImposedSlipRatesYoffeFactory<Config>>(drParameters);
  case FrictionLawType::ImposedSlipRatesGaussian:
    return std::make_unique<ImposedSlipRatesGaussianFactory<Config>>(drParameters);
  case FrictionLawType::LinearSlipWeakening:
    return std::make_unique<LinearSlipWeakeningFactory<Config>>(drParameters);
  case FrictionLawType::LinearSlipWeakeningBimaterial:
    return std::make_unique<LinearSlipWeakeningBimaterialFactory<Config>>(drParameters);
  case FrictionLawType::RateAndStateAgingLaw:
    return std::make_unique<RateAndStateAgingFactory<Config>>(drParameters);
  case FrictionLawType::RateAndStateSlipLaw:
    return std::make_unique<RateAndStateSlipFactory<Config>>(drParameters);
  case FrictionLawType::RateAndStateVelocityWeakening:
    logError() << "friction law 7 currently disabled";
    return std::unique_ptr<AbstractFactory<Config>>(nullptr);
  case FrictionLawType::RateAndStateAgingNucleation:
    logError() << "friction law 101 currently disabled";
    return std::unique_ptr<AbstractFactory<Config>>(nullptr);
  case FrictionLawType::RateAndStateFastVelocityWeakening:
    return std::make_unique<RateAndStateFastVelocityWeakeningFactory<Config>>(drParameters);
  default:
    logError() << "unknown friction law";
    return nullptr;
  }
}

template <typename Config>
DynamicRuptureTuple<Config> NoFaultFactory<Config>::produce() {
  return {std::make_unique<seissol::initializers::DynamicRupture<Config>>(),
          std::make_unique<initializers::NoFaultInitializer<Config>>(drParameters),
          std::make_unique<friction_law::NoFault<Config>>(drParameters.get()),
          std::make_unique<output::OutputManager>(std::make_unique<output::NoFault>())};
}

template <typename Config>
DynamicRuptureTuple<Config> LinearSlipWeakeningFactory<Config>::produce() {
  return {
      std::make_unique<seissol::initializers::LTSLinearSlipWeakening<Config>>(),
      std::make_unique<initializers::LinearSlipWeakeningInitializer<Config>>(drParameters),
      std::make_unique<
          friction_law_impl::LinearSlipWeakeningLaw<Config,
                                                    friction_law_impl::NoSpecialization<Config>>>(
          drParameters.get()),
      std::make_unique<output::OutputManager>(std::make_unique<output::LinearSlipWeakening>())};
}

template <typename Config>
DynamicRuptureTuple<Config> RateAndStateAgingFactory<Config>::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializers::LTSRateAndState<Config>>(),
            std::make_unique<initializers::RateAndStateInitializer<Config>>(drParameters),
            std::make_unique<
                friction_law::AgingLaw<Config, friction_law::ThermalPressurization<Config>>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>())};
  } else {
    return {std::make_unique<seissol::initializers::LTSRateAndState<Config>>(),
            std::make_unique<initializers::RateAndStateInitializer<Config>>(drParameters),
            std::make_unique<friction_law_impl::AgingLaw<Config, friction_law_impl::NoTP<Config>>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>())};
  }
}

template <typename Config>
DynamicRuptureTuple<Config> RateAndStateSlipFactory<Config>::produce() {
  if (drParameters->isThermalPressureOn) {
    return {std::make_unique<seissol::initializers::LTSRateAndState<Config>>(),
            std::make_unique<initializers::RateAndStateInitializer<Config>>(drParameters),
            std::make_unique<
                friction_law::SlipLaw<Config, friction_law::ThermalPressurization<Config>>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(
                std::make_unique<output::RateAndStateThermalPressurization>())};
  } else {
    return {std::make_unique<seissol::initializers::LTSRateAndState<Config>>(),
            std::make_unique<initializers::RateAndStateInitializer<Config>>(drParameters),
            std::make_unique<friction_law_impl::SlipLaw<Config, friction_law_impl::NoTP<Config>>>(
                drParameters.get()),
            std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>())};
  }
}

template <typename Config>
DynamicRuptureTuple<Config> LinearSlipWeakeningBimaterialFactory<Config>::produce() {
  using Specialization = friction_law_impl::BiMaterialFault<Config>;
  using FrictionLawType = friction_law_impl::LinearSlipWeakeningLaw<Config, Specialization>;

  return {std::make_unique<seissol::initializers::LTSLinearSlipWeakeningBimaterial<Config>>(),
          std::make_unique<initializers::LinearSlipWeakeningBimaterialInitializer<Config>>(
              drParameters),
          std::make_unique<FrictionLawType>(drParameters.get()),
          std::make_unique<output::OutputManager>(
              std::make_unique<output::LinearSlipWeakeningBimaterial>())};
}

template <typename Config>
DynamicRuptureTuple<Config> ImposedSlipRatesYoffeFactory<Config>::produce() {
  return {std::make_unique<seissol::initializers::LTSImposedSlipRatesYoffe<Config>>(),
          std::make_unique<initializers::ImposedSlipRatesYoffeInitializer<Config>>(drParameters),
          std::make_unique<friction_law::ImposedSlipRates<Config, friction_law::YoffeSTF<Config>>>(
              drParameters.get()),
          std::make_unique<output::OutputManager>(std::make_unique<output::ImposedSlipRates>())};
}

template <typename Config>
DynamicRuptureTuple<Config> ImposedSlipRatesGaussianFactory<Config>::produce() {
  return {
      std::make_unique<seissol::initializers::LTSImposedSlipRatesGaussian<Config>>(),
      std::make_unique<initializers::ImposedSlipRatesGaussianInitializer<Config>>(drParameters),
      std::make_unique<friction_law::ImposedSlipRates<Config, friction_law::GaussianSTF<Config>>>(
          drParameters.get()),
      std::make_unique<output::OutputManager>(std::make_unique<output::ImposedSlipRates>())};
}

template <typename Config>
DynamicRuptureTuple<Config> RateAndStateFastVelocityWeakeningFactory<Config>::produce() {
  if (drParameters->isThermalPressureOn) {
    return {
        std::make_unique<seissol::initializers::LTSRateAndStateThermalPressurization<Config>>(),
        std::make_unique<initializers::RateAndStateThermalPressurizationInitializer<Config>>(
            drParameters),
        std::make_unique<
            friction_law::FastVelocityWeakeningLaw<Config,
                                                   friction_law::ThermalPressurization<Config>>>(
            drParameters.get()),
        std::make_unique<output::OutputManager>(
            std::make_unique<output::RateAndStateThermalPressurization>())};
  } else {
    return {
        std::make_unique<seissol::initializers::LTSRateAndStateFastVelocityWeakening<Config>>(),
        std::make_unique<initializers::RateAndStateFastVelocityInitializer<Config>>(drParameters),
        std::make_unique<
            friction_law_impl::FastVelocityWeakeningLaw<Config, friction_law_impl::NoTP<Config>>>(
            drParameters.get()),
        std::make_unique<output::OutputManager>(std::make_unique<output::RateAndState>())};
  }
}
} // namespace seissol::dr::factory
