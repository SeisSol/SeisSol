#ifndef SEISSOL_FACTORY_H
#define SEISSOL_FACTORY_H

#include <iostream>
#include <stdexcept>
#include <tuple>

#include "Output.h"
#include "DynamicRupture/FrictionLaws/FrictionLaws.h"
#include "DynamicRupture/Initializers/Initializers.h"
#include "Initializer/DynamicRupture.h"

namespace seissol::dr::factory {
struct products {
  std::shared_ptr<seissol::initializers::DynamicRupture> ltsTree;
  std::shared_ptr<seissol::dr::initializers::BaseDRInitializer> initializer;
  std::shared_ptr<seissol::dr::friction_law::BaseFrictionLaw> frictionLaw;
  std::shared_ptr<seissol::dr::output::OutputBase> output;
};
class AbstractFactory;
struct NoFaultFactory;
struct LinearSlipWeakeningFactory;
struct RateAndStateAgingFactory;
struct RateAndStateSlipFactory;
struct LinearSlipWeakeningBimaterialFactory;
struct LinearSlipWeakeningForcedRuptureTimeFactory;
struct ImposedSlipRatesFactory;
struct RateAndStateFastVelocityWeakeningFactory;
struct RateAndStateThermalPressurisationFactory;

std::shared_ptr<AbstractFactory> getFactory(dr::DRParameters& dynRupParameter);
} // namespace seissol::dr::factory

class seissol::dr::factory::AbstractFactory {
  protected:
  dr::DRParameters& drParameters;

  public:
  AbstractFactory(dr::DRParameters& drParameters) : drParameters(drParameters){};
  virtual ~AbstractFactory() {}
  virtual products produce() = 0;
};

class seissol::dr::factory::NoFaultFactory : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class seissol::dr::factory::LinearSlipWeakeningFactory
    : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class seissol::dr::factory::RateAndStateAgingFactory
    : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce();
};

class seissol::dr::factory::RateAndStateSlipFactory : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class seissol::dr::factory::LinearSlipWeakeningBimaterialFactory
    : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class seissol::dr::factory::LinearSlipWeakeningForcedRuptureTimeFactory
    : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class seissol::dr::factory::ImposedSlipRatesFactory : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class seissol::dr::factory::RateAndStateFastVelocityWeakeningFactory
    : public seissol::dr::factory::AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class seissol::dr::factory::RateAndStateThermalPressurisationFactory
    : public seissol::dr::factory::AbstractFactory {
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

std::shared_ptr<seissol::dr::factory::AbstractFactory>
    seissol::dr::factory::getFactory(dr::DRParameters& dynRupParameter);

#endif // SEISSOL_FACTORY_H
