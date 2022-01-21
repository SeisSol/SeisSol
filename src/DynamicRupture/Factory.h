#ifndef SEISSOL_FACTORY_H
#define SEISSOL_FACTORY_H

#include <iostream>
#include <stdexcept>
#include <tuple>

#include "Output/Output.hpp"
#include "FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Initializers/Initializers.h"
#include "Initializer/DynamicRupture.h"

namespace seissol::dr::factory {
struct products {
  std::unique_ptr<seissol::initializers::DynamicRupture> ltsTree;
  std::unique_ptr<seissol::dr::initializers::BaseDRInitializer> initializer;
  std::unique_ptr<seissol::dr::friction_law::FrictionSolver> frictionLaw;
  std::unique_ptr<seissol::dr::output::Base> output;
};

class AbstractFactory {
  protected:
  dr::DRParameters& drParameters;

  public:
  AbstractFactory(dr::DRParameters& drParameters) : drParameters(drParameters){};
  virtual ~AbstractFactory() {}
  virtual products produce() = 0;
};

class NoFaultFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class LinearSlipWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class RateAndStateAgingFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce();
};

class RateAndStateSlipFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class LinearSlipWeakeningBimaterialFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class LinearSlipWeakeningForcedRuptureTimeFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class ImposedSlipRatesFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class RateAndStateFastVelocityWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

class RateAndStateThermalPressurisationFactory : public AbstractFactory {
  using AbstractFactory::AbstractFactory;
  virtual products produce() override;
};

std::unique_ptr<seissol::dr::factory::AbstractFactory>
    getFactory(dr::DRParameters& dynRupParameter);

} // namespace seissol::dr::factory
#endif // SEISSOL_FACTORY_H
