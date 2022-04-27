#ifndef SEISSOL_FACTORY_H
#define SEISSOL_FACTORY_H

#include <iostream>
#include <stdexcept>
#include <tuple>

#include "DynamicRupture/Initializers/Initializers.h"
#include "FrictionLaws/FrictionSolver.h"
#include "Initializer/DynamicRupture.h"
#include "Output/Output.hpp"

namespace seissol::dr::factory {
struct Products {
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
  virtual Products produce() = 0;
};

class NoFaultFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

class LinearSlipWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

class RateAndStateAgingFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce();
};

class RateAndStateSlipFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

class LinearSlipWeakeningBimaterialFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

class LinearSlipWeakeningForcedRuptureTimeFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

class ImposedSlipRatesYoffeFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

class ImposedSlipRatesGaussianFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

class RateAndStateFastVelocityWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  virtual Products produce() override;
};

std::unique_ptr<seissol::dr::factory::AbstractFactory>
    getFactory(dr::DRParameters& dynRupParameter);

} // namespace seissol::dr::factory
#endif // SEISSOL_FACTORY_H
