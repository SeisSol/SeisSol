#ifndef SEISSOL_FACTORY_H
#define SEISSOL_FACTORY_H

#include <stdexcept>
#include <tuple>

#include "DynamicRupture/Initializers/Initializers.h"
#include "FrictionLaws/FrictionSolver.h"
#include "Initializer/DynamicRupture.h"
#include "Output/Output.hpp"

namespace seissol {
class SeisSol;
namespace dr::factory {
/**
 * This struct stores all ingredients, needed for Dynamic Rupture:
 * ltsTree: holds all the data, like parameters (e.g. friction coefficients) or results (e.g.
 * sliprate) initializer: reads the parameter file and stores default values in ltsTree frictionLaw:
 * evaluates friction during the actual simulation output: manages, when we need to write which
 * quantity to which output
 */
struct DynamicRuptureTuple {
  std::unique_ptr<seissol::initializers::DynamicRupture> ltsTree;
  std::unique_ptr<seissol::dr::initializers::BaseDRInitializer> initializer;
  std::unique_ptr<seissol::dr::friction_law::FrictionSolver> frictionLaw;
  std::unique_ptr<seissol::dr::output::OutputManager> output;
};

class AbstractFactory {
  protected:
  std::shared_ptr<dr::DRParameters> drParameters;
  seissol::SeisSol& seissolInstance;

  public:
  AbstractFactory(std::shared_ptr<dr::DRParameters> drParameters, seissol::SeisSol& seissolInstance)
      : drParameters(drParameters), seissolInstance(seissolInstance){};
  virtual ~AbstractFactory() = default;
  virtual DynamicRuptureTuple produce() = 0;
};

class NoFaultFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

class LinearSlipWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

class RateAndStateAgingFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

class RateAndStateSlipFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

class LinearSlipWeakeningBimaterialFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

class ImposedSlipRatesYoffeFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

class ImposedSlipRatesGaussianFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

class RateAndStateFastVelocityWeakeningFactory : public AbstractFactory {
  public:
  using AbstractFactory::AbstractFactory;
  DynamicRuptureTuple produce() override;
};

std::unique_ptr<seissol::dr::factory::AbstractFactory>
    getFactory(std::shared_ptr<dr::DRParameters> dynRupParameter,
               seissol::SeisSol& seissolInstance);

} // namespace dr::factory
} // namespace seissol
#endif // SEISSOL_FACTORY_H
