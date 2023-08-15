#ifndef SEISSOL_FACTORY_H
#define SEISSOL_FACTORY_H

#include <stdexcept>
#include <tuple>

#include "DynamicRupture/Initializers/Initializers.h"
#include "FrictionLaws/FrictionSolver.h"
#include "Initializer/DynamicRupture.h"
#include "Output/Output.hpp"

namespace seissol::dr::factory {
/**
 * This struct stores all ingredients, needed for Dynamic Rupture:
 * ltsTree: holds all the data, like parameters (e.g. friction coefficients) or results (e.g.
 * sliprate) initializer: reads the parameter file and stores default values in ltsTree frictionLaw:
 * evaluates friction during the actual simulation output: manages, when we need to write which
 * quantity to which output
 */
template <typename Config>
struct DynamicRuptureTuple {
  std::unique_ptr<seissol::initializers::DynamicRupture<Config>> ltsTree;
  std::unique_ptr<seissol::dr::initializers::BaseDRInitializer<Config>> initializer;
  std::unique_ptr<seissol::dr::friction_law::FrictionSolver<Config>> frictionLaw;
  std::unique_ptr<seissol::dr::output::OutputManager> output;
};

template <typename Config>
class AbstractFactory {
  protected:
  std::shared_ptr<dr::DRParameters> drParameters;

  public:
  AbstractFactory<Config>(std::shared_ptr<dr::DRParameters> drParameters)
      : drParameters(drParameters){};
  virtual ~AbstractFactory<Config>() = default;
  virtual DynamicRuptureTuple<Config> produce() = 0;
};

template <typename Config>
std::unique_ptr<seissol::dr::factory::AbstractFactory<Config>>
    getFactory(std::shared_ptr<dr::DRParameters> dynRupParameter);

} // namespace seissol::dr::factory
#endif // SEISSOL_FACTORY_H
