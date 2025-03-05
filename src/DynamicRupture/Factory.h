// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_FACTORY_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_FACTORY_H_

#include <stdexcept>
#include <tuple>
#include <utility>

#include "DynamicRupture/Initializer/Initializers.h"
#include "FrictionLaws/FrictionSolver.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Output/Output.h"

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
  std::unique_ptr<seissol::initializer::DynamicRupture> ltsTree;
  std::unique_ptr<seissol::dr::initializer::BaseDRInitializer> initializer;
  std::unique_ptr<seissol::dr::friction_law::FrictionSolver> frictionLaw;
  std::unique_ptr<seissol::dr::friction_law::FrictionSolver> frictionLawDevice;
  std::unique_ptr<seissol::dr::output::OutputManager> output;
};

class AbstractFactory {
  protected:
  std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters;
  seissol::SeisSol& seissolInstance;

  public:
  AbstractFactory(std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters,
                  seissol::SeisSol& seissolInstance)
      : drParameters(std::move(drParameters)), seissolInstance(seissolInstance) {};
  virtual ~AbstractFactory() = default;
  virtual DynamicRuptureTuple produce() = 0;
};

std::unique_ptr<seissol::dr::factory::AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
               seissol::SeisSol& seissolInstance);

} // namespace dr::factory
} // namespace seissol

#endif // SEISSOL_SRC_DYNAMICRUPTURE_FACTORY_H_
