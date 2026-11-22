// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FACTORY_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FACTORY_H_

#include "Core/Data/DataReader.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Features/DynamicRupture/Output/Output.h"
#include "Features/DynamicRupture/Setup/Initializer/Initializers.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"

#include <stdexcept>
#include <tuple>
#include <utility>

namespace seissol::dr::factory {
/**
 * This struct stores all ingredients, needed for Dynamic Rupture:
 * storage: holds all the data, like parameters (e.g. friction coefficients) or results (e.g.
 * sliprate) initializer: reads the parameter file and stores default values in storage frictionLaw:
 * evaluates friction during the actual simulation output: manages, when we need to write which
 * quantity to which output
 */
struct DynamicRuptureTuple {
  std::unique_ptr<seissol::DynamicRupture> storage;
  std::unique_ptr<seissol::dr::initializer::BaseDRInitializer> initializer;
  seissol::dr::friction_law::FrictionSolverFactory frictionLaw;
  seissol::dr::friction_law::FrictionSolverFactory frictionLawDevice;
  // what the output of the fault writes for the friction law (see output::OutputManager)
  std::unique_ptr<seissol::dr::output::ReceiverOutput> output;
};

class AbstractFactory {
  protected:
  std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters_;
  initializer::FaultMesh faultMesh_;
  data::ModelLoader models_;

  public:
  /// `models` opens the fault parameter files.
  AbstractFactory(std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters,
                  const initializer::FaultMesh& faultMesh,
                  data::ModelLoader models)
      : drParameters_(std::move(drParameters)), faultMesh_(faultMesh),
        models_(std::move(models)) {};
  virtual ~AbstractFactory() = default;
  virtual DynamicRuptureTuple produce() = 0;
};

std::unique_ptr<seissol::dr::factory::AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
               const initializer::FaultMesh& faultMesh,
               const data::ModelLoader& models);

} // namespace seissol::dr::factory

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FACTORY_H_
