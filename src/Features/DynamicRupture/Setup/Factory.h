// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FACTORY_H_
#define SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FACTORY_H_

#include "Core/Data/DataReader.h"
#include "Features/DynamicRupture/Output/Output.h"
#include "Features/DynamicRupture/Setup/Initializer/Initializers.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"
#include "Features/DynamicRupture/Solver/Storage.h"

#include <memory>
#include <optional>
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

/**
 * What the solver of a friction law needs: the storage of the faces, with the variables of the
 * friction law, and its friction solvers on the host and the device.
 */
struct FrictionLawSolver {
  std::unique_ptr<seissol::DynamicRupture> storage;
  seissol::dr::friction_law::FrictionSolverFactory frictionLaw;
  seissol::dr::friction_law::FrictionSolverFactory frictionLawDevice;
};

/**
 * The parts of the friction law of a run (see DynamicRuptureTuple). Its solver and its initializer
 * need only the parameters of the run; only the initializer of a factory with the mesh and the
 * fault files reads the parameters of the faults from these (see BaseDRInitializer).
 */
class AbstractFactory {
  protected:
  std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters_;
  std::optional<initializer::FaultMesh> faultMesh_;
  data::ModelLoader models_;

  /// The initializer `T`: with the mesh and the fault files if the factory has them.
  template <typename T>
  [[nodiscard]] std::unique_ptr<initializer::BaseDRInitializer> makeInitializer() const {
    if (faultMesh_.has_value()) {
      return std::make_unique<T>(drParameters_, *faultMesh_, models_);
    }
    return std::make_unique<T>(drParameters_);
  }

  public:
  /// `models` opens the fault parameter files.
  AbstractFactory(std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters,
                  const initializer::FaultMesh& faultMesh,
                  data::ModelLoader models)
      : drParameters_(std::move(drParameters)), faultMesh_(faultMesh),
        models_(std::move(models)) {};

  /// without the mesh and the fault files
  explicit AbstractFactory(
      std::shared_ptr<seissol::initializer::parameters::DRParameters> drParameters)
      : drParameters_(std::move(drParameters)) {};

  virtual ~AbstractFactory() = default;

  /// the storage of the faces and the friction solvers
  virtual FrictionLawSolver solver() = 0;

  /// the initializer of the friction law
  virtual std::unique_ptr<initializer::BaseDRInitializer> initializer() = 0;

  /// what the output of the fault writes for the friction law
  virtual std::unique_ptr<seissol::dr::output::ReceiverOutput> output() = 0;

  /// all parts together
  DynamicRuptureTuple produce();
};

std::unique_ptr<seissol::dr::factory::AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters,
               const initializer::FaultMesh& faultMesh,
               const data::ModelLoader& models);

/// The factory of the friction law of `drParameters` without the mesh and the fault files: for its
/// solver, and for the initial state derived from parameters in the storage (e.g. on the
/// synthetic faults of the proxy).
std::unique_ptr<seissol::dr::factory::AbstractFactory>
    getFactory(const std::shared_ptr<seissol::initializer::parameters::DRParameters>& drParameters);

} // namespace seissol::dr::factory

#endif // SEISSOL_SRC_FEATURES_DYNAMICRUPTURE_SETUP_FACTORY_H_
