// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_
#define SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_

#include "Geometry/MeshReader.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Model/MaterialType.h"
#include "Model/ScriptParameters.h"
#include "Modules/Module.h"
#include "Reader/Scripting/DataReader.h"
#include "Reader/Scripting/DataTable.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/TimeCluster.h"

#include <algorithm>
#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol {
class SeisSol;
namespace physics {

bool isAnisotropicReflectionTypeSupported(
    seissol::initializer::parameters::ReflectionType reflectionType);
double getSwaveScaledLambda(double lambda, double mu, double velocityScalingFactor);
double
    getElasticTimeStepScalingFactor(seissol::initializer::parameters::ReflectionType reflectionType,
                                    double velocityScalingFactor);

/// The script of a time mirror of the reflection type `reflectionType` for a material of the type
/// `materialType`: what the reflection type does to the material, from its parameters before the
/// mirror (rho0, mu0, ...) and the scaling factor `n` of the mirror. Nothing for a material the
/// reflection type is not defined for.
std::optional<std::string>
    mirrorScript(model::MaterialType materialType,
                 seissol::initializer::parameters::ReflectionType reflectionType);

/// The names a time mirror script reads, for a material with the parameters `parameters`: x, y,
/// z (the center of a cell), n (the scaling factor of the mirror) and each parameter before the
/// mirror (<parameter>0).
std::vector<std::string> mirrorInputs(const std::vector<std::string>& parameters);

/// Checks that the time mirror script `reader`, which `name` names in messages, reads only what a
/// mirror has (mirrorInputs) and gives only parameters of a material of the type `MaterialT` that
/// a material file gives (MaterialT::ParameterList): anything else is an error.
template <typename MaterialT>
void checkMirrorScript(reader::scripting::DataReader& reader, const std::string& name) {
  const auto inputs = mirrorInputs(model::scriptParameterNames<MaterialT>());
  for (const auto& input : reader.inputVars()) {
    if (std::find(inputs.begin(), inputs.end(), input) == inputs.end()) {
      logError() << "The time mirror script" << name << "reads" << input
                 << ", which a mirror does not have; it reads" << model::parameterList(inputs);
    }
  }
  std::vector<std::string> parameters;
  parameters.reserve(MaterialT::ParameterList.size());
  for (const auto& parameter : MaterialT::ParameterList) {
    parameters.push_back(parameter.first);
  }
  for (const auto& output : reader.outputVars()) {
    if (std::find(parameters.begin(), parameters.end(), output) == parameters.end()) {
      logError() << "The time mirror script" << name << "gives" << output
                 << ", which is no parameter of the material that a mirror can change; it gives"
                 << model::parameterList(parameters);
    }
  }
}

/// The materials `materials` after the time mirror script `reader` with the scaling factor `n`,
/// in the cells with the centers `centers`: what the script does not give stays as it was. Each
/// material is then initialized anew (the attenuation of a material with relaxation is fitted to
/// the new parameters). A material that is not admissible -- rho not positive, mu or lambda
/// negative, or any of them not a number -- is an error. `name` names the script in messages.
template <typename MaterialT>
void mirrorMaterials(reader::scripting::DataReader& reader,
                     const std::string& name,
                     double n,
                     const std::vector<MaterialT*>& materials,
                     const std::vector<std::array<double, 3>>& centers,
                     const initializer::parameters::ModelParameters& modelParameters) {
  using reader::scripting::Direction;
  const std::size_t count = materials.size();
  if (count == 0) {
    return;
  }
  // the parameters before the mirror, one column each; those of the material file come first
  const auto names = model::scriptParameterNames<MaterialT>();
  std::vector<std::vector<double>> before(names.size(), std::vector<double>(count));
  std::vector<double> values(names.size());
  for (std::size_t i = 0; i < count; ++i) {
    materials[i]->scriptParameters(values.data());
    for (std::size_t k = 0; k < names.size(); ++k) {
      before[k][i] = values[k];
    }
  }
  const auto& list = MaterialT::ParameterList;
  std::vector<std::vector<double>> after(before.begin(),
                                         before.begin() + static_cast<std::ptrdiff_t>(list.size()));

  reader::scripting::DataTable table(count);
  table.bindViewConst<double>("x", Direction::In, centers.data()->data(), 3, 0);
  table.bindViewConst<double>("y", Direction::In, centers.data()->data(), 3, 1);
  table.bindViewConst<double>("z", Direction::In, centers.data()->data(), 3, 2);
  table.bindConstant<double>("n", n);
  for (std::size_t k = 0; k < names.size(); ++k) {
    table.bindViewConst<double>(names[k] + "0", Direction::In, before[k].data());
  }
  for (const auto& output : reader.outputVars()) {
    for (std::size_t k = 0; k < list.size(); ++k) {
      if (list[k].first == output) {
        table.bindView<double>(output, Direction::Out, after[k].data());
      }
    }
  }
  reader.call(table);

  const auto indexOf = [&](const char* parameter) -> std::optional<std::size_t> {
    for (std::size_t k = 0; k < list.size(); ++k) {
      if (list[k].first == parameter) {
        return k;
      }
    }
    return std::nullopt;
  };
  const auto rho = indexOf("rho");
  const auto mu = indexOf("mu");
  const auto lambda = indexOf("lambda");
  for (std::size_t i = 0; i < count; ++i) {
    // written such that NaN is not admissible either
    if ((rho.has_value() && !(after[*rho][i] > 0.0)) ||
        (mu.has_value() && !(after[*mu][i] >= 0.0)) ||
        (lambda.has_value() && !(after[*lambda][i] >= 0.0))) {
      logError() << "The time mirror script" << name << "gives a material that is not admissible at"
                 << centers[i][0] << centers[i][1] << centers[i][2] << ".";
    }
    for (std::size_t k = 0; k < list.size(); ++k) {
      (*materials[i]).*(list[k].second) = after[k][i];
    }
    materials[i]->initialize(modelParameters);
  }
}

class InstantaneousTimeMirrorManager : public Module {

  seissol::SeisSol& seissolInstance_;
  bool isEnabled_{false};
  double velocityScalingFactor_{1.0};
  double triggerTime_{};

  seissol::geometry::MeshReader* meshReader_{nullptr};
  LTS::Storage* ltsStorage_{nullptr};
  const initializer::ClusterLayout* clusterLayout_{nullptr};

  std::vector<seissol::time_stepping::AbstractTimeCluster*> clusters_;

  public:
  explicit InstantaneousTimeMirrorManager(seissol::SeisSol& seissolInstance);
  ~InstantaneousTimeMirrorManager() override;

  // delete all other constructors due to it being a module.

  InstantaneousTimeMirrorManager(const InstantaneousTimeMirrorManager&) = delete;
  auto operator=(const InstantaneousTimeMirrorManager&) = delete;
  InstantaneousTimeMirrorManager(InstantaneousTimeMirrorManager&&) = delete;
  auto operator=(InstantaneousTimeMirrorManager&&) = delete;

  void init(double velocityScalingFactor,
            double triggerTime,
            seissol::geometry::MeshReader* meshReader,
            LTS::Storage& ltsStorage,
            const initializer::ClusterLayout* clusterLayout);

  void setClusterVector(const std::vector<seissol::time_stepping::AbstractTimeCluster*>& clusters);

  void syncPoint(double currentTime) override;

  private:
  void scaleClusterTimes(double scalingFactor);
  /// The materials of all cells, also of the ghost cells, after the mirror, by its script.
  void updateVelocities();
  void updateTimeSteps();

  /// Per configuration, the script of the mirror -- the material script of the parameters, else
  /// that of the reflection type -- and its name in messages.
  std::vector<std::unique_ptr<reader::scripting::DataReader>> scripts_;
  std::vector<std::string> scriptNames_;
};

void initializeTimeMirrorManagers(double scalingFactor,
                                  double triggerTime,
                                  seissol::geometry::MeshReader* meshReader,
                                  LTS::Storage& ltsStorage,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  seissol::SeisSol& seissolInstance,
                                  const initializer::ClusterLayout* clusterLayout);
} // namespace physics
} // namespace seissol

#endif // SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_
