// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InstantaneousTimeMirrorManager.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Expr/SderivFrontend.h"
#include "Geometry/MeshReader.h"
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Model/CellLocalMatrices.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/MaterialType.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "Reader/Scripting/CompiledReader.h"
#include "Reader/Scripting/ReaderBuilder.h"
#include "SeisSol.h"

#include <array>
#include <cstddef>
#include <memory>
#include <optional>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::physics {

bool isAnisotropicReflectionTypeSupported(
    seissol::initializer::parameters::ReflectionType reflectionType) {
  return reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves;
}

double getSwaveScaledLambda(double lambda, double mu, double velocityScalingFactor) {
  return (lambda + 2.0 * mu) / velocityScalingFactor - 2.0 * velocityScalingFactor * mu;
}

double
    getElasticTimeStepScalingFactor(seissol::initializer::parameters::ReflectionType reflectionType,
                                    double velocityScalingFactor) {
  if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves ||
      reflectionType == seissol::initializer::parameters::ReflectionType::Pwave) {
    return 1.0 / velocityScalingFactor;
  }
  if (reflectionType == seissol::initializer::parameters::ReflectionType::Swave) {
    return velocityScalingFactor;
  }
  return 1.0;
}

std::optional<std::string>
    mirrorScript(model::MaterialType materialType,
                 seissol::initializer::parameters::ReflectionType reflectionType) {
  using seissol::initializer::parameters::ReflectionType;
  const std::string header = "# the time mirror of the reflection type " +
                             std::to_string(static_cast<int>(reflectionType)) +
                             ", written by SeisSol\n";
  if (materialType == model::MaterialType::Elastic) {
    switch (reflectionType) {
    case ReflectionType::BothWaves:
      return header + "out def mu = mu0 * n * n\nout def lambda = lambda0 * n * n\n";
    case ReflectionType::BothWavesVelocity:
      return header +
             "out def rho = rho0 * n\nout def mu = mu0 * n\nout def lambda = lambda0 * n\n";
    case ReflectionType::Pwave:
      return header + "out def lambda = lambda0 * n * n\n";
    case ReflectionType::Swave:
      return header + "out def rho = n * rho0\nout def mu = n * mu0\n"
                      "out def lambda = (lambda0 + 2.0 * mu0) / n - 2.0 * n * mu0\n";
    }
  }
  if (materialType == model::MaterialType::Anisotropic &&
      isAnisotropicReflectionTypeSupported(reflectionType)) {
    // for anisotropic materials, you could scale down density or scale up all the
    // direction-dependent coefficients; the density is simpler
    return header + "out def rho = rho0 / (n * n)\n";
  }
  return std::nullopt;
}

std::vector<std::string> mirrorInputs(const std::vector<std::string>& parameters) {
  std::vector<std::string> inputs{"x", "y", "z", "n"};
  for (const auto& parameter : parameters) {
    inputs.push_back(parameter + "0");
  }
  return inputs;
}

InstantaneousTimeMirrorManager::InstantaneousTimeMirrorManager(seissol::SeisSol& seissolInstance)
    : seissolInstance_(seissolInstance) {}

InstantaneousTimeMirrorManager::~InstantaneousTimeMirrorManager() = default;

void InstantaneousTimeMirrorManager::init(double velocityScalingFactor,
                                          double triggerTime,
                                          seissol::geometry::MeshReader* meshReader,
                                          LTS::Storage& ltsStorage,
                                          const initializer::ClusterLayout* clusterLayout) {

  const auto itmParameters = seissolInstance_.parameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;
  const auto& path = itmParameters.itmMaterialScript;

  // the script of each configuration of the run, checked now: the material script, else the one
  // SeisSol writes for the reflection type
  scripts_.clear();
  scripts_.resize(builtConfigCount());
  scriptNames_.assign(builtConfigCount(), path);
  for (const auto config : ltsStorage.configs()) {
    dispatchConfig(config, [&](auto cfg) {
      using MaterialT = model::MaterialOf<decltype(cfg)>;
      if (path.empty()) {
        const auto source = mirrorScript(MaterialT::Type, reflectionType);
        if (!source.has_value()) {
          if (MaterialT::Type == model::MaterialType::Anisotropic) {
            logError() << "Anisotropic materials cannot have Pwave, Swave, and BothWavesVelocity.";
          } else {
            logError() << "ITM material update is not implemented for this material type; an "
                          "ITM material script gives the material after the mirror.";
          }
          return;
        }
        scriptNames_[config] = "of the reflection type";
        scripts_[config] = std::make_unique<reader::scripting::CompiledReader>(
            expr::compileSderivModule(*source), nullptr);
      } else {
        scripts_[config] = reader::scripting::buildReader(path, {"x", "y", "z"});
      }
      checkMirrorScript<MaterialT>(*scripts_[config], scriptNames_[config]);
    });
  }

  isEnabled_ = true; // This is to sync just before and after the ITM. This does not toggle the ITM.
                     // Need this by default as true for it to work.
  this->velocityScalingFactor_ = velocityScalingFactor;
  this->triggerTime_ = triggerTime;
  this->meshReader_ = meshReader;
  this->ltsStorage_ = &ltsStorage;
  this->clusterLayout_ = clusterLayout;
  setSyncInterval(triggerTime);
  Modules::registerHook(*this, ModuleHook::SynchronizationPoint);
}

void InstantaneousTimeMirrorManager::syncPoint(double currentTime) {
  Module::syncPoint(currentTime);

  logInfo() << "InstantaneousTimeMirrorManager: Factor " << velocityScalingFactor_;
  if (!isEnabled_) {
    logInfo() << "InstantaneousTimeMirrorManager: Skipping syncing at " << currentTime
              << "as it is disabled";
    return;
  }

  logInfo() << "InstantaneousTimeMirrorManager Syncing at " << currentTime;

  logInfo() << "Scaling velocitites by factor of " << velocityScalingFactor_;
  updateVelocities();

  logInfo() << "Updating CellLocalMatrices";
  initializer::initializeCellLocalMatrices(
      *meshReader_, *ltsStorage_, *clusterLayout_, seissolInstance_.parameters().model);

#ifdef ACL_DEVICE
  void* stream = device::DeviceInstance::instance().api().getDefaultStream();
  ltsStorage_->varSynchronizeTo<LTS::LocalIntegration>(
      seissol::initializer::AllocationPlace::Device, stream);
  ltsStorage_->varSynchronizeTo<LTS::NeighboringIntegration>(
      seissol::initializer::AllocationPlace::Device, stream);
  ltsStorage_->varSynchronizeTo<LTS::MaterialParameters>(
      seissol::initializer::AllocationPlace::Device, stream);
  device::DeviceInstance::instance().api().syncDefaultStreamWithHost();
#endif

  logInfo() << "Updating TimeSteps by a factor of " << 1 / velocityScalingFactor_;
  updateTimeSteps();

  logInfo() << "Finished flipping.";
  isEnabled_ = false;
}

void InstantaneousTimeMirrorManager::updateVelocities() {
  const auto& elements = meshReader_->getElements();
  const auto& vertices = meshReader_->getVertices();
  for (const auto config : ltsStorage_->configs()) {
    dispatchConfig(config, [&](auto cfg) {
      using MaterialT = model::MaterialOf<decltype(cfg)>;
      // the cells of the configuration, gathered from its layers, so that the script is bound once;
      // also those of the ghost layers, whose material the fluxes across the boundaries of the
      // rank see
      std::vector<MaterialT*> materials;
      std::vector<std::array<double, Cell::Dim>> centers;
      for (auto& layer : ltsStorage_->leaves()) {
        if (layer.getIdentifier().config != config) {
          continue;
        }
        const bool ghost = layer.getIdentifier().halo == HaloType::Ghost;
        auto* layerMaterials = layer.var<LTS::MaterialData>(cfg);
        const auto* secondary = layer.var<LTS::SecondaryInformation>();
        for (std::size_t cell = 0; cell < layer.size(); ++cell) {
          materials.push_back(&layerMaterials[cell]);
          std::array<double, Cell::Dim> center{};
          if (ghost) {
            // the position of a ghost cell is in the ghost layer, by the rank that owns it
            const auto& linear = meshReader_->linearGhostlayer()[secondary[cell].meshId];
            const auto& metadata =
                meshReader_->getGhostlayerMetadata().at(linear.rank)[linear.inRankIndices[0]];
            for (const auto& vertex : metadata.vertices) {
              for (std::size_t d = 0; d < Cell::Dim; ++d) {
                center[d] += vertex[d] / Cell::NumVertices;
              }
            }
          } else {
            for (const auto vertex : elements[secondary[cell].meshId].vertices) {
              for (std::size_t d = 0; d < Cell::Dim; ++d) {
                center[d] += vertices[vertex].coords[d] / Cell::NumVertices;
              }
            }
          }
          centers.push_back(center);
        }
      }
      mirrorMaterials<MaterialT>(*scripts_[config],
                                 scriptNames_[config],
                                 velocityScalingFactor_,
                                 materials,
                                 centers,
                                 seissolInstance_.parameters().model);
    });
  }
}

void InstantaneousTimeMirrorManager::updateTimeSteps() {
  const auto itmParameters = seissolInstance_.parameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;

  const double timeStepScaling =
      getElasticTimeStepScalingFactor(reflectionType, velocityScalingFactor_);

  if (timeStepScaling != 1.0) {
    scaleClusterTimes(timeStepScaling);
  }
}

void InstantaneousTimeMirrorManager::scaleClusterTimes(double scalingFactor) {
  for (auto& cluster : clusters_) {
    cluster->setClusterTimes(cluster->getClusterTimes() * scalingFactor);
    auto* neighborClusters = cluster->getNeighborClusters();
    for (auto& neighborCluster : *neighborClusters) {
      neighborCluster.ct.setTimeStepSize(neighborCluster.ct.getTimeStepSize() * scalingFactor);
    }
  }
}

void InstantaneousTimeMirrorManager::setClusterVector(
    const std::vector<seissol::time_stepping::AbstractTimeCluster*>& clusters) {
  this->clusters_ = clusters;
}

void initializeTimeMirrorManagers(double scalingFactor,
                                  double triggerTime,
                                  seissol::geometry::MeshReader* meshReader,
                                  LTS::Storage& ltsStorage,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  seissol::SeisSol& seissolInstance,
                                  const initializer::ClusterLayout* clusterLayout) {
  increaseManager.init(scalingFactor, triggerTime, meshReader, ltsStorage, clusterLayout);
  auto itmParameters = seissolInstance.parameters().model.itmParameters;
  const double eps = itmParameters.itmDuration;

  decreaseManager.init(1 / scalingFactor, triggerTime + eps, meshReader, ltsStorage, clusterLayout);
};
} // namespace seissol::physics
