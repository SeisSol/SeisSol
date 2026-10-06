// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "InstantaneousTimeMirrorManager.h"

#include "Common/ConfigDispatch.h"
#include "Common/Constants.h"
#include "Initializer/Model/CellLocalMatrices.h"
#include "Initializer/Parameters/ModelParameters.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/CommonDatastructures.h"
#include "Modules/Module.h"
#include "Modules/Modules.h"
#include "Reader/Scripting/DataTable.h"
#include "Reader/Scripting/ReaderBuilder.h"
#include "SeisSol.h"

#include <array>
#include <cstddef>
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

namespace {
void checkSupported(const model::Material& material,
                    initializer::parameters::ReflectionType reflectionType) {
  if (material.getMaterialType() != model::MaterialType::Elastic &&
      material.getMaterialType() != model::MaterialType::Anisotropic) {
    logError() << "ITM material update is not implemented for this material type.";
  }
  if (material.getMaterialType() == model::MaterialType::Anisotropic &&
      !isAnisotropicReflectionTypeSupported(reflectionType)) {
    logError() << "Anisotropic materials cannot have Pwave, Swave, and BothWavesVelocity.";
  }
}
} // namespace

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

  // check over all cells (cheap; though it can be reduced to at most one per layer)
  const bool scripted = !itmParameters.itmMaterialScript.empty();
  for (auto& layer : ltsStorage.leaves()) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      const auto* materials = layer.var<LTS::MaterialData>(cfg);
      for (std::size_t i = 0; i < layer.size(); ++i) {
        checkSupported(materials[i], reflectionType);
        if (scripted && materials[i].getMaterialType() != model::MaterialType::Elastic) {
          logError() << "An ITM material script gives an elastic material (rho, mu, lambda); "
                        "the mesh has cells of another material.";
        }
      }
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
  device::DeviceInstance::instance().api().syncDefaultStreamWithHost();
#endif

  logInfo() << "Updating TimeSteps by a factor of " << 1 / velocityScalingFactor_;
  updateTimeSteps();

  logInfo() << "Finished flipping.";
  isEnabled_ = false;
}

void InstantaneousTimeMirrorManager::updateMaterialsByScript(const std::string& path) {
  const auto& elements = meshReader_->getElements();
  const auto& vertices = meshReader_->getVertices();
  for (const auto config : ltsStorage_->configs()) {
    dispatchConfig(config, [&](auto cfg) {
      // the cells of the configuration, gathered from its layers, so that the script is bound once
      std::vector<model::Material*> materials;
      std::vector<std::array<double, Cell::Dim>> centers;
      for (auto& layer : ltsStorage_->leaves(Ghost)) {
        if (layer.getIdentifier().config != config) {
          continue;
        }
        auto* layerMaterials = layer.var<LTS::MaterialData>(cfg);
        const auto* secondary = layer.var<LTS::SecondaryInformation>();
        for (std::size_t cell = 0; cell < layer.size(); ++cell) {
          materials.push_back(&layerMaterials[cell]);
          std::array<double, Cell::Dim> center{};
          for (const auto vertex : elements[secondary[cell].meshId].vertices) {
            for (std::size_t d = 0; d < Cell::Dim; ++d) {
              center[d] += vertices[vertex].coords[d] / Cell::NumVertices;
            }
          }
          centers.push_back(center);
        }
      }
      if (materials.empty()) {
        return;
      }

      const std::size_t count = materials.size();
      std::vector<double> rho(count);
      std::vector<double> mu(count);
      std::vector<double> lambda(count);
      for (std::size_t i = 0; i < count; ++i) {
        rho[i] = materials[i]->getDensity();
        mu[i] = materials[i]->getMuBar();
        lambda[i] = materials[i]->getLambdaBar();
      }
      // what the script does not give stays as it is
      auto newRho = rho;
      auto newMu = mu;
      auto newLambda = lambda;

      using reader::scripting::Direction;
      reader::scripting::DataTable table(count);
      table.bindViewConst<double>("x", Direction::In, centers.data()->data(), Cell::Dim, 0);
      table.bindViewConst<double>("y", Direction::In, centers.data()->data(), Cell::Dim, 1);
      table.bindViewConst<double>("z", Direction::In, centers.data()->data(), Cell::Dim, 2);
      // the material before the mirror; the script gives the one after it
      table.bindViewConst<double>("rho0", Direction::In, rho.data());
      table.bindViewConst<double>("mu0", Direction::In, mu.data());
      table.bindViewConst<double>("lambda0", Direction::In, lambda.data());
      table.bindConstant<double>("n", velocityScalingFactor_);

      const auto reader = reader::scripting::buildReader(path, {"x", "y", "z"});
      for (const auto& name : reader->outputVars()) {
        if (name == "rho") {
          table.bindView<double>(name, Direction::Out, newRho.data());
        } else if (name == "mu") {
          table.bindView<double>(name, Direction::Out, newMu.data());
        } else if (name == "lambda") {
          table.bindView<double>(name, Direction::Out, newLambda.data());
        } else {
          logError() << "The ITM material script" << path << "gives" << name
                     << "; it gives rho, mu and lambda.";
        }
      }
      reader->call(table);

      for (std::size_t i = 0; i < count; ++i) {
        if (newLambda[i] < 0.0 || newMu[i] < 0.0 || newRho[i] <= 0.0) {
          logError() << "The ITM material script" << path
                     << "gives a material that is not admissible at" << centers[i][0]
                     << centers[i][1] << centers[i][2] << ".";
        }
        materials[i]->setDensity(newRho[i]);
        materials[i]->setLameParameters(newMu[i], newLambda[i]);
      }
    });
  }
}

void InstantaneousTimeMirrorManager::updateVelocities() {
  const auto itmParameters = seissolInstance_.parameters().model.itmParameters;
  const auto reflectionType = itmParameters.itmReflectionType;

  if (!itmParameters.itmMaterialScript.empty()) {
    updateMaterialsByScript(itmParameters.itmMaterialScript);
    return;
  }

  const auto updateMaterial = [&](model::Material& material) {
    if (material.getMaterialType() == model::MaterialType::Elastic) {
      const auto rho = material.getDensity();
      const auto lambda = material.getLambdaBar();
      const auto mu = material.getMuBar();

      if (reflectionType == seissol::initializer::parameters::ReflectionType::BothWaves) {
        material.setLameParameters(mu * velocityScalingFactor_ * velocityScalingFactor_,
                                   lambda * velocityScalingFactor_ * velocityScalingFactor_);
      } else if (reflectionType ==
                 seissol::initializer::parameters::ReflectionType::BothWavesVelocity) {
        material.setDensity(rho * velocityScalingFactor_);
        material.setLameParameters(mu * velocityScalingFactor_, lambda * velocityScalingFactor_);
      } else if (reflectionType == seissol::initializer::parameters::ReflectionType::Pwave) {
        material.setLameParameters(mu, lambda * velocityScalingFactor_ * velocityScalingFactor_);
      } else if (reflectionType == seissol::initializer::parameters::ReflectionType::Swave) {
        const auto newLambda = getSwaveScaledLambda(lambda, mu, velocityScalingFactor_);
        if (newLambda < 0.0) {
          logError() << "New lambda is negative. This is not allowed. Please adjust your scaling "
                        "factor.";
        }
        material.setDensity(velocityScalingFactor_ * rho);
        material.setLameParameters(velocityScalingFactor_ * mu, newLambda);
      } else {
        logError() << "Unknown reflection type; material cannot be updated.";
      }
    } else if (material.getMaterialType() == model::MaterialType::Anisotropic) {
      // for anisotropic materials, you could scale down density
      // or scale up all the direction-dependent coefficients.
      // we scale density for code simplicity
      material.setDensity(material.getDensity() /
                          (velocityScalingFactor_ * velocityScalingFactor_));
    }
  };

  for (auto& layer : ltsStorage_->leaves(Ghost)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      auto* materials = layer.var<LTS::MaterialData>(cfg);

#pragma omp parallel for schedule(static)
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        updateMaterial(materials[cell]);
      }
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
