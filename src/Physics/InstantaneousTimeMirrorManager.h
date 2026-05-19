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
#include "Modules/Module.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include "Solver/TimeStepping/AbstractTimeCluster.h"
#include "Solver/TimeStepping/TimeCluster.h"

#include <vector>

namespace seissol {
class SeisSol;
namespace ITM {

bool isAnisotropicReflectionTypeSupported(
    seissol::initializer::parameters::ReflectionType reflectionType);
double getSwaveScaledLambda(double lambda, double mu, double velocityScalingFactor);
double
    getElasticTimeStepScalingFactor(seissol::initializer::parameters::ReflectionType reflectionType,
                                    double velocityScalingFactor);

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
  explicit InstantaneousTimeMirrorManager(seissol::SeisSol& seissolInstance)
      : seissolInstance_(seissolInstance) {};

  void init(double velocityScalingFactor,
            double triggerTime,
            seissol::geometry::MeshReader* meshReader,
            LTS::Storage& ltsStorage,
            const initializer::ClusterLayout* clusterLayout);

  void setClusterVector(const std::vector<seissol::time_stepping::AbstractTimeCluster*>& clusters);

  void syncPoint(double currentTime) override;

  private:
  void scaleClusterTimes(double scalingFactor);
  void updateVelocities();
  void updateTimeSteps();

  template <typename MaterialType>
  void updateVelocitiesForMaterialType();

  template <typename MaterialType>
  void updateTimeStepsForMaterialType();
};

void initializeTimeMirrorManagers(double scalingFactor,
                                  double triggerTime,
                                  seissol::geometry::MeshReader* meshReader,
                                  LTS::Storage& ltsStorage,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  seissol::SeisSol& seissolInstance,
                                  const initializer::ClusterLayout* clusterLayout);
} // namespace ITM
} // namespace seissol

#endif // SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_
