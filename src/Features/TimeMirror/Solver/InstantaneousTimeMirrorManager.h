// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_FEATURES_TIMEMIRROR_SOLVER_INSTANTANEOUSTIMEMIRRORMANAGER_H_
#define SEISSOL_SRC_FEATURES_TIMEMIRROR_SOLVER_INSTANTANEOUSTIMEMIRRORMANAGER_H_

#include "Core/Geometry/MeshReader.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/LTSTree.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Modules/Module.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Compute/CellCluster.h"
#include "Core/TimeStepping/Halo/GhostCluster.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "Features/TimeMirror/ITMParameters.h"

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

class InstantaneousTimeMirrorManager : public Module {

  seissol::SeisSol& seissolInstance_;
  bool isEnabled_{false};
  double velocityScalingFactor_{1.0};
  double triggerTime_{};

  seissol::geometry::MeshReader* meshReader_{nullptr};
  LTS::Storage* ltsStorage_{nullptr};
  const initializer::ClusterLayout* clusterLayout_{nullptr};

  std::vector<seissol::solver::AbstractTimeCluster*> clusters_;

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

  void setClusterVector(const std::vector<seissol::solver::AbstractTimeCluster*>& clusters);

  void syncPoint(double currentTime) override;

  private:
  void scaleClusterTimes(double scalingFactor);
  void updateVelocities();
  void updateTimeSteps();
};

void initializeTimeMirrorManagers(double scalingFactor,
                                  double triggerTime,
                                  seissol::geometry::MeshReader* meshReader,
                                  LTS::Storage& ltsStorage,
                                  InstantaneousTimeMirrorManager& increaseManager,
                                  InstantaneousTimeMirrorManager& decreaseManager,
                                  const seissol::SeisSol& seissolInstance,
                                  const initializer::ClusterLayout* clusterLayout);
} // namespace physics
} // namespace seissol

#endif // SEISSOL_SRC_FEATURES_TIMEMIRROR_SOLVER_INSTANTANEOUSTIMEMIRRORMANAGER_H_
