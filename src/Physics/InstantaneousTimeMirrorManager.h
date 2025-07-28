// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_
#define SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_

#include "Geometry/MeshReader.h"
#include "Initializer/Typedefs.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/LTSTree.h"
#include "Modules/Module.h"
#include "Solver/TimeStepping/AbstractGhostTimeCluster.h"
#include "Solver/TimeStepping/TimeCluster.h"
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Solver/TimeStepping/AbstractTimeCluster.h>

namespace seissol {
class SeisSol;
namespace ITM {

class InstantaneousTimeMirrorManager : Module {
  seissol::SeisSol& seissolInstance;
  bool isEnabled{false};
  double velocityScalingFactor{1.0};
  double timeStepScalingFactor{1.0};
  double triggerTime{};

  seissol::geometry::MeshReader* meshReader{};
  LTS::Storage* ltsTree{};
  const initializer::ClusterLayout* clusterLayout;

  std::vector<seissol::time_stepping::AbstractTimeCluster*> clusters;

  public:
  explicit InstantaneousTimeMirrorManager(seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance) {};

  void init(
      double velocityScalingFactor,
      double triggerTime,
      seissol::geometry::MeshReader* meshReader,
      LTS::Storage* ltsTree,
      const initializer::ClusterLayout* clusterLayout); // An empty timestepping is added. Need to
                                                        // discuss what exactly is to be sent here

  void setClusterVector(const std::vector<seissol::time_stepping::AbstractTimeCluster*>& clusters);

  void syncPoint(double currentTime) override;

  private:
  void updateVelocities();
  void updateTimeSteps();
};

void initializeTimeMirrorManagers(
    double scalingFactor,
    double triggerTime,
    seissol::geometry::MeshReader* meshReader,
    LTS::Storage* ltsTree,
    InstantaneousTimeMirrorManager& increaseManager,
    InstantaneousTimeMirrorManager& decreaseManager,
    seissol::SeisSol& seissolInstance,
    const initializer::ClusterLayout* clusterLayout); // An empty timestepping is added. Need to
                                                      // discuss what exactly is to be sent here
} // namespace ITM
} // namespace seissol

#endif // SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_
