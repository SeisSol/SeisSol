// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
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
#include "Memory/Tree/Lut.h"
#include "Modules/Module.h"
#include "Solver/time_stepping/AbstractGhostTimeCluster.h"
#include "Solver/time_stepping/TimeCluster.h"

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
  initializer::LTSTree* ltsTree{};
  initializer::LTS* lts{};
  initializer::Lut* ltsLut{};
  const TimeStepping* timestepping{};

  std::vector<std::unique_ptr<seissol::time_stepping::TimeCluster>>* timeClusters{};
  std::vector<std::unique_ptr<seissol::time_stepping::AbstractGhostTimeCluster>>*
      ghostTimeClusters{};

  public:
  InstantaneousTimeMirrorManager(seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance) {};

  void init(double velocityScalingFactor,
            double triggerTime,
            seissol::geometry::MeshReader* meshReader,
            initializer::LTSTree* ltsTree,
            initializer::LTS* lts,
            initializer::Lut* ltsLut,
            const TimeStepping* timestepping); // An empty timestepping is added. Need to discuss
                                               // what exactly is to be sent here

  void setTimeClusterVector(
      std::vector<std::unique_ptr<seissol::time_stepping::TimeCluster>>* clusters);

  void setGhostClusterVector(
      std::vector<std::unique_ptr<seissol::time_stepping::AbstractGhostTimeCluster>>* clusters);

  void syncPoint(double currentTime) override;

  private:
  void updateVelocities();
  void updateTimeSteps();
};

void initializeTimeMirrorManagers(
    double scalingFactor,
    double triggerTime,
    seissol::geometry::MeshReader* meshReader,
    initializer::LTSTree* ltsTree,
    initializer::LTS* lts,
    initializer::Lut* ltsLut,
    InstantaneousTimeMirrorManager& increaseManager,
    InstantaneousTimeMirrorManager& decreaseManager,
    seissol::SeisSol& seissolInstance,
    const TimeStepping* timestepping); // An empty timestepping is added. Need to discuss what
                                       // exactly is to be sent here
} // namespace ITM
} // namespace seissol

#endif // SEISSOL_SRC_PHYSICS_INSTANTANEOUSTIMEMIRRORMANAGER_H_
