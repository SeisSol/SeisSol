#ifndef SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H
#define SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H

#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include "Modules/Module.h"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/typedefs.hpp"
#include "Solver/time_stepping/TimeCluster.h"
#include "Solver/time_stepping/AbstractGhostTimeCluster.h"

namespace seissol::ITM {

class InstantaneousTimeMirrorManager : Module {
  bool isEnabled;
  double velocityScalingFactor{};
  double triggerTime{};

  seissol::geometry::MeshReader* meshReader{};
  initializers::LTSTree* ltsTree{};
  initializers::LTS* lts{};
  initializers::Lut* ltsLut{};
  const TimeStepping* timestepping{};

  std::vector<std::unique_ptr<seissol::time_stepping::TimeCluster>>* timeClusters;
  std::vector<std::unique_ptr<seissol::time_stepping::AbstractGhostTimeCluster>>* ghostTimeClusters;

  public:
  InstantaneousTimeMirrorManager() : isEnabled(false){};

  void init(double velocityScalingFactor,
            double triggerTime,
            seissol::geometry::MeshReader* meshReader,
            initializers::LTSTree* ltsTree,
            initializers::LTS* lts,
            initializers::Lut* ltsLut,
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
    initializers::LTSTree* ltsTree,
    initializers::LTS* lts,
    initializers::Lut* ltsLut,
    InstantaneousTimeMirrorManager& increaseManager,
    InstantaneousTimeMirrorManager& decreaseManager,
    const TimeStepping* timestepping); // An empty timestepping is added. Need to discuss what
                                       // exactly is to be sent here
} // namespace seissol::ITM

#endif // SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H
