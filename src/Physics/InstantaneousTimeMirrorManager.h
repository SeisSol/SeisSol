#ifndef SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H
#define SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H

#include "Geometry/MeshReader.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Lut.h"
#include "Initializer/Typedefs.h"
#include "Modules/Module.h"
#include "Solver/time_stepping/AbstractGhostTimeCluster.h"
#include "Solver/time_stepping/TimeCluster.h"

namespace seissol {
class SeisSol;
namespace ITM {

class InstantaneousTimeMirrorManager : Module {
  seissol::SeisSol& seissolInstance;
  bool isEnabled{false};
  double velocityScalingFactor{};
  double triggerTime{};

  seissol::geometry::MeshReader* meshReader{};
  initializer::LTSTree* ltsTree{};
  initializer::LTS* lts{};
  initializer::Lut* ltsLut{};
  const TimeStepping* timestepping{};

  std::vector<std::unique_ptr<seissol::time_stepping::TimeCluster>>* timeClusters;
  std::vector<std::unique_ptr<seissol::time_stepping::AbstractGhostTimeCluster>>* ghostTimeClusters;

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

#endif // SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H
