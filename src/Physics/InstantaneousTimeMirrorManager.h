#ifndef SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H
#define SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H

#include "Geometry/MeshReader.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Lut.h"
#include "Initializer/Typedefs.h"
#include "Modules/Module.h"
#include <Solver/Clustering/AbstractTimeCluster.h>

namespace seissol {
class SeisSol;
namespace ITM {

class InstantaneousTimeMirrorManager : Module {
  seissol::SeisSol& seissolInstance;
  bool isEnabled;
  double velocityScalingFactor{};
  double triggerTime{};

  seissol::geometry::MeshReader* meshReader{};
  initializer::LTSTree* ltsTree{};
  initializer::LTS* lts{};
  initializer::Lut* ltsLut{};
  const TimeStepping* timestepping{};

  std::vector<std::unique_ptr<seissol::solver::clustering::AbstractTimeCluster>>* timeClusters;

  public:
  InstantaneousTimeMirrorManager(seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance), isEnabled(false) {};

  void init(double velocityScalingFactor,
            double triggerTime,
            seissol::geometry::MeshReader* meshReader,
            initializer::LTSTree* ltsTree,
            initializer::LTS* lts,
            initializer::Lut* ltsLut,
            const TimeStepping* timestepping); // An empty timestepping is added. Need to discuss
                                               // what exactly is to be sent here

  void setTimeClusterVector(
      std::vector<std::unique_ptr<seissol::solver::clustering::AbstractTimeCluster>>* timeClusters);

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
