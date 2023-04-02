#ifndef SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H
#define SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H

#include <Geometry/MeshReader.h>
#include <Initializer/LTS.h>
#include "Modules/Module.h"
#include "Initializer/tree/LTSTree.hpp"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/typedefs.hpp"

namespace seissol {

class InstantaneousTimeMirrorManager : Module {
  bool isEnabled;
  double velocityScalingFactor{};
  double triggerTime{};

  MeshReader* meshReader{};
  initializers::LTSTree* ltsTree{};
  initializers::LTS* lts{};
  initializers::Lut* ltsLut{};
  TimeStepping* timestepping;
  
public:
  InstantaneousTimeMirrorManager() : isEnabled(false) {};

  void init(double velocityScalingFactor,
            double triggerTime,
            MeshReader* meshReader,
            initializers::LTSTree* ltsTree,
            initializers::LTS* lts,
            initializers::Lut* ltsLut,
            TimeStepping* timestepping);

  void syncPoint(double currentTime) override;

private:
  void updateVelocities();
  
};

void initializeTimeMirrorManagers(
    double scalingFactor, double triggerTime,
    MeshReader* meshReader,
    initializers::LTSTree* ltsTree,
    initializers::LTS* lts,
    initializers::Lut* ltsLut,
    InstantaneousTimeMirrorManager& increaseManager,
    InstantaneousTimeMirrorManager& decreaseManager,
    TimeStepping* timestepping);
} // namespace SeisSol


#endif //SEISSOL_INSTANTANEOUSTIMEMIRRORMANAGER_H
