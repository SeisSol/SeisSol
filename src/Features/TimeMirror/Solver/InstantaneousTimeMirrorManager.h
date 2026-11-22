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
#include "Core/Run/Context.h"
#include "Core/TimeStepping/Actor/AbstractTimeCluster.h"
#include "Core/TimeStepping/Compute/CellCluster.h"
#include "Core/TimeStepping/Halo/GhostCluster.h"
#include "Core/TimeStepping/Layout/ClusterLayout.h"
#include "Features/TimeMirror/ITMParameters.h"
#include "WavePropagation/WavePropagationParameters.h"

#include <vector>

namespace seissol::physics {

bool isAnisotropicReflectionTypeSupported(
    seissol::initializer::parameters::ReflectionType reflectionType);
double getSwaveScaledLambda(double lambda, double mu, double velocityScalingFactor);
double
    getElasticTimeStepScalingFactor(seissol::initializer::parameters::ReflectionType reflectionType,
                                    double velocityScalingFactor);

class InstantaneousTimeMirrorManager : public Module {

  RunContext& run_;
  initializer::parameters::ITMParameters itm_;
  initializer::parameters::WavePropagationParameters wavePropagation_;
  bool isEnabled_{false};
  double velocityScalingFactor_{1.0};
  double triggerTime_{};

  seissol::geometry::MeshReader* meshReader_{nullptr};
  LTS::Storage* ltsStorage_{nullptr};
  const initializer::ClusterLayout* clusterLayout_{nullptr};

  std::vector<seissol::solver::AbstractTimeCluster*> clusters_;

  public:
  InstantaneousTimeMirrorManager(
      RunContext& run,
      initializer::parameters::ITMParameters itm,
      initializer::parameters::WavePropagationParameters wavePropagation);
  ~InstantaneousTimeMirrorManager() override;

  // delete all other constructors due to it being a module.

  InstantaneousTimeMirrorManager(const InstantaneousTimeMirrorManager&) = delete;
  auto operator=(const InstantaneousTimeMirrorManager&) = delete;
  InstantaneousTimeMirrorManager(InstantaneousTimeMirrorManager&&) = delete;
  auto operator=(InstantaneousTimeMirrorManager&&) = delete;

  /// Mirrors at the time `triggerTime`, by the factor `velocityScalingFactor`; scales the time
  /// steps of the clusters `clusters` (all clusters of the time manager, also the ghost clusters).
  void init(double velocityScalingFactor,
            double triggerTime,
            seissol::geometry::MeshReader* meshReader,
            LTS::Storage& ltsStorage,
            const initializer::ClusterLayout* clusterLayout,
            const std::vector<seissol::solver::AbstractTimeCluster*>& clusters);

  void syncPoint(double currentTime) override;

  private:
  void scaleClusterTimes(double scalingFactor);
  void updateVelocities();
  void updateTimeSteps();
};

void initializeTimeMirrorManagers(
    double scalingFactor,
    double triggerTime,
    seissol::geometry::MeshReader* meshReader,
    LTS::Storage& ltsStorage,
    InstantaneousTimeMirrorManager& increaseManager,
    InstantaneousTimeMirrorManager& decreaseManager,
    const initializer::parameters::ITMParameters& itm,
    const initializer::ClusterLayout* clusterLayout,
    const std::vector<seissol::solver::AbstractTimeCluster*>& clusters);
} // namespace seissol::physics

#endif // SEISSOL_SRC_FEATURES_TIMEMIRROR_SOLVER_INSTANTANEOUSTIMEMIRRORMANAGER_H_
