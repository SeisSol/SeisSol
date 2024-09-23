// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Ravil Dorozhinskii (ravil.dorozhinskii AT tum.de)
 *
 */

#ifndef SEISSOL_SRC_SOLVER_PIPELINE_DRPIPELINE_H_
#define SEISSOL_SRC_SOLVER_PIPELINE_DRPIPELINE_H_

#include "Solver/Pipeline/GenericPipeline.h"
#include "Solver/Pipeline/DrTuner.h"
#include "generated_code/tensor.h"
#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::time_stepping {
  class TimeCluster;
}

namespace seissol::dr::pipeline {
  using DrPipeline = seissol::GenericPipeline<3, 1024, DrPipelineTuner>;

  struct DrContext {
    using QInterpolatedPtrT = real (*)[ConvergenceOrder][tensor::QInterpolated::size()];
    using imposedStatePlusT = real (*)[tensor::QInterpolated::size()];

    real* QInterpolatedPlusOnDevice{nullptr};
    real* QInterpolatedMinusOnDevice{nullptr};
    QInterpolatedPtrT QInterpolatedPlusOnHost{nullptr};
    QInterpolatedPtrT QInterpolatedMinusOnHost{nullptr};
    imposedStatePlusT imposedStatePlusOnHost{nullptr};
    imposedStatePlusT imposedStateMinusOnHost{nullptr};
    DRFaceInformation* faceInformation{nullptr};
    real (*devImposedStatePlus)[tensor::QInterpolated::size()]{nullptr};
    real (*devImposedStateMinus)[tensor::QInterpolated::size()]{nullptr};
    model::IsotropicWaveSpeeds* waveSpeedsPlus{nullptr};
    model::IsotropicWaveSpeeds* waveSpeedsMinus{nullptr};
  };

  struct DrBaseCallBack : public DrPipeline::PipelineCallBack {
    explicit DrBaseCallBack(DrContext userContext, time_stepping::TimeCluster *cluster)
        : context(userContext), cluster(cluster) {}
  protected:
    DrContext context{};
    time_stepping::TimeCluster *cluster{nullptr};
#ifdef ACL_DEVICE
    device::DeviceInstance &device = device::DeviceInstance::getInstance();
#endif
  };
}



#endif // SEISSOL_SRC_SOLVER_PIPELINE_DRPIPELINE_H_

