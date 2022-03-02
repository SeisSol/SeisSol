/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Ravil Dorozhinskii (ravil.dorozhinskii AT tum.de)
 *
 * @section LICENSE
 * Copyright (c) 2020-2021, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * An implementation of a pipeline.
 *
 * Specification of Dynamic Rupture pipeline.
 **/

#ifndef DR_PIPELINE_H
#define DR_PIPELINE_H

#include <Solver/Pipeline/GenericPipeline.h>
#include <Solver/Pipeline/DrTuner.h>
#include <generated_code/tensor.h>
#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol::time_stepping {
  class TimeCluster;
}

namespace seissol::dr::pipeline {
  using DrPipeline = seissol::GenericPipeline<3, 1024, DrPipelineTuner>;

  struct DrContext {
    using QInterpolatedPtrT = real (*)[CONVERGENCE_ORDER][tensor::QInterpolated::size()];
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


#endif //DR_PIPELINE_H
