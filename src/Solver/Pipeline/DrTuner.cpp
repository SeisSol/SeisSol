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
 * A custom tuner of DR pipeline which is based on the golden bisection method
 *
 * Note, this can be deprecated once friction solvers are adapted for GPU computing
 **/

#include "DrTuner.h"
#include <cmath>
#include <iostream>

namespace seissol::dr::pipeline {

  DrPipelineTuner::DrPipelineTuner() {
    assert(minBatchSize > 1.0 && "min batch size must be at least 1");

    invPhi = 0.5 * (std::sqrt(5.0) - 1);
    invPhiSquared = invPhi * invPhi;

    stepSize = maxBatchSize - minBatchSize;
    leftPoint = minBatchSize + invPhiSquared * stepSize;
    rightPoint = minBatchSize + invPhi * stepSize;

    currBatchSize = rightPoint;
    action = Action::BeginRecordingRightEvaluation;
  }

  /**
   * Implements a golden-section search to find a optimal batch size.
   *
   * The objective function is given as million DR cells updates per second which is
   * supposed to get maximized. Note, that a callee evaluates the function. Values are returned back
   * during the next call. Actions are used to handle the logic. The method sets `SkipAction` which
   * results in fixing a batch size once a convergence achieved.
   *
   * @param  stageTiming average CPU time (in seconds) step on each stage for a batch processing.
   **/
  void DrPipelineTuner::tune(const std::array<double, NumStages>& stageTiming) {
    constexpr size_t ComputeStageId{1};
    double currPerformance = 1e6 * currBatchSize / (stageTiming[ComputeStageId] + 1e-12);

    switch (action) {
      case Action::SkipAction: {
        return;
      }
      case Action::BeginRecordingRightEvaluation : {
        rightValue = currPerformance;

        // set next action to take (i.e. to evaluate and record left value)
        action = Action::BeginRecordingLeftEvaluation;
        currBatchSize = leftPoint;
        return;
      }
      case Action::BeginRecordingLeftEvaluation : {
        leftValue = currPerformance;
        break;
      }
      case Action::RecordRightEvaluation : {
        rightValue = currPerformance;
        break;
      }
      case Action::RecordLeftEvaluation : {
        leftValue = currPerformance;
        break;
      }
      default: {
        break;
      }
    }

    const auto diff = std::abs(leftPoint - rightPoint);
    if (eps > diff) {
      // convergence achieved
      action = Action::SkipAction;
      isConverged = true;
      return;
    }

    if (leftValue > rightValue) {
      maxBatchSize = rightPoint;
      rightPoint = leftPoint;
      rightValue = leftValue;
      stepSize = invPhi * stepSize;
      leftPoint = minBatchSize + invPhiSquared * stepSize;
      action = Action::RecordLeftEvaluation;
      currBatchSize = leftPoint;
    }
    else {
      minBatchSize = leftPoint;
      leftPoint = rightPoint;
      leftValue = rightValue;
      stepSize = invPhi * stepSize;
      rightPoint = minBatchSize + invPhi * stepSize;
      action = Action::RecordRightEvaluation;
      currBatchSize = rightPoint;
    }
  }
}
