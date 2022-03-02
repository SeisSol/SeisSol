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
 * A custom pipeline tuner of DR pipeline which is based on the golden bisection method
 *
 * Note, this can be deprecated once friction solvers are adapted for GPU computing
 **/

#ifndef DR_TUNER_H
#define DR_TUNER_H

#include <Solver/Pipeline/GenericPipeline.h>


namespace seissol::dr::pipeline {
  class DrPipelineTuner: public PipelineTuner<3, 1024> {
  public:
    DrPipelineTuner();
    ~DrPipelineTuner() override = default;
    void tune(const std::array<double, NumStages>& stageTiming) override;
    [[nodiscard]] bool isTunerConverged() const {return isConverged;}
    [[nodiscard]] double getMaxBatchSize() const {return maxBatchSize;}
    [[nodiscard]] double getMinBatchSize() const {return minBatchSize;}
  private:
    enum class Action {
      BeginRecordingLeftEvaluation,
      BeginRecordingRightEvaluation,
      RecordLeftEvaluation,
      RecordRightEvaluation,
      SkipAction,
    };
    Action action{Action::BeginRecordingRightEvaluation};

    double maxBatchSize{DefaultBatchSize};
    double minBatchSize{0.1 * DefaultBatchSize};

    double stepSize{};
    double leftPoint{};
    double rightPoint{};
    double leftValue{0.0};
    double rightValue{0.0};
    double invPhi{};
    double invPhiSquared{};
    double eps{5e-2};
    bool isConverged{false};
  };
}

#endif //DR_TUNER_H
