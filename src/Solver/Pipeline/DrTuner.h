// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Ravil Dorozhinskii (ravil.dorozhinskii AT tum.de)
 *
 */

#ifndef SEISSOL_SRC_SOLVER_PIPELINE_DRTUNER_H_
#define SEISSOL_SRC_SOLVER_PIPELINE_DRTUNER_H_

#include "Solver/Pipeline/GenericPipeline.h"


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


#endif // SEISSOL_SRC_SOLVER_PIPELINE_DRTUNER_H_

