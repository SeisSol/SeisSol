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

#ifndef SEISSOL_SRC_SOLVER_PIPELINE_GENERICPIPELINE_H_
#define SEISSOL_SRC_SOLVER_PIPELINE_GENERICPIPELINE_H_

#include "Monitoring/Stopwatch.h"
#include <limits>
#include <array>
#include <cassert>


namespace seissol {

  template<unsigned NumStagesP, unsigned DefaultBatchSizeP>
  class PipelineTuner {
  public:
    virtual ~PipelineTuner() = default;
    constexpr static decltype(NumStagesP) NumStages{NumStagesP};
    constexpr static decltype(DefaultBatchSizeP) DefaultBatchSize{DefaultBatchSizeP};

    virtual void tune(const std::array<double, NumStages>& stageTiming) {
      /**no default implementation provided**/
    };
    size_t getBatchSize() {
      const auto batchSize = static_cast<size_t>(currBatchSize);
      assert(batchSize != 0 && "error: autotuner derived zero-length batched size");
      return batchSize;
    }

  protected:
    double currBatchSize{static_cast<double>(DefaultBatchSize)};
  };


  template<unsigned NumStagesP, unsigned DefaultBatchSizeP, typename TunerT = PipelineTuner<NumStagesP, DefaultBatchSizeP>>
  class GenericPipeline {
  public:
    struct PipelineCallBack {
      virtual ~PipelineCallBack() = default;
      virtual void operator()(size_t begin, size_t batchSize, size_t callCounter) = 0;
      virtual void finalize() = 0;
    };

    GenericPipeline(bool resetAfterRun = true) : resetAfterRun(resetAfterRun) {
      using BaseTunerT = PipelineTuner<NumStagesP, DefaultBatchSizeP>;
      static_assert(std::is_base_of_v<BaseTunerT, TunerT>,
                    "ConcreteTunerT must be derived from PipelineTuner");
    }
    ~GenericPipeline() = default;
    constexpr static decltype(NumStagesP) NumStages{NumStagesP};
    constexpr static decltype(NumStagesP) TailSize{NumStagesP - 1};
    constexpr static decltype(DefaultBatchSizeP) DefaultBatchSize{DefaultBatchSizeP};

    void registerCallBack(unsigned id, PipelineCallBack* callBack) {
      assert(id < NumStages);
      callBacks[id] = callBack;
    }

    void run(size_t size) {
      init(size);
      fill();
      run();
      drain();
      clean();

      for (auto& time: stageTiming)
        time /= numTotalIterations;
      tuner.tune(stageTiming);
    }

  private:
    struct Range {
      Range() = default;
      Range(const Range& other) = default;
      explicit Range(size_t step, size_t limit) :  begin(0), end(step), stepSize(step), limit(limit) {
        clamp();
      }
      size_t begin{0};
      size_t end{0};
      size_t stepSize{0};
      size_t limit{std::numeric_limits<size_t>::max()};
      size_t currentIteration{0};
      size_t size() {
        assert(end > begin);
        return (end - begin);
      }
      Range& operator++(int) {
        this->begin += stepSize;
        this->end += stepSize;
        ++currentIteration;
        clamp();
        return *this;
      }

      void clamp() {
        begin = (begin > limit) ? limit : begin;
        end = (end > limit) ? limit : end;
      }
    };

    void init(size_t size) {
      const auto currBatchSize = tuner.getBatchSize();
      for (auto& range: ranges) {
        range = Range(currBatchSize, size);
      }
      numTotalIterations = (size + currBatchSize - 1) / currBatchSize;
      numFullPipeIterations = (numTotalIterations > TailSize) ? numTotalIterations - TailSize : 0;
      stageTiming = decltype(stageTiming){};
    }

    void fill() {
      for (unsigned i = 0; i < TailSize; ++i) {
        for (unsigned stage = 0; stage < (i + 1); ++stage) {
          if (ranges[stage].currentIteration < numTotalIterations) {
            execute(stage, ranges[stage]);
            ranges[stage]++;
            checkToFinalize(stage, ranges[stage].currentIteration);
          }
        }
      }
    }

    void run() {
      // reduce numFullPipeIterations by 1 to handle a corner case in the `drain` stage
      const size_t length = (numFullPipeIterations > 0) ? numFullPipeIterations - 1 : 0;
      for (size_t i = 0; i < length; ++i) {
        for (unsigned stage = 0; stage < NumStages; ++stage) {
          execute(stage, ranges[stage]);
          ranges[stage]++;
        }
      }
    }

    void drain() {
      for (unsigned i = 0; i < TailSize + 1; ++i) {
        for (unsigned stage = 0; stage < NumStages; ++stage) {
          if (ranges[stage].currentIteration < numTotalIterations) {
            execute(stage, ranges[stage]);
            ranges[stage]++;
            checkToFinalize(stage, ranges[stage].currentIteration);
          }
        }
      }
    }

    void clean() {
      if (resetAfterRun) {
        for (auto& callBack: callBacks) {
          callBack = nullptr;
        }
      }
    }

    void execute(unsigned stage, Range range) {
      assert(callBacks[stage] != nullptr && "call back has not been setup");
      stopwatch.start();
      (*callBacks[stage])(range.begin, range.size(), range.currentIteration);
      stageTiming[stage] += stopwatch.stop();
    }

    void checkToFinalize(unsigned stage, size_t currentIteration) {
      assert(callBacks[stage] != nullptr && "call back has not been setup");
      if (currentIteration == numTotalIterations) {
        stopwatch.start();
        callBacks[stage]->finalize();
        stageTiming[stage] += stopwatch.stop();
      }
    }

    size_t numTotalIterations{0};
    size_t numFullPipeIterations{0};
    std::array<Range, NumStages> ranges;
    std::array<double, NumStages> stageTiming{};
    std::array<PipelineCallBack*, NumStages> callBacks{nullptr};
    bool resetAfterRun{true};
    TunerT tuner;
    Stopwatch stopwatch{};
  };
}


#endif // SEISSOL_SRC_SOLVER_PIPELINE_GENERICPIPELINE_H_

