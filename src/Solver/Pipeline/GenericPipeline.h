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
 * A user sets up a num. of stages, default batch size and callbacks. A instance of
 * a pipeline is going to automatically handle corner cases, namely: filling
 * and draining pipeline. A user defines callbacks as instance of classes derived
 * from PipelineCallBack inner class. This approach allows one to customize a
 * pipeline for her/his needs.
 **/

#ifndef GENERIC_PIPELINE_H
#define GENERIC_PIPELINE_H

#include <Monitoring/Stopwatch.h>
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

#endif //GENERIC_PIPELINE_H
