#include <cxxtest/TestSuite.h>
#include "Solver/Pipeline/DrTuner.h"
#include <array>


namespace seissol {
  namespace unit_test {
    class DrTunerTest;
  }
}


class seissol::unit_test::DrTunerTest : public CxxTest::TestSuite {
public:
  void testGoesToLeft() {
    dr::pipeline::DrPipelineTuner tuner;
    // this is going to results in: Performance ~ 1 / batchSize
    auto squareFunction = [] (size_t x) {return static_cast<double>(x * x);};

    while (!tuner.isTunerConverged()) {
      batchSize = tuner.getBatchSize();
      timing[ComputeStageId] = squareFunction(batchSize);
      tuner.tune(timing);
    }
    TS_ASSERT_DELTA(batchSize, tuner.getMinBatchSize(), eps);
  }

  void testGoesToRight() {
    dr::pipeline::DrPipelineTuner tuner;
    auto hyperbolicTime = [](size_t x) {return 1.0 / (static_cast<double>(x + 1.0));};

    while (!tuner.isTunerConverged()) {
      batchSize = tuner.getBatchSize();
      timing[ComputeStageId] = hyperbolicTime(batchSize);
      tuner.tune(timing);
    }
    TS_ASSERT_DELTA(batchSize, tuner.getMaxBatchSize(), eps);
  }

  void testMaxWithinRange() {
    dr::pipeline::DrPipelineTuner tuner;
    const auto midPoint = 0.5 * (tuner.getMaxBatchSize() + tuner.getMinBatchSize());

    auto hatFunction = [midPoint] (size_t x) {
      return std::abs(midPoint - x);
    };

    while (!tuner.isTunerConverged()) {
      batchSize = tuner.getBatchSize();
      timing[ComputeStageId] = hatFunction(batchSize);
      tuner.tune(timing);
    }
    TS_ASSERT_DELTA(batchSize, midPoint, eps);
  }

private:
  constexpr static size_t ComputeStageId{1};
  std::array<double, 3> timing{};
  constexpr static double eps{2.0};
  size_t batchSize{0};
};
