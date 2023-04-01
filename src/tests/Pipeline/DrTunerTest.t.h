#include "Solver/Pipeline/DrTuner.h"
#include <array>

namespace seissol::unit_test {

TEST_CASE("Dr tuner") {
  constexpr static size_t ComputeStageId{1};
  std::array<double, 3> timing{};
  constexpr static double eps{2.0};
  size_t batchSize{0};
  dr::pipeline::DrPipelineTuner tuner;

  SUBCASE("Goes to left") {
    // this is going to results in: Performance ~ 1 / batchSize
    auto squareFunction = [](size_t x) { return static_cast<double>(x * x); };

    while (!tuner.isTunerConverged()) {
      batchSize = tuner.getBatchSize();
      timing[ComputeStageId] = squareFunction(batchSize);
      tuner.tune(timing);
    }
    REQUIRE(batchSize == AbsApprox(tuner.getMinBatchSize()).epsilon(eps));
  }

  SUBCASE("Goes to right") {
    auto hyperbolicTime = [](size_t x) { return 1.0 / (static_cast<double>(x + 1.0)); };

    while (!tuner.isTunerConverged()) {
      batchSize = tuner.getBatchSize();
      timing[ComputeStageId] = hyperbolicTime(batchSize);
      tuner.tune(timing);
    }
    REQUIRE(batchSize == AbsApprox(tuner.getMinBatchSize()).epsilon(eps));
  }

  SUBCASE("Max is withing range") {
    const auto midPoint = 0.5 * (tuner.getMaxBatchSize() + tuner.getMinBatchSize());

    auto hatFunction = [midPoint](size_t x) { return std::abs(midPoint - x); };

    while (!tuner.isTunerConverged()) {
      batchSize = tuner.getBatchSize();
      timing[ComputeStageId] = hatFunction(batchSize);
      tuner.tune(timing);
    }
    REQUIRE(batchSize == AbsApprox(midPoint).epsilon(eps));
  }
}
} // namespace seissol::unit_test
