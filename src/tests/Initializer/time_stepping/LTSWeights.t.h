#include "Geometry/PUMLReader.h"
#include "Initializer/time_stepping/LtsWeights/WeightsModels.h"
#include <memory>

TEST_CASE("LTS Weights") {
// PUMLReader is only available with MPI
#ifdef USE_MPI
  std::cout.setstate(std::ios_base::failbit);
  using namespace seissol::initializers::time_stepping;
  LtsWeightsConfig config{"Testing/material.yaml", 2, 1, 1, 1 };

  auto ltsWeights = std::make_unique<ExponentialWeights>(config);
  seissol::PUMLReader pumlReader("Testing/mesh.h5", 5000.0, "", ltsWeights.get());
  std::cout.clear();

  std::array<unsigned, 24> expectedWeights = {
      2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1
  };
  for (int i = 0; i < 24; i++) {
    REQUIRE(ltsWeights->vertexWeights()[i] == expectedWeights[i]);
  }
#endif
}
