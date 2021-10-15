#include <cxxtest/TestSuite.h>

#include "Geometry/PUMLReader.h"
#include "Initializer/time_stepping/LtsWeights/WeightsModels.h"
#include <memory>

namespace seissol {
  namespace unit_test {
    class LTSWeightsTestSuite;
  }
}

class seissol::unit_test::LTSWeightsTestSuite : public CxxTest::TestSuite
{
  public:
    const double epsilon = std::numeric_limits<double>::epsilon();

    void testTimeStepsState()
    {
// PUMLReader is only available with MPI
#ifdef USE_MPI
      std::cout.setstate(std::ios_base::failbit);
      using namespace seissol::initializers::time_stepping;
      LtsWeightsConfig config{"Testing/material.yaml", 2, 1, 1, 1 };

      auto ltsWeights = std::make_unique<ExponentialWeights>(config);
      PUMLReader pumlReader("Testing/mesh.h5", 5000.0, "", ltsWeights.get());
      std::cout.clear();

      std::array<unsigned, 24> expectedWeights = {
        2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1
      };
      for (int i = 0; i < 24; i++) {
        TS_ASSERT_EQUALS(ltsWeights->vertexWeights()[i], expectedWeights[i]);
      }
#endif
    }

};
