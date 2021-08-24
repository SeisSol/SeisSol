#include <cxxtest/TestSuite.h>

#include "Geometry/PUMLReader.h"
#include "Initializer/time_stepping/LtsWeights.h"

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
      std::cout.setstate(std::ios_base::failbit);
      seissol::initializers::time_stepping::LtsWeights ltsWeights("Testing/material.yaml", 2, 1, 1, 1);
      PUMLReader pumlReader("Testing/mesh.h5", "", &ltsWeights);
      std::cout.clear();

      std::array<unsigned, 24> expectedWeights = {
        2, 2, 1, 1, 1, 1, 1, 2, 1, 1, 2, 2, 2, 1, 1, 1, 1, 2, 2, 2, 1, 1, 2, 1
      };
      for (int i = 0; i < 24; i++) {
        TS_ASSERT_EQUALS(ltsWeights.vertexWeights()[i], expectedWeights[i]);
      }
    }

};
