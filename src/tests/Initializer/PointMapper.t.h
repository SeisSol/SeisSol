#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

#include "tests/Geometry/MockReader.h"
#include "Initializer/PointMapper.h"

namespace unit_tests {
  class PointMapperTestSuite;
}

class unit_tests::PointMapperTestSuite: public CxxTest::TestSuite {
  public:
    //We do all tests in double precision
    const double epsilon = std::numeric_limits<double>::epsilon();
    std::array<Eigen::Vector3d, 4> vertices;
    
    PointMapperTestSuite() {
      std::srand(321);
      vertices = {{
        Eigen::Vector3d(0.0, 0.0, 0.0),
        Eigen::Vector3d(1.0, 0.0, 0.0),
        Eigen::Vector3d(0.0, 1.0, 0.0),
        Eigen::Vector3d(0.0, 0.0, 1.0)}};
    }

    void testFindMeshId() {
      const seissol::MockReader mockReader(vertices);

      const Eigen::Vector3d points[3] = {
        Eigen::Vector3d((double)std::rand()/RAND_MAX,(double)std::rand()/RAND_MAX,(double)std::rand()/RAND_MAX),
        Eigen::Vector3d((double)std::rand()/RAND_MAX,(double)std::rand()/RAND_MAX,(double)std::rand()/RAND_MAX),
        0.25*(vertices[0] + vertices[1] + vertices[2] + vertices[3])
      };
      short contained[3] = {0, 0, 0};
      unsigned meshId[3] = {
        std::numeric_limits<unsigned>::max(),
        std::numeric_limits<unsigned>::max(),
        std::numeric_limits<unsigned>::max()
      };
      seissol::initializers::findMeshIds( points, mockReader, 3, contained, meshId);

      std::array<short, 3> expectedContained = {0, 0, 1};
      std::array<unsigned, 3> expectedMeshId = {
        std::numeric_limits<unsigned>::max(),
        std::numeric_limits<unsigned>::max(),
        0
      };

      for (int i = 0; i < 3; i++) {
        TS_ASSERT_EQUALS(contained[i], expectedContained[i]);
        TS_ASSERT_EQUALS(meshId[i], expectedMeshId[i]);
      }
  };
};
