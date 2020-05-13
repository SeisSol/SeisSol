#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>
#include <glm/vec3.hpp>

#include <Numerical_aux/Transformation.h>

namespace seissol {
  namespace unit_test {
    class TransformationTestSuite;
  }
}

class seissol::unit_test::TransformationTestSuite : public CxxTest::TestSuite
{
public:
    //We do all tests in double precision
    const real epsilon = std::numeric_limits<double>::epsilon();
    std::array<Eigen::Vector3d, 4> vertices;
    
    TransformationTestSuite() {
      std::srand(9);
      vertices = {{
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX),
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX),
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX),
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX)}};
    }
	void testTetrahedronGlobalToRefernce()
	{

          auto center = 0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3]);
          glm::dvec3 xyz;
          xyz.x = center[0];
          xyz.y = center[1];
          xyz.z = center[2];

          auto res = seissol::transformations::tetrahedronGlobalToReference(  
              vertices[0].data(),
              vertices[1].data(),
              vertices[2].data(),
              vertices[3].data(),
              xyz );
          TS_ASSERT_DELTA(res.x, 0.25, epsilon);
          TS_ASSERT_DELTA(res.y, 0.25, epsilon);
          TS_ASSERT_DELTA(res.z, 0.25, epsilon);
        }
};
