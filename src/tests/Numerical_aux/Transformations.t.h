#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

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

          auto res = seissol::transformations::tetrahedronGlobalToReference(  
              vertices[0].data(),
              vertices[1].data(),
              vertices[2].data(),
              vertices[3].data(),
              center );
          TS_ASSERT_DELTA(res(0), 0.25, epsilon*10);
          TS_ASSERT_DELTA(res(1), 0.25, epsilon*10);
          TS_ASSERT_DELTA(res(2), 0.25, epsilon*10);
        }
};
