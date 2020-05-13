#include <array>

#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

#include "MockReader.h"
#include "Geometry/refinement/MeshRefiner.h"
#include "Geometry/refinement/RefinerUtils.h"


namespace seissol {
  namespace unit_test {
    class RefinerTestSuite;
  }
}

class seissol::unit_test::RefinerTestSuite : public CxxTest::TestSuite
{
  public:
    //We do all tests in double precision
    const real epsilon = std::numeric_limits<double>::epsilon();
    std::array<Eigen::Vector3d, 4> vertices;
    
    RefinerTestSuite() {
      std::srand(0);
      vertices = {{
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX),
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX),
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX),
        Eigen::Vector3d((double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX, (double)std::rand()/RAND_MAX)}};
    }

    void testDivideBy4() {
      const MockReader mockReader(vertices);
      const std::array<Eigen::Vector3d, 5> exptectedVerticesDivideBy8 = {
        vertices[0], vertices[1], vertices[2], vertices[3],
        0.25*(vertices[0] + vertices[1] + vertices[2] + vertices[3]) 
      };

      const std::array<Eigen::Vector4i, 4> exptectedCellsDivideBy8 {
        Eigen::Vector4i(0, 1, 2, 4),
        Eigen::Vector4i(0, 1, 3, 4),
        Eigen::Vector4i(0, 2, 3, 4),
        Eigen::Vector4i(1, 2, 3, 4)
      };

      seissol::refinement::DivideTetrahedronBy4<double> refineBy4;
      seissol::refinement::MeshRefiner<double>meshRefiner(mockReader, refineBy4);
      TS_ASSERT_EQUALS(meshRefiner.getNumCells(), 4);
      TS_ASSERT_EQUALS(meshRefiner.getNumVertices(), 5);
      for (unsigned i = 0; i < meshRefiner.getNumVertices(); i++) {
        assertPoint(&meshRefiner.getVertexData()[3*i], exptectedVerticesDivideBy8[i]);
      }
      for (unsigned i = 0; i < meshRefiner.getNumCells(); i++) {
        assertCell(&meshRefiner.getCellData()[4*i], exptectedCellsDivideBy8[i]);
      }
    }

    void testDivideBy8() {
      const MockReader mockReader(vertices);
      const std::array<Eigen::Vector3d, 10> exptectedVerticesDivideBy8 = {
        vertices[0], vertices[1], vertices[2], vertices[3],
        0.5*(vertices[0] + vertices[1]), 
        0.5*(vertices[0] + vertices[2]), 
        0.5*(vertices[0] + vertices[3]), 
        0.5*(vertices[1] + vertices[2]), 
        0.5*(vertices[1] + vertices[3]), 
        0.5*(vertices[2] + vertices[3])
      };

      const std::array<Eigen::Vector4i, 8> exptectedCellsDivideBy8 {
        Eigen::Vector4i(0, 4, 5, 6),
        Eigen::Vector4i(1, 4, 7, 8),
        Eigen::Vector4i(2, 5, 7, 9),
        Eigen::Vector4i(3, 6, 8, 9),
        Eigen::Vector4i(4, 5, 6, 8),
        Eigen::Vector4i(4, 5, 7, 8),
        Eigen::Vector4i(5, 6, 8, 9),
        Eigen::Vector4i(5, 7, 8, 9)
      };

      seissol::refinement::DivideTetrahedronBy8<double> refineBy8;
      seissol::refinement::MeshRefiner<double>meshRefiner(mockReader, refineBy8);
      TS_ASSERT_EQUALS(meshRefiner.getNumCells(), 8);
      TS_ASSERT_EQUALS(meshRefiner.getNumVertices(), 10);
      for (unsigned i = 0; i < meshRefiner.getNumVertices(); i++) {
        assertPoint(&meshRefiner.getVertexData()[3*i], exptectedVerticesDivideBy8[i]);
      }
      for (unsigned i = 0; i < meshRefiner.getNumCells(); i++) {
        assertCell(&meshRefiner.getCellData()[4*i], exptectedCellsDivideBy8[i]);
      }
    }

    void testDivideBy32() {
      const MockReader mockReader(vertices);
      const std::array<Eigen::Vector3d, 18> exptectedVerticesDivideBy32 = {
        vertices[0], vertices[1], vertices[2], vertices[3],
        0.5*(vertices[0] + vertices[1]), 
        0.5*(vertices[0] + vertices[2]), 
        0.5*(vertices[0] + vertices[3]), 
        0.5*(vertices[1] + vertices[2]), 
        0.5*(vertices[1] + vertices[3]), 
        0.5*(vertices[2] + vertices[3]), 

        0.625*vertices[0] + 0.125*vertices[1] + 0.125*vertices[2] + 0.125*vertices[3], 
        0.125*vertices[0] + 0.625*vertices[1] + 0.125*vertices[2] + 0.125*vertices[3], 
        0.125*vertices[0] + 0.125*vertices[1] + 0.625*vertices[2] + 0.125*vertices[3], 
        0.125*vertices[0] + 0.125*vertices[1] + 0.125*vertices[2] + 0.625*vertices[3], 

        0.375*vertices[0] + 0.25 *vertices[1] + 0.125*vertices[2] + 0.25 *vertices[3], 
        0.25 *vertices[0] + 0.375*vertices[1] + 0.25 *vertices[2] + 0.125*vertices[3], 
        0.25 *vertices[0] + 0.125*vertices[1] + 0.25 *vertices[2] + 0.375*vertices[3],
        0.125*vertices[0] + 0.25 *vertices[1] + 0.375*vertices[2] + 0.25 *vertices[3]
      };

      const std::array<Eigen::Vector4i, 32> exptectedCellsDivideBy32 {
        Eigen::Vector4i(0,4,5,10),
        Eigen::Vector4i(0,4,6,10),
        Eigen::Vector4i(0,5,6,10),
        Eigen::Vector4i(4,5,6,10),
        Eigen::Vector4i(1,4,7,11),
        Eigen::Vector4i(1,4,8,11),
        Eigen::Vector4i(1,7,8,11),
        Eigen::Vector4i(4,7,8,11),
        Eigen::Vector4i(2,5,7,12),
        Eigen::Vector4i(2,5,9,12),
        Eigen::Vector4i(2,7,9,12),
        Eigen::Vector4i(5,7,9,12),
        Eigen::Vector4i(3,6,8,13),
        Eigen::Vector4i(3,6,9,13),
        Eigen::Vector4i(3,8,9,13),
        Eigen::Vector4i(6,8,9,13),
        Eigen::Vector4i(4,5,6,14),
        Eigen::Vector4i(4,5,8,14),
        Eigen::Vector4i(4,6,8,14),
        Eigen::Vector4i(5,6,8,14),
        Eigen::Vector4i(4,5,7,15),
        Eigen::Vector4i(4,5,8,15),
        Eigen::Vector4i(4,7,8,15),
        Eigen::Vector4i(5,7,8,15),
        Eigen::Vector4i(5,6,8,16),
        Eigen::Vector4i(5,6,9,16),
        Eigen::Vector4i(5,8,9,16),
        Eigen::Vector4i(6,8,9,16),
        Eigen::Vector4i(5,7,8,17),
        Eigen::Vector4i(5,7,9,17),
        Eigen::Vector4i(5,8,9,17),
        Eigen::Vector4i(7,8,9,17)
      };

      seissol::refinement::DivideTetrahedronBy32<double> refineBy32;
      seissol::refinement::MeshRefiner<double>meshRefiner(mockReader, refineBy32);
      TS_ASSERT_EQUALS(meshRefiner.getNumCells(), 32);
      TS_ASSERT_EQUALS(meshRefiner.getNumVertices(), 18);
      for (unsigned i = 0; i < meshRefiner.getNumVertices(); i++) {
        assertPoint(&meshRefiner.getVertexData()[3*i], exptectedVerticesDivideBy32[i]);
      }
      for (unsigned i = 0; i < meshRefiner.getNumCells(); i++) {
        assertCell(&meshRefiner.getCellData()[4*i], exptectedCellsDivideBy32[i]);
      }
    }

    void assertPoint(const double* a, const Eigen::Vector3d& b) {
      for (int i = 0; i < 3; i++) {
        TS_ASSERT_DELTA(a[i], b[i], epsilon);
      }
    }

    void assertCell(const unsigned int* a, const Eigen::Vector4i& b) {
      for (int i = 0; i < 4; i++) {
        TS_ASSERT_EQUALS(a[i], b[i]);
      }
    }
}; 
