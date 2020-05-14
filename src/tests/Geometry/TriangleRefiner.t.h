#include <array>
#include <iostream>
#include <iomanip>

#include <cxxtest/TestSuite.h>
#include <Eigen/Dense>

#include "Geometry/refinement/TriangleRefiner.h"


namespace seissol {
  namespace unit_test {
    class TriangleRefinerTestSuite;
  }
}

class seissol::unit_test::TriangleRefinerTestSuite : public CxxTest::TestSuite
{
  public:
    //We do all tests in double precision
    const double epsilon = std::numeric_limits<double>::epsilon();

    void testDivideBy1() {
      std::array<std::array<Eigen::Vector2d,3>,4> expetectedTriangles = {{
        {Eigen::Vector2d(0.0,0.0), Eigen::Vector2d(0.5,0.0), Eigen::Vector2d(0.0,0.5)},
        {Eigen::Vector2d(0.5,0.0), Eigen::Vector2d(1.0,0.0), Eigen::Vector2d(0.5,0.5)},
        {Eigen::Vector2d(0.0,0.5), Eigen::Vector2d(0.5,0.5), Eigen::Vector2d(0.0,1.0)},
        {Eigen::Vector2d(0.5,0.5), Eigen::Vector2d(0.0,0.5), Eigen::Vector2d(0.5,0.0)}
      }};
      const double area = 0.25; 

      seissol::refinement::TriangleRefiner tr;
      tr.refine(1);
      
      for (unsigned i = 0; i < 4; i ++){
        assertTriangle(tr.subTris.at(i), expetectedTriangles[i], area);
      }
    }

    void testDivideBy3() {
      std::array<std::array<Eigen::Vector2d,3>,64> expetectedTriangles = {{
        {Eigen::Vector2d(    0,    0), Eigen::Vector2d(0.125,    0), Eigen::Vector2d(    0,0.125)},
        {Eigen::Vector2d(0.125,    0), Eigen::Vector2d( 0.25,    0), Eigen::Vector2d(0.125,0.125)},
        {Eigen::Vector2d(    0,0.125), Eigen::Vector2d(0.125,0.125), Eigen::Vector2d(    0, 0.25)},
        {Eigen::Vector2d(0.125,0.125), Eigen::Vector2d(    0,0.125), Eigen::Vector2d(0.125,    0)},
        {Eigen::Vector2d( 0.25,    0), Eigen::Vector2d(0.375,    0), Eigen::Vector2d( 0.25,0.125)},
        {Eigen::Vector2d(0.375,    0), Eigen::Vector2d(  0.5,    0), Eigen::Vector2d(0.375,0.125)},
        {Eigen::Vector2d( 0.25,0.125), Eigen::Vector2d(0.375,0.125), Eigen::Vector2d( 0.25, 0.25)},
        {Eigen::Vector2d(0.375,0.125), Eigen::Vector2d( 0.25,0.125), Eigen::Vector2d(0.375,    0)},
        {Eigen::Vector2d(    0, 0.25), Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d(    0,0.375)},
        {Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d( 0.25, 0.25), Eigen::Vector2d(0.125,0.375)},
        {Eigen::Vector2d(    0,0.375), Eigen::Vector2d(0.125,0.375), Eigen::Vector2d(    0,  0.5)},
        {Eigen::Vector2d(0.125,0.375), Eigen::Vector2d(    0,0.375), Eigen::Vector2d(0.125, 0.25)},
        {Eigen::Vector2d( 0.25, 0.25), Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d( 0.25,0.125)},
        {Eigen::Vector2d(0.125, 0.25), Eigen::Vector2d(    0, 0.25), Eigen::Vector2d(0.125,0.125)},
        {Eigen::Vector2d( 0.25,0.125), Eigen::Vector2d(0.125,0.125), Eigen::Vector2d( 0.25,    0)},
        {Eigen::Vector2d(0.125,0.125), Eigen::Vector2d( 0.25,0.125), Eigen::Vector2d(0.125, 0.25)},
        {Eigen::Vector2d(  0.5,    0), Eigen::Vector2d(0.625,    0), Eigen::Vector2d(  0.5,0.125)},
        {Eigen::Vector2d(0.625,    0), Eigen::Vector2d( 0.75,    0), Eigen::Vector2d(0.625,0.125)},
        {Eigen::Vector2d(  0.5,0.125), Eigen::Vector2d(0.625,0.125), Eigen::Vector2d(  0.5, 0.25)},
        {Eigen::Vector2d(0.625,0.125), Eigen::Vector2d(  0.5,0.125), Eigen::Vector2d(0.625,    0)},
        {Eigen::Vector2d( 0.75,    0), Eigen::Vector2d(0.875,    0), Eigen::Vector2d( 0.75,0.125)},
        {Eigen::Vector2d(0.875,    0), Eigen::Vector2d(    1,    0), Eigen::Vector2d(0.875,0.125)},
        {Eigen::Vector2d( 0.75,0.125), Eigen::Vector2d(0.875,0.125), Eigen::Vector2d( 0.75, 0.25)},
        {Eigen::Vector2d(0.875,0.125), Eigen::Vector2d( 0.75,0.125), Eigen::Vector2d(0.875,    0)},
        {Eigen::Vector2d(  0.5, 0.25), Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d(  0.5,0.375)},
        {Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d( 0.75, 0.25), Eigen::Vector2d(0.625,0.375)},
        {Eigen::Vector2d(  0.5,0.375), Eigen::Vector2d(0.625,0.375), Eigen::Vector2d(  0.5,  0.5)},
        {Eigen::Vector2d(0.625,0.375), Eigen::Vector2d(  0.5,0.375), Eigen::Vector2d(0.625, 0.25)},
        {Eigen::Vector2d( 0.75, 0.25), Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d( 0.75,0.125)},
        {Eigen::Vector2d(0.625, 0.25), Eigen::Vector2d(  0.5, 0.25), Eigen::Vector2d(0.625,0.125)},
        {Eigen::Vector2d( 0.75,0.125), Eigen::Vector2d(0.625,0.125), Eigen::Vector2d( 0.75,    0)},
        {Eigen::Vector2d(0.625,0.125), Eigen::Vector2d( 0.75,0.125), Eigen::Vector2d(0.625, 0.25)},
        {Eigen::Vector2d(    0,  0.5), Eigen::Vector2d(0.125,  0.5), Eigen::Vector2d(    0,0.625)},
        {Eigen::Vector2d(0.125,  0.5), Eigen::Vector2d( 0.25,  0.5), Eigen::Vector2d(0.125,0.625)},
        {Eigen::Vector2d(    0,0.625), Eigen::Vector2d(0.125,0.625), Eigen::Vector2d(    0, 0.75)},
        {Eigen::Vector2d(0.125,0.625), Eigen::Vector2d(    0,0.625), Eigen::Vector2d(0.125,  0.5)},
        {Eigen::Vector2d( 0.25,  0.5), Eigen::Vector2d(0.375,  0.5), Eigen::Vector2d( 0.25,0.625)},
        {Eigen::Vector2d(0.375,  0.5), Eigen::Vector2d(  0.5,  0.5), Eigen::Vector2d(0.375,0.625)},
        {Eigen::Vector2d( 0.25,0.625), Eigen::Vector2d(0.375,0.625), Eigen::Vector2d( 0.25, 0.75)},
        {Eigen::Vector2d(0.375,0.625), Eigen::Vector2d( 0.25,0.625), Eigen::Vector2d(0.375,  0.5)},
        {Eigen::Vector2d(    0, 0.75), Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d(    0,0.875)},
        {Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d( 0.25, 0.75), Eigen::Vector2d(0.125,0.875)},
        {Eigen::Vector2d(    0,0.875), Eigen::Vector2d(0.125,0.875), Eigen::Vector2d(    0,    1)},
        {Eigen::Vector2d(0.125,0.875), Eigen::Vector2d(    0,0.875), Eigen::Vector2d(0.125, 0.75)},
        {Eigen::Vector2d( 0.25, 0.75), Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d( 0.25,0.625)},
        {Eigen::Vector2d(0.125, 0.75), Eigen::Vector2d(    0, 0.75), Eigen::Vector2d(0.125,0.625)},
        {Eigen::Vector2d( 0.25,0.625), Eigen::Vector2d(0.125,0.625), Eigen::Vector2d( 0.25,  0.5)},
        {Eigen::Vector2d(0.125,0.625), Eigen::Vector2d( 0.25,0.625), Eigen::Vector2d(0.125, 0.75)},
        {Eigen::Vector2d(  0.5,  0.5), Eigen::Vector2d(0.375,  0.5), Eigen::Vector2d(  0.5,0.375)},
        {Eigen::Vector2d(0.375,  0.5), Eigen::Vector2d( 0.25,  0.5), Eigen::Vector2d(0.375,0.375)},
        {Eigen::Vector2d(  0.5,0.375), Eigen::Vector2d(0.375,0.375), Eigen::Vector2d(  0.5, 0.25)},
        {Eigen::Vector2d(0.375,0.375), Eigen::Vector2d(  0.5,0.375), Eigen::Vector2d(0.375,  0.5)},
        {Eigen::Vector2d( 0.25,  0.5), Eigen::Vector2d(0.125,  0.5), Eigen::Vector2d( 0.25,0.375)},
        {Eigen::Vector2d(0.125,  0.5), Eigen::Vector2d(    0,  0.5), Eigen::Vector2d(0.125,0.375)},
        {Eigen::Vector2d( 0.25,0.375), Eigen::Vector2d(0.125,0.375), Eigen::Vector2d( 0.25, 0.25)},
        {Eigen::Vector2d(0.125,0.375), Eigen::Vector2d( 0.25,0.375), Eigen::Vector2d(0.125,  0.5)},
        {Eigen::Vector2d(  0.5, 0.25), Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d(  0.5,0.125)},
        {Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d( 0.25, 0.25), Eigen::Vector2d(0.375,0.125)},
        {Eigen::Vector2d(  0.5,0.125), Eigen::Vector2d(0.375,0.125), Eigen::Vector2d(  0.5,    0)},
        {Eigen::Vector2d(0.375,0.125), Eigen::Vector2d(  0.5,0.125), Eigen::Vector2d(0.375, 0.25)},
        {Eigen::Vector2d( 0.25, 0.25), Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d( 0.25,0.375)},
        {Eigen::Vector2d(0.375, 0.25), Eigen::Vector2d(  0.5, 0.25), Eigen::Vector2d(0.375,0.375)},
        {Eigen::Vector2d( 0.25,0.375), Eigen::Vector2d(0.375,0.375), Eigen::Vector2d( 0.25,  0.5)},
        {Eigen::Vector2d(0.375,0.375), Eigen::Vector2d( 0.25,0.375), Eigen::Vector2d(0.375, 0.25)}
      }};
      const double area = 0.015625; 

      seissol::refinement::TriangleRefiner tr;
      tr.refine(3);
      
      for (unsigned i = 0; i < 64; i ++){
        assertTriangle(tr.subTris.at(i), expetectedTriangles[i], area);
      }
    }
     
    void assertTriangle(seissol::refinement::Triangle& a, 
        std::array<Eigen::Vector2d, 3>& b, double area) {
      TS_ASSERT_DELTA(a.area, area, epsilon);
      for (int i = 0; i < 3; i++) {
        TS_ASSERT_DELTA(a.x[i].x, b[i][0], epsilon);
        TS_ASSERT_DELTA(a.x[i].y, b[i][1], epsilon);
      }
    }

}; 
