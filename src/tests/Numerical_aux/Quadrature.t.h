#include <Numerical_aux/Quadrature.h>

namespace seissol::unit_test {

TEST_CASE("Test generation of Gauss Jacobi quadrature points") {
  constexpr auto epsilon = std::numeric_limits<double>::epsilon() * 10;
  double points[5];
  double weights[5];
  seissol::quadrature::GaussJacobi(points, weights, 5, 1, 3);
  // Compare to Maple reference solution
  CHECK(points[0] == AbsApprox(0.86698568210542769702).epsilon(epsilon));
  CHECK(points[1] == AbsApprox(0.57652877512667440772).epsilon(epsilon));
  CHECK(points[2] == AbsApprox(0.17976783188823737401).epsilon(epsilon));
  CHECK(points[3] == AbsApprox(-.25499675973326581341).epsilon(epsilon));
  CHECK(points[4] == AbsApprox(-.65399981510135937963).epsilon(epsilon));
  CHECK(weights[0] == AbsApprox(0.18915446768616357329).epsilon(epsilon));
  CHECK(weights[1] == AbsApprox(0.58714974961811369751).epsilon(epsilon));
  CHECK(weights[2] == AbsApprox(0.57657004957734461768).epsilon(epsilon));
  CHECK(weights[3] == AbsApprox(0.22255926867518051648).epsilon(epsilon));
  CHECK(weights[4] == AbsApprox(0.024566464443197594119).epsilon(epsilon));
}

TEST_CASE("test triangle quadrature") {
  constexpr auto epsilon = std::numeric_limits<double>::epsilon() * 10;

  double points[4][2];
  double weights[4];
  seissol::quadrature::TriangleQuadrature(points, weights, 2);
  // Compare to Maple reference solution
  CHECK(points[0][0] == AbsApprox(0.64494897427831780982).epsilon(epsilon));
  CHECK(points[1][0] == AbsApprox(0.64494897427831780982).epsilon(epsilon));
  CHECK(points[2][0] == AbsApprox(0.15505102572168219018).epsilon(epsilon));
  CHECK(points[3][0] == AbsApprox(0.15505102572168219018).epsilon(epsilon));
  CHECK(points[0][1] == AbsApprox(0.28001991549907407200).epsilon(epsilon));
  CHECK(points[1][1] == AbsApprox(0.075031110222608118175).epsilon(epsilon));
  CHECK(points[2][1] == AbsApprox(0.66639024601470138669).epsilon(epsilon));
  CHECK(points[3][1] == AbsApprox(0.17855872826361642311).epsilon(epsilon));
  CHECK(weights[0] == AbsApprox(0.090979309128011415315).epsilon(epsilon));
  CHECK(weights[1] == AbsApprox(0.090979309128011415315).epsilon(epsilon));
  CHECK(weights[2] == AbsApprox(0.15902069087198858472).epsilon(epsilon));
  CHECK(weights[3] == AbsApprox(0.15902069087198858472).epsilon(epsilon));
}

} // namespace seissol::unit_test