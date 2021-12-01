#include <Eigen/Dense>
#include <Numerical_aux/Transformation.h>

namespace seissol::unit_test {

TEST_CASE("Test tetrahedron global to reference") {
  // We do all tests in double precision
  constexpr real epsilon = 10 * std::numeric_limits<double>::epsilon();

  std::srand(9);
  const auto vertices =
      std::array<Eigen::Vector3d, 4>{{Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX),
                                      Eigen::Vector3d((double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX,
                                                      (double)std::rand() / RAND_MAX)}};
  const auto center = 0.25 * (vertices[0] + vertices[1] + vertices[2] + vertices[3]);

  const auto res = seissol::transformations::tetrahedronGlobalToReference(
      vertices[0].data(), vertices[1].data(), vertices[2].data(), vertices[3].data(), center);
  REQUIRE(res(0) == AbsApprox(0.25).epsilon(epsilon));
  REQUIRE(res(1) == AbsApprox(0.25).epsilon(epsilon));
  REQUIRE(res(2) == AbsApprox(0.25).epsilon(epsilon));
}

} // namespace seissol::unit_test