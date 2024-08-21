#include "Kernels/precision.hpp"
#include "Numerical/Eigenvalues.h"

namespace seissol::unit_test {

template <unsigned dim>
void testResidual(std::array<std::complex<double>, dim * dim>& M,
                  seissol::eigenvalues::Eigenpair<std::complex<double>, dim>& eigenpair) {
  constexpr auto epsilon = std::numeric_limits<double>::epsilon() * 10;

  // compute M*R
  std::array<std::complex<double>, dim * dim> M_R{};
  for (unsigned int i = 0; i < dim; i++) {
    for (unsigned int j = 0; j < dim; j++) {
      for (unsigned int k = 0; k < dim; k++) {
        M_R[i + dim * j] += M[i + dim * k] * eigenpair.vectors[k + dim * j];
      }
    }
  }
  // compute R*L
  std::array<std::complex<double>, dim * dim> R_L{};
  for (unsigned int i = 0; i < dim; i++) {
    for (unsigned int j = 0; j < dim; j++) {
      R_L[i + dim * j] = eigenpair.vectors[i + dim * j] * eigenpair.values[j];
    }
  }
  // compare residual
  for (unsigned i = 0; i < dim * dim; i++) {
    REQUIRE(std::abs(M_R[i] - R_L[i]) == AbsApprox(0.0).epsilon(epsilon));
  }
}

TEST_CASE("Eigenvalues are correctly computed") {
  SUBCASE("Eigenvalues 3") {
    const unsigned dim = 3;

    std::array<std::array<std::complex<double>, dim * dim>, 2> matrices = {{
        {2.0, -3.0, -3.0, -2.0, 1.0, -2.0, -1.0, 1.0, 4.0},
        {-2.0, 3.0, 3.0, 2.0, -1.0, 2.0, 1.0, -1.0, -4.0},
    }};
    for (auto& M : matrices) {
      seissol::eigenvalues::Eigenpair<std::complex<double>, dim> eigenpair{};
      seissol::eigenvalues::computeEigenvaluesWithEigen3(M, eigenpair);
      testResidual<dim>(M, eigenpair);
#ifdef USE_POROELASTIC
      seissol::eigenvalues::computeEigenvaluesWithLapack(M, eigenpair);
      testResidual<dim>(M, eigenpair);
#endif
    }
  }
}

} // namespace seissol::unit_test