#include "Kernels/precision.hpp"
#include "Numerical/Eigenvalues.h"

namespace seissol::unit_test {

template <unsigned dim>
void testResidual(std::array<std::complex<double>, dim * dim>& m,
                  seissol::eigenvalues::Eigenpair<std::complex<double>, dim>& eigenpair) {
  constexpr auto Epsilon = std::numeric_limits<double>::epsilon() * 10;

  // compute M*R
  std::array<std::complex<double>, dim * dim> mR{};
  for (unsigned int i = 0; i < dim; i++) {
    for (unsigned int j = 0; j < dim; j++) {
      for (unsigned int k = 0; k < dim; k++) {
        mR[i + dim * j] += m[i + dim * k] * eigenpair.vectors[k + dim * j];
      }
    }
  }
  // compute R*L
  std::array<std::complex<double>, dim * dim> rL{};
  for (unsigned int i = 0; i < dim; i++) {
    for (unsigned int j = 0; j < dim; j++) {
      rL[i + dim * j] = eigenpair.vectors[i + dim * j] * eigenpair.values[j];
    }
  }
  // compare residual
  for (unsigned i = 0; i < dim * dim; i++) {
    REQUIRE(std::abs(mR[i] - rL[i]) == AbsApprox(0.0).epsilon(Epsilon));
  }
}

TEST_CASE("Eigenvalues are correctly computed") {
  SUBCASE("Eigenvalues 3") {
    const unsigned dim = 3;

    std::array<std::array<std::complex<double>, dim * dim>, 2> matrices = {{
        {2.0, -3.0, -3.0, -2.0, 1.0, -2.0, -1.0, 1.0, 4.0},
        {-2.0, 3.0, 3.0, 2.0, -1.0, 2.0, 1.0, -1.0, -4.0},
    }};
    for (auto& m : matrices) {
      seissol::eigenvalues::Eigenpair<std::complex<double>, dim> eigenpair{};
      seissol::eigenvalues::computeEigenvaluesWithEigen3(m, eigenpair);
      testResidual<dim>(m, eigenpair);
#ifdef USE_POROELASTIC
      seissol::eigenvalues::computeEigenvaluesWithLapack(M, eigenpair);
      testResidual<dim>(M, eigenpair);
#endif
    }
  }
}

} // namespace seissol::unit_test