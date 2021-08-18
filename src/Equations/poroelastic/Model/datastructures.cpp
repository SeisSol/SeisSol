#include "datastructures.hpp"
#include "PoroelasticSetup.h"
#include "Numerical_aux/Eigenvalues.h"


#ifdef USE_POROELASTIC
double seissol::model::PoroElasticMaterial::getPWaveSpeed() const {
  eigenvalues::Eigenpair<std::complex<double>, 13> eigendecomposition;
  std::array<std::complex<double>, 169> AT_values{};
  auto AT = yateto::DenseTensorView<2,std::complex<double>>(AT_values.data(), {13, 13});
  seissol::model::getTransposedCoefficientMatrix(*this, 0, AT);
  seissol::eigenvalues::computeEigenvaluesWithArmadillo(AT_values, eigendecomposition);
  double max_ev = std::numeric_limits<double>::lowest();
  for (int i = 0; i < 13; i++) {
    max_ev = eigendecomposition.values.at(i).real() > max_ev ? eigendecomposition.values.at(i).real() : max_ev;
  }
  return max_ev;
}
#endif
