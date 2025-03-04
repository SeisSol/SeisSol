// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Datastructures.h"
#include <cmath>

// Also if we don't use poroelastic materials, the material has to be properly defined
// such that e.g. seissol::Interoperability::initializeModel work without ifdefs
#ifdef USE_POROELASTIC
#include "Numerical/Eigenvalues.h"
#include "PoroelasticSetup.h"
double seissol::model::PoroElasticMaterial::getPWaveSpeed() const {
  eigenvalues::Eigenpair<std::complex<double>, 13> eigendecomposition;
  std::array<std::complex<double>, 169> AT_values{};
  auto AT = yateto::DenseTensorView<2, std::complex<double>>(AT_values.data(), {13, 13});
  seissol::model::getTransposedCoefficientMatrix(*this, 0, AT);
  seissol::eigenvalues::computeEigenvalues(AT_values, eigendecomposition);
  double maxEv = std::numeric_limits<double>::lowest();
  for (int i = 0; i < 13; i++) {
    maxEv = eigendecomposition.values.at(i).real() > maxEv ? eigendecomposition.values.at(i).real()
                                                           : maxEv;
  }
  return maxEv;
}
#else
// Return an estimate which neglects fluid effects
double seissol::model::PoroElasticMaterial::getPWaveSpeed() const {
  return std::sqrt(lambda + 2 * mu) / rho;
}
#endif
