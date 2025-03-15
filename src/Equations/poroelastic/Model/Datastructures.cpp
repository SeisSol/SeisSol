// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Datastructures.h"
#include <Equations/Datastructures.h>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

// Also if we don't use poroelastic materials, the material has to be properly defined
// such that e.g. seissol::Interoperability::initializeModel work without ifdefs
#include "Numerical/Eigenvalues.h"
double seissol::model::PoroElasticMaterial::getPWaveSpeed() const {
  eigenvalues::Eigenpair<std::complex<double>, 13> eigendecomposition;
  std::array<std::complex<double>, 169> atValues{};
  auto at = yateto::DenseTensorView<2, std::complex<double>>(atValues.data(), {13, 13});
#ifdef USE_POROELASTIC
  seissol::model::getTransposedCoefficientMatrix(*this, 0, at);
#endif
  seissol::eigenvalues::computeEigenvalues(atValues, eigendecomposition);
  double maxEv = std::numeric_limits<double>::lowest();
  for (int i = 0; i < 13; i++) {
    maxEv = eigendecomposition.values.at(i).real() > maxEv ? eigendecomposition.values.at(i).real()
                                                           : maxEv;
  }
  return maxEv;
}
