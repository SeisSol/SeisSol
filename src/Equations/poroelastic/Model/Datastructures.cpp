// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Datastructures.h"
#include <Equations/Datastructures.h>
#include <Equations/Setup.h> // IWYU pragma: keep
#include <Model/Common.h>
#include <Model/CommonDatastructures.h>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

// Also if we don't use poroelastic materials, the material has to be properly defined
// such that e.g. seissol::Interoperability::initializeModel work without ifdefs
#include "Numerical/Eigenvalues.h"
double seissol::model::PoroElasticMaterial::getPWaveSpeed() const {
  eigenvalues::Eigenpair<std::complex<double>, NumQuantities> eigendecomposition;
  std::array<std::complex<double>, NumQuantities * NumQuantities> atValues{};
  auto at = yateto::DenseTensorView<2, std::complex<double>>(atValues.data(),
                                                             {NumQuantities, NumQuantities});

  // TODO: remove this if constexpr guard (needs multi-equation build support)
  if constexpr (seissol::model::MaterialT::Type == seissol::model::MaterialType::Poroelastic) {
    seissol::model::getTransposedCoefficientMatrix(*this, 0, at);
  }

  seissol::eigenvalues::computeEigenvalues(atValues, eigendecomposition);
  double maxEv = std::numeric_limits<double>::lowest();
  for (int i = 0; i < NumQuantities; i++) {
    maxEv = eigendecomposition.values.at(i).real() > maxEv ? eigendecomposition.values.at(i).real()
                                                           : maxEv;
  }
  return maxEv;
}
