// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Datastructures.h"

#include "Equations/Setup.h" // IWYU pragma: keep

#include <array>
#include <complex>
#include <cstddef>
#include <limits>

// Also if we don't use poroelastic materials, the material has to be properly defined
// such that e.g. seissol::Interoperability::initializeModel work without ifdefs
#include "Numerical/Eigenvalues.h"
double seissol::model::PoroElasticMaterial::getPWaveSpeed() const {
  eigenvalues::Eigenpair<std::complex<double>, NumQuantities> eigendecomposition;
  std::array<std::complex<double>, NumQuantities * NumQuantities> atValues{};
  auto at = yateto::DenseTensorView<2, std::complex<double>>(atValues.data(),
                                                             {NumQuantities, NumQuantities});

  // the setup of the material exists in builds with a configuration of it, which the space-time
  // predictor advances
#ifdef SEISSOL_KERNELS_STP
  seissol::model::MaterialSetup<PoroElasticMaterial>::getTransposedCoefficientMatrix(*this, 0, at);
#endif

  // LAPACK where it is linked, as for every eigenvalue decomposition of this material
  seissol::eigenvalues::computeEigenvalues<eigenvalues::LapackIfLinked>(atValues,
                                                                        eigendecomposition);
  double maxEv = std::numeric_limits<double>::lowest();
  for (std::size_t i = 0; i < NumQuantities; i++) {
    maxEv = eigendecomposition.values.at(i).real() > maxEv ? eigendecomposition.values.at(i).real()
                                                           : maxEv;
  }
  return maxEv;
}
