// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_INTEGRATIONDATA_H_
#define SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_INTEGRATIONDATA_H_

#include "Common/Real.h"

#include <array>

namespace seissol::model {

template <typename Cfg>
struct AnisotropicEnergyData {
  std::array<Real<Cfg>, 36> matS{};
};

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_ANISOTROPIC_MODEL_INTEGRATIONDATA_H_
