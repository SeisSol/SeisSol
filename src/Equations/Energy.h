// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_ENERGY_H_

namespace seissol::model {

template <typename MaterialT>
struct EnergyCompute;

} // namespace seissol::model

// IWYU pragma: begin_exports

// Gather all Setup Headers here
#include "Equations/acoustic/Model/Energy.h"
#include "Equations/anisotropic/Model/Energy.h"
#include "Equations/elastic/Model/Energy.h"
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/Energy.h"
#endif
#ifdef USE_VISCOELASTIC
#include "Equations/viscoelastic2/Model/Energy.h"
#endif
#ifdef USE_VISCOELASTIC2
#include "Equations/viscoelastic2/Model/Energy.h"
#endif

#endif // SEISSOL_SRC_EQUATIONS_ENERGY_H_
