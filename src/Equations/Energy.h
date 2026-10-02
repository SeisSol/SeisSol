// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_ENERGY_H_
#define SEISSOL_SRC_EQUATIONS_ENERGY_H_

#include "Equations/EnergyBase.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

// IWYU pragma: begin_exports

// Gather all Energy headers here.
// Note: these specializations are keyed on the *material*, not on the solver
// variant. A header that needs generated tensors of one solver chooses by the
// solver of the configuration itself.
#include "Equations/acoustic/Model/Energy.h"
#include "Equations/anisotropic/Model/Energy.h"
#include "Equations/elastic/Model/Energy.h"
#include "Equations/poroelastic/Model/Energy.h"
#include "Equations/viscoacoustic/Model/Energy.h"
#include "Equations/viscoelastic/Model/Energy.h"

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_EQUATIONS_ENERGY_H_
