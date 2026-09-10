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
// variant. ViscoElasticMaterial is defined in viscoelastic2/ and shared by both
// USE_VISCOELASTIC and USE_VISCOELASTIC2 (mirroring Datastructures.h), so there
// is no viscoelastic/Model/Energy.h. The guards are only needed where the header
// pulls in generated tensors that exist for that build alone.
#include "Equations/acoustic/Model/Energy.h"
#include "Equations/anisotropic/Model/Energy.h"
#include "Equations/elastic/Model/Energy.h"
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/Energy.h"
#endif
#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
#include "Equations/viscoelastic2/Model/Energy.h"
#endif

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_EQUATIONS_ENERGY_H_
