// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_WAVEPROPAGATION_PHYSICS_ENERGY_H_
#define SEISSOL_SRC_WAVEPROPAGATION_PHYSICS_ENERGY_H_

#include "WavePropagation/Physics/EnergyBase.h"

#include <array>
#include <cstddef>
#include <cstdint>
#include <string_view>

// IWYU pragma: begin_exports

// Gather all Energy headers here.
// Note: these specializations are keyed on the *material*, not on the solver
// variant. A header that needs generated tensors of one solver chooses by the
// solver of the configuration itself.
#include "WavePropagation/Physics/acoustic/Energy.h"
#include "WavePropagation/Physics/anisotropic/Energy.h"
#include "WavePropagation/Physics/elastic/Energy.h"
#include "WavePropagation/Physics/poroelastic/Energy.h"
#include "WavePropagation/Physics/viscoacoustic/Energy.h"
#include "WavePropagation/Physics/viscoelastic/Energy.h"

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_WAVEPROPAGATION_PHYSICS_ENERGY_H_
