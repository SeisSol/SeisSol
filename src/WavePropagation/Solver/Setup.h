// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SETUP_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SETUP_H_

// IWYU pragma: begin_exports

// Gather all Setup Headers here
#ifdef SEISSOL_KERNELS_LINEARCK
#include "WavePropagation/Solver/LinearCK/Setup.h"
#endif
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
#include "WavePropagation/Solver/LinearCKAnelastic/Setup.h"
#endif
#ifdef SEISSOL_KERNELS_STP
#include "WavePropagation/Solver/STP/Setup.h"
#endif

#include "WavePropagation/Physics/acoustic/Setup.h"
#include "WavePropagation/Physics/anisotropic/Setup.h"
#include "WavePropagation/Physics/elastic/Setup.h"
#include "WavePropagation/Physics/poroelastic/Setup.h"
#include "WavePropagation/Physics/viscoacoustic/Setup.h"
#include "WavePropagation/Physics/viscoelastic/Setup.h"

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_SETUP_H_
