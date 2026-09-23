// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_SETUP_H_
#define SEISSOL_SRC_EQUATIONS_SETUP_H_

// IWYU pragma: begin_exports

// Gather all Setup Headers here
#ifdef SEISSOL_KERNELS_LINEARCK
#include "Kernels/LinearCK/Setup.h"
#endif
#ifdef SEISSOL_KERNELS_LINEARCKANELASTIC
#include "Kernels/LinearCKAnelastic/Setup.h"
#endif
#ifdef SEISSOL_KERNELS_STP
#include "Kernels/STP/Setup.h"
#endif

#include "Equations/acoustic/Model/Setup.h"
#include "Equations/anisotropic/Model/Setup.h"
#include "Equations/elastic/Model/Setup.h"
#include "Equations/poroelastic/Model/Setup.h"
#include "Equations/viscoacoustic/Model/Setup.h"
#include "Equations/viscoelastic/Model/Setup.h"

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_EQUATIONS_SETUP_H_
