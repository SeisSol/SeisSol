// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_EQUATIONS_IMPEDANCE_H_
#define SEISSOL_SRC_EQUATIONS_IMPEDANCE_H_

#include "Equations/ImpedanceBase.h"

// IWYU pragma: begin_exports

// Gather all Impedance headers here.
// Unlike Setup.h and Energy.h, no guards are needed: the admittance of a material only depends on
// its parameters, not on the generated code of the build. The acoustic material does not support
// dynamic rupture and has no specialization.
#include "Equations/anisotropic/Model/Impedance.h"
#include "Equations/elastic/Model/Impedance.h"
#include "Equations/poroelastic/Model/Impedance.h"
#include "Equations/viscoelastic2/Model/Impedance.h"

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_EQUATIONS_IMPEDANCE_H_
