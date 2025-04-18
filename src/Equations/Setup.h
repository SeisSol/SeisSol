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
#include "Equations/acoustic/Model/AcousticSetup.h"
#include "Equations/anisotropic/Model/AnisotropicSetup.h"
#include "Equations/elastic/Model/ElasticSetup.h"
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/PoroelasticSetup.h"
#endif
#ifdef USE_VISCOELASTIC
#include "Equations/viscoelastic/Model/ViscoelasticSetup.h"
#endif
#ifdef USE_VISCOELASTIC2
#include "Equations/viscoelastic2/Model/ViscoelasticSetup.h"
#endif

// IWYU pragma: end_exports

#endif // SEISSOL_SRC_EQUATIONS_SETUP_H_
