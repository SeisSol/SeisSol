// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_EQUATIONS_SETUP_H_
#define SEISSOL_SRC_EQUATIONS_SETUP_H_

// Gather all Setup Headers here
#include "Equations/anisotropic/Model/AnisotropicSetup.h" // IWYU pragma: keep
#include "Equations/elastic/Model/ElasticSetup.h"         // IWYU pragma: keep
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/PoroelasticSetup.h" // IWYU pragma: keep
#endif
#ifdef USE_VISCOELASTIC
#include "Equations/viscoelastic/Model/ViscoelasticSetup.h" // IWYU pragma: keep
#endif
#ifdef USE_VISCOELASTIC2
#include "Equations/viscoelastic2/Model/ViscoelasticSetup.h" // IWYU pragma: keep
#endif

#endif // SEISSOL_SRC_EQUATIONS_SETUP_H_
