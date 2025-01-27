// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_DATASTRUCTURES_H_

// Gather all datastructure Headers here
#include "Equations/acoustic/Model/Datastructures.h"      // IWYU pragma: keep
#include "Equations/anisotropic/Model/Datastructures.h"   // IWYU pragma: keep
#include "Equations/elastic/Model/Datastructures.h"       // IWYU pragma: keep
#include "Equations/poroelastic/Model/Datastructures.h"   // IWYU pragma: keep
#include "Equations/viscoelastic2/Model/Datastructures.h" // IWYU pragma: keep

#include "Equations/acoustic/Model/IntegrationData.h"    // IWYU pragma: keep
#include "Equations/anisotropic/Model/IntegrationData.h" // IWYU pragma: keep
#include "Equations/elastic/Model/IntegrationData.h"     // IWYU pragma: keep
#ifdef USE_POROELASTIC
#include "Equations/poroelastic/Model/IntegrationData.h" // IWYU pragma: keep
#endif
#ifdef USE_VISCOELASTIC
#include "Equations/viscoelastic/Model/IntegrationData.h" // IWYU pragma: keep
#endif
#ifdef USE_VISCOELASTIC2
#include "Equations/viscoelastic2/Model/IntegrationData.h" // IWYU pragma: keep
#endif

namespace seissol::model {
#if defined(USE_ANISOTROPIC)
using MaterialT = AnisotropicMaterial;
#elif defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
using MaterialT = ViscoElasticMaterial;
#elif defined(USE_ELASTIC)
using MaterialT = ElasticMaterial;
#elif defined(USE_ACOUSTIC)
using MaterialT = AcousticMaterial;
#elif defined(USE_POROELASTIC)
using MaterialT = PoroElasticMaterial;
#else
#error "Material class unknown."
#endif
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DATASTRUCTURES_H_
