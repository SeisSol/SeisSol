// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_DATASTRUCTURES_H_

#include "Config.h"
#include "Model/CommonDatastructures.h"

// IWYU pragma: begin_exports

// Gather all datastructure Headers here
#include "Equations/acoustic/Model/Datastructures.h"
#include "Equations/anisotropic/Model/Datastructures.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/IntegrationData.h"
#include "Equations/viscoelastic/Model/IntegrationData.h"
#include "Equations/viscoelastic2/Model/Datastructures.h"
#include "Equations/viscoelastic2/Model/IntegrationData.h"

// IWYU pragma: end_exports

namespace seissol::model {
template <MaterialType Type, std::size_t Mechanisms>
struct MaterialTypeSelector;

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Elastic, Mechanisms> {
  using Type = ElasticMaterial;
};

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Anisotropic, Mechanisms> {
  using Type = AnisotropicMaterial;
};

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Viscoelastic, Mechanisms> {
  using Type = ViscoElasticMaterialParametrized<Mechanisms>;
};

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Acoustic, Mechanisms> {
  using Type = AcousticMaterial;
};

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Poroelastic, Mechanisms> {
  using Type = PoroElasticMaterial;
};

template <typename Config>
using MaterialTT =
    typename MaterialTypeSelector<Config::MaterialType, Config::RelaxationMechanisms>::Type;

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DATASTRUCTURES_H_
