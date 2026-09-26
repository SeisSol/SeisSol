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
#include "Equations/anisotropic/Model/IntegrationData.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/viscoacoustic/Model/Datastructures.h"
#include "Equations/viscoelastic/Model/Datastructures.h"

// IWYU pragma: end_exports

namespace seissol::model {
template <MaterialType Type>
struct MaterialTypeSelector;

template <>
struct MaterialTypeSelector<MaterialType::Elastic> {
  using Type = ElasticMaterial;
};

template <>
struct MaterialTypeSelector<MaterialType::Anisotropic> {
  using Type = AnisotropicMaterial;
};

template <>
struct MaterialTypeSelector<MaterialType::Viscoelastic> {
  using Type = ViscoElasticMaterial<Config::RelaxationMechanisms>;
};

template <>
struct MaterialTypeSelector<MaterialType::Viscoacoustic> {
  using Type = ViscoAcousticMaterial<Config::RelaxationMechanisms>;
};

template <>
struct MaterialTypeSelector<MaterialType::Acoustic> {
  using Type = AcousticMaterial;
};

template <>
struct MaterialTypeSelector<MaterialType::Poroelastic> {
  using Type = PoroElasticMaterial;
};

using MaterialT = MaterialTypeSelector<Config::MaterialType>::Type;

} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_DATASTRUCTURES_H_
