// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_PHYSICS_DATASTRUCTURES_H_
#define SEISSOL_SRC_WAVEPROPAGATION_PHYSICS_DATASTRUCTURES_H_

#include "Config.h"
#include "WavePropagation/Physics/CommonDatastructures.h"

// IWYU pragma: begin_exports

// Gather all datastructure Headers here
#include "WavePropagation/Physics/acoustic/Datastructures.h"
#include "WavePropagation/Physics/anisotropic/Datastructures.h"
#include "WavePropagation/Physics/anisotropic/IntegrationData.h"
#include "WavePropagation/Physics/elastic/Datastructures.h"
#include "WavePropagation/Physics/poroelastic/Datastructures.h"
#include "WavePropagation/Physics/viscoacoustic/Datastructures.h"
#include "WavePropagation/Physics/viscoelastic/Datastructures.h"

// IWYU pragma: end_exports

#include <cstddef>

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
  using Type = ViscoElasticMaterial<Mechanisms>;
};

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Viscoacoustic, Mechanisms> {
  using Type = ViscoAcousticMaterial<Mechanisms>;
};

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Acoustic, Mechanisms> {
  using Type = AcousticMaterial;
};

template <std::size_t Mechanisms>
struct MaterialTypeSelector<MaterialType::Poroelastic, Mechanisms> {
  using Type = PoroElasticMaterial;
};

/// The material a configuration computes with.
template <typename Cfg>
using MaterialOf = MaterialTypeSelector<Cfg::MaterialType, Cfg::RelaxationMechanisms>::Type;

} // namespace seissol::model

#endif // SEISSOL_SRC_WAVEPROPAGATION_PHYSICS_DATASTRUCTURES_H_
