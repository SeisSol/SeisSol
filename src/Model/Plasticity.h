// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_PLASTICITY_H_
#define SEISSOL_SRC_MODEL_PLASTICITY_H_

#include "Alignment.h"
#include "Common/Real.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Model/CommonDatastructures.h"
#include "Model/PlasticityQuantities.h"
#include "Solver/MultipleSimulations.h"

#include <array>
#include <cmath>
#include <cstddef>
#include <string>

namespace seissol::model {

// plasticity information per cell
template <typename Cfg>
struct PlasticityData {
  static constexpr auto PointCount = tensor::vNodes<Cfg>::Shape[0];

  // initial loading (stress tensor)
  alignas(Alignment) Real<Cfg> initialLoading[tensor::initialLoading<Cfg>::size()]{};
  alignas(Alignment) Real<Cfg> cohesionTimesCosAngularFriction[tensor::meanStress<Cfg>::size()]{};
  alignas(Alignment) Real<Cfg> sinAngularFriction[tensor::meanStress<Cfg>::size()]{};

  // depends only on the material (i.e. only relevant for #1297 or multi-fused-material)
  Real<Cfg> mufactor{};

  PlasticityData(const std::array<const Plasticity*, Cfg::NumSimulations>& plasticity,
                 const Material* material,
                 bool pointwise) {
    auto initialLoadingV = init::initialLoading<Cfg>::view::create(initialLoading);
    initialLoadingV.setZero();

    auto cohesionTimesCosAngularFrictionV =
        init::meanStress<Cfg>::view::create(cohesionTimesCosAngularFriction);
    cohesionTimesCosAngularFrictionV.setZero();

    auto sinAngularFrictionV = init::meanStress<Cfg>::view::create(sinAngularFriction);
    sinAngularFrictionV.setZero();

    using Multisim = multisim::MultisimHelperWrapper<Cfg>;
    for (std::size_t s = 0; s < Cfg::NumSimulations; ++s) {
      auto initialLoadingVS = Multisim::simtensor(initialLoadingV, s);
      auto cohesionTimesCosAngularFrictionVS =
          Multisim::simtensor(cohesionTimesCosAngularFrictionV, s);
      auto sinAngularFrictionVS = Multisim::simtensor(sinAngularFrictionV, s);

      for (std::size_t i = 0; i < PointCount; ++i) {
        const auto ii = pointwise ? i : 0;
        initialLoadingVS(i, 0) = plasticity[s][ii].sXX;
        initialLoadingVS(i, 1) = plasticity[s][ii].sYY;
        initialLoadingVS(i, 2) = plasticity[s][ii].sZZ;
        initialLoadingVS(i, 3) = plasticity[s][ii].sXY;
        initialLoadingVS(i, 4) = plasticity[s][ii].sYZ;
        initialLoadingVS(i, 5) = plasticity[s][ii].sXZ;

        const double angularFriction = std::atan(plasticity[s][ii].bulkFriction);

        cohesionTimesCosAngularFrictionVS(i) =
            plasticity[s][ii].plastCo * std::cos(angularFriction);
        sinAngularFrictionVS(i) = std::sin(angularFriction);
      }
    }

    const auto mubar = material->getMuBar();
    mufactor = 1.0 / (2.0 * mubar);
  }

  static constexpr std::size_t NumQuantities = PlasticityQuantityCount;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Parameters = 9;
  static const inline std::string Text = "plasticity";
  static inline const std::array<std::string, NumQuantities> Quantities = PlasticityQuantities;
};

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_PLASTICITY_H_
