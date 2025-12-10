// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_PLASTICITY_H_
#define SEISSOL_SRC_MODEL_PLASTICITY_H_

#include "Alignment.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Model/CommonDatastructures.h"
#include "Solver/MultipleSimulations.h"

#include <cmath>
#include <cstddef>
#include <string>

namespace seissol::model {

// plasticity information per cell
struct PlasticityData {
  static constexpr auto PointCount = tensor::vNodes::Shape[0];

  // initial loading (stress tensor)
  alignas(Alignment) real initialLoading[tensor::QStressNodal::size()]{};
  alignas(Alignment) real cohesionTimesCosAngularFriction[tensor::meanStress::size()]{};
  alignas(Alignment) real sinAngularFriction[tensor::meanStress::size()]{};

  // depends on the material only (i.e. #1297 or multi-fused-material relevant only)
  real mufactor{};

  PlasticityData(const std::array<const Plasticity*, seissol::multisim::NumSimulations>& plasticity,
                 const Material* material) {
    auto initialLoadingV = init::QStressNodal::view::create(initialLoading);
    initialLoadingV.setZero();

    auto cohesionTimesCosAngularFrictionV =
        init::meanStress::view::create(cohesionTimesCosAngularFriction);
    cohesionTimesCosAngularFrictionV.setZero();

    auto sinAngularFrictionV = init::meanStress::view::create(sinAngularFriction);
    sinAngularFrictionV.setZero();

    for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
      auto initialLoadingVS = multisim::simtensor(initialLoadingV, s);
      auto cohesionTimesCosAngularFrictionVS =
          multisim::simtensor(cohesionTimesCosAngularFrictionV, s);
      auto sinAngularFrictionVS = multisim::simtensor(sinAngularFrictionV, s);

      for (std::size_t i = 0; i < PointCount; ++i) {
        initialLoadingVS(i, 0) = plasticity[s][i].sXX;
        initialLoadingVS(i, 1) = plasticity[s][i].sYY;
        initialLoadingVS(i, 2) = plasticity[s][i].sZZ;
        initialLoadingVS(i, 3) = plasticity[s][i].sXY;
        initialLoadingVS(i, 4) = plasticity[s][i].sYZ;
        initialLoadingVS(i, 5) = plasticity[s][i].sXZ;

        const double angularFriction = std::atan(plasticity[s][i].bulkFriction);

        cohesionTimesCosAngularFrictionVS(i) = plasticity[s][i].plastCo * std::cos(angularFriction);
        sinAngularFrictionVS(i) = std::sin(angularFriction);
      }
    }

    const auto mubar = material->getMuBar();
    mufactor = 1.0 / (2.0 * mubar);
  }

  static constexpr std::size_t NumQuantities = 7;
  static constexpr std::size_t NumberPerMechanism = 0;
  static constexpr std::size_t Parameters = 9;
  static const inline std::string Text = "plasticity";
  static inline const std::array<std::string, NumQuantities> Quantities = {
      "ep_xx", "ep_yy", "ep_zz", "ep_xy", "ep_yz", "ep_xz", "eta"};
};

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_PLASTICITY_H_
