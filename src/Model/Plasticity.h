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
  alignas(Alignment) real initialLoading[tensor::initialLoading::size()]{};
  alignas(Alignment) real cohesionTimesCosAngularFriction[tensor::meanStress::size()]{};
  alignas(Alignment) real sinAngularFriction[tensor::meanStress::size()]{};

  // the one quantity here that comes from the material rather than from the
  // plasticity parameters, and so the one that varies within a cell once the
  // material does
  alignas(Alignment) real mufactor[tensor::meanStress::size()]{};

  PlasticityData(const std::array<const Plasticity*, seissol::multisim::NumSimulations>& plasticity,
                 const std::array<double, PointCount>& muBar,
                 bool pointwise) {
    auto initialLoadingV = init::initialLoading::view::create(initialLoading);
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

    auto mufactorV = init::meanStress::view::create(mufactor);
    mufactorV.setZero();
    for (std::size_t s = 0; s < multisim::NumSimulations; ++s) {
      auto mufactorVS = multisim::simtensor(mufactorV, s);
      for (std::size_t i = 0; i < PointCount; ++i) {
        mufactorVS(i) = 1.0 / (2.0 * muBar[i]);
      }
    }
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
