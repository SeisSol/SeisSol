// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_PLASTICITY_H_
#define SEISSOL_SRC_MODEL_PLASTICITY_H_

#include "Model/CommonDatastructures.h"
#include <Alignment.h>
#include <Kernels/Precision.h>
#include <Solver/MultipleSimulations.h>
#include <cmath>
#include <equation-elastic-6-double/init.h>
#include <generated_code/init.h>
#include <generated_code/tensor.h>
#include <string>

namespace seissol::model {

// plasticity information per cell
struct PlasticityData {
  // initial loading (stress tensor)
  alignas(Alignment) real initialLoading[tensor::QStressNodal::size()];
  alignas(Alignment) real cohesionTimesCosAngularFriction[tensor::meanStress::size()];
  alignas(Alignment) real sinAngularFriction[tensor::meanStress::size()];

  // depends on the material only (i.e. #1297 or multi-fused-material relevant only)
  real mufactor;

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
      for (std::size_t i = 0; i < tensor::meanStress::size(); ++i) {
        initialLoadingV(i, 0) = plasticity[s][i].sXX;
        initialLoadingV(i, 1) = plasticity[s][i].sYY;
        initialLoadingV(i, 2) = plasticity[s][i].sZZ;
        initialLoadingV(i, 3) = plasticity[s][i].sXY;
        initialLoadingV(i, 4) = plasticity[s][i].sYZ;
        initialLoadingV(i, 5) = plasticity[s][i].sXZ;

        const double angularFriction = std::atan(plasticity[s][i].bulkFriction);

        cohesionTimesCosAngularFrictionV(i) = plasticity[s][i].plastCo * std::cos(angularFriction);
        sinAngularFrictionV(i) = std::sin(angularFriction);
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
