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
#include <cmath>
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

  // TODO: make vary as well? (at least for #1297)
  real mufactor;

  PlasticityData(const Plasticity* plasticity, const Material* material) {
    auto initialLoadingV = init::QStressNodal::view::create(initialLoading);
    initialLoadingV.setZero();

    for (std::size_t i = 0; i < tensor::meanStress::size(); ++i) {
      initialLoadingV(i, 0) = plasticity[i].sXX;
      initialLoadingV(i, 1) = plasticity[i].sYY;
      initialLoadingV(i, 2) = plasticity[i].sZZ;
      initialLoadingV(i, 3) = plasticity[i].sXY;
      initialLoadingV(i, 4) = plasticity[i].sYZ;
      initialLoadingV(i, 5) = plasticity[i].sXZ;

      const double angularFriction = std::atan(plasticity[i].bulkFriction);

      cohesionTimesCosAngularFriction[i] = plasticity[i].plastCo * std::cos(angularFriction);
      sinAngularFriction[i] = std::sin(angularFriction);
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
