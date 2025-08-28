// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_PLASTICITY_H_
#define SEISSOL_SRC_MODEL_PLASTICITY_H_

#include "Model/CommonDatastructures.h"
#include <Kernels/Precision.h>
#include <cmath>
#include <string>

namespace seissol::model {
// plasticity information per cell

struct PlasticityData {
  // initial loading (stress tensor)
  real initialLoading[MULTIPLE_SIMULATIONS][6];
  real cohesionTimesCosAngularFriction[MULTIPLE_SIMULATIONS];
  real sinAngularFriction[MULTIPLE_SIMULATIONS];
  real mufactor; // Only dependent on mu which is to be constant per simulation

  PlasticityData(const std::array<Plasticity, MULTIPLE_SIMULATIONS>& plasticity, const Material* material) {
    for(std::size_t i = 0; i < MULTIPLE_SIMULATIONS; ++i) {
      initialLoading[i][0] = plasticity[i].sXX;
      initialLoading[i][1] = plasticity[i].sYY;
      initialLoading[i][2] = plasticity[i].sZZ;
      initialLoading[i][3] = plasticity[i].sXY;
      initialLoading[i][4] = plasticity[i].sYZ;
      initialLoading[i][5] = plasticity[i].sXZ;

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
