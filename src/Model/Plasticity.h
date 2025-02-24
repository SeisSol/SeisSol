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
  real initialLoading[6];
  real cohesionTimesCosAngularFriction;
  real sinAngularFriction;
  real mufactor;

  PlasticityData(const Plasticity& plasticity, const Material* material) {
    initialLoading[0] = plasticity.sXX;
    initialLoading[1] = plasticity.sYY;
    initialLoading[2] = plasticity.sZZ;
    initialLoading[3] = plasticity.sXY;
    initialLoading[4] = plasticity.sYZ;
    initialLoading[5] = plasticity.sXZ;

    const double angularFriction = std::atan(plasticity.bulkFriction);

    cohesionTimesCosAngularFriction = plasticity.plastCo * std::cos(angularFriction);
    sinAngularFriction = std::sin(angularFriction);

    mufactor = 1.0 / (2.0 * material->getMuBar());
  }

  static constexpr std::size_t NumQuantities = 7;
  static constexpr std::size_t NumberPerMechanism = 0;
  static const inline std::string Text = "plasticity";
  static inline const std::array<std::string, NumQuantities> Quantities = {
      "ep_xx", "ep_yy", "ep_zz", "ep_xy", "ep_yz", "ep_xz", "eta"};
};

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_PLASTICITY_H_
