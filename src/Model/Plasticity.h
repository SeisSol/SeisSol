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
#include <Solver/MultipleSimulations.h>
#include <cmath>
#include <cstddef>
#include <string>

namespace seissol::model {
// plasticity information per cell. In case of multiple simulations, it contains data of all
// simulations

struct PlasticityData {
  // initial loading (stress tensor)
  real initialLoading[6 * seissol::multisim::NumSimulations];
  real cohesionTimesCosAngularFriction[seissol::multisim::NumSimulations];
  real sinAngularFriction[seissol::multisim::NumSimulations];
  real mufactor; // Only dependent on mu which is to be constant for all simulations

  PlasticityData(const std::array<Plasticity, seissol::multisim::NumSimulations>& plasticity,
                 const Material* material) {
    for (std::size_t i = 0; i < seissol::multisim::NumSimulations; ++i) {
      // interleave these so that the kernel does not need any modifications
      initialLoading[static_cast<std::size_t>(0 * seissol::multisim::NumSimulations) + i] =
          plasticity[i].sXX;
      initialLoading[static_cast<std::size_t>(1 * seissol::multisim::NumSimulations) + i] =
          plasticity[i].sYY;
      initialLoading[static_cast<std::size_t>(2 * seissol::multisim::NumSimulations) + i] =
          plasticity[i].sZZ;
      initialLoading[static_cast<std::size_t>(3 * seissol::multisim::NumSimulations) + i] =
          plasticity[i].sXY;
      initialLoading[static_cast<std::size_t>(4 * seissol::multisim::NumSimulations) + i] =
          plasticity[i].sYZ;
      initialLoading[static_cast<std::size_t>(5 * seissol::multisim::NumSimulations) + i] =
          plasticity[i].sXZ;

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
