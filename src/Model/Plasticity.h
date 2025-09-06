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
#include <string>

namespace seissol::model {
// plasticity information per cell. In case of multiple simulations, it contains data of all
// simulations

template <typename Cfg>
struct PlasticityData {
  using real = Real<Cfg>;
  // initial loading (stress tensor)
  real initialLoading[6 * seissol::multisim::NumSimulations<Cfg>];
  real cohesionTimesCosAngularFriction[seissol::multisim::NumSimulations<Cfg>];
  real sinAngularFriction[seissol::multisim::NumSimulations<Cfg>];
  real mufactor; // Only dependent on mu which is to be constant for all simulations

  PlasticityData(const std::array<Plasticity, seissol::multisim::NumSimulations<Cfg>>& plasticity,
                 const Material* material) {
    for (std::size_t i = 0; i < seissol::multisim::NumSimulations<Cfg>; ++i) {
      // interleave these so that the kernel does not need any modifications
      initialLoading[0 * seissol::multisim::NumSimulations<Cfg> + i] = plasticity[i].sXX;
      initialLoading[1 * seissol::multisim::NumSimulations<Cfg> + i] = plasticity[i].sYY;
      initialLoading[2 * seissol::multisim::NumSimulations<Cfg> + i] = plasticity[i].sZZ;
      initialLoading[3 * seissol::multisim::NumSimulations<Cfg> + i] = plasticity[i].sXY;
      initialLoading[4 * seissol::multisim::NumSimulations<Cfg> + i] = plasticity[i].sYZ;
      initialLoading[5 * seissol::multisim::NumSimulations<Cfg> + i] = plasticity[i].sXZ;

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
