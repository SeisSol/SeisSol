#ifndef SEISSOL_MODEL_PLASTICITY_HPP_
#define SEISSOL_MODEL_PLASTICITY_HPP_

#include "Model/common_datastructures.hpp"
#include <Kernels/precision.hpp>
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

  static constexpr std::size_t NumberOfQuantities = 7;
  static constexpr std::size_t NumberPerMechanism = 0;
  static const inline std::string Text = "plasticity";
  static inline const std::array<std::string, NumberOfQuantities> Quantities = {
      "ep_xx", "ep_yy", "ep_zz", "ep_xy", "ep_yz", "ep_xz", "eta"};
};

} // namespace seissol::model

#endif // SEISSOL_MODEL_PLASTICITY_HPP_
