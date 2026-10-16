// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_MODEL_WEIGHTADJUSTED_H_
#define SEISSOL_SRC_INITIALIZER_MODEL_WEIGHTADJUSTED_H_

#include "Model/CommonDatastructures.h"

#include <algorithm>
#include <limits>

namespace seissol::initializer {

/**
 * Which cells are weight-adjusted (see WeightAdjustable), and how much shorter their time step is.
 *
 * A cell is where its impedances spread by more than the ratio the run asks for
 * (ModelParameters::weightAdjustedRatio) over the points its material is asked at: the samples the
 * operator is formed at and the points of its faces the flux is integrated at. Lifting the fluxes
 * with the material at the face grows where the material jumps inside a cell, the more the higher
 * the contrast. In an acoustic prototype it did so at order 4 from an impedance ratio of 4 on, at
 * order 5 at 10, at orders 3 and 6 not up to 10, and nowhere for a transition over a tenth of a
 * cell; where it starts depends on where the jump lies against the points.
 *
 * A weight-adjusted cell pays with its time step. Its faces lift the dissipation of the Riemann
 * problem without material, and the inverse weight at the samples scales it with the stiffest
 * material in the cell where the face may see the softest. On a periodic cube of 1296 cells at
 * order 4 and CFL 0.5, a slab whose jumps cut the cells stayed stable up to a time step between
 * 1.5 / r and 2 / r times the one its wave speeds permit, for acoustic jumps by impedance ratios r
 * of 3, 5, 10 and 20 and an elastic one by 10; with 2 / (1 + r), under which the spectral radius
 * of the prototype stayed, it grows at 10 and 20. The time step is shortened by 1 / r. An elastic
 * medium carries the same for its shear waves, which are slower: what they ask for counts with the
 * ratio of the speeds.
 */
class ImpedanceSpread {
  public:
  /// One more point of the cell.
  void add(const seissol::model::Material& material) {
    const double density = material.getDensity();
    const double pWave = material.getPWaveSpeed();
    const double sWave = material.getSWaveSpeed();
    include(density * pWave, minimumP_, maximumP_);
    pWaveSpeed_ = std::max(pWaveSpeed_, pWave);
    // a point without shear has no shear impedance to compare
    if (sWave > 0.0) {
      include(density * sWave, minimumS_, maximumS_);
      sWaveSpeed_ = std::max(sWaveSpeed_, sWave);
    }
  }

  /// The largest ratio between two P impedances of the cell, or between two S impedances of its
  /// points that carry shear. One where the material does not vary.
  [[nodiscard]] double ratio() const { return std::max(ratioP(), ratioS()); }

  /// What the time step of a weight-adjusted cell is shortened by against the one its largest wave
  /// speed permits: 1 / r for the P waves, and the same for the S waves weighed with their speed
  /// against the largest P wave speed.
  [[nodiscard]] double timestepFactor() const {
    const double pWaves = ratioP();
    const double sWaves = pWaveSpeed_ > 0.0 ? sWaveSpeed_ / pWaveSpeed_ * ratioS() : 0.0;
    return 1.0 / std::max(pWaves, sWaves);
  }

  private:
  static void include(double impedance, double& minimum, double& maximum) {
    minimum = std::min(minimum, impedance);
    maximum = std::max(maximum, impedance);
  }

  [[nodiscard]] double ratioP() const { return maximumP_ > 0.0 ? maximumP_ / minimumP_ : 1.0; }
  [[nodiscard]] double ratioS() const { return maximumS_ > 0.0 ? maximumS_ / minimumS_ : 1.0; }

  double minimumP_{std::numeric_limits<double>::max()};
  double maximumP_{0.0};
  double minimumS_{std::numeric_limits<double>::max()};
  double maximumS_{0.0};
  double pWaveSpeed_{0.0};
  double sWaveSpeed_{0.0};
};

/// Whether a cell whose impedances spread by `ratio` is weight-adjusted, for the ratio the run asks
/// for; none is where that is zero.
inline bool isWeightAdjusted(double ratio, double threshold) {
  return threshold > 0.0 && ratio > threshold;
}

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_MODEL_WEIGHTADJUSTED_H_
