// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_WAVEPROPAGATIONPARAMETERS_H_
#define SEISSOL_SRC_WAVEPROPAGATION_WAVEPROPAGATIONPARAMETERS_H_

#include <string>

namespace seissol::initializer::parameters {

enum class NumericalFlux { Godunov, Rusanov };

inline std::string fluxToString(NumericalFlux flux) {
  if (flux == NumericalFlux::Godunov) {
    return "Godunov flux";
  }
  if (flux == NumericalFlux::Rusanov) {
    return "Rusanov flux";
  }
  return "(unknown flux)";
}

struct WavePropagationParameters {
  NumericalFlux flux{NumericalFlux::Godunov};
  // the flux at the faces of the cells next to a fault (one of the two cells has a fault face)
  NumericalFlux fluxNearFault{NumericalFlux::Godunov};
  // the attenuation of the viscoelastic and viscoacoustic materials: its central frequency, and
  // the ratio of the highest to the lowest frequency of the band it holds for
  double freqCentral{};
  double freqRatio{1.0};
  // the material of a cell averaged over the cell, rather than taken at its barycenter
  bool useCellHomogenizedMaterial{true};
  // of the free surface with gravity
  double gravitationalAcceleration{};
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_WAVEPROPAGATION_WAVEPROPAGATIONPARAMETERS_H_
