// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_SCENARIO_PLANEWAVE_H_
#define SEISSOL_SRC_PHYSICS_SCENARIO_PLANEWAVE_H_

#include "Common/ConfigRegistry.h"
#include "Initializer/Typedefs.h"
#include "Physics/InitialField.h"

#include <Eigen/Core>
#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

namespace seissol::physics {

/// One mode of a plane wave: the real part of exp(i (omega t - k.(x - o) + phi)) times
/// `coefficients`, the eigenvector of the plane wave operator that belongs to the eigenvalue
/// omega, times an amplitude.
struct PlaneWaveMode {
  Eigen::Vector3d k;
  std::complex<double> omega;
  // one per component of the state
  std::vector<std::complex<double>> coefficients;
};

/// The modes `modes` of the plane waves with the wave vector `k` in the material `materialData`
/// of the configuration `config`, with the amplitudes `amplitudes`. A mode is an index into the
/// eigenpairs of the plane wave operator, which are sorted by the real part of the eigenvalue.
std::vector<PlaneWaveMode> planeWaveModes(const CellMaterialData& materialData,
                                          ConfigId config,
                                          const Eigen::Vector3d& k,
                                          const std::vector<int>& modes,
                                          const std::vector<std::complex<double>>& amplitudes);

/// The modes a plane wave imposes in a material unless told otherwise, with their amplitudes.
std::pair<std::vector<int>, std::vector<std::complex<double>>>
    defaultPlaneWaveModes(const CellMaterialData& materialData, ConfigId config);

/// A sum of plane waves, each a sum of modes. All share the phase `phase`, an sderiv expression.
/// With an origin, the modes travel from it, and a mode is cut to one period: it is there where
/// its argument lies in (-pi/2, 3pi/2) (the travelling wave).
struct PlaneWaveSum {
  std::vector<std::vector<PlaneWaveMode>> waves;
  std::string phase{"0.0"};
  std::optional<Eigen::Vector3d> origin;
};

/// The sderiv module of a sum of plane waves, giving the components `names` of the state in the
/// order of the coefficients of the modes. It computes in real arithmetic: a mode with the
/// eigenvalue omega and the coefficient a + ib of a component adds
/// exp(-Im(omega) t) (a cos(theta) - b sin(theta)), theta = Re(omega) t - k.(x - o) + phi.
std::string planeWaveModule(const PlaneWaveSum& sum, const std::vector<std::string>& names);

/// The planar wave of the scenario Planarwave, one field per fused simulation (shifted in phase
/// by 2 pi sim / simulations), in the material `materialData` of the configuration `config`.
std::vector<std::unique_ptr<InitialField>> planarwave(const CellMaterialData& materialData,
                                                      ConfigId config);

/// Three planar waves along the axes, superimposed (SuperimposedPlanarwave).
std::vector<std::unique_ptr<InitialField>>
    superimposedPlanarwave(const CellMaterialData& materialData, ConfigId config);

/// A part of a planar wave travelling in one direction (Travelling).
std::vector<std::unique_ptr<InitialField>>
    travellingWave(const CellMaterialData& materialData,
                   ConfigId config,
                   const TravellingWaveParameters& parameters);

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_SCENARIO_PLANEWAVE_H_
