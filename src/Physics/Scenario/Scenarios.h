// SPDX-FileCopyrightText: 2019 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_PHYSICS_SCENARIO_SCENARIOS_H_
#define SEISSOL_SRC_PHYSICS_SCENARIO_SCENARIOS_H_

#include "Common/ConfigRegistry.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Initializer/Typedefs.h"
#include "Physics/InitialField.h"

#include <memory>
#include <string>
#include <vector>

// The scenarios with a closed form, which SeisSol writes as sderiv modules at setup (see
// ScriptedScenario.h): the constants they take from the material and the parameters are numbers in
// the source. The planar waves are in PlaneWave.h.

namespace seissol::physics {

/// No wave: the state is zero.
std::vector<std::unique_ptr<InitialField>> zeroField(ConfigId config);

/// A Scholte wave along the interface z = 0 between an elastic material below and an acoustic one
/// (mu = 0) above. A cell is acoustic if its mu is zero.
std::vector<std::unique_ptr<InitialField>> scholteWave(ConfigId config);

/// A plane wave refracted at the interface z = 0 between an elastic material below and an
/// acoustic one (mu = 0) above, as Snell's law says.
std::vector<std::unique_ptr<InitialField>> snellsLaw(ConfigId config);

/// A standing wave of the mode `mode` (0, the gravity mode, 1 or 2) in an ocean of an acoustic
/// material (mu = 0) with the gravitational acceleration 9.81e-3, from
///
/// Abrahams, L. S., Krenz, L., Dunham, E. M., & Gabriel, A. A. (2019, December).
/// Verification of a 3D fully-coupled earthquake and tsunami model.
/// In AGU Fall Meeting Abstracts (Vol. 2019, pp. NH43F-1000).
///
/// a 3D extension of the 2D scenario in
///
/// Lotto, G. C., & Dunham, E. M. (2015).
/// High-order finite difference modeling of tsunami generation in a compressible ocean from
/// offshore earthquakes. Computational Geosciences, 19(2), 327-340.
///
/// It reads the density of each cell; a cell with shear is an error.
std::vector<std::unique_ptr<InitialField>>
    ocean(int mode, double gravitationalAcceleration, ConfigId config);

/// A gaussian perturbation of the pressure of a poroelastic material, with the origin, the
/// magnitude and the width of `parameters`.
std::vector<std::unique_ptr<InitialField>>
    pressureInjection(const initializer::parameters::InitializationParameters& parameters,
                      ConfigId config);

/// An acoustic wave travelling along x, through the time mirror of `parameters`, in the material
/// `materialData` (rho and the speed sqrt(lambda / rho) of the wave before the mirror).
std::vector<std::unique_ptr<InitialField>>
    acousticTravellingWaveITM(const CellMaterialData& materialData,
                              const AcousticTravellingWaveParametersITM& parameters,
                              ConfigId config);

} // namespace seissol::physics

#endif // SEISSOL_SRC_PHYSICS_SCENARIO_SCENARIOS_H_
