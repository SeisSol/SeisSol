// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SOLVER_PLANEWAVEOPERATOR_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SOLVER_PLANEWAVEOPERATOR_H_

#include "Core/Common/ConfigRegistry.h"
#include "WavePropagation/Physics/CommonDatastructures.h"

#include <Eigen/Dense>
#include <complex>
#include <vector>

namespace seissol::solver {

/// The plane wave operator of the material `material` of the configuration `config` for the wave
/// vector `k`, as the solver of the configuration lays the quantities out: a square matrix over
/// the quantities of the material, column by column. The plane waves of the initial conditions
/// take it so (cf. physics::PlaneWaveOperator), without depending on the solvers.
std::vector<std::complex<double>>
    planeWaveOperator(ConfigId config, const model::Material& material, const Eigen::Vector3d& k);

} // namespace seissol::solver

#endif // SEISSOL_SRC_WAVEPROPAGATION_SOLVER_PLANEWAVEOPERATOR_H_
