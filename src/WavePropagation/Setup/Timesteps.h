// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SETUP_TIMESTEPS_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SETUP_TIMESTEPS_H_

#include "Core/Data/DataReader.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/TimeStepping/Layout/GlobalTimestep.h"

namespace seissol::initializer {

namespace parameters {
struct SeisSolParameters;
} // namespace parameters

/// The time steps of the cells by the CFL condition, with the materials of the material files of
/// the run (which `models` opens) and the order of the configuration of each cell.
GlobalTimestep computeTimesteps(const CellToVertexArray& cellToVertex,
                                const parameters::SeisSolParameters& seissolParams,
                                const data::ModelLoader& models);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_WAVEPROPAGATION_SETUP_TIMESTEPS_H_
