// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_SETUP_TIMESTEPS_H_
#define SEISSOL_SRC_WAVEPROPAGATION_SETUP_TIMESTEPS_H_

#include "Core/Common/RunConfigs.h"
#include "Core/Data/DataReader.h"
#include "Core/Geometry/CellToVertexArray.h"
#include "Core/TimeStepping/Layout/GlobalTimestep.h"
#include "Core/TimeStepping/LtsParameters.h"
#include "WavePropagation/WavePropagationParameters.h"

#include <string>
#include <vector>

namespace seissol::initializer {

/// The time steps of the cells by the CFL condition of `timeStepping`, with the materials of the
/// files `materialFileNames` of the simulations of the run (which `models` opens), sampled as
/// `wavePropagation` says, and the order of the configuration of each cell.
GlobalTimestep computeTimesteps(const CellToVertexArray& cellToVertex,
                                const parameters::RunConfigs& runConfigs,
                                const parameters::WavePropagationParameters& wavePropagation,
                                const parameters::TimeSteppingParameters& timeStepping,
                                const std::vector<std::string>& materialFileNames,
                                const data::ModelLoader& models);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_WAVEPROPAGATION_SETUP_TIMESTEPS_H_
