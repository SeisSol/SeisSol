// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_CHECKPOINT_H_
#define SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_CHECKPOINT_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/IO/Instance/Checkpoint/CheckpointManager.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "WavePropagation/Output/SurfaceStorage.h"

namespace seissol::solver {

/// Registers what the checkpoints hold of the cells of `storage`, in the configuration `config`,
/// after their tree "lts": their unknowns, the anelastic ones if the configuration has them, and
/// their time integrals if the outputs integrate them.
void registerCellCheckpoint(io::instance::checkpoint::CheckpointManager& manager,
                            LTS::Storage& storage,
                            ConfigId config);

/// Registers what the checkpoints hold of the faces of the free surface of `storage`, after their
/// tree "surface": their displacements.
void registerSurfaceCheckpoint(io::instance::checkpoint::CheckpointManager& manager,
                               SurfaceLTS::Storage& storage);

} // namespace seissol::solver

#endif // SEISSOL_SRC_WAVEPROPAGATION_OUTPUT_CHECKPOINT_H_
