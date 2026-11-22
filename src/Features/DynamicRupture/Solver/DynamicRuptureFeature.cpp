// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRuptureFeature.h"

#include "Core/Geometry/MeshReader.h"
#include "Core/IO/Instance/Checkpoint/CheckpointManager.h"
#include "Core/Kernels/GlobalDataStore.h"
#include "Core/Memory/Descriptor/DynamicRupture.h"
#include "Core/Memory/Descriptor/LTS.h"
#include "Core/Memory/Tree/Layer.h"
#include "Core/Run/Context.h"
#include "Features/DynamicRupture/Output/OutputManager.h"
#include "Features/DynamicRupture/Solver/FrictionLaws/FrictionSolver.h"

#include <cstddef>
#include <vector>

namespace seissol::dr {

DynamicRuptureFeature::DynamicRuptureFeature(
    RunContext& run,
    DynamicRupture::Storage& storage,
    const DynamicRupture& descriptor,
    const geometry::MeshReader& mesh,
    std::size_t copyCount,
    GlobalDataStore& globalData,
    const friction_law::FrictionSolverFactory& frictionLaw,
    const friction_law::FrictionSolverFactory& frictionLawDevice,
    output::OutputManager* outputManager,
    double endTime)
    : storage_(storage), descriptor_(descriptor), mesh_(mesh), copyCount_(copyCount),
      clusters_(run, storage, globalData, frictionLaw, frictionLawDevice, outputManager, endTime) {}

void DynamicRuptureFeature::registerCheckpoint(
    io::instance::checkpoint::CheckpointManager& checkpoint, LTS::Storage& /*cells*/) {
  std::vector<std::size_t> faceIdentifiers(storage_.size(initializer::LayerMask(Ghost)));
  const auto* faceInformation = storage_.var<DynamicRupture::FaceInformation>();
  const auto& fault = mesh_.getFault();

#pragma omp parallel for schedule(static)
  for (std::size_t i = 0; i < faceIdentifiers.size(); ++i) {
    const auto& face = fault[faceInformation[i].meshFace];
    // take the positive cell and side as fault face identifier
    // (should result in roughly twice as large numbers as when indexing all faces; cf.
    // handshake theorem)
    faceIdentifiers[i] = (face.globalId * 4 + face.side) * copyCount_ + faceInformation[i].copy;
  }
  checkpoint.registerTree("dynrup", storage_, faceIdentifiers);
  descriptor_.registerCheckpointVariables(checkpoint, storage_);
}

} // namespace seissol::dr
