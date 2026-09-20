// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "FaceTypeCheck.h"

#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/Solver.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <sstream>
#include <utils/logger.h>

namespace seissol::initializer::internal {

namespace {

std::uint32_t faceTypeBit(FaceType faceType) {
  return std::uint32_t{1} << static_cast<std::uint8_t>(faceType);
}

std::uint32_t collectFaceTypes(LTS::Storage& storage) {
  std::uint32_t present = 0;

  const LayerMask ghostMask(Ghost);
  for (auto& layer : storage.leaves(ghostMask)) {
    const auto* cellInformation = layer.var<LTS::CellInformation>();
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        present |= faceTypeBit(cellInformation[cell].faceTypes[face]);
      }
    }
  }

  // A rank without a face of some type must not conclude that the type is absent; all ranks
  // have to reach the same verdict, or the run deadlocks instead of aborting.
  MPI_Allreduce(MPI_IN_PLACE, &present, 1, MPI_UINT32_T, MPI_BOR, Mpi::mpi.comm());

  return present;
}

} // namespace

void checkFaceTypeSupport(LTS::Storage& storage) {
  const auto present = collectFaceTypes(storage);

  std::stringstream problems;
  std::size_t problemCount = 0;

  for (const auto faceType : FaceTypes) {
    if ((present & faceTypeBit(faceType)) == 0) {
      continue;
    }

    const auto material = model::faceTypeSupport<model::MaterialT>(faceType);
    const auto solver = kernels::Solver::implementsFaceType(faceType);

    if (!material.supported) {
      ++problemCount;
      problems << "\n  " << faceTypeName(faceType) << ": not defined for material "
               << model::MaterialT::Text << " (" << material.reason << ")";
    } else if (!solver.supported) {
      ++problemCount;
      problems << "\n  " << faceTypeName(faceType) << ": not implemented by the solver used for "
               << model::MaterialT::Text << " (" << solver.reason << ")";
    }
  }

  if (problemCount > 0) {
    logError() << "The mesh contains" << problemCount
               << "boundary condition(s) that this build cannot handle:" << problems.str();
  }
}

} // namespace seissol::initializer::internal
