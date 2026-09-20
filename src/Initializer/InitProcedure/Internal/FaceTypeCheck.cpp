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
#include "Physics/Scenario/Registry.h"

#include <array>
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

struct FaceTypeCensus {
  std::uint32_t present{0};
  std::uint32_t cellRejected{0};
};

FaceTypeCensus collectFaceTypes(LTS::Storage& storage) {
  FaceTypeCensus census;

  const LayerMask ghostMask(Ghost);
  for (auto& layer : storage.leaves(ghostMask)) {
    const auto* cellInformation = layer.var<LTS::CellInformation>();
    const auto* materialData = layer.var<LTS::MaterialData>();
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        const auto faceType = cellInformation[cell].faceTypes[face];
        census.present |= faceTypeBit(faceType);
        if (!model::faceTypeCellAdmissible<model::MaterialT>(faceType, materialData[cell])) {
          census.cellRejected |= faceTypeBit(faceType);
        }
      }
    }
  }

  // A rank without a face of some type must not conclude that the type is absent; all ranks
  // have to reach the same verdict, or the run deadlocks instead of aborting.
  std::array<std::uint32_t, 2> reduced{census.present, census.cellRejected};
  MPI_Allreduce(
      MPI_IN_PLACE, reduced.data(), reduced.size(), MPI_UINT32_T, MPI_BOR, Mpi::mpi.comm());
  census.present = reduced[0];
  census.cellRejected = reduced[1];

  return census;
}

} // namespace

void checkFaceTypeSupport(LTS::Storage& storage, parameters::InitializationType scenarioType) {
  const auto census = collectFaceTypes(storage);

  std::stringstream problems;
  std::size_t problemCount = 0;

  for (const auto faceType : FaceTypes) {
    if ((census.present & faceTypeBit(faceType)) == 0) {
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
    } else if ((census.cellRejected & faceTypeBit(faceType)) != 0) {
      ++problemCount;
      const auto requirement = model::faceTypeCellRequirement<model::MaterialT>(faceType);
      problems << "\n  " << faceTypeName(faceType) << ": present on a cell that does not qualify ("
               << requirement.reason << ")";
    } else if (faceType == FaceType::Analytical) {
      const auto scenario = physics::scenario::analyticalBoundaryAvailability(scenarioType);
      if (!scenario.available) {
        ++problemCount;
        problems << "\n  " << faceTypeName(faceType) << ": the configured scenario "
                 << physics::scenario::name(scenarioType) << " cannot serve it, "
                 << scenario.reason;
      }
    }
  }

  if (problemCount > 0) {
    logError() << "The mesh contains" << problemCount
               << "boundary condition(s) that this build cannot handle:" << problems.str();
  }
}

} // namespace seissol::initializer::internal
