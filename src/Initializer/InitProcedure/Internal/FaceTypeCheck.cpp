// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "FaceTypeCheck.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Equations/Setup.h" // IWYU pragma: keep
#include "Initializer/BasicTypedefs.h"
#include "Initializer/Parameters/InitializationParameters.h"
#include "Kernels/SolverSelector.h"
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
#include <vector>

namespace seissol::initializer::internal {

namespace {

std::uint32_t faceTypeBit(FaceType faceType) {
  return std::uint32_t{1} << static_cast<std::uint8_t>(faceType);
}

struct FaceTypeCensus {
  std::uint32_t present{0};
  std::uint32_t cellRejected{0};
};

/// The face types of the cells of each configuration, by the id of the configuration.
std::vector<FaceTypeCensus> collectFaceTypes(LTS::Storage& storage) {
  std::vector<FaceTypeCensus> census(builtConfigCount());

  const LayerMask ghostMask(Ghost);
  for (auto& layer : storage.leaves(ghostMask)) {
    dispatchConfig(layer.getIdentifier().config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      auto& configCensus = census[configIdOf<Cfg>()];
      const auto* cellInformation = layer.var<LTS::CellInformation>();
      const auto* materialData = layer.var<LTS::MaterialData>(Cfg());
      for (std::size_t cell = 0; cell < layer.size(); ++cell) {
        for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
          const auto faceType = cellInformation[cell].faceTypes[face];
          configCensus.present |= faceTypeBit(faceType);
          // NOLINTNEXTLINE
          if (!model::faceTypeCellAdmissible<model::MaterialOf<Cfg>>(faceType,
                                                                     materialData[cell])) {
            configCensus.cellRejected |= faceTypeBit(faceType);
          }
        }
      }
    });
  }

  // A rank without a face of some type must not conclude that the type is absent; all ranks
  // have to reach the same verdict, or the run deadlocks instead of aborting.
  std::vector<std::uint32_t> reduced;
  for (const auto& configCensus : census) {
    reduced.push_back(configCensus.present);
    reduced.push_back(configCensus.cellRejected);
  }
  MPI_Allreduce(MPI_IN_PLACE,
                reduced.data(),
                static_cast<int>(reduced.size()),
                Mpi::castToMpiType<std::uint32_t>(),
                MPI_BOR,
                Mpi::mpi.comm());
  for (std::size_t config = 0; config < census.size(); ++config) {
    census[config].present = reduced[2 * config];
    census[config].cellRejected = reduced[2 * config + 1];
  }

  return census;
}

} // namespace

void checkFaceTypeSupport(LTS::Storage& storage, parameters::InitializationType scenarioType) {
  const auto census = collectFaceTypes(storage);

  std::stringstream problems;
  std::size_t problemCount = 0;

  // the cells of a configuration have to support the face types they have
  for (ConfigId config = 0; config < census.size(); ++config) {
    dispatchConfig(config, [&](auto cfg) {
      using Cfg = decltype(cfg);
      using MaterialT = model::MaterialOf<Cfg>;
      const auto& configCensus = census[config];
      for (const auto faceType : FaceTypes) {
        if ((configCensus.present & faceTypeBit(faceType)) == 0) {
          continue;
        }

        const auto material = model::faceTypeSupport<MaterialT>(faceType);
        const auto solver = kernels::SolverOf<Cfg>::implementsFaceType(faceType);

        if (!material.supported) {
          ++problemCount;
          problems << "\n  " << faceTypeName(faceType) << ": not defined for material "
                   << MaterialT::Text << " (" << material.reason << ")";
        } else if (!solver.supported) {
          ++problemCount;
          problems << "\n  " << faceTypeName(faceType)
                   << ": not implemented by the solver used for " << MaterialT::Text << " ("
                   << solver.reason << ")";
        } else if ((configCensus.cellRejected & faceTypeBit(faceType)) != 0) {
          ++problemCount;
          const auto requirement = model::faceTypeCellRequirement<MaterialT>(faceType);
          problems << "\n  " << faceTypeName(faceType)
                   << ": present on a cell that does not qualify (" << requirement.reason << ")";
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
    });
  }

  if (problemCount > 0) {
    logError() << "The mesh contains" << problemCount
               << "boundary condition(s) that this build cannot handle:" << problems.str();
  }
}

} // namespace seissol::initializer::internal
