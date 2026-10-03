// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ConfigBoundaryCheck.h"

#include "Common/ConfigDispatch.h"
#include "Common/ConfigRegistry.h"
#include "Common/ConfigValue.h"
#include "Common/Constants.h"
#include "Equations/Datastructures.h"
#include "Initializer/BasicTypedefs.h"
#include "Kernels/Common.h"
#include "Memory/Descriptor/LTS.h"
#include "Memory/Tree/Layer.h"
#include "Model/Common.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <cstdint>
#include <mpi.h>
#include <sstream>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer::internal {

namespace {

constexpr std::uint32_t AnyFace = 1;
constexpr std::uint32_t DynamicRuptureFace = 2;

/// The faces between the cells of each pair of configurations, as a bit set of the kinds above,
/// indexed by the configuration of the cell times the number of configurations plus the one of
/// the neighbor; the same on all ranks.
std::vector<std::uint32_t> collectConfigBoundaries(LTS::Storage& storage) {
  const auto configCount = builtConfigCount();
  std::vector<std::uint32_t> boundaries(configCount * configCount);

  for (auto& layer : storage.leaves(Ghost)) {
    const auto config = layer.getIdentifier().config;
    const auto* cellInformation = layer.var<LTS::CellInformation>();
    for (std::size_t cell = 0; cell < layer.size(); ++cell) {
      for (std::size_t face = 0; face < Cell::NumFaces; ++face) {
        // without a neighbor, the face has no configuration on the other side
        const auto neighborConfig = cellInformation[cell].neighborConfigIds[face];
        if (neighborConfig != config && neighborConfig < configCount) {
          auto& boundary = boundaries[config * configCount + neighborConfig];
          boundary |= AnyFace;
          if (cellInformation[cell].faceTypes[face] == FaceType::DynamicRupture) {
            boundary |= DynamicRuptureFace;
          }
        }
      }
    }
  }

  MPI_Allreduce(MPI_IN_PLACE,
                boundaries.data(),
                static_cast<int>(boundaries.size()),
                Mpi::castToMpiType<std::uint32_t>(),
                MPI_BOR,
                Mpi::mpi.comm());
  return boundaries;
}

/// Whether the cells of the configurations `first` and `second` can be face neighbors by their
/// materials.
bool materialsCanNeighbor(ConfigId first, ConfigId second) {
  return dispatchConfig(first, [&](auto firstCfg) {
    return dispatchConfig(second, [&](auto secondCfg) {
      return model::CanNeighbor<model::MaterialOf<decltype(firstCfg)>,
                                model::MaterialOf<decltype(secondCfg)>>;
    });
  });
}

} // namespace

void checkConfigBoundaries(LTS::Storage& storage) {
  const auto boundaries = collectConfigBoundaries(storage);
  const auto configCount = builtConfigCount();

  std::stringstream problems;
  std::size_t problemCount = 0;
  for (ConfigId first = 0; first < configCount; ++first) {
    for (ConfigId second = first + 1; second < configCount; ++second) {
      const auto boundary =
          boundaries[first * configCount + second] | boundaries[second * configCount + first];
      if ((boundary & AnyFace) == 0) {
        continue;
      }

      const auto& firstValue = configValue(first);
      const auto& secondValue = configValue(second);
      const auto pair = "between " + configName(firstValue) + " and " + configName(secondValue);
      if (!materialsCanNeighbor(first, second)) {
        ++problemCount;
        problems << "\n  faces " << pair
                 << ": the materials pose the Riemann problem at their faces differently";
      }
      if (firstValue.numSimulations != secondValue.numSimulations) {
        ++problemCount;
        problems << "\n  faces " << pair << ": they fuse a different number of simulations";
      }
      if ((boundary & DynamicRuptureFace) != 0) {
        ++problemCount;
        problems << "\n  dynamic rupture faces " << pair << ": not supported";
      }
      if (isDeviceOn()) {
        ++problemCount;
        problems << "\n  faces " << pair << ": not supported on GPUs yet";
      }
    }
  }

  if (problemCount > 0) {
    logError() << "Some faces between cells of different configurations cannot be computed:"
               << problems.str();
  }
}

} // namespace seissol::initializer::internal
