// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_MEMORY_MEMORYCONTAINER_H_
#define SEISSOL_SRC_MEMORY_MEMORYCONTAINER_H_

#include <Config.h>
#include <Initializer/BasicTypedefs.h>
#include <Initializer/Typedefs.h>
#include <Memory/Descriptor/Boundary.h>
#include <Memory/Descriptor/DynamicRupture.h>
#include <Memory/Descriptor/LTS.h>
#include <Memory/MemoryAllocator.h>
#include <Memory/Tree/Backmap.h>
#include <Memory/Tree/Colormap.h>
#include <Memory/Tree/LTSTree.h>

namespace seissol::memory {

struct MemoryContainer {
  static std::vector<std::size_t> clusters(std::size_t maxCluster) {
    std::vector<std::size_t> clusterMap(maxCluster + 1);
    std::iota(clusterMap.begin(), clusterMap.end(), 0);
    return clusterMap;
  }

  MemoryContainer(std::size_t maxCluster, const std::vector<ConfigVariant>& configs)
      : colorMap(initializer::TraitLayer<ConfigVariant>(configs),
                 initializer::EnumLayer<std::size_t>(clusters(maxCluster)),
                 initializer::EnumLayer<HaloType>(
                     {HaloType::Ghost, HaloType::Copy, HaloType::Interior})) {}

  memory::ManagedAllocator allocator;

  seissol::initializer::LTSColorMap colorMap;

  seissol::initializer::LTSTree volume;
  seissol::initializer::LTSTree dynrup;
  seissol::initializer::LTSTree boundary;

  seissol::initializer::LTS wpdesc;
  std::shared_ptr<seissol::initializer::DynamicRupture> drdesc;
  seissol::initializer::Boundary bnddesc;

  seissol::initializer::StorageBackmap<4, seissol::initializer::LTSTree> clusterBackmap;
  seissol::initializer::StorageBackmap<1, seissol::initializer::LTSTree> ghostClusterBackmap;
  seissol::initializer::StorageBackmap<1, seissol::initializer::LTSTree> dynrupBackmap;

  CompoundGlobalData globalDataStorage;
};

} // namespace seissol::memory
#endif // SEISSOL_SRC_MEMORY_MEMORYCONTAINER_H_
