// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
#define SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_

#include <Geometry/MeshReader.h>
#include <Initializer/BasicTypedefs.h>
#include <mpi.h>
#include <stdexcept>
#include <vector>

namespace seissol::initializer::internal {

struct TransferRegion {
  int rank;
  std::size_t start;
  std::size_t size;
  int tag;
};

struct ClusterLayout {
  std::vector<std::size_t> interior;
  std::vector<std::size_t> copy;
  std::vector<std::size_t> ghost;
  std::vector<TransferRegion> copyRegions;
  std::vector<TransferRegion> ghostRegions;

  [[nodiscard]] const std::vector<std::size_t>& cells(HaloType halo) const {
    switch (halo) {
    case seissol::HaloType::Interior:
      return interior;
    case seissol::HaloType::Copy:
      return copy;
    case seissol::HaloType::Ghost:
      return ghost;
    default:
      throw std::runtime_error("Invalid LayerType at setup");
    }
  }

  [[nodiscard]] bool empty() const { return interior.empty() && copy.empty() && ghost.empty(); }

  [[nodiscard]] bool interiorOnly() const {
    return !interior.empty() && copy.empty() && ghost.empty();
  }
};

std::vector<ClusterLayout> layoutCells(const std::vector<std::size_t>& color,
                                       const std::vector<std::size_t>& ghostColor,
                                       std::size_t maxColors,
                                       const geometry::MeshReader& meshReader);

} // namespace seissol::initializer::internal
#endif // SEISSOL_SRC_INITIALIZER_INITPROCEDURE_INTERNAL_MESHLAYOUT_H_
