// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Memory/Tree/LTSTree.h>
#include <cstddef>
#include <mpi.h>
#include <vector>
namespace seissol::initializer {
    
struct RemoteCellRegion {
  std::size_t count;
  int rank;
};

struct HaloStructure {
  std::vector<std::vector<RemoteCellRegion>> ghost;
  std::vector<std::vector<RemoteCellRegion>> copy;
};

template <typename T>
void haloCommunication(const HaloStructure& comm,
                       seissol::initializer::Variable<T>& var,
                       seissol::initializer::LTSTree& tree,
                       MPI_Datatype datatype) {
  haloCommunication(comm, var.handle, tree, datatype);
}

void haloCommunication(const HaloStructure& comm,
                       unsigned varIndex,
                       seissol::initializer::LTSTree& tree,
                       MPI_Datatype datatype);

} // namespace seissol::initializer
