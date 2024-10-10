// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#pragma once

#include <Initializer/Tree/LTSTree.h>
#include <Initializer/Tree/Layer.h>
#include <Solver/Clustering/Communication/NeighborCluster.h>
#include <mpi.h>
namespace seissol::solver::clustering::communication {

template <typename T>
void haloCommunication(const HaloCommunication& comm,
                       seissol::initializer::Variable<T>& var,
                       seissol::initializer::LTSTree& tree,
                       MPI_Datatype datatype) {
  haloCommunication(comm, var.handle, tree, datatype);
}

void haloCommunication(const HaloCommunication& comm,
                       unsigned varIndex,
                       seissol::initializer::LTSTree& tree,
                       MPI_Datatype datatype);

} // namespace seissol::solver::clustering::communication
