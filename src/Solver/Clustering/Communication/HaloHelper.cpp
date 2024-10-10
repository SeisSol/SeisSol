// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "HaloHelper.h"

#include <Initializer/Tree/LTSTree.h>
#include <Initializer/Tree/Layer.h>
#include <Parallel/MPI.h>
#include <mpi.h>
namespace seissol::solver::clustering::communication {

void haloCommunication(const HaloCommunication& comm,
                       unsigned varIndex,
                       seissol::initializer::LTSTree& tree,
                       MPI_Datatype datatype) {
  int datatypesize = 0;
  MPI_Type_size(datatype, &datatypesize);
  const std::size_t elemsize = datatypesize;

  std::vector<MPI_Request> requests;
  for (auto& layer : tree.leaves(Interior)) {
    int cluster = layer.getClusterId();
    auto* data = reinterpret_cast<char*>(layer.varUntyped(varIndex));
    if (layer.getLayerType() == Ghost) {
      for (const auto& remote : comm.ghost.at(cluster)) {
        MPI_Request request;
        MPI_Irecv(data,
                  remote.size,
                  datatype,
                  remote.rank,
                  remote.tag,
                  seissol::MPI::mpi.comm(),
                  &request);
        requests.emplace_back(request);
        data += elemsize * remote.size;
      }
    } else if (layer.getLayerType() == Copy) {
      for (const auto& remote : comm.copy.at(cluster)) {
        MPI_Request request;
        MPI_Isend(data,
                  remote.size,
                  datatype,
                  remote.rank,
                  remote.tag,
                  seissol::MPI::mpi.comm(),
                  &request);
        requests.emplace_back(request);
        data += elemsize * remote.size;
      }
    }
  }
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

} // namespace seissol::solver::clustering::communication
