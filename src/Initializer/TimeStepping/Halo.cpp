// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Halo.h"

#include <Initializer/BasicTypedefs.h>
#include <Memory/Tree/LTSTree.h>
#include <Memory/Tree/Layer.h>
#include <Parallel/MPI.h>
#include <mpi.h>
namespace seissol::initializer {

void haloCommunication(const MeshLayout& comm,
                       std::size_t varIndex,
                       LTS::Storage& storage,
                       MPI_Datatype datatype) {
  std::vector<MPI_Request> requests;
  for (auto& layer : storage.leaves(Interior)) {
    auto* data = reinterpret_cast<uint8_t*>(layer.varUntyped(varIndex));
    const auto elemsize = storage.info(varIndex).bytesLayer(layer.getIdentifier());

    if (layer.getIdentifier().halo == HaloType::Ghost) {
      for (const auto& remote : comm[layer.id()].regions) {
        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Irecv(data,
                  remote.count,
                  datatype,
                  remote.rank,
                  remote.tag,
                  seissol::MPI::mpi.comm(),
                  &request);
        requests.emplace_back(request);
        data += elemsize * remote.count;
      }
    } else if (layer.getIdentifier().halo == HaloType::Copy) {
      for (const auto& remote : comm[layer.id()].regions) {
        MPI_Request request = MPI_REQUEST_NULL;
        MPI_Isend(data,
                  remote.count,
                  datatype,
                  remote.rank,
                  remote.tag,
                  seissol::MPI::mpi.comm(),
                  &request);
        requests.emplace_back(request);
        data += elemsize * remote.count;
      }
    }
  }
  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);
}

} // namespace seissol::initializer
