// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_HALO_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_HALO_H_

#include <Memory/Descriptor/LTS.h>
#include <cstddef>
#include <mpi.h>
#include <vector>

namespace seissol::initializer {
    
struct RemoteCellRegion {
  std::size_t tag{0};
  std::size_t localId{0};
  std::size_t remoteId{0};
  std::size_t count{0};
  int rank{-1};

  RemoteCellRegion(std::size_t tag, std::size_t localId, std::size_t remoteId, std::size_t count, int rank)
    :tag(tag), localId(localId), remoteId(remoteId), count(count), rank(rank) {}
};

struct HaloStructure {
  std::vector<std::vector<RemoteCellRegion>> ghost;
  std::vector<std::vector<RemoteCellRegion>> copy;
};

template <typename HandleT>
void haloCommunication(const HaloStructure& comm,
                       const HandleT& var,
                       LTS::Storage& storage,
                       MPI_Datatype datatype) {
  haloCommunication(comm, storage.info(var).index, storage, datatype);
}

template <typename StorageT>
void haloCommunication(const HaloStructure& comm,
                       LTS::Storage& storage,
                       MPI_Datatype datatype) {
  haloCommunication(comm, storage.info<StorageT>().index, storage, datatype);
}

void haloCommunication(const HaloStructure& comm,
                       std::size_t varIndex,
                       LTS::Storage& storage,
                       MPI_Datatype datatype);

} // namespace seissol::initializer
#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_HALO_H_
