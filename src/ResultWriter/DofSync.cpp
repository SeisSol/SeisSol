// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DofSync.h"

#include "Modules/Modules.h"
#include "Parallel/MPI.h"

#include <limits>

namespace seissol::writer {

void DofSync::setup(const initializer::MeshLayout& layout, LTS::Storage* storage) {
  layout_ = layout;
  storage_ = storage;
}

void DofSync::syncDofs(double time) {
  if (time_ == time) {
    return;
  }

  time_ = time;

  logInfo() << "Synchronizing halo dofs for output.";

  MPI_Barrier(seissol::Mpi::mpi.comm());

  // TODO: unify with the haloCommunication functions

  const auto varIndexCopy = storage_->info<LTS::Dofs>().index;
  const auto varIndexGhost = storage_->info<LTS::DofsHalo>().index;

  std::vector<MPI_Request> requests;
  for (auto& layer : storage_->leaves(Interior)) {
    const auto varIndex =
        layer.getIdentifier().halo == HaloType::Ghost ? varIndexGhost : varIndexCopy;
    auto* data = reinterpret_cast<uint8_t*>(layer.varUntyped(varIndex));
    const auto elemsize = storage_->info(varIndex).bytesLayer(layer.getIdentifier());

    if (layer.getIdentifier().halo == HaloType::Ghost) {
      for (const auto& remote : layout_[layer.id()].regions) {
        auto& request = requests.emplace_back(MPI_REQUEST_NULL);
        MPI_Irecv(data,
                  remote.count * elemsize,
                  MPI_BYTE,
                  remote.rank,
                  remote.tag,
                  seissol::Mpi::mpi.comm(),
                  &request);
        data += elemsize * remote.count;
      }
    } else if (layer.getIdentifier().halo == HaloType::Copy) {
      for (const auto& remote : layout_[layer.id()].regions) {
        auto& request = requests.emplace_back(MPI_REQUEST_NULL);
        MPI_Isend(data,
                  remote.count * elemsize,
                  MPI_BYTE,
                  remote.rank,
                  remote.tag,
                  seissol::Mpi::mpi.comm(),
                  &request);
        data += elemsize * remote.count;
      }
    }
  }

  MPI_Waitall(requests.size(), requests.data(), MPI_STATUSES_IGNORE);

  logInfo() << "Synchronizing halo dofs for output. Done.";
}

} // namespace seissol::writer
