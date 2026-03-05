// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_RESULTWRITER_PARALLELHDF5RECEIVERWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_PARALLELHDF5RECEIVERWRITER_H_

#include <Eigen/Dense>
#include <array>
#include <cstdint>
#include <hdf5.h>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>

namespace seissol::writer {

class ParallelHdf5ReceiverWriter {
  public:
  ParallelHdf5ReceiverWriter(MPI_Comm comm,
                             const std::string& filename,
                             hsize_t totalReceivers,
                             hsize_t numVariables);

  void writeChunk(hsize_t timeOffset,
                  hsize_t timeCount,
                  const std::vector<std::uint64_t>& pointIds,
                  const std::vector<double>& data);

  void writeCoordinates(const std::vector<Eigen::Vector3d>& points);
  void flush();

  ~ParallelHdf5ReceiverWriter();

  private:
  static constexpr int Rank = 3;

  MPI_Comm comm_;
  hid_t fileId_;
  hid_t dsetId_;
  hid_t filespaceId_;

  std::array<hsize_t, Rank> dims_;
};

} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_PARALLELHDF5RECEIVERWRITER_H_
