// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_PARALLEL_RUNTIME_STREAMMPI_H_
#define SEISSOL_SRC_PARALLEL_RUNTIME_STREAMMPI_H_

#include <memory>
#include <mpi.h>
#include <vector>
namespace seissol {

struct MPIStreamData;

class MPIStream {
  public:
  explicit MPIStream(void* stream);

  void send(void* data,
            std::size_t count,
            MPI_Datatype datatype,
            int other,
            int tag,
            MPI_Request request);
  void recv(void* data,
            std::size_t count,
            MPI_Datatype datatype,
            int other,
            int tag,
            MPI_Request request);

  void commit();
  void complete();

  private:
  void init();

  void* stream_;
  MPI_Comm comm_;
  std::shared_ptr<MPIStreamData> data_;
  std::shared_ptr<std::vector<MPI_Request>> requests_;
};

} // namespace seissol
#endif // SEISSOL_SRC_PARALLEL_RUNTIME_STREAMMPI_H_
