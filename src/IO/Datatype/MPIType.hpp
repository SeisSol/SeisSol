// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_DATATYPE_MPITYPE_HPP_
#define SEISSOL_SRC_IO_DATATYPE_MPITYPE_HPP_

#include "Datatype.hpp"
#include <memory>
#include <mpi.h>

namespace seissol::io::datatype {
MPI_Datatype convertToMPI(std::shared_ptr<Datatype> datatype, bool autocommit = true);
} // namespace seissol::io::datatype

#endif // SEISSOL_SRC_IO_DATATYPE_MPITYPE_HPP_
