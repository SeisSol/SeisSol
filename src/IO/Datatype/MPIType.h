// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_DATATYPE_MPITYPE_H_
#define SEISSOL_SRC_IO_DATATYPE_MPITYPE_H_

#include "Datatype.h"
#include <memory>
#include <mpi.h>

namespace seissol::io::datatype {
MPI_Datatype convertToMPI(const std::shared_ptr<Datatype>& datatype, bool autocommit = true);
} // namespace seissol::io::datatype

#endif // SEISSOL_SRC_IO_DATATYPE_MPITYPE_H_
