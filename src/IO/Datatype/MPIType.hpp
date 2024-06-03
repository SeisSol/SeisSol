#pragma once

#include "Datatype.hpp"
#include <memory>
#include <mpi.h>

namespace seissol::io::datatype {
MPI_Datatype convertToMPI(std::shared_ptr<Datatype> datatype, bool autocommit = true);
} // namespace seissol::io::datatype
