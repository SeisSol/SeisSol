#pragma once

#include "Datatype.hpp"
#include <hdf5.h>
#include <memory>

namespace seissol::io::datatype {
hid_t convertToHdf5(std::shared_ptr<Datatype> datatype);
} // namespace seissol::io::datatype
