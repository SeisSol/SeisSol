// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_DATATYPE_HDF5TYPE_HPP_
#define SEISSOL_SRC_IO_DATATYPE_HDF5TYPE_HPP_

#include "Datatype.hpp"
#include <hdf5.h>
#include <memory>

namespace seissol::io::datatype {
hid_t convertToHdf5(std::shared_ptr<Datatype> datatype);
} // namespace seissol::io::datatype

#endif // SEISSOL_SRC_IO_DATATYPE_HDF5TYPE_HPP_
