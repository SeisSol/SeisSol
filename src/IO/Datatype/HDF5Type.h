// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_DATATYPE_HDF5TYPE_H_
#define SEISSOL_SRC_IO_DATATYPE_HDF5TYPE_H_

#include "Datatype.h"
#include <hdf5.h>
#include <memory>

namespace seissol::io::datatype {
hid_t convertToHdf5(const std::shared_ptr<Datatype>& datatype);
} // namespace seissol::io::datatype

#endif // SEISSOL_SRC_IO_DATATYPE_HDF5TYPE_H_
