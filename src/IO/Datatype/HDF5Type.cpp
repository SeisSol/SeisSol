// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "HDF5Type.h"

#include "Datatype.h"
#include <H5Ipublic.h>
#include <H5Tpublic.h>
#include <H5public.h>
#include <H5version.h>
#include <memory>
#include <vector>

#include "utils/logger.h"

namespace {
hid_t _eh(hid_t data) {
  if (data < 0) {
    logError() << "HDF5 error concerning datatype creation:" << data;
  }
  return data;
}

hid_t convertOpaque(const seissol::io::datatype::Datatype& datatype) {
  return _eh(H5Tcreate(H5T_OPAQUE, datatype.size()));
}

hid_t convertString(const seissol::io::datatype::Datatype& datatype) {
  const hid_t handle = _eh(H5Tcopy(H5T_C_S1));
  _eh(H5Tset_size(handle, datatype.size()));
  return handle;
}

hid_t convertArray(const seissol::io::datatype::ArrayDatatype& datatype) {
  std::vector<hsize_t> h5dims(datatype.dimensions().begin(), datatype.dimensions().end());
  return _eh(H5Tarray_create(
      seissol::io::datatype::convertToHdf5(datatype.base()), h5dims.size(), h5dims.data()));
}

hid_t convertStruct(const seissol::io::datatype::StructDatatype& datatype) {
  const hid_t handle = _eh(H5Tcreate(H5T_COMPOUND, datatype.size()));
  for (const auto& member : datatype.members()) {
    _eh(H5Tinsert(handle,
                  member.name.c_str(),
                  member.offset,
                  seissol::io::datatype::convertToHdf5(member.datatype)));
  }
  return handle;
}

hid_t convertInteger(const seissol::io::datatype::IntegerDatatype& datatype) {
  if (datatype.size() == 1) {
    return datatype.sign() ? H5T_STD_I8LE : H5T_STD_U8LE;
  } else if (datatype.size() == 2) {
    return datatype.sign() ? H5T_STD_I16LE : H5T_STD_U16LE;
  } else if (datatype.size() == 4) {
    return datatype.sign() ? H5T_STD_I32LE : H5T_STD_U32LE;
  } else if (datatype.size() == 8) {
    return datatype.sign() ? H5T_STD_I64LE : H5T_STD_U64LE;
  } else {
    const hid_t copy = _eh(H5Tcopy(H5T_STD_I32LE));
    _eh(H5Tset_sign(copy, datatype.sign() ? H5T_SGN_2 : H5T_SGN_NONE));
    _eh(H5Tset_size(copy, datatype.size()));
    return copy;
  }
}
} // namespace

namespace seissol::io::datatype {
hid_t convertToHdf5(const std::shared_ptr<Datatype>& datatype) {
  if (dynamic_cast<const ArrayDatatype*>(datatype.get()) != nullptr) {
    return convertArray(dynamic_cast<const ArrayDatatype&>(*datatype));
  } else if (dynamic_cast<const StructDatatype*>(datatype.get()) != nullptr) {
    return convertStruct(dynamic_cast<const StructDatatype&>(*datatype));
  } else if (dynamic_cast<const IntegerDatatype*>(datatype.get()) != nullptr) {
    return convertInteger(dynamic_cast<const IntegerDatatype&>(*datatype));
  } else if (dynamic_cast<const F32Datatype*>(datatype.get()) != nullptr) {
    return H5T_NATIVE_FLOAT;
  } else if (dynamic_cast<const F64Datatype*>(datatype.get()) != nullptr) {
    return H5T_NATIVE_DOUBLE;
  } else if (dynamic_cast<const F80Datatype*>(datatype.get()) != nullptr) {
    return H5T_NATIVE_LDOUBLE;
  } else if (dynamic_cast<const OpaqueDatatype*>(datatype.get()) != nullptr) {
    return convertOpaque(dynamic_cast<const OpaqueDatatype&>(*datatype));
  } else if (dynamic_cast<const StringDatatype*>(datatype.get()) != nullptr) {
    return convertString(dynamic_cast<const StringDatatype&>(*datatype));
  }
  return H5T_NATIVE_INT;
}
} // namespace seissol::io::datatype
