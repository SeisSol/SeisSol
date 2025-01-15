// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MPIType.h"
#include <IO/Datatype/Datatype.h>
#include <cstddef>
#include <iosfwd>
#include <memory>
#include <mpi.h>
#include <string>
#include <unordered_map>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {
// (needs to stay non-const)
// NOLINTNEXTLINE
static std::unordered_map<std::string, MPI_Datatype> autocommitRegistry;
} // namespace

namespace {
MPI_Datatype convertOpaque(const seissol::io::datatype::Datatype& datatype) {
  MPI_Datatype type = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(datatype.size(), MPI_BYTE, &type);
  return type;
}

MPI_Datatype convertString(const seissol::io::datatype::Datatype& datatype) {
  MPI_Datatype type = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(datatype.size(), MPI_CHAR, &type);
  return type;
}

MPI_Datatype convertArray(const seissol::io::datatype::ArrayDatatype& datatype) {
  std::size_t total = 1;
  for (const std::size_t dim : datatype.dimensions()) {
    total *= dim;
  }
  MPI_Datatype type = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(total, seissol::io::datatype::convertToMPI(datatype.base(), false), &type);
  return type;
}

MPI_Datatype convertStruct(const seissol::io::datatype::StructDatatype& datatype) {
  MPI_Datatype type = MPI_DATATYPE_NULL;
  std::vector<MPI_Datatype> subtypes;
  std::vector<MPI_Aint> suboffsets;
  std::vector<int> subsizes;
  for (const auto& member : datatype.members()) {
    subsizes.push_back(1);
    suboffsets.push_back(member.offset);
    subtypes.push_back(seissol::io::datatype::convertToMPI(member.datatype, false));
  }
  MPI_Type_create_struct(
      subtypes.size(), subsizes.data(), suboffsets.data(), subtypes.data(), &type);
  return type;
}

MPI_Datatype convertInteger(const seissol::io::datatype::IntegerDatatype& datatype) {
  if (datatype.size() == 1) {
    return datatype.sign() ? MPI_INT8_T : MPI_UINT8_T;
  } else if (datatype.size() == 2) {
    return datatype.sign() ? MPI_INT16_T : MPI_UINT16_T;
  } else if (datatype.size() == 4) {
    return datatype.sign() ? MPI_INT32_T : MPI_UINT32_T;
  } else if (datatype.size() == 8) {
    return datatype.sign() ? MPI_INT64_T : MPI_UINT64_T;
  } else {
    return MPI_BYTE; // TODO
  }
}
} // namespace

namespace seissol::io::datatype {
MPI_Datatype convertToMPI(const std::shared_ptr<Datatype>& datatype, bool autocommit) {
  std::string serialized;
  if (autocommit) {
    std::stringstream sstr;
    {
      YAML::Emitter output(sstr);
      output << datatype->serialize();
    }
    serialized = sstr.str();
  }
  if (autocommit) {
    if (autocommitRegistry.find(serialized) != autocommitRegistry.end()) {
      return autocommitRegistry.at(serialized);
    }
  }
  MPI_Datatype type = MPI_DATATYPE_NULL;
  bool needsCommit = false;
  if (dynamic_cast<const ArrayDatatype*>(datatype.get()) != nullptr) {
    type = convertArray(dynamic_cast<const ArrayDatatype&>(*datatype));
    needsCommit = true;
  } else if (dynamic_cast<const StructDatatype*>(datatype.get()) != nullptr) {
    type = convertStruct(dynamic_cast<const StructDatatype&>(*datatype));
    needsCommit = true;
  } else if (dynamic_cast<const IntegerDatatype*>(datatype.get()) != nullptr) {
    type = convertInteger(dynamic_cast<const IntegerDatatype&>(*datatype));
    needsCommit = false;
  } else if (dynamic_cast<const F32Datatype*>(datatype.get()) != nullptr) {
    type = MPI_FLOAT;
    needsCommit = false;
  } else if (dynamic_cast<const F64Datatype*>(datatype.get()) != nullptr) {
    type = MPI_DOUBLE;
    needsCommit = false;
  } else if (dynamic_cast<const F80Datatype*>(datatype.get()) != nullptr) {
    type = MPI_LONG_DOUBLE;
    needsCommit = false;
  } else if (dynamic_cast<const OpaqueDatatype*>(datatype.get()) != nullptr) {
    type = convertOpaque(dynamic_cast<const OpaqueDatatype&>(*datatype));
    needsCommit = true;
  } else if (dynamic_cast<const StringDatatype*>(datatype.get()) != nullptr) {
    type = convertString(dynamic_cast<const StringDatatype&>(*datatype));
    needsCommit = true;
  }
  if (needsCommit && autocommit) {
    MPI_Type_commit(&type);
    autocommitRegistry[serialized] = type;
  }
  return type;
}
} // namespace seissol::io::datatype
