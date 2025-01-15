// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "TableWriter.h"

#include <IO/Datatype/Datatype.h>
#include <algorithm>
#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::io::instance::point {

TableWriter::TableWriter() = default;

void TableWriter::addQuantity(const TableQuantity& quantity) { quantities.push_back(quantity); }

void TableWriter::addCellRaw(const void* data, std::size_t size) {
  const auto position = rowstorage.size();
  rowstorage.resize(position + size);
  std::copy_n(reinterpret_cast<const char*>(data), size, rowstorage.begin() + position);
  inrowPos += 1;
  if (inrowPos == quantities.size()) {
    rowCount += 1;
    inrowPos = 0;
  }
}

std::shared_ptr<datatype::Datatype> TableWriter::getRowDatatype() const {
  std::vector<datatype::StructDatatype::MemberInfo> memberInfo(quantities.size());
  std::size_t offset = 0;
  for (std::size_t i = 0; i < memberInfo.size(); ++i) {
    memberInfo[i] =
        datatype::StructDatatype::MemberInfo{quantities[i].name, offset, quantities[i].datatype};
    offset += quantities[i].datatype->size();
  }
  return std::make_shared<datatype::StructDatatype>(memberInfo);
}

void TableWriter::resetStorage() {
  inrowPos = 0;
  rowCount = 0;
  rowstorage.clear();
}

} // namespace seissol::io::instance::point
