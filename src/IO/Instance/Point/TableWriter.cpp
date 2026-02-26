// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "TableWriter.h"

#include "IO/Datatype/Datatype.h"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <vector>

namespace seissol::io::instance::point {

TableWriter::TableWriter() = default;

void TableWriter::addQuantity(const TableQuantity& quantity) { quantities_.push_back(quantity); }

void TableWriter::addCellRaw(const void* data, std::size_t size) {
  const auto position = rowstorage_.size();
  rowstorage_.resize(position + size);
  std::copy_n(reinterpret_cast<const char*>(data), size, rowstorage_.begin() + position);
  inrowPos_ += 1;
  if (inrowPos_ == quantities_.size()) {
    rowCount_ += 1;
    inrowPos_ = 0;
  }
}

std::shared_ptr<datatype::Datatype> TableWriter::getRowDatatype() const {
  std::vector<datatype::StructDatatype::MemberInfo> memberInfo(quantities_.size());
  std::size_t offset = 0;
  for (std::size_t i = 0; i < memberInfo.size(); ++i) {
    memberInfo[i] =
        datatype::StructDatatype::MemberInfo{quantities_[i].name, offset, quantities_[i].datatype};
    offset += quantities_[i].datatype->size();
  }
  return std::make_shared<datatype::StructDatatype>(memberInfo);
}

void TableWriter::resetStorage() {
  inrowPos_ = 0;
  rowCount_ = 0;
  rowstorage_.clear();
}

} // namespace seissol::io::instance::point
