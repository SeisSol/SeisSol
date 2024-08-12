#pragma once

#include <IO/Datatype/Datatype.hpp>
#include <IO/Writer/Writer.hpp>
#include <functional>
#include <memory>
#include <string>

#include "utils/logger.h"

namespace seissol::io::instance::point {

struct TableQuantity {
  std::string name;
  std::shared_ptr<datatype::Datatype> datatype;
};

class TableWriter {
  public:
  TableWriter() = default;

  void addQuantity(const TableQuantity& quantity) { quantities.push_back(quantity); }

  void addCellRaw(const void* data, std::size_t size) {
    const auto position = rowstorage.size();
    rowstorage.resize(position + size);
    std::copy_n(reinterpret_cast<const char*>(data), size, rowstorage.begin() + position);
  }

  template <typename T>
  void addCell(const T& data) {
    addCellRaw(&data, sizeof(T));
  }

  std::shared_ptr<datatype::Datatype> getRowDatatype() const {
    std::vector<datatype::StructDatatype::MemberInfo> memberInfo(quantities.size());
    std::size_t offset = 0;
    for (std::size_t i = 0; i < memberInfo.size(); ++i) {
      memberInfo[i] =
          datatype::StructDatatype::MemberInfo{quantities[i].name, offset, quantities[i].datatype};
      offset += quantities[i].datatype->size();
    }
    return std::make_shared<datatype::StructDatatype>(memberInfo);
  }

  virtual std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() = 0;

  protected:
  std::string name;
  std::vector<TableQuantity> quantities;
  std::vector<char> rowstorage;
};

} // namespace seissol::io::instance::point
