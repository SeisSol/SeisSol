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
  TableWriter();

  void addQuantity(const TableQuantity& quantity);

  void addCellRaw(const void* data, std::size_t size);

  template <typename T>
  void addCell(const T& data) {
    addCellRaw(&data, sizeof(T));
  }

  std::shared_ptr<datatype::Datatype> getRowDatatype() const;

  void resetStorage();

  virtual std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() = 0;

  protected:
  std::string name;
  std::vector<TableQuantity> quantities;
  std::vector<char> rowstorage;
  std::size_t rowCount{0};
  std::size_t inrowPos{0};
};

} // namespace seissol::io::instance::point
