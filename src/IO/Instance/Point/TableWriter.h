// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_POINT_TABLEWRITER_H_
#define SEISSOL_SRC_IO_INSTANCE_POINT_TABLEWRITER_H_

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Writer/Writer.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::instance::point {

struct TableQuantity {
  std::string name;
  std::shared_ptr<datatype::Datatype> datatype;
};

class TableWriter {
  public:
  virtual ~TableWriter() = default;
  explicit TableWriter(std::string name);

  //! @brief What the file this writes is named after.
  [[nodiscard]] const std::string& name() const;

  void addQuantity(const TableQuantity& quantity);

  //! @brief Declares a column holding a value of type @p T.
  template <typename T>
  void addColumn(const std::string& name) {
    addQuantity(TableQuantity{name, datatype::inferDatatype<T>()});
  }

  /**
   * @brief Declares a column holding text, of at most @p maxLength bytes.
   *
   * A row has the same size wherever it is written, so text is held in a field of its own length
   * rather than one that follows the value.
   */
  void addTextColumn(const std::string& name, std::size_t maxLength);

  void addCellRaw(const void* data, std::size_t size);

  template <typename T>
  void addCell(const T& data) {
    addCellRaw(&data, sizeof(T));
  }

  //! @brief Adds the next cell as text, cut or padded to the length its column was declared with.
  void addText(const std::string& value);

  [[nodiscard]] std::shared_ptr<datatype::Datatype> getRowDatatype() const;

  void resetStorage();

  virtual std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() = 0;

  protected:
  std::string name_;
  std::vector<TableQuantity> quantities_;
  std::vector<char> rowstorage_;
  std::size_t rowCount_{0};
  std::size_t inrowPos_{0};
};

} // namespace seissol::io::instance::point

#endif // SEISSOL_SRC_IO_INSTANCE_POINT_TABLEWRITER_H_
