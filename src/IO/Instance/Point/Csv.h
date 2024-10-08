// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_IO_INSTANCE_POINT_CSV_H_
#define SEISSOL_SRC_IO_INSTANCE_POINT_CSV_H_

#include <IO/Datatype/Datatype.h>
#include <IO/Instance/Point/TableWriter.h>
#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Data.h>
#include <memory>
#include <sstream>
#include <string>

namespace seissol::io::instance::point {

class Csv : public TableWriter {
  public:
  Csv();

  std::string quote(const std::string& str) const;

  std::ostringstream& quote(std::ostringstream& stream, const std::string& str) const;

  std::string header() const;

  std::string rows() const;

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() override;

  private:
  std::string rowcache;
  char delimiterStr{';'};
  char quoteStr{'\"'};
  char newlineStr{'\n'};
};

} // namespace seissol::io::instance::point

#endif // SEISSOL_SRC_IO_INSTANCE_POINT_CSV_H_
