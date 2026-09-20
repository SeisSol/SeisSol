// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Csv.h"

#include "IO/Datatype/Datatype.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Writer.h"

#include <algorithm>
#include <cstddef>
#include <fstream>
#include <functional>
#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::instance::point {

namespace {
//! @brief Whether a column holds text, i.e. something a delimiter could turn up in.
bool isText(const datatype::Datatype& datatype) {
  return dynamic_cast<const datatype::StringDatatype*>(&datatype) != nullptr;
}
} // namespace

Csv::Csv(std::string name, CsvFormat format) : TableWriter(std::move(name)), format_(format) {}

std::ostringstream&
    Csv::field(std::ostringstream& stream, const std::string& value, bool text) const {
  const bool quoted =
      format_.quoting == CsvQuoting::All || (format_.quoting == CsvQuoting::Text && text);
  if (!quoted) {
    stream << value;
    return stream;
  }
  stream << format_.quote;
  for (const char character : value) {
    stream << character;
    if (character == format_.quote) {
      stream << format_.quote;
    }
  }
  stream << format_.quote;
  return stream;
}

std::string Csv::header() const {
  std::ostringstream stream;
  for (std::size_t i = 0; i < quantities_.size(); ++i) {
    if (i > 0) {
      stream << format_.delimiter;
    }
    // a column name is text whatever the column holds
    field(stream, quantities_[i].name, true);
  }
  stream << format_.newline;
  return stream.str();
}

std::string Csv::rows() const {
  std::ostringstream stream;
  const char* data = rowstorage_.data();

  for (std::size_t row = 0; row < rowCount_; ++row) {
    for (std::size_t i = 0; i < quantities_.size(); ++i) {
      if (i > 0) {
        stream << format_.delimiter;
      }
      field(stream, quantities_[i].datatype->toStringRaw(data), isText(*quantities_[i].datatype));
      data += quantities_[i].datatype->size();
    }
    stream << format_.newline;
  }
  return stream.str();
}

std::function<writer::Writer(const std::string&, std::size_t, double)> Csv::makeWriter() {
  return [this](const std::string& prefix, std::size_t counter, double /*time*/) -> writer::Writer {
    const auto filename = prefix + "-" + name() + ".csv";
    auto writer = writer::Writer();

    if (counter == 0) {
      // the same on every rank, so it goes in once
      writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
          filename, writer::WriteInline::createString(header()), 0, true));
    }
    this->rowcache_ = rows();
    this->resetStorage();

    // the rows of the ranks land behind one another, in rank order
    writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
        filename, writer::WriteBuffer::create(rowcache_.c_str(), rowcache_.size()), 0, true));
    return writer;
  };
}

std::size_t CsvTable::column(const std::string& name) const {
  const auto position = std::find(header.begin(), header.end(), name);
  if (position == header.end()) {
    logError() << "The column" << name << "is not in the table.";
  }
  return static_cast<std::size_t>(position - header.begin());
}

CsvTable parseCsv(const std::string& content, const CsvFormat& format) {
  CsvTable table;
  std::vector<std::string> row;
  std::string value;
  bool quoted = false;
  bool started = false;

  const auto endField = [&]() {
    row.push_back(value);
    value.clear();
    started = false;
  };
  const auto endRow = [&]() {
    endField();
    if (table.header.empty()) {
      table.header = std::move(row);
    } else {
      table.rows.push_back(std::move(row));
    }
    row = {};
  };

  for (std::size_t position = 0; position < content.size(); ++position) {
    const char character = content[position];
    if (quoted) {
      if (character != format.quote) {
        value.push_back(character);
      } else if (position + 1 < content.size() && content[position + 1] == format.quote) {
        // a quote inside a quoted value is written twice
        value.push_back(format.quote);
        ++position;
      } else {
        quoted = false;
      }
    } else if (character == format.quote && !started) {
      quoted = true;
      started = true;
    } else if (character == format.delimiter) {
      endField();
    } else if (character == format.newline) {
      endRow();
    } else {
      value.push_back(character);
      started = true;
    }
  }
  // a file that does not end in a newline still has that last row
  if (!row.empty() || started || !value.empty()) {
    endRow();
  }
  return table;
}

CsvTable readCsv(const std::string& path, const CsvFormat& format) {
  std::ifstream stream(path);
  if (!stream.good()) {
    logError() << "Could not read the table" << path << ".";
  }
  std::ostringstream buffer;
  buffer << stream.rdbuf();
  return parseCsv(buffer.str(), format);
}

} // namespace seissol::io::instance::point
