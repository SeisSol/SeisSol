// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Csv.h"

#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Writer.h>
#include <cstddef>
#include <functional>
#include <memory>
#include <sstream>
#include <string>

namespace seissol::io::instance::point {

Csv::Csv() = default;

std::string Csv::quote(const std::string& str) const { return quoteStr + str + quoteStr; }

std::ostringstream& Csv::quote(std::ostringstream& stream, const std::string& str) const {
  stream << quoteStr;
  for (std::size_t i = 0; i < str.size(); ++i) {
    stream << str[i];
    if (str[i] == quoteStr) {
      stream << quoteStr;
    }
  }
  stream << quoteStr;
  return stream;
}

std::string Csv::header() const {
  std::ostringstream stream;
  quote(stream, quantities[0].name);
  for (std::size_t i = 1; i < quantities.size(); ++i) {
    stream << delimiterStr;
    quote(stream, quantities[i].name);
  }
  stream << newlineStr;
  return stream.str();
}

std::string Csv::rows() const {
  std::ostringstream stream;
  const char* data = reinterpret_cast<const char*>(rowstorage.data());

  for (std::size_t j = 0; j < rowCount; ++j) {
    quote(stream, quantities[0].datatype->toStringRaw(data));
    data += quantities[0].datatype->size();
    for (std::size_t i = 1; i < quantities.size(); ++i) {
      stream << delimiterStr;
      quote(stream, quantities[i].datatype->toStringRaw(data));
      data += quantities[i].datatype->size();
    }
    stream << newlineStr;
  }
  return stream.str();
}

std::function<writer::Writer(const std::string&, std::size_t, double)> Csv::makeWriter() {
  return [this](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
    const auto filename = prefix + "-" + name + ".csv";
    auto writer = writer::Writer();

    if (counter == 0) {
      writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
          filename, writer::WriteInline::createString(header())));
    }
    this->rowcache = rows();
    this->resetStorage();

    writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
        filename, writer::WriteBuffer::create(rowcache.c_str(), rowcache.size())));
    return writer;
  };
}

} // namespace seissol::io::instance::point
