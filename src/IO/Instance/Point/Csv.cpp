// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Csv.h"

#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Writer.h"

#include <cstddef>
#include <functional>
#include <memory>
#include <sstream>
#include <string>

namespace seissol::io::instance::point {

Csv::Csv() = default;

std::string Csv::quote(const std::string& str) const { return quoteStr_ + str + quoteStr_; }

std::ostringstream& Csv::quote(std::ostringstream& stream, const std::string& str) const {
  stream << quoteStr_;
  for (std::size_t i = 0; i < str.size(); ++i) {
    stream << str[i];
    if (str[i] == quoteStr_) {
      stream << quoteStr_;
    }
  }
  stream << quoteStr_;
  return stream;
}

std::string Csv::header() const {
  std::ostringstream stream;
  quote(stream, quantities_[0].name);
  for (std::size_t i = 1; i < quantities_.size(); ++i) {
    stream << delimiterStr_;
    quote(stream, quantities_[i].name);
  }
  stream << newlineStr_;
  return stream.str();
}

std::string Csv::rows() const {
  std::ostringstream stream;
  const char* data = reinterpret_cast<const char*>(rowstorage_.data());

  for (std::size_t j = 0; j < rowCount_; ++j) {
    quote(stream, quantities_[0].datatype->toStringRaw(data));
    data += quantities_[0].datatype->size();
    for (std::size_t i = 1; i < quantities_.size(); ++i) {
      stream << delimiterStr_;
      quote(stream, quantities_[i].datatype->toStringRaw(data));
      data += quantities_[i].datatype->size();
    }
    stream << newlineStr_;
  }
  return stream.str();
}

std::function<writer::Writer(const std::string&, std::size_t, double)> Csv::makeWriter() {
  return [this](const std::string& prefix, std::size_t counter, double /*time*/) -> writer::Writer {
    const auto filename = prefix + "-" + name_ + ".csv";
    auto writer = writer::Writer();

    if (counter == 0) {
      writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
          filename, writer::WriteInline::createString(header())));
    }
    this->rowcache_ = rows();
    this->resetStorage();

    writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
        filename, writer::WriteBuffer::create(rowcache_.c_str(), rowcache_.size())));
    return writer;
  };
}

} // namespace seissol::io::instance::point
