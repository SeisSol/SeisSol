#pragma once

#include <IO/Datatype/Datatype.hpp>
#include <IO/Instance/Point/TableWriter.hpp>
#include <IO/Writer/Instructions/Binary.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <memory>
#include <sstream>
#include <string>

namespace seissol::io::instance::point {

class Csv : public TableWriter {
  public:
  Csv() = default;

  std::string quote(const std::string& str) const { return quoteStr + str + quoteStr; }

  std::ostringstream& quote(std::ostringstream& stream, const std::string& str) const {
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

  std::string header() const {
    std::ostringstream stream;
    quote(stream, quantities[0].name);
    for (std::size_t i = 1; i < quantities.size(); ++i) {
      stream << delimiterStr;
      quote(stream, quantities[i].name);
    }
    stream << newlineStr;
    return stream.str();
  }

  std::string rows() const {
    std::ostringstream stream;
    const char* data = reinterpret_cast<const char*>(rowstorage.data());

    while (data < rowstorage.data() + rowstorage.size()) {
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

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() override {
    return [this](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
      const auto filename = prefix + "-" + name + ".csv";
      auto writer = writer::Writer();

      if (counter == 0) {
        writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
            filename, writer::WriteInline::createString(header())));
      }
      this->rowcache = rows();
      this->rowstorage.clear();
      writer.addInstruction(std::make_shared<writer::instructions::BinaryWrite>(
          filename, writer::WriteBuffer::create(rowcache.c_str(), rowcache.size())));
      return writer;
    };
  }

  private:
  std::string rowcache;
  char delimiterStr{';'};
  char quoteStr{'\"'};
  char newlineStr{'\n'};
};

} // namespace seissol::io::instance::point
