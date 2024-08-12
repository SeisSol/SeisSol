#pragma once

#include <IO/Datatype/Datatype.hpp>
#include <IO/Instance/Point/TableWriter.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <IO/Writer/Instructions/Hdf5.hpp>
#include <memory>
#include <string>
namespace seissol::io::instance::point {

// reference: https://www.pytables.org/usersguide/file_format.html

class Pytables : public TableWriter {
  public:
  Pytables() = default;

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() override {
    return [this](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
      const auto filename = prefix + "-" + name + ".h5";
      auto writer = writer::Writer();

      this->rowstorageCopy = this->rowstorage;
      this->rowstorage.clear();

      writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {}),
          "receiverdata",
          writer::WriteBuffer::create(rowstorageCopy.data(), rowstorageCopy.size()),
          getRowDatatype()));

      // first write
      if (counter == 0) {
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}),
            "CLASS",
            writer::WriteInline::createString("GROUP")));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}),
            "PYTABLES_FORMAT_VERSION",
            writer::WriteInline::createString("2.0")));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}),
            "VERSION",
            writer::WriteInline::createString("1.0")));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}),
            "TITLE",
            writer::WriteInline::createString(name)));

        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}, "receiverdata"),
            "CLASS",
            writer::WriteInline::createString("TABLE")));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}, "receiverdata"),
            "VERSION",
            writer::WriteInline::createString("2.0")));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}, "receiverdata"),
            "TITLE",
            writer::WriteInline::createString("receiverdata")));
      }
      return writer;
    };
  }

  private:
  std::vector<char> rowstorageCopy;
};

} // namespace seissol::io::instance::point
