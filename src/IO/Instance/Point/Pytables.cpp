// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Pytables.h"

#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <IO/Writer/Writer.h>
#include <cstddef>
#include <functional>
#include <memory>
#include <string>
namespace seissol::io::instance::point {

// reference: https://www.pytables.org/usersguide/file_format.html

Pytables::Pytables() = default;

std::function<writer::Writer(const std::string&, std::size_t, double)> Pytables::makeWriter() {
  return [this](const std::string& prefix, std::size_t counter, double time) -> writer::Writer {
    const auto filename = prefix + "-" + name + ".h5";
    auto writer = writer::Writer();

    this->rowstorageCopy = this->rowstorage;
    const auto rowCount = this->rowCount;
    this->resetStorage();

    auto rowDatatype = getRowDatatype();
    writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
        writer::instructions::Hdf5Location(filename, {}),
        "receiverdata",
        writer::WriteBuffer::create(rowstorageCopy.data(), rowCount, {}, rowDatatype),
        rowDatatype));

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
          writer::WriteInline::createString("2.6")));
      writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
          writer::instructions::Hdf5Location(filename, {}, "receiverdata"),
          "TITLE",
          writer::WriteInline::createString("receiverdata")));
      for (std::size_t i = 0; i < quantities.size(); ++i) {
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}, "receiverdata"),
            std::string("FIELD_") + std::to_string(i) + std::string("_NAME"),
            writer::WriteInline::createString(quantities[i].name)));
        // TODO: fix
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {}, "receiverdata"),
            std::string("FIELD_") + std::to_string(i) + std::string("_FILL"),
            writer::WriteInline::create(0ULL, quantities[i].datatype)));
      }
    }
    // TODO: fix
    /*writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
          writer::instructions::Hdf5Location(filename, {}, "receiverdata"),
          "NROWS",
          writer::WriteInline::create(0)));*/
    return writer;
  };
}

} // namespace seissol::io::instance::point
