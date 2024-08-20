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
  Pytables();

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() override;

  private:
  std::vector<char> rowstorageCopy;
};

} // namespace seissol::io::instance::point
