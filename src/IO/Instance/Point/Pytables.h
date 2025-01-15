// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_POINT_PYTABLES_H_
#define SEISSOL_SRC_IO_INSTANCE_POINT_PYTABLES_H_

#include <IO/Datatype/Datatype.h>
#include <IO/Instance/Point/TableWriter.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <memory>
#include <string>
namespace seissol::io::instance::point {

// reference: https://www.pytables.org/usersguide/file_format.html

class Pytables : public TableWriter {
  public:
  ~Pytables() override = default;
  Pytables();

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() override;

  private:
  std::vector<char> rowstorageCopy;
};

} // namespace seissol::io::instance::point

#endif // SEISSOL_SRC_IO_INSTANCE_POINT_PYTABLES_H_
