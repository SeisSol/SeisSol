// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_FILE_FILE_H_
#define SEISSOL_SRC_IO_WRITER_FILE_FILE_H_

#include <string>
namespace seissol::io::writer::file {
class File {
  public:
  std::string filename() { return name; }

  private:
  std::string name;
};
} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_FILE_H_
