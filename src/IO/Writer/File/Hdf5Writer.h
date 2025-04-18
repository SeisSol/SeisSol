// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_FILE_HDF5WRITER_H_
#define SEISSOL_SRC_IO_WRITER_FILE_HDF5WRITER_H_
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <unordered_map>
#include <vector>

namespace seissol::io::writer::file {

class Hdf5File {
  public:
  Hdf5File(MPI_Comm comm);
  void openFile(const std::string& name);
  void openGroup(const std::string& name);
  void openDataset(const std::string& name);
  void writeAttribute(const async::ExecInfo& info,
                      const std::string& name,
                      const std::shared_ptr<DataSource>& source);
  void writeData(const async::ExecInfo& info,
                 const std::string& name,
                 const std::shared_ptr<DataSource>& source,
                 const std::shared_ptr<datatype::Datatype>& targetType,
                 int compress);
  void closeDataset();
  void closeGroup();
  void closeFile();

  private:
  hid_t file;
  std::stack<hid_t> handles; // TODO: have something more sophisticated than a single stack
  MPI_Comm comm;
};

class Hdf5Writer {
  public:
  Hdf5Writer(MPI_Comm comm);

  void writeAttribute(const async::ExecInfo& info, const instructions::Hdf5AttributeWrite& write);

  void writeData(const async::ExecInfo& info, const instructions::Hdf5DataWrite& write);

  void finalize();

  private:
  std::unordered_map<std::string, Hdf5File> openFiles;
  MPI_Comm comm;
};
} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_HDF5WRITER_H_
