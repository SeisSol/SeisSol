// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_FILE_HDF5WRITER_H_
#define SEISSOL_SRC_IO_WRITER_FILE_HDF5WRITER_H_
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"

#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace seissol::io::writer::file {

/**
 * @brief The files the outputs of this run have written, remembered from one write to the next.
 *
 * A run owns its output files. The first time it writes one, the file is created anew, whatever an
 * earlier run left under that name; later writes of the same run open it again to append. Only a
 * run that resumes from a checkpoint opens a file it has not written itself, since it continues
 * what the run before it wrote.
 */
struct RunFiles {
  bool resumed{false};
  std::unordered_set<std::string> written;
};

class Hdf5File {
  public:
  explicit Hdf5File(MPI_Comm comm);
  /**
   * @brief Opens @p name, or creates it if it does not exist yet.
   *
   * With @p fresh , the file is created in any case, replacing what is there.
   */
  void openFile(const std::string& name, bool fresh = false);
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
  void writeLinkExternal(const std::string& name,
                         const std::string& targetFile,
                         const std::string& targetPath);
  void closeDataset();
  void closeGroup();
  void closeFile();

  private:
  hid_t file_{-1};
  std::stack<hid_t> handles_; // TODO: have something more sophisticated than a single stack
  MPI_Comm comm_{MPI_COMM_NULL};
};

class Hdf5Writer {
  public:
  /**
   * @brief A writer for one write plan.
   *
   * @p runFiles is what the run has written so far, and is updated by the writes of this plan.
   * Without it, every file that exists is opened and appended to.
   */
  explicit Hdf5Writer(MPI_Comm comm, RunFiles* runFiles = nullptr);

  void writeAttribute(const async::ExecInfo& info, const instructions::Hdf5AttributeWrite& write);

  void writeData(const async::ExecInfo& info, const instructions::Hdf5DataWrite& write);

  void writeLinkExternal(const async::ExecInfo& info,
                         const instructions::Hdf5LinkExternalWrite& write);

  void finalize();

  private:
  //! The file @p name , opened on first use in this plan.
  Hdf5File file(const std::string& name);

  std::unordered_map<std::string, Hdf5File> openFiles_;
  MPI_Comm comm_{MPI_COMM_NULL};
  RunFiles* runFiles_{nullptr};
};
} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_HDF5WRITER_H_
