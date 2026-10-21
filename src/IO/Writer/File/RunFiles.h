// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_FILE_RUNFILES_H_
#define SEISSOL_SRC_IO_WRITER_FILE_RUNFILES_H_

#include "Common/Filesystem.h"

#include <ctime>
#include <string>
#include <unordered_set>
#include <utils/timeutils.h>

namespace seissol::io::writer::file {

/**
 * @brief The files the outputs of this run have written, remembered from one write to the next.
 *
 * A run owns its output files. The first time it writes one, the file is created anew; later
 * writes of the same run open it again to append. What an earlier run left under that name is
 * replaced by a fresh run. A run that resumes from a checkpoint keeps it instead, under a backup
 * name, since it holds the output up to the checkpoint that the resumed run does not write again.
 */
struct RunFiles {
  bool resumed{false};
  std::unordered_set<std::string> written;

  /**
   * @brief Records that this run writes @p name , and says whether that is the first time.
   */
  bool firstWrite(const std::string& name) { return written.insert(name).second; }
};

/**
 * @brief Moves an existing file @p name out of the way, to a name carrying the time of the move.
 *
 * Meant to be called on one rank only, before the file is created anew.
 */
inline void backUpFile(const std::string& name) {
  const auto path = seissol::filesystem::path(name);
  if (!seissol::directoryExists(seissol::filesystem::directory_entry(path))) {
    return;
  }
  // one stamp for all files of this process, so that the backups of one restart belong together
  static const std::string Stamp =
      utils::TimeUtils::timeAsString("%Y-%m-%d_%H-%M-%S", std::time(nullptr));
  auto backup = path;
  backup.replace_filename(path.stem().string() + ".bak_" + Stamp + path.extension().string());
  seissol::filesystem::rename(path, backup);
}

} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_RUNFILES_H_
