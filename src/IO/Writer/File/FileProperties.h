// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_FILE_FILEPROPERTIES_H_
#define SEISSOL_SRC_IO_WRITER_FILE_FILEPROPERTIES_H_

#include <cstddef>
#include <mpi.h>
#include <string>
#include <utils/env.h>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::io::writer::file {

//! The metadata of one file is gathered into blocks of this size before it is written.
constexpr std::size_t MetaBlockSize = 1024 * 1024;

/**
 * @brief The value of an output tuning knob, from the environment.
 *
 * SEISSOL_IO_ is the name to use; XDMFWRITER_ is the prefix the previous writer used and is
 * accepted so that existing job scripts keep working.
 */
inline std::optional<std::string> outputEnv(const std::string& name) {
  auto value = utils::Env("SEISSOL_IO_").getOptional<std::string>(name);
  if (!value.has_value()) {
    value = utils::Env("XDMFWRITER_").getOptional<std::string>(name);
  }
  return value;
}

/**
 * @brief How far the bulk data of an output should be aligned in the file, or zero.
 *
 * Bulk writes that straddle a stripe boundary make more than one storage target take part in a
 * single write, which serialises them. What the right value is depends on the file system, so
 * there is nothing sensible to default to.
 */
inline std::size_t outputAlignment() {
  const auto value = outputEnv("ALIGNMENT");
  return value.has_value() ? utils::StringUtils::parse<std::size_t>(value.value()) : 0;
}

/**
 * @brief The MPI-IO hints for the output, from the environment.
 *
 * Which ones help depends entirely on the file system, so they are not something SeisSol can
 * pick; the names are the ones ROMIO understands.
 */
inline MPI_Info outputMpioHints() {
  static const std::vector<std::string> Hints = {"ind_rd_buffer_size",
                                                 "ind_wr_buffer_size",
                                                 "romio_ds_read",
                                                 "romio_ds_write",
                                                 "cb_buffer_size",
                                                 "cb_nodes",
                                                 "romio_cb_read",
                                                 "romio_cb_write",
                                                 "striping_factor",
                                                 "striping_unit"};

  static MPI_Info info = MPI_INFO_NULL;
  if (info == MPI_INFO_NULL) {
    MPI_Info_create(&info);
    for (const auto& hint : Hints) {
      auto name = hint;
      utils::StringUtils::toUpper(name);
      const auto value = outputEnv("MPIO_" + name);
      if (value.has_value()) {
        logInfo() << "Output: MPI-IO hint" << hint << "=" << value.value();
        MPI_Info_set(info, hint.c_str(), value.value().c_str());
      }
    }
  }
  return info;
}

} // namespace seissol::io::writer::file

#endif // SEISSOL_SRC_IO_WRITER_FILE_FILEPROPERTIES_H_
