// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_POINT_CSV_H_
#define SEISSOL_SRC_IO_INSTANCE_POINT_CSV_H_

#include "IO/Instance/Point/TableWriter.h"

#include <cstddef>
#include <functional>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::io::instance::point {

//! @brief Which values are written between quotes.
enum class CsvQuoting : std::uint8_t {
  //! Every one of them.
  All,
  //! Only the ones that are text, which is where a delimiter could turn up.
  Text,
  //! None, which leaves it to the writer of the values that none of them needs it.
  None
};

/**
 * @brief The punctuation of a CSV file.
 *
 * The same settings read a file back, so a reader takes this rather than guessing.
 */
struct CsvFormat {
  char delimiter{','};
  char quote{'"'};
  char newline{'\n'};
  CsvQuoting quoting{CsvQuoting::Text};
};

class Csv : public TableWriter {
  public:
  ~Csv() override = default;
  explicit Csv(std::string name, CsvFormat format = {});

  [[nodiscard]] std::string header() const;

  [[nodiscard]] std::string rows() const;

  std::function<writer::Writer(const std::string&, std::size_t, double)> makeWriter() override;

  /**
   * @brief Writes the header and the rows collected so far to @p path, from the calling rank.
   *
   * For the tables a run leaves next to its output, which one rank assembles in full and which
   * nothing in the simulation waits for. What a simulation produces while it runs goes through
   * makeWriter instead, so that the ranks write their own rows and no one gathers.
   */
  void writeFile(const std::string& path);

  private:
  void field(std::ostringstream& stream, const std::string& value, bool text) const;

  std::string rowcache_;
  CsvFormat format_;
};

//! @brief A CSV file as it was read: the column names, and the rows as the text they hold.
struct CsvTable {
  std::vector<std::string> header;
  std::vector<std::vector<std::string>> rows;

  //! @brief The position of the column named @p name; fails if there is none.
  [[nodiscard]] std::size_t column(const std::string& name) const;
};

/**
 * @brief Reads a CSV file written with @p format.
 *
 * The values come back as the text they are written as rather than converted, since what a column
 * holds is not in the file. Enough to compare two files or to pick a column out of one.
 */
CsvTable readCsv(const std::string& path, const CsvFormat& format = {});

//! @brief As readCsv, on text that is already in hand.
CsvTable parseCsv(const std::string& content, const CsvFormat& format = {});

} // namespace seissol::io::instance::point

#endif // SEISSOL_SRC_IO_INSTANCE_POINT_CSV_H_
