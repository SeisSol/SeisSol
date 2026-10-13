// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_RESULTWRITER_MINISEISSOLWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_MINISEISSOLWRITER_H_

#include <optional>
#include <string>
#include <vector>

namespace seissol::writer {
/**
 * The timings mini SeisSol measures to weight the ranks in the partitioning.
 *
 * They are measured while the mesh is read, which is before the output directory is created;
 * hence they are recorded first and written once the directory exists.
 */
class MiniSeisSolWriter {
  public:
  //! Keeps this rank's measurement until it can be written.
  void record(double elapsedTime, double weight);

  //! Writes the recorded measurements of all ranks next to @p outputPrefix, if there are any.
  //! Collective: every rank has to call it.
  void write(const std::string& outputPrefix) const;

  private:
  struct Measurement {
    double elapsedTime;
    double weight;
  };
  std::optional<Measurement> measurement_;
};
} // namespace seissol::writer

#endif // SEISSOL_SRC_RESULTWRITER_MINISEISSOLWRITER_H_
