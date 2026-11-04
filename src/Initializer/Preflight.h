// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PREFLIGHT_H_
#define SEISSOL_SRC_INITIALIZER_PREFLIGHT_H_

#include <optional>
#include <string>
#include <vector>

namespace seissol::initializer {

/**
 * @brief Checks a run before it starts, as `SeisSol --check` does.
 *
 * Reads the parameter file with the @p overrides (each as section.field=value) as a run reads
 * them, and checks what a run finds out only later: that the files the parameters name exist,
 * that the YAML files among them parse, that the mesh opens, that the output can be written where
 * the parameters put it, and that the @p checkpoint , if one is given, exists. It reports the size
 * of the mesh as well. Runs on one process, without MPI and devices, so that it can run where a
 * job is prepared, e.g. on a login node.
 *
 * An error in the parameters stops the check as it would stop a run; what the check itself finds,
 * it collects. If @p reportPath is not empty, the outcome is written there as JSON, also when an
 * error stops the check.
 *
 * @return The exit code: 0 if the run can start, 1 otherwise.
 */
int preflight(const std::string& parameterFile,
              const std::vector<std::string>& overrides,
              const std::optional<std::string>& checkpoint,
              const std::string& reportPath);

/**
 * @brief Takes --check, or --check=<report>, out of the command line.
 *
 * @return Nothing if it is not there; otherwise where the report goes, empty if nowhere.
 */
std::optional<std::string> extractCheckOption(int& argc, char** argv);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_PREFLIGHT_H_
