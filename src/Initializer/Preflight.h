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

//! @brief How far `SeisSol --check` goes.
enum class CheckDepth {
  //! the parameters, the files they name, the head of the mesh, the output directory
  Files,
  //! also the whole mesh: its boundary conditions, groups and cells, and where the receivers, the
  //! on-fault receivers and the point sources lie in it
  Mesh
};

//! @brief The options of `SeisSol --check`.
struct CheckOptions {
  //! where the report goes; empty if nowhere
  std::string reportPath;
  CheckDepth depth{CheckDepth::Mesh};
};

/**
 * @brief Checks a run before it starts, as `SeisSol --check` does.
 *
 * Reads the parameter file with the @p overrides (each as section.field=value) as a run reads
 * them, and checks what a run finds out only later: that the files the parameters name exist,
 * that the YAML files among them parse, that the mesh opens, that the output can be written where
 * the parameters put it, and that the @p checkpoint , if one is given, exists. It reports the size
 * of the mesh as well. With CheckDepth::Mesh, it reads the whole mesh once (see
 * meshscan::scanMesh): it checks the boundary conditions of the faces against the face map and
 * what this build supports, the cells as the run checks them, and that the receivers, the on-fault
 * receivers and the point sources lie where the run can put them.
 *
 * Runs on one process, without MPI and devices, so that it can run where a job is prepared, e.g.
 * on a login node; it uses the threads of OpenMP.
 *
 * An error in the parameters stops the check as it would stop a run; what the check itself finds,
 * it collects: as problems, which keep the run from starting or going as intended, and as
 * warnings, which do not. If the report path is not empty, the outcome is written there as JSON,
 * also when an error stops the check.
 *
 * @return The exit code: 0 if the run can start, 1 otherwise.
 */
int preflight(const std::string& parameterFile,
              const std::vector<std::string>& overrides,
              const std::optional<std::string>& checkpoint,
              const CheckOptions& options);

/**
 * @brief Takes --check, or --check=<report>, and --check-depth=<files|mesh> out of the command
 * line; --check-depth alone asks for the check as well.
 *
 * @return Nothing if neither is there; otherwise the options of the check.
 */
std::optional<CheckOptions> extractCheckOptions(int& argc, char** argv);

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_PREFLIGHT_H_
