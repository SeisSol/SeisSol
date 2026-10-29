// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_PARAMETEROVERRIDES_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_PARAMETEROVERRIDES_H_

#include "Initializer/Parameters/ParameterReader.h"

#include <optional>
#include <string>
#include <string_view>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::initializer::parameters {

/**
 * @brief A value given on the command line in place of the one in the parameter file.
 *
 * Written as section.field=value, e.g. output.outputfile=runs/a/data or abortcriteria.endtime=2.
 * The value is read as if it stood in the parameter file, quotes and Fortran exponents included.
 */
struct ParameterOverride {
  //! The section and the field the value goes to, in lower case, as the parameter file is read.
  std::vector<std::string> path;
  //! The value, as it would be written in the parameter file.
  std::string value;

  //! @brief The field the value goes to, as section.field.
  [[nodiscard]] std::string name() const;
};

/**
 * @brief Reads @p text as an override, or nothing if it is not of the form section.field=value.
 */
std::optional<ParameterOverride> parseOverride(std::string_view text);

/**
 * @brief Puts the values of @p overrides into @p parameters, the parameter file as it was read.
 *
 * A later override of the same field wins over an earlier one; a section that the file does not
 * have is added.
 */
void applyOverrides(const YAML::Node& parameters, const std::vector<ParameterOverride>& overrides);

/**
 * @brief The overrides of fields that reading the parameters with @p reader did not look at.
 *
 * Those would have no effect: a misspelt name is the usual reason. Unlike a field the parameter
 * file has but SeisSol does not know, which only warrants a warning, one given on the command line
 * is an error, since whoever gave it expects it to change the run.
 */
std::vector<ParameterOverride> unknownOverrides(const ParameterReader& reader,
                                                const std::vector<ParameterOverride>& overrides);

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_PARAMETEROVERRIDES_H_
