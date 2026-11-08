// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_ENVIRONMENT_H_
#define SEISSOL_SRC_COMMON_ENVIRONMENT_H_

namespace seissol {

/**
 * @brief Sets the environment variable @p name to @p value, whatever it held before.
 *
 * On Windows, an empty @p value removes the variable instead, as the C runtime there does.
 *
 * The environment belongs to the whole process: no other thread may read it meanwhile.
 */
void setEnvironmentVariable(const char* name, const char* value) noexcept;

/**
 * @brief Removes the environment variable @p name, if there is one.
 *
 * The environment belongs to the whole process: no other thread may read it meanwhile.
 */
void unsetEnvironmentVariable(const char* name) noexcept;

} // namespace seissol

#endif // SEISSOL_SRC_COMMON_ENVIRONMENT_H_
