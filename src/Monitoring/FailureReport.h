// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MONITORING_FAILUREREPORT_H_
#define SEISSOL_SRC_MONITORING_FAILUREREPORT_H_

#include "Common/Json.h"

#include <string>
#include <string_view>
#include <utils/logger.h>

/**
 * @brief The file that a failure rewrites to say that it failed, e.g. the run summary.
 *
 * Its text is set ahead, in two parts, between which the failure goes as a JSON object: so that a
 * signal handler can write it, which may neither allocate nor take a lock. Only the first failure
 * of the process writes it.
 */
namespace seissol::monitoring {

/**
 * @brief Sets the failure report: @p before , the failure and @p after replace @p path .
 *
 * Set it again whenever the text changes; the one set before is replaced.
 */
void setFailureReport(const std::string& path, std::string before, std::string after);

/**
 * @brief Stands for the failure in a report given as JSON, as the value of a member; the control
 * characters keep it apart from any other text.
 */
inline constexpr std::string_view FailurePlaceholder = "\1failure\1";

/**
 * @brief Sets the failure report to @p report , as JSON text with a newline, with the failure in
 * place of the string FailurePlaceholder.
 */
void setFailureReport(const std::string& path, const json::Value& report);

//! @brief Sets no failure report, e.g. once the run has completed.
void clearFailureReport();

/**
 * @brief Writes the failure report with @p failure , a JSON object, unless a failure has written it
 * already or none is set. Async-signal-safe.
 */
void writeFailureReport(std::string_view failure);

//! @brief A fatal error as a JSON object: its message, its source in the code, and its time.
std::string failureOf(const utils::LogRecord& record);

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_FAILUREREPORT_H_
