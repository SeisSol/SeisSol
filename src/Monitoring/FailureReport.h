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
 * @brief The file that a failure or a termination rewrites, e.g. the run summary.
 *
 * Its text is set ahead, in two parts, between which the detail of the failure or termination goes
 * as a JSON object: so that a signal handler can write it, which may neither allocate nor take a
 * lock. Only the first of either writes it; both kinds may be set at the same time, for whichever
 * comes.
 */
namespace seissol::monitoring {

/**
 * @brief Stands for the detail in a report given as JSON, as the value of a member; the control
 * characters keep it apart from any other text.
 */
inline constexpr std::string_view FailurePlaceholder = "\1failure\1";

/**
 * @brief Sets the report a failure writes to @p report , as JSON text with a newline, with the
 * failure in place of the string FailurePlaceholder. Set it again whenever the text changes.
 */
void setFailureReport(const std::string& path, const json::Value& report);

/**
 * @brief Sets the report a termination (a request to stop) writes, like setFailureReport(). A
 * termination and a failure may each have a report; the first to come writes it.
 */
void setTerminationReport(const std::string& path, const json::Value& report);

//! @brief Sets no failure and no termination report, e.g. once the run has completed.
void clearFailureReports();

/**
 * @brief Writes the failure report with @p failure , a JSON object, unless a failure or a
 * termination has written it already, or no failure report is set. Async-signal-safe.
 */
void writeFailureReport(std::string_view failure);

/**
 * @brief Writes the termination report with @p termination , a JSON object, unless a failure or a
 * termination has written it already, or no termination report is set. Async-signal-safe.
 */
void writeTerminationReport(std::string_view termination);

//! @brief A fatal error as a JSON object: its message, its source in the code, and its time.
std::string failureOf(const utils::LogRecord& record);

} // namespace seissol::monitoring

#endif // SEISSOL_SRC_MONITORING_FAILUREREPORT_H_
