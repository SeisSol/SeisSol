// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_LOGGING_H_
#define SEISSOL_SRC_COMMON_LOGGING_H_

#include <utils/logger.h>

/**
 * @brief The settings of the logger, for code outside of SeisSol, e.g. bindings to Python.
 *
 * The logger is header-only, so a library that compiles it and does not export its symbols keeps a
 * copy of its settings of its own: what such a library sets through utils::Logger does not reach
 * the messages of SeisSol. These functions are compiled into SeisSol, so that what they set does,
 * for the libraries compiled into SeisSol as well. Code of SeisSol itself may use utils::Logger
 * directly.
 */
namespace seissol::logging {

/**
 * @brief Registers @p callback, which gets every message of this process from @p minimum on.
 *
 * @return The ID to remove the callback with.
 */
int addCallback(utils::Logger::Callback callback, utils::LogLevel minimum);

//! @brief Removes the callback registered with @p id.
void removeCallback(int id);

//! @brief Sets whether messages go to the console, e.g. not when another framework shows them.
void enableConsole(bool enable);

//! @brief Sets the least severe level of the messages that are handed on.
void setLevel(utils::LogLevel minimum);

/**
 * @brief Sets what ends the program after a fatal error.
 *
 * @return The handler set before.
 */
utils::Logger::AbortHandler setAbortHandler(utils::Logger::AbortHandler handler);

} // namespace seissol::logging

#endif // SEISSOL_SRC_COMMON_LOGGING_H_
