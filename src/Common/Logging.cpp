// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Logging.h"

#include <utility>
#include <utils/logger.h>

namespace seissol::logging {

int addCallback(utils::Logger::Callback callback, utils::LogLevel minimum) {
  return utils::Logger::addCallback(std::move(callback), minimum);
}

void removeCallback(int id) { utils::Logger::removeCallback(id); }

void enableConsole(bool enable) { utils::Logger::enableConsole(enable); }

void setLevel(utils::LogLevel minimum) { utils::Logger::setLevel(minimum); }

utils::Logger::AbortHandler setAbortHandler(utils::Logger::AbortHandler handler) {
  return utils::Logger::setAbortHandler(std::move(handler));
}

} // namespace seissol::logging
