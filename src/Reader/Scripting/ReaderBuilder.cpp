// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#include "ReaderBuilder.h"

#include "Reader/Scripting/EasiReader.h"
#include "Reader/Scripting/LuaReader.h"

#include <fstream>
#include <memory>
#include <sstream>
#include <utils/logger.h>
#include <utils/stringutils.h>

namespace seissol::reader::scripting {

std::unique_ptr<DataReader> buildReader(const std::string& path,
                                        const std::vector<std::string>& defaultInArgs) {
  const auto parts = utils::StringUtils::split(path, ':');
  if (parts.size() == 1) {
    return std::make_unique<EasiReader>(path, defaultInArgs);
  } else {
    if (parts[0] == "easi") {
      return std::make_unique<EasiReader>(path, defaultInArgs);
    } else if (parts[0] == "lua") {
      std::ifstream file(path);
      std::stringstream code;
      code << file.rdbuf();
      return std::make_unique<LuaReader>(code.str());
    } else {
      logError() << "";
    }
  }

  return nullptr;
}

} // namespace seissol::reader::scripting
