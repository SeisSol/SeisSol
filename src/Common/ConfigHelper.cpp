// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigHelper.h"

#include <Config.h>
#include <Equations/Datastructures.h>
#include <array>
#include <sstream>
#include <string>
#include <variant>

namespace {
template <typename Cfg>
std::string configToString() {
  std::stringstream str;
  str << seissol::model::MaterialTT<Cfg>::Text << "-p" << Cfg::ConvergenceOrder << "-"
      << stringRealType(Cfg::Precision);
  if (Cfg::NumSimulations > 1) {
    str << "-s" << Cfg::NumSimulations;
  }
  return str.str();
}
} // namespace

namespace seissol {

const std::array<ConfigVariant, std::variant_size_v<ConfigVariant>> ConfigVariantList{
#define SEISSOL_CONFIGITER(cfg) cfg(),
#include "ConfigInclude.h"
};

const std::array<std::string, std::variant_size_v<ConfigVariant>> ConfigString{
#define SEISSOL_CONFIGITER(cfg) configToString<cfg>(),
#include "ConfigInclude.h"
};

} // namespace seissol
