// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ConfigHelper.h"

#include "Common/Real.h"
#include "Common/Typedefs.h"
#include "Config.h"
#include "Equations/Datastructures.h"

#include <sstream>
#include <string>

namespace seissol {

namespace {
template <typename Cfg>
std::string configToString() {
  std::stringstream str;
  str << seissol::model::MaterialT::Text << "-o" << Cfg::ConvergenceOrder << "-"
      << stringRealType(Cfg::Precision);
  if (Cfg::NumSimulations > 1) {
    str << "-s" << Cfg::NumSimulations;
  }
  return str.str();
}

template <typename Cfg>
std::string configToDescriptor() {
  std::stringstream str;

  str << "Material: \"" << seissol::model::MaterialT::Text << "\"\n";
  str << "Convergence order: " << Cfg::ConvergenceOrder << "\n";
  str << "Precision: " << stringRealType(Cfg::Precision) << "\n";
  str << "Simulation count: " << Cfg::NumSimulations << "\n";
  str << "Dynamic rupture quadrature rule: \"";
  switch (Cfg::DRQuadRule) {
  case DRQuadRuleType::Dunavant:
    str << "Dunavant";
    break;
  case DRQuadRuleType::Stroud:
    str << "Stroud";
    break;
  case DRQuadRuleType::WitherdenVincent:
    str << "Witherden-Vincent";
    break;
  default:
    str << "Unknown (consult the developers for more info)";
    break;
  }
  str << "\"\n";

  return str.str();
}

} // namespace

const std::string ConfigString = configToString<Config>();
const std::string ConfigDescriptor = configToDescriptor<Config>();

} // namespace seissol
