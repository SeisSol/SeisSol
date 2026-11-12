// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Model/ScriptTypes.h"

#include "Equations/Datastructures.h"
#include "Expr/SderivFrontend.h"
#include "Model/CommonDatastructures.h"

#include <mutex>
#include <string>
#include <vector>

namespace seissol::model {

namespace {

template <typename MaterialT>
void registerMaterial() {
  std::vector<std::string> fields;
  fields.reserve(MaterialT::ParameterList.size());
  for (const auto& parameter : MaterialT::ParameterList) {
    fields.push_back(parameter.first);
  }
  expr::registerRecordType(scriptTypeName(MaterialT::Text), fields);
}

} // namespace

std::string scriptTypeName(const std::string& text) { return text.substr(0, text.find('-')); }

void registerScriptTypes() {
  static std::once_flag once;
  std::call_once(once, [] {
    registerMaterial<ElasticMaterial>();
    registerMaterial<AcousticMaterial>();
    registerMaterial<AnisotropicMaterial>();
    registerMaterial<PoroElasticMaterial>();
    // the parameters a material file gives do not depend on the number of mechanisms
    registerMaterial<ViscoElasticMaterial<1>>();
    registerMaterial<ViscoAcousticMaterial<1>>();
    registerMaterial<Plasticity>();
  });
}

} // namespace seissol::model
