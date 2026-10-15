// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_SCRIPTPARAMETERS_H_
#define SEISSOL_SRC_MODEL_SCRIPTPARAMETERS_H_

#include <array>
#include <cstddef>
#include <optional>
#include <string>
#include <string_view>
#include <tuple>
#include <type_traits>
#include <vector>

namespace seissol::model {

/// The parameters a material derives for scripts besides those a material file gives: none,
/// unless the material says otherwise. They keep `mu` and `lambda`, which scripts read of every
/// material before, readable where the material has none of its own.
template <typename MaterialT>
struct DerivedScriptParameters {
  static constexpr std::array<std::string_view, 0> Names{};
  static void values(const MaterialT& /*material*/, double* /*values*/) {}
};

/// How many parameters of the material `MaterialT` scripts read.
template <typename MaterialT>
inline constexpr std::size_t ScriptParameterCount =
    std::tuple_size_v<std::remove_cv_t<decltype(MaterialT::ParameterList)>> +
    DerivedScriptParameters<MaterialT>::Names.size();

/// The names of the parameters of the material `MaterialT` as scripts read them, in order: those
/// a material file gives (MaterialT::ParameterList), then those the material derives
/// (DerivedScriptParameters).
template <typename MaterialT>
std::vector<std::string> scriptParameterNames() {
  std::vector<std::string> names;
  names.reserve(ScriptParameterCount<MaterialT>);
  for (const auto& parameter : MaterialT::ParameterList) {
    names.push_back(parameter.first);
  }
  for (const auto& name : DerivedScriptParameters<MaterialT>::Names) {
    names.emplace_back(name);
  }
  return names;
}

/// The parameters of `material` as scripts read them, in the order of scriptParameterNames, into
/// `values`: what its members hold. A material whose members do not hold what the material file
/// gave (the moduli of one with relaxation, after the fit) corrects that afterwards.
template <typename MaterialT>
void fillScriptParameters(const MaterialT& material, double* values) {
  std::size_t i = 0;
  for (const auto& parameter : MaterialT::ParameterList) {
    values[i++] = material.*(parameter.second);
  }
  DerivedScriptParameters<MaterialT>::values(material, values + i);
}

/// The position of the parameter `name` among the parameters of the material `MaterialT` as
/// scripts read them, or nothing if it has none of that name.
template <typename MaterialT>
std::optional<std::size_t> scriptParameterIndex(std::string_view name) {
  const auto names = scriptParameterNames<MaterialT>();
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (names[i] == name) {
      return i;
    }
  }
  return std::nullopt;
}

/// The position of `name` among `names`, or nothing.
inline std::optional<std::size_t> parameterIndexOf(const std::vector<std::string>& names,
                                                   std::string_view name) {
  for (std::size_t i = 0; i < names.size(); ++i) {
    if (names[i] == name) {
      return i;
    }
  }
  return std::nullopt;
}

/// The names `names` as a list for a message: "rho, mu, lambda".
inline std::string parameterList(const std::vector<std::string>& names) {
  std::string list;
  for (std::size_t i = 0; i < names.size(); ++i) {
    list.append(i == 0 ? "" : ", ").append(names[i]);
  }
  return list;
}

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_SCRIPTPARAMETERS_H_
