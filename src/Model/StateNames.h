// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_MODEL_STATENAMES_H_
#define SEISSOL_SRC_MODEL_STATENAMES_H_

#include "Model/Quantities.h"

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace seissol::model {

/// The names of the components of a group of quantities called `name`: the name itself for a
/// scalar, the name and 1, 2, 3 for a vector, the name and _xx, _yy, _zz, _xy, _yz, _xz for a
/// symmetric tensor (in Voigt order).
inline void appendComponentNames(std::vector<std::string>& names,
                                 const std::string& name,
                                 QuantityKind kind) {
  switch (kind) {
  case QuantityKind::Scalar:
    names.push_back(name);
    break;
  case QuantityKind::Vector:
    for (const char* component : {"1", "2", "3"}) {
      names.push_back(name + component);
    }
    break;
  case QuantityKind::SymTensor2:
    for (const char* component : {"_xx", "_yy", "_zz", "_xy", "_yz", "_xz"}) {
      names.push_back(name + component);
    }
    break;
  }
}

/// The names of all components of the state of the material `MaterialT`, as scripts read and give
/// them, in the order of an initial condition: the quantities, then the memory variables of each
/// relaxation mechanism -- theta1_xx, ..., theta1_xz, theta2_xx, ... for a stress, and theta1,
/// theta2, ... for a pressure.
template <typename MaterialT>
std::vector<std::string> stateNames() {
  std::vector<std::string> names(MaterialT::Quantities.begin(), MaterialT::Quantities.end());
  if constexpr (MaterialT::Mechanisms > 0) {
    for (std::size_t mechanism = 0; mechanism < MaterialT::Mechanisms; ++mechanism) {
      for (const auto& group : MaterialT::MechanismGroups) {
        appendComponentNames(
            names, std::string(group.name) + std::to_string(mechanism + 1), group.kind);
      }
    }
  }
  return names;
}

} // namespace seissol::model

#endif // SEISSOL_SRC_MODEL_STATENAMES_H_
