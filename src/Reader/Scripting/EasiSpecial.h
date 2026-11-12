// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_READER_SCRIPTING_EASISPECIAL_H_
#define SEISSOL_SRC_READER_SCRIPTING_EASISPECIAL_H_

// The special components of easi -- !STRESS_STR_DIP_SLIP_AM, !OptimalStress and
// !AndersonianStress, the stress under which a fault of a given orientation and friction is
// optimally oriented -- as closed formulas of their inputs. A program computes them as easi's C++
// does (easi's src/component/Special.cpp, OptimalStress.cpp and AndersonianStress.cpp): the
// formulas are easi's, their names changed, so the operations and their order are easi's. What the
// constants fix the translation computes itself, with the same operations; a branch on an input of
// the point becomes a select.

#include "Expr/Ir.h"

#include <map>
#include <string>
#include <vector>

namespace seissol::reader::scripting {

/// A special component of easi: its tag, and the names of its inputs and of its outputs (the
/// members of its `in` and `out`).
struct EasiSpecial {
  std::string tag;
  std::vector<std::string> in;
  std::vector<std::string> out;
};

/// The special component of the tag `tag`, or null for another tag.
[[nodiscard]] const EasiSpecial* easiSpecialOf(const std::string& tag);

/// The outputs of `special` by name, at `in`: a node of `arena` for every input, by name.
[[nodiscard]] std::map<std::string, expr::NodeId> easiSpecialValues(
    expr::Arena& arena, const EasiSpecial& special, const std::map<std::string, expr::NodeId>& in);

} // namespace seissol::reader::scripting

#endif // SEISSOL_SRC_READER_SCRIPTING_EASISPECIAL_H_
