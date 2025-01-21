// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_INTEGERMASKPARSER_H_
#define SEISSOL_SRC_COMMON_INTEGERMASKPARSER_H_

#include <optional>
#include <string>
#include <vector>

namespace seissol {
class IntegerMaskParser {
  public:
  using MaskType = std::vector<std::vector<int>>;
  static auto parse(const std::string& mask) -> std::optional<MaskType>;

  private:
  using OptionalIntVectorType = std::optional<std::vector<int>>;
  static auto parseIntRange(const std::string& str) -> OptionalIntVectorType;
  static auto parseIntList(const std::string& str) -> OptionalIntVectorType;
};
} // namespace seissol

#endif // SEISSOL_SRC_COMMON_INTEGERMASKPARSER_H_
