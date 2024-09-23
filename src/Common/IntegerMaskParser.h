// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#ifndef SEISSOL_SRC_COMMON_INTEGERMASKPARSER_H_
#define SEISSOL_SRC_COMMON_INTEGERMASKPARSER_H_

#include <string>
#include <optional>
#include <vector>


namespace seissol {
class IntegerMaskParser {
  public:
  using MaskType = std::vector<std::vector<int>>;
  static std::optional<MaskType> parse(std::string mask);

  private:
  using OptionalIntVectorType = std::optional<std::vector<int>>;
  static OptionalIntVectorType parseIntRange(const std::string& str);
  static OptionalIntVectorType parseIntList(const std::string& str);
};
} // namespace seissol

#endif // SEISSOL_SRC_COMMON_INTEGERMASKPARSER_H_

