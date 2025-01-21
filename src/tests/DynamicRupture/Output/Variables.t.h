// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "Memory/MemoryAllocator.h"

namespace seissol::unit_test {

TEST_CASE("DR Variables") {
  SUBCASE("GeneralVariablesCount") {
    output::DrVarsT drVars;

    unsigned variableCounter = 0;
    auto countVariables = [&variableCounter](auto& var, int) { ++variableCounter; };

    misc::forEach(drVars, countVariables);
    REQUIRE(variableCounter == 12);
  }

  SUBCASE("TotalVariablesCount") {
    output::DrVarsT drVars;

    std::array<bool, std::tuple_size<output::DrVarsT>::value> mask{};
    for (size_t i = 0; i < std::tuple_size<output::DrVarsT>::value; ++i) {
      mask[i] = true;
    }

    auto assignMask = [&mask](auto& var, int index) { var.isActive = mask[index]; };

    misc::forEach(drVars, assignMask);

    unsigned variableCounter = 0;
    auto countVariables = [&variableCounter](auto& var, int) {
      if (var.isActive) {
        variableCounter += var.dim();
      }
    };

    misc::forEach(drVars, countVariables);
    REQUIRE(variableCounter == 20);
  }

  SUBCASE("AllocationDeallocationOfVariables") {
    output::DrVarsT drVars;

    std::array<bool, std::tuple_size<output::DrVarsT>::value> mask{};
    for (size_t i = 0; i < std::tuple_size<output::DrVarsT>::value; ++i) {
      mask[i] = true;
    }

    auto assignMask = [&mask](auto& var, int index) { var.isActive = mask[index]; };

    misc::forEach(drVars, assignMask);

    const seissol::memory::ManagedAllocator allocator;
    const unsigned numElements = 1024;
    auto allocateVariables = [numElements](auto& var, int) {
      var.maxCacheLevel = 3;
      var.allocateData(numElements);
    };
    misc::forEach(drVars, allocateVariables);

    real assignValue = 0.0;
    auto initVariables = [assignValue](auto& var, int) {
      if (var.isActive) {
        for (size_t dim = 0; dim < var.data.size(); ++dim) {
          for (size_t level = 0; level < var.maxCacheLevel; ++level) {
            for (size_t i = 0; i < var.size; ++i) {
              var(dim, level, i) = assignValue;
            }
          }
        }
      }
    };
    misc::forEach(drVars, initVariables);

    auto deallocateVariables = [](auto& var, int) { var.releaseData(); };
    misc::forEach(drVars, deallocateVariables);
  }
}
} // namespace seissol::unit_test
