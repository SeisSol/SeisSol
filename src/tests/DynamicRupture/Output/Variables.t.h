#ifndef SEISSOL_DRVARIABLES_T_H
#define SEISSOL_DRVARIABLES_T_H

#include <cxxtest/TestSuite.h>
#include <DynamicRupture/Output/DataTypes.hpp>
#include <DynamicRupture/Math.h>
#include <Initializer/MemoryAllocator.h>

namespace seissol::unit_test::dr {
class Variables;
}

using namespace seissol;
using namespace seissol::dr;

class seissol::unit_test::dr::Variables : public CxxTest::TestSuite {
  public:
  void testGeneralVariablesCount() {
    output::DrVarsT drVars;

    unsigned variableCounter = 0;
    auto countVariables = [&variableCounter](auto& var, int) { ++variableCounter; };

    aux::forEach(drVars, countVariables);
    TS_ASSERT_EQUALS(variableCounter, 12);
  }

  void testTotalVariablesCount() {
    output::DrVarsT drVars;

    std::array<bool, std::tuple_size<output::DrVarsT>::value> mask{};
    for (size_t i = 0; i < std::tuple_size<output::DrVarsT>::value; ++i)
      mask[i] = true;

    auto assignMask = [&mask](auto& var, int index) { var.isActive = mask[index]; };

    aux::forEach(drVars, assignMask);

    unsigned variableCounter = 0;
    auto countVariables = [&variableCounter](auto& var, int) {
      if (var.isActive) {
        variableCounter += var.dim();
      }
    };

    aux::forEach(drVars, countVariables);
    TS_ASSERT_EQUALS(variableCounter, 20);
  }

  void testAllocationDeallocationOfVariables() {
    output::DrVarsT drVars;

    std::array<bool, std::tuple_size<output::DrVarsT>::value> mask{};
    for (size_t i = 0; i < std::tuple_size<output::DrVarsT>::value; ++i)
      mask[i] = true;

    auto assignMask = [&mask](auto& var, int index) { var.isActive = mask[index]; };

    aux::forEach(drVars, assignMask);

    seissol::memory::ManagedAllocator allocator;
    const unsigned numElements = 1024;
    auto allocateVariables = [numElements, &allocator](auto& var, int) {
      var.maxCacheLevel = 3;
      var.allocateData(numElements);
    };
    aux::forEach(drVars, allocateVariables);

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
    aux::forEach(drVars, initVariables);

    auto deallocateVariables = [](auto& var, int) { var.releaseData(); };
    aux::forEach(drVars, deallocateVariables);
  }
};

#endif // SEISSOL_DRVARIABLES_T_H
