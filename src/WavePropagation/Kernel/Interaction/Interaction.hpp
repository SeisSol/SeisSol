#pragma once

#include <Initializer/typedefs.hpp>
#include <Kernels/common.hpp>
#include <type_traits>
#include <variant>
#include "Common/configs.hpp"
#include "Common/cellconfigconv.hpp"

namespace seissol::waveprop::kernel {
bool canPassThrough(int source, int target);

template <typename TargetConfig>
struct TimeIntegrator {
  static typename TargetConfig::RealT* timeIntegrateFace(
      int source, void* sourcePtr, void* tempmem, uint16_t ltsSetup, double start, double end);
  static void timeIntegrateCell(typename TargetConfig::RealT* target[4],
                                void* sourcePtr[4],
                                void* tempmem[4],
                                const CellLocalInformation& info,
                                double start,
                                double end);
}
} // namespace seissol::waveprop::kernel
