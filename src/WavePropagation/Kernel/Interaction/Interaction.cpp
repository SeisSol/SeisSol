#include "Interaction.hpp"

#include <Kernels/common.hpp>
#include <type_traits>
#include <variant>
#include "Common/configs.hpp"
#include "Common/cellconfigconv.hpp"

namespace {
// currently a kernel replacement

template <std::size_t Quantities,
          std::size_t ThisOrder,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          std::size_t Offset>
static inline void taylorSeriesInner(TargetRealT* target,
                                     SourceRealT* source,
                                     TargetRealT start,
                                     TargetRealT end,
                                     TargetRealT startCoeff,
                                     TargetRealT endCoeff) {
  constexpr std::size_t SourceStride =
      seissol::kernels::NumberOfAlignedBasisFunctions(SourceOrder - ThisOrder);
  constexpr std::size_t TargetStride = seissol::kernels::NumberOfAlignedBasisFunctions(TargetOrder);

  TargetRealT coeff = endCoeff - startCoeff;
  constexpr std::size_t BasisFunctionsSize = std::min(SourceStride, TargetStride);
#pragma omp for simd collapse(2)
  for (std::size_t j = 0; j < Quantities; ++j) {
    for (std::size_t k = 0; k < BasisFunctionsSize; ++k) {
      target[TargetStride * j + k] +=
          coeff * static_cast<TargetRealT>((source + Offset)[SourceStride * j + k]);
    }
  }
  if constexpr (ThisOrder <
                std::min(SourceOrder, TargetOrder)) { // TODO(David): are we sure about this? Or is
                                                      // ThisOrder < SourceOrder enough?
    TargetRealT newStartCoeff = startCoeff * start / static_cast<TargetRealT>(ThisOrder + 1);
    TargetRealT newEndCoeff = endCoeff * end / static_cast<TargetRealT>(ThisOrder + 1);
    taylorSeriesInner<Quantities,
                      ThisOrder + 1,
                      SourceRealT,
                      TargetRealT,
                      SourceOrder,
                      Offset + SourceStride * Quantities>(
        target, source, start, end, newStartCoeff, newEndCoeff);
  }
}

template <std::size_t SourceQuantities,
          std::size_t TargetQuantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder>
static void
    taylorSeries(TargetRealT* target, SourceRealT* source, TargetRealT start, TargetRealT end) {
  constexpr std::size_t Quantities = std::min(SourceQuantities, TargetQuantities);
  constexpr std::size_t TargetStride = seissol::kernels::NumberOfAlignedBasisFunctions(TargetOrder);

  std::memset(target, 0, TargetStride * TargetQuantities * sizeof(TargetRealT));
  taylorSeriesInner<Quantities, 0, SourceRealT, TargetRealT, SourceOrder, TargetOrder, 0>(
      target, source, start, end, start, end);
}
} // namespace

namespace seissol::waveprop::kernel {
template <typename SourceMaterialT, typename TargetMaterialT>
struct SameDofs {
  constexpr static Result = false;
};

template <typename MaterialT>
struct SameDofs<MaterialT, MaterialT> {
  constexpr static Result = true;
};

template <typename SourceMaterialT, typename TargetMaterialT>
struct SubsetDofs {
  constexpr static Result = SameDofs<SourceMaterialT, TargetMaterialT>;
};

bool canPassThrough(int source, int target) {
  if (source == target) {
    return true;
  } else {
    return std::visit(
        [&](const auto& sourceConfig) {
          return std::visit(
              [&](const auto& targetConfig) {
                using SourceConfig = std::decay_t<decltype(sourceConfig)>;
                using TargetConfig = std::decay_t<decltype(targetConfig)>;

                using SourceMaterialT = typename SourceConfig::MaterialT;
                using TargetMaterialT = typename TargetConfig::MaterialT;

                using SourceRealT = typename SourceConfig::RealT;
                using TargetRealT = typename TargetConfig::RealT;

                // (plasticity doesn't matter)
                // also, the SameDofs relation is symmetric, hence check
                return (SameDofs<SourceMaterialT, TargetMaterialT>::Result ||
                        SameDofs<TargetMaterialT, SourceMaterialT>::Result) &&
                       std::is_same_v<SourceRealT, TargetRealT> &&
                       SourceConfig::ConvergenceOrder == TargetConfig::ConvergenceOrder;
              },
              ConfigInstances[target]);
        },
        ConfigInstances[source]);
  }
}

template <typename TargetConfig>
typename TargetConfig::RealT* TimeIntegrator<TargetConfig>::timeIntegrateFace(
    int source, void* sourcePtr, void* tempmem, uint16_t ltsSetup, double start, double end) {
  if (((ltsSetup >> 8) & 1) != 0) {
    return sourcePtr;
  } else {
    std::visit(
        [&](auto sourceConfig) {
          using SourceConfig = std::decay_t<decltype(sourceConfig)>;

          using SourceMaterialT = typename SourceConfig::MaterialT;
          using TargetMaterialT = typename TargetConfig::MaterialT;

          using SourceRealT = typename SourceConfig::RealT;
          using TargetRealT = typename TargetConfig::RealT;

          static_assert(SubsetDofs<SourceMaterialT, TargetMaterialT>::Result,
                        "Only subset DOFS mode supported right now.");
          taylorSeries<SourceMaterialT::NumberOfQuantities,
                       TargetMaterialT::NumberOfQuantities,
                       SourceRealT,
                       TargetRealT,
                       SourceConfig::ConvergenceOrder,
                       TargetConfig::ConvergenceOrder>(tempmem, sourcePtr, start, end);
        },
        ConfigInstances[source]);
    return tempmem;
  }
}

template <typename TargetConfig>
void TimeIntegrator<TargetConfig>::timeIntegrateCell(typename TargetConfig::RealT* target[4],
                                                     void* sourcePtr[4],
                                                     void* tempmem[4],
                                                     const CellLocalInformation& info,
                                                     double start,
                                                     double end) {
  for (int face = 0; face < 4; ++face) {
    if (info.faceTypes[face] != FaceType::dynamicRupture &&
        info.faceTypes[face] != FaceType::outflow) {
      target[face] = timeIntegrateFace(
          info.neighborConfigIds[face], sourcePtr[face], tempmem[face], info.ltsSetup, start, end);
    }
  }
}

} // namespace seissol::waveprop::kernel

namespace seissol::_definitions {
const seissol::DeclareForAllConfigs<seissol::waveprop::kernel::TimeIntegrator> declTimeIntegrator;
} // namespace seissol::_definitions
