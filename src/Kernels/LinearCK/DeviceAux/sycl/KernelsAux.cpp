// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/Datastructures.h"
#include "Kernels/Common.h"
#include "Kernels/Precision.h"
#include <cstdio>
#include <init.h>
#include <sycl/sycl.hpp>
#include <tensor.h>
#include <yateto.h>

#include <Solver/MultipleSimulations.h>

#ifdef DEVICE_EXPERIMENTAL_EXPLICIT_KERNELS
namespace {
constexpr std::size_t Blocksize = 128;

template <typename SourceRealT>
constexpr std::size_t RestFunctions(std::size_t SourceOrder, std::size_t ThisOrder) {
  std::size_t total = 0;
  for (std::size_t j = ThisOrder; j < SourceOrder; ++j) {
    total +=
        seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder - ThisOrder);
  }
  return total;
}

template <bool Integral,
          std::size_t Quantities,
          std::size_t ThisOrder,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder,
          std::size_t Offset,
          std::size_t SharedOffset>
static void taylorSumInner(sycl::nd_item<1>& item,
                           TargetRealT* const __restrict__ target,
                           const SourceRealT* const __restrict__ source,
                           TargetRealT start,
                           TargetRealT end,
                           TargetRealT startCoeff,
                           TargetRealT endCoeff,
                           SourceRealT* const __restrict__ shmem,
                           TargetRealT reg[Quantities]) {
  constexpr std::size_t MemorySize =
      seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder) * Quantities;
  constexpr std::size_t SourceStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder - ThisOrder);
  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);
  constexpr std::size_t RestMemSize =
      SharedOffset > 0 ? 0 : RestFunctions<SourceRealT>(SourceOrder, ThisOrder) * Quantities;
  constexpr std::size_t SourceMemSize = SourceStride * Quantities;
  constexpr bool UseShared = MemorySize >= RestMemSize;
  constexpr std::size_t LoadSize = UseShared ? RestMemSize : SourceMemSize;
  constexpr TargetRealT DivisionCoefficient =
      Integral ? static_cast<TargetRealT>(ThisOrder + 2) : static_cast<TargetRealT>(ThisOrder + 1);

  static_assert(seissol::tensor::dQ::size(ThisOrder) == SourceMemSize,
                "Tensor size mismatch in explicit kernel.");

  if constexpr (LoadSize > 0) {
    item.barrier();
    constexpr std::size_t Rounds = LoadSize / Blocksize;
    constexpr std::size_t Rest = LoadSize % Blocksize;
    if constexpr (Rounds > 0) {
#pragma unroll
      for (std::size_t j = 0; j < Rounds; ++j) {
        shmem[j * Blocksize + item.get_local_id(0)] =
            source[Offset + j * Blocksize + item.get_local_id(0)];
      }
    }
    if constexpr (Rest > 0) {
      if (item.get_local_id(0) < Rest) {
        shmem[Rounds * Blocksize + item.get_local_id(0)] =
            source[Offset + Rounds * Blocksize + item.get_local_id(0)];
      }
    }
    item.barrier();
  }

  const TargetRealT coeff = endCoeff - startCoeff;
  constexpr std::size_t BasisFunctionsSize = std::min(SourceStride, TargetStride);
  if (item.get_local_id(0) < BasisFunctionsSize) {
#pragma unroll
    for (std::size_t j = 0; j < Quantities; ++j) {
      // FIXME: non-optimal warp utilization (as before... But now it's in registers)
      reg[j] += coeff * static_cast<TargetRealT>(
                            shmem[SharedOffset + SourceStride * j + item.get_local_id(0)]);
    }
  }
  // TODO(David): are we sure about this? Or is ThisOrder + 1 < SourceOrder enough?
  if constexpr (ThisOrder + 1 < std::min(SourceOrder, TargetOrder)) {
    constexpr std::size_t SharedPosition = UseShared ? SharedOffset + SourceMemSize : 0;
    const TargetRealT newStartCoeff = startCoeff * start / DivisionCoefficient;
    const TargetRealT newEndCoeff = endCoeff * end / DivisionCoefficient;
    taylorSumInner<Integral,
                   Quantities,
                   ThisOrder + 1,
                   SourceRealT,
                   TargetRealT,
                   SourceOrder,
                   TargetOrder,
                   Offset + LoadSize,
                   SharedPosition>(
        item, target, source, start, end, newStartCoeff, newEndCoeff, shmem, reg);
  }
}

template <bool Integral,
          std::size_t SourceQuantities,
          std::size_t TargetQuantities,
          typename SourceRealT,
          typename TargetRealT,
          std::size_t SourceOrder,
          std::size_t TargetOrder>
static void taylorSumInternal(std::size_t count,
                              TargetRealT** targetBatch,
                              const SourceRealT** sourceBatch,
                              TargetRealT start,
                              TargetRealT end,
                              void* stream) {
  constexpr std::size_t Quantities = std::min(SourceQuantities, TargetQuantities);
  constexpr std::size_t TargetStride =
      seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

  sycl::nd_range rng{{count * Blocksize}, {Blocksize}};

  auto queue = reinterpret_cast<sycl::queue*>(stream);

  queue->submit([&](sycl::handler& cgh) {
    sycl::local_accessor<real> shmem(
        sycl::range<1>(
            seissol::kernels::getNumberOfAlignedBasisFunctions<SourceRealT>(SourceOrder) *
            Quantities),
        cgh);

    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      int batchId = item.get_group().get_group_id(0);

      TargetRealT reg[Quantities] = {0};
      const SourceRealT* const __restrict__ source =
          const_cast<const SourceRealT*>(sourceBatch[batchId]);
      TargetRealT* const __restrict__ target = targetBatch[batchId];
      const TargetRealT startCoeff = Integral ? start : 0;
      const TargetRealT endCoeff = Integral ? end : 1;

      taylorSumInner<Integral,
                     Quantities,
                     0,
                     SourceRealT,
                     TargetRealT,
                     SourceOrder,
                     TargetOrder,
                     0,
                     0>(
          item, target, source, start, end, startCoeff, endCoeff, shmem.get_pointer(), reg);

      constexpr std::size_t TargetStride =
          seissol::kernels::getNumberOfAlignedBasisFunctions<TargetRealT>(TargetOrder);

      if (item.get_local_id(0) < TargetStride) {
#pragma unroll
        for (std::size_t j = 0; j < Quantities; ++j) {
          target[TargetStride * j + item.get_local_id(0)] = reg[j];
        }
      }
    });
  });
}
} // namespace

namespace seissol::kernels::time::aux {
void taylorSum(bool integral,
               std::size_t count,
               real** target,
               const real** source,
               real start,
               real end,
               void* stream) {
  if (integral) {
    taylorSumInternal<true,
                      seissol::model::MaterialT::NumQuantities,
                      seissol::model::MaterialT::NumQuantities,
                      real,
                      real,
                      ConvergenceOrder,
                      ConvergenceOrder>(count, target, source, start, end, stream);
  } else {
    taylorSumInternal<false,
                      seissol::model::MaterialT::NumQuantities,
                      seissol::model::MaterialT::NumQuantities,
                      real,
                      real,
                      ConvergenceOrder,
                      ConvergenceOrder>(count, target, source, start, end, stream);
  }
}
} // namespace seissol::kernels::time::aux
#endif

namespace {
using namespace seissol::multisim;

template <typename Tensor>
constexpr size_t leadDim() {
  if constexpr (MultisimEnabled) {
    return Tensor::Stop[1] - Tensor::Start[1];
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

template <typename Tensor>
constexpr size_t linearDim() {
  if constexpr (MultisimEnabled) {
    return (Tensor::Stop[1] - Tensor::Start[1]) * (Tensor::Stop[0] - Tensor::Start[0]);
  } else {
    return Tensor::Stop[0] - Tensor::Start[0];
  }
}

auto getrange(std::size_t size, std::size_t numElements) {
  if constexpr (MultisimEnabled) {
    return sycl::nd_range<1>({numElements * NumSimulations * size}, {NumSimulations * size});
  } else {
    return sycl::nd_range<1>({numElements * size}, {size});
  }
}
} // namespace

namespace seissol::kernels::local_flux::aux::details {

void launchFreeSurfaceGravity(real** dofsFaceBoundaryNodalPtrs,
                              real** displacementDataPtrs,
                              double* rhos,
                              double g,
                              size_t numElements,
                              void* deviceStream) {

  auto queue = reinterpret_cast<sycl::queue*>(deviceStream);
  const size_t workGroupSize = leadDim<seissol::nodal::init::nodes2D>();
  auto rng = getrange(workGroupSize, numElements);

  queue->parallel_for(rng, [=](sycl::nd_item<1> item) {
    const int tid = item.get_local_id(0);
    const int elementId = item.get_group().get_group_id(0);
    if (elementId < numElements) {
      const double rho = rhos[elementId];
      real* elementBoundaryDofs = dofsFaceBoundaryNodalPtrs[elementId];
      real* elementDisplacement = displacementDataPtrs[elementId];

      constexpr auto numNodes = linearDim<seissol::nodal::init::nodes2D>();
      if (tid < numNodes) {
        constexpr auto ldINodal = linearDim<seissol::init::INodal>();

        const auto pressureAtBnd = static_cast<real>(-1.0) * rho * g * elementDisplacement[tid];

#pragma unroll
        for (int component{0}; component < 3; ++component) {
          elementBoundaryDofs[tid + component * ldINodal] =
              2.0 * pressureAtBnd - elementBoundaryDofs[tid + component * ldINodal];
        }
      }
    }
  });
}

void launchEasiBoundary(real** dofsFaceBoundaryNodalPtrs,
                        real** easiBoundaryMapPtrs,
                        real** easiBoundaryConstantPtrs,
                        size_t numElements,
                        void* deviceStream) {

  auto queue = reinterpret_cast<sycl::queue*>(deviceStream);
  const size_t workGroupSize = leadDim<seissol::init::INodal>();
  auto rng = getrange(workGroupSize, numElements);

  constexpr auto ldINodalDim = linearDim<seissol::init::INodal>();
  constexpr auto INodalDim0 = seissol::tensor::INodal::Shape[multisim::BasisFunctionDimension + 0];
  constexpr auto INodalDim1 = seissol::tensor::INodal::Shape[multisim::BasisFunctionDimension + 1];

  constexpr auto ldConstantDim = linearDim<seissol::init::easiBoundaryConstant>();
  constexpr auto ConstantDim0 =
      seissol::tensor::easiBoundaryConstant::Shape[multisim::BasisFunctionDimension + 0];
  constexpr auto ConstantDim1 =
      seissol::tensor::easiBoundaryConstant::Shape[multisim::BasisFunctionDimension + 1];

  constexpr auto ldMapDim = yateto::leadDim<seissol::init::easiBoundaryMap>();
  constexpr auto MapDim0 = seissol::tensor::easiBoundaryMap::Shape[0];
  constexpr auto MapDim1 = seissol::tensor::easiBoundaryMap::Shape[1];
  constexpr auto MapDim2 = seissol::tensor::easiBoundaryMap::Shape[2];

  static_assert(INodalDim1 == ConstantDim0, "supposed to be equal");
  static_assert(INodalDim1 == MapDim0, "supposed to be equal");

  queue->submit([&](sycl::handler& cgh) {
    sycl::local_accessor<real, 3> resultTerm(
        sycl::range<3>(INodalDim1, INodalDim0, multisim::NumSimulations), cgh);
    sycl::local_accessor<real, 3> rightTerm(
        sycl::range<3>(ConstantDim0, ConstantDim1, multisim::NumSimulations), cgh);
    sycl::local_accessor<real, 2> leftTerm(sycl::range<2>(MapDim0, MapDim2), cgh);

    cgh.parallel_for(rng, [=](sycl::nd_item<1> item) {
      const int tid = item.get_local_id(0);
      const int elementId = item.get_group().get_group_id(0);
      const int sid = tid / multisim::NumSimulations;

      if (elementId < numElements) {
        real* dofsFaceBoundaryNodal = dofsFaceBoundaryNodalPtrs[elementId];
        real* easiBoundaryMap = easiBoundaryMapPtrs[elementId];
        auto easiBoundaryConstant = easiBoundaryConstantPtrs[elementId];

        for (int i = tid; i < (ldConstantDim * ConstantDim1); i += item.get_local_range(0)) {
          const auto sim = i % multisim::NumSimulations;
          const auto subsim = i / multisim::NumSimulations;
          const auto quantity = subsim % ConstantDim0;
          const auto quadpoint = subsim / ConstantDim0;
          rightTerm[quantity][quadpoint][sim] = easiBoundaryConstant[i];
        }
        item.barrier();

        for (int i = 0; i < INodalDim1; ++i) {
          if (tid < INodalDim0) {
            resultTerm[i][tid][sid] = 0.0;
          }
        }
        item.barrier();

        for (int quantity = 0; quantity < MapDim1; ++quantity) {
          for (int quadpoint = 0; quadpoint < MapDim2; ++quadpoint) {
            if (tid < MapDim0) {
              leftTerm[tid][quadpoint] =
                  easiBoundaryMap[tid + ldMapDim * (quantity + quadpoint * MapDim1)];
            }
          }
          item.barrier();

          if (tid < MapDim2) {
            const real col = dofsFaceBoundaryNodal[tid + quantity * ldINodalDim];
            for (int quantity2 = 0; quantity2 < MapDim0; ++quantity2) {
              resultTerm[quantity2][tid][sid] += leftTerm[quantity2][tid] * col;
            }
          }
          item.barrier();
        }

        if (tid < INodalDim0) {
          for (int quantity2 = 0; quantity2 < INodalDim1; ++quantity2) {
            dofsFaceBoundaryNodal[tid + quantity2 * ldINodalDim] =
                resultTerm[quantity2][tid][sid] + rightTerm[quantity2][tid][sid];
          }
        }
      }
    });
  });
}
} // namespace seissol::kernels::local_flux::aux::details

namespace seissol::kernels::time::aux {
void extractRotationMatrices(real** displacementToFaceNormalPtrs,
                             real** displacementToGlobalDataPtrs,
                             real** TPtrs,
                             real** TinvPtrs,
                             size_t numElements,
                             void* deviceStream) {
  auto queue = reinterpret_cast<sycl::queue*>(deviceStream);

  constexpr size_t workGroupSize = 9;
  sycl::nd_range<1> rng{{numElements * workGroupSize}, {workGroupSize}};

  queue->parallel_for(rng, [=](sycl::nd_item<1> item) {
    const int elementId = item.get_group().get_group_id(0);
    if (elementId < numElements) {
      real* displacementToFaceNormal = displacementToFaceNormalPtrs[elementId];
      real* displacementToGlobalData = displacementToGlobalDataPtrs[elementId];
      auto* T = TPtrs[elementId];
      auto* Tinv = TinvPtrs[elementId];

      constexpr auto ldTinv = yateto::leadDim<seissol::init::Tinv>();
      constexpr auto ldT = yateto::leadDim<seissol::init::T>();
      constexpr auto ldDisplacement = yateto::leadDim<seissol::init::displacementRotationMatrix>();

      const int i = item.get_local_id(0) % 3;
      const int j = item.get_local_id(0) / 3;

      displacementToFaceNormal[i + j * ldDisplacement] = Tinv[(i + 6) + (j + 6) * ldTinv];
      displacementToGlobalData[i + j * ldDisplacement] = T[(i + 6) + (j + 6) * ldT];
    }
  });
}

void initializeTaylorSeriesForGravitationalBoundary(real** prevCoefficientsPtrs,
                                                    real** integratedDisplacementNodalPtrs,
                                                    real** rotatedFaceDisplacementPtrs,
                                                    double deltaTInt,
                                                    size_t numElements,
                                                    void* deviceStream) {

  auto queue = reinterpret_cast<sycl::queue*>(deviceStream);
  const size_t workGroupSize = leadDim<seissol::nodal::init::nodes2D>();
  auto rng = getrange(workGroupSize, numElements);

  queue->parallel_for(rng, [=](sycl::nd_item<1> item) {
    const int elementId = item.get_group().get_group_id(0);
    if (elementId < numElements) {
      auto* prevCoefficients = prevCoefficientsPtrs[elementId];
      auto* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
      const auto* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];

      assert(linearDim<seissol::nodal::init::nodes2D>() <=
             linearDim<seissol::init::rotatedFaceDisplacement>());

      const int tid = item.get_local_id(0);
      constexpr auto num2dNodes = linearDim<seissol::nodal::init::nodes2D>();
      if (tid < num2dNodes) {
        prevCoefficients[tid] = rotatedFaceDisplacement[tid];
        integratedDisplacementNodal[tid] = deltaTInt * rotatedFaceDisplacement[tid];
      }
    }
  });
}

void computeInvAcousticImpedance(
    double* invImpedances, double* rhos, double* lambdas, size_t numElements, void* deviceStream) {
  auto queue = reinterpret_cast<sycl::queue*>(deviceStream);
  sycl::nd_range rng({numElements}, {1024});

  queue->parallel_for(rng, [=](sycl::nd_item<1> item) {
    size_t index = item.get_global_id(0);
    if (index < numElements) {
      invImpedances[index] = 1.0 / sycl::sqrt(lambdas[index] * rhos[index]);
    }
  });
}

void updateRotatedFaceDisplacement(real** rotatedFaceDisplacementPtrs,
                                   real** prevCoefficientsPtrs,
                                   real** integratedDisplacementNodalPtrs,
                                   real** dofsFaceNodalPtrs,
                                   double* invImpedances,
                                   double* rhos,
                                   double g,
                                   double factorEvaluated,
                                   double factorInt,
                                   size_t numElements,
                                   void* deviceStream) {

  auto queue = reinterpret_cast<sycl::queue*>(deviceStream);
  const size_t workGroupSize = leadDim<seissol::nodal::init::nodes2D>();
  auto rng = getrange(workGroupSize, numElements);

  queue->parallel_for(rng, [=](sycl::nd_item<1> item) {
    const int elementId = item.get_group().get_group_id(0);
    if (elementId < numElements) {
      constexpr int pIdx = 0;
      constexpr int uIdx = 6;
      constexpr auto num2dNodes = linearDim<seissol::nodal::init::nodes2D>();

      const int tid = item.get_local_id(0);
      if (tid < num2dNodes) {

        real* dofsFaceNodal = dofsFaceNodalPtrs[elementId];
        constexpr auto ldINodal = linearDim<seissol::init::INodal>();

        const auto uInside = dofsFaceNodal[tid + (uIdx + 0) * ldINodal];
        const auto vInside = dofsFaceNodal[tid + (uIdx + 1) * ldINodal];
        const auto wInside = dofsFaceNodal[tid + (uIdx + 2) * ldINodal];
        const auto pressureInside = dofsFaceNodal[tid + pIdx * ldINodal];

        real* prevCoefficients = prevCoefficientsPtrs[elementId];
        const auto rho = rhos[elementId];
        const auto invImpedance = invImpedances[elementId];

        const double curCoeff =
            uInside - invImpedance * (rho * g * prevCoefficients[tid] + pressureInside);
        prevCoefficients[tid] = curCoeff;

        constexpr auto ldFaceDisplacement = linearDim<seissol::init::faceDisplacement>();
        static_assert(num2dNodes <= ldFaceDisplacement, "");

        real* rotatedFaceDisplacement = rotatedFaceDisplacementPtrs[elementId];
        rotatedFaceDisplacement[tid + 0 * ldFaceDisplacement] += factorEvaluated * curCoeff;
        rotatedFaceDisplacement[tid + 1 * ldFaceDisplacement] += factorEvaluated * vInside;
        rotatedFaceDisplacement[tid + 2 * ldFaceDisplacement] += factorEvaluated * wInside;

        constexpr auto ldIntegratedFaceDisplacement =
            linearDim<seissol::init::averageNormalDisplacement>();
        static_assert(num2dNodes <= ldIntegratedFaceDisplacement, "");

        real* integratedDisplacementNodal = integratedDisplacementNodalPtrs[elementId];
        integratedDisplacementNodal[tid] += factorInt * curCoeff;
      }
    }
  });
}
} // namespace seissol::kernels::time::aux
