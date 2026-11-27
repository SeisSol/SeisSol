// SPDX-FileCopyrightText: 2017 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Stephanie Wollherr

#include "Plasticity.h"

#include "Alignment.h"
#include "Config.h"
#include "Core/Common/Marker.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Memory/Batch/ConditionalTable.h"
#include "Core/Memory/Typedefs.h"
#include "Core/Monitoring/Metric.h"
#include "Core/Parallel/Runtime/Stream.h"
#include "Features/Plasticity/Physics/DruckerPrager.h"
#include "Features/Plasticity/Physics/Plasticity.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <algorithm>
#include <cassert>
#include <cstddef>
#include <cstdint>
#include <utility>

#ifdef ACL_DEVICE
#include "Core/Memory/Batch/ConditionalKey.h"
#include "Core/Memory/Batch/EncodedConstants.h"
#include "DeviceAux/PlasticityAux.h"
using namespace device;
#endif

#ifndef ACL_DEVICE
#include <utils/logger.h>
#endif

namespace seissol::kernels {
template <typename Cfg>
std::size_t
    Plasticity<Cfg>::computePlasticity(real oneMinusIntegratingFactor,
                                       const GlobalData<Cfg>* global,
                                       const seissol::model::PlasticityData<Cfg>* plasticityData,
                                       real degreesOfFreedom[tensor::Q<Cfg>::size()],
                                       real* pstrain) {
  // a material without the stresses of a solid has nothing that yields
  if constexpr (!model::MaterialOf<Cfg>::SupportsPlasticity) {
    return 0;
  }

  assert(reinterpret_cast<uintptr_t>(degreesOfFreedom) % Vectorsize == 0);

  alignas(Alignment) real qStressNodal[tensor::QStressNodal<Cfg>::size()]{};

  // the nodes of all simulations of the cell; the stress components of a node lie this far apart
  constexpr std::size_t NumNodes = init::QStressNodal<Cfg>::Stop[multisim::BasisDim<Cfg>] -
                                   init::QStressNodal<Cfg>::Start[multisim::BasisDim<Cfg>];
  constexpr auto NumTotalPoints = NumNodes * Cfg::NumSimulations;
  static_assert(tensor::yieldFactor<Cfg>::size() == NumTotalPoints &&
                    tensor::meanStress<Cfg>::size() == NumTotalPoints,
                "the yield factor and the mean stress have one value per node");
  static_assert(tensor::QStressNodal<Cfg>::size() ==
                    NumTotalPoints * model::plasticity::StressComponents,
                "the nodal stresses hold the six stress components of every node");

  alignas(Alignment) real yieldFactor[NumTotalPoints]{};

  /* Convert modal to nodal and add sigma0.
   * Stores s_{ij} := sigma_{ij} + sigma0_{ij} for every node.
   * sigma0 is constant
   * also stores the previous DOFs before adding the new initial loading
   */

  kernel::plConvertToNodal<Cfg> m2nKrnl;
  m2nKrnl.bindGlobals(*global);
  m2nKrnl.QStress = degreesOfFreedom;
  m2nKrnl.QStressNodal = qStressNodal;
  m2nKrnl.initialLoading = plasticityData->initialLoading;
  m2nKrnl.execute();

  // the yield criterion at every node (see model::plasticity::yieldOfNode); the nodal stresses
  // keep their deviatoric part for the relaxation
  std::int32_t adjust = 0;

#pragma omp simd reduction(max : adjust)
  for (std::size_t ip = 0; ip < NumTotalPoints; ++ip) {
    real stress[model::plasticity::StressComponents];
    for (int x = 0; x < model::plasticity::StressComponents; ++x) {
      stress[x] = qStressNodal[ip + NumTotalPoints * x];
    }
    const auto node =
        model::plasticity::yieldOfNode(stress,
                                       plasticityData->cohesionTimesCosAngularFriction[ip],
                                       plasticityData->sinAngularFriction[ip],
                                       oneMinusIntegratingFactor);
    for (int x = 0; x < model::plasticity::StressComponents; ++x) {
      qStressNodal[ip + NumTotalPoints * x] = node.deviator[x];
    }
    yieldFactor[ip] = node.factor;
    // Combine with the previous value: under the reduction, each SIMD lane keeps only its own
    // last value, so a plain assignment would discard all but the last chunk of nodes.
    adjust = std::max(adjust, static_cast<std::int32_t>(node.yields ? 1 : 0));
  }

  if (adjust != 0) {
    real* __restrict qEtaNodal = &pstrain[tensor::QStressNodal<Cfg>::size()];

    // the relaxation of every node (see model::plasticity::relaxNode); the nodal stresses hold
    // the update of the stress, yield s_{ij}, afterwards, which is zero where a node does not
    // yield
#pragma omp simd
    for (std::size_t qp = 0; qp < NumTotalPoints; ++qp) {
      model::plasticity::NodeYield<real> node{};
      for (int x = 0; x < model::plasticity::StressComponents; ++x) {
        node.deviator[x] = qStressNodal[qp + NumTotalPoints * x];
      }
      node.factor = yieldFactor[qp];

      real stressUpdate[model::plasticity::StressComponents];
      real pstrainUpdate[model::plasticity::StressComponents];
      qEtaNodal[qp] +=
          model::plasticity::relaxNode(node, plasticityData->mufactor, stressUpdate, pstrainUpdate);

      for (int x = 0; x < model::plasticity::StressComponents; ++x) {
        const auto q = qp + NumTotalPoints * x;
        pstrain[q] += pstrainUpdate[x];
        qStressNodal[q] = stressUpdate[x];
      }
    }

    // added to the DOFs once, for the reason given in Local::computeIntegral
    alignas(Alignment) real stressUpdate[tensor::QStress<Cfg>::size()]{};
    kernel::plConvertToModal<Cfg> adjKrnl;
    adjKrnl.QStress = stressUpdate;
    adjKrnl.bindGlobals(*global);
    adjKrnl.QStressNodal = qStressNodal;
    adjKrnl.execute();
#pragma omp simd
    for (std::size_t i = 0; i < tensor::QStress<Cfg>::size(); ++i) {
      degreesOfFreedom[i] += stressUpdate[i];
    }

    return 1;
  }
  return 0;
}

template <typename Cfg>
void Plasticity<Cfg>::computePlasticityBatched(
    SEISSOL_GPU_PARAM real timeStepWidth,
    SEISSOL_GPU_PARAM real tV,
    SEISSOL_GPU_PARAM const GlobalData<Cfg>* global,
    SEISSOL_GPU_PARAM recording::ConditionalTable& table,
    SEISSOL_GPU_PARAM seissol::model::PlasticityData<Cfg>* plasticityData,
    SEISSOL_GPU_PARAM std::size_t* yieldCounter,
    SEISSOL_GPU_PARAM unsigned* isAdjustableVector,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // a material without the stresses of a solid has nothing that yields
  if constexpr (!model::MaterialOf<Cfg>::SupportsPlasticity) {
    return;
  }

  using namespace seissol::recording;

  // the number of nodes is that of the plasticity point set and need not be the number of basis
  // functions; buffers and kernels take it from QStressNodal
  static_assert(tensor::QStressNodal<Cfg>::Shape[multisim::BasisDim<Cfg>] ==
                    tensor::v<Cfg>::Shape[0],
                "nodal stresses and Vandermonde matrix must have the same number of nodes");
  static_assert(tensor::Q<Cfg>::Shape[multisim::BasisDim<Cfg>] == tensor::v<Cfg>::Shape[1],
                "modal dofs and Vandermonde matrix must have the same number of basis functions");

  const ConditionalKey key(*KernelNames::Plasticity);
  auto* defaultStream = runtime.stream();

  if (table.contains(key)) {
    const auto oneMinusIntegratingFactor = computeRelaxTime(tV, timeStepWidth);

    auto& entry = table[key];
    const size_t numElements = entry.get<real*>(inner_keys::Wp::Id::Dofs)->getSize();

    // Convert modal to nodal
    real** modalStressTensors = entry.get<real*>(inner_keys::Wp::Id::Dofs)->getDeviceDataPtr();
    real** nodalStressTensors =
        entry.get<real*>(inner_keys::Wp::Id::NodalStressTensor)->getDeviceDataPtr();

    static_assert(kernel::gpu_plConvertToNodal<Cfg>::TmpMaxMemRequiredInBytes == 0);
    real** initLoad = entry.get<real*>(inner_keys::Wp::Id::InitialLoad)->getDeviceDataPtr();
    kernel::gpu_plConvertToNodal<Cfg> m2nKrnl;
    m2nKrnl.bindGlobals(*global);
    m2nKrnl.QStress = const_cast<const real**>(modalStressTensors);
    m2nKrnl.QStressNodal = nodalStressTensors;
    m2nKrnl.initialLoading = const_cast<const real**>(initLoad);
    m2nKrnl.streamPtr = defaultStream;
    m2nKrnl.numElements = numElements;
    m2nKrnl.execute();

    real** pstrains = entry.get<real*>(inner_keys::Wp::Id::Pstrains)->getDeviceDataPtr();
    // the batch holds only the plastic cells; plasticityData is indexed by the layer cell
    const unsigned* cellIndices =
        table[key].get<unsigned>(inner_keys::Indices::Id::Cells)->getDeviceDataPtr();

    device::aux::plasticity::plasticityNonlinear(nodalStressTensors,
                                                 pstrains,
                                                 isAdjustableVector,
                                                 yieldCounter,
                                                 plasticityData,
                                                 cellIndices,
                                                 oneMinusIntegratingFactor,
                                                 numElements,
                                                 defaultStream);

    kernel::gpu_plConvertToModal<Cfg> n2mKrnl;
    n2mKrnl.bindGlobals(*global);
    n2mKrnl.QStressNodal = const_cast<const real**>(nodalStressTensors);
    n2mKrnl.QStress = modalStressTensors;
    n2mKrnl.streamPtr = defaultStream;
    n2mKrnl.flags = isAdjustableVector;
    n2mKrnl.numElements = numElements;
    n2mKrnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif // ACL_DEVICE
}

template <typename Cfg>
std::pair<PerformanceEstimate, PerformanceEstimate> Plasticity<Cfg>::metrics() {
  // reset flops
  PerformanceEstimate check;
  PerformanceEstimate yield;

  // flops from checking, i.e. outside if (adjust) {}
  check += PerformanceEstimate::fromKernel<kernel::plConvertToNodal<Cfg>>();

  // compute mean stress
  check += PerformanceEstimate::fromKernel<kernel::plComputeMean<Cfg>>();

  // subtract mean stress
  check += PerformanceEstimate::fromKernel<kernel::plSubtractMean<Cfg>>();

  // compute second invariant
  check += PerformanceEstimate::fromKernel<kernel::plComputeSecondInvariant<Cfg>>();

  // compute taulim (1 add, 1 mul, max NOT counted)
  check.nonzeroFlop += static_cast<std::uint64_t>(2 * tensor::meanStress<Cfg>::size());
  check.hardwareFlop += static_cast<std::uint64_t>(2 * tensor::meanStress<Cfg>::size());

  // check for yield (NOT counted, as it would require counting the number of yielding points)

  // flops from plastic yielding, i.e. inside if (adjust) {}
  yield += PerformanceEstimate::fromKernel<kernel::plConvertToModal<Cfg>>();

  // manually counted
  yield.nonzeroFlop += static_cast<std::uint64_t>(tensor::QStressNodal<Cfg>::size() * 6);
  yield.hardwareFlop += static_cast<std::uint64_t>(tensor::QStressNodal<Cfg>::size() * 6);
  yield.nonzeroFlop += static_cast<std::uint64_t>(tensor::QEtaNodal<Cfg>::size() * 3);
  yield.hardwareFlop += static_cast<std::uint64_t>(tensor::QEtaNodal<Cfg>::size() * 3);

  return {check, yield};
}
#define SEISSOL_INSTANTIATE(Cfg) template class Plasticity<Cfg>;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels
