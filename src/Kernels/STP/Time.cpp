// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Time.h"

#include "Common/Marker.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Common.h"
#include "Kernels/MemoryOps.h"
#include "Kernels/STP/Setup.h"
#include "Monitoring/Metric.h"

#include <Eigen/Dense>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <stdint.h>
#include <yateto.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

namespace seissol::kernels::solver::stp {

template <typename Cfg>
void Spacetime<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {
  krnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

template <typename Cfg>
void Spacetime<Cfg>::executeSTP(double timeStepWidth,
                                LTS::Ref<Cfg>& data,
                                real* timeIntegrated,
                                real* stp)

{
  assert((reinterpret_cast<uintptr_t>(stp)) % Vectorsize == 0);
  std::fill(stp, stp + tensor::spaceTimePredictor<Cfg>::size(), 0);
  kernel::spaceTimePredictor<Cfg> krnl = krnlPrototype_;

  // libxsmm can not generate GEMMs with alpha!=1. As a workaround we multiply the
  // star matrices with dt before we execute the kernel.
  real A_values[init::star<Cfg>::size(0)];
  real B_values[init::star<Cfg>::size(1)];
  real C_values[init::star<Cfg>::size(2)];
  for (std::size_t i = 0; i < init::star<Cfg>::size(0); i++) {
    A_values[i] = timeStepWidth * data.template get<LTS::LocalIntegration>().starMatrices[0][i];
    B_values[i] = timeStepWidth * data.template get<LTS::LocalIntegration>().starMatrices[1][i];
    C_values[i] = timeStepWidth * data.template get<LTS::LocalIntegration>().starMatrices[2][i];
  }
  krnl.star(0) = A_values;
  krnl.star(1) = B_values;
  krnl.star(2) = C_values;

  for (std::size_t i = 0; i < model::MaterialOf<Cfg>::StiffSourceRows.size(); ++i) {
    krnl.G(i) = data.template get<LTS::LocalIntegration>().specific.G[i] * timeStepWidth;
  }

  krnl.Q = const_cast<real*>(data.template get<LTS::Dofs>());
  krnl.I = timeIntegrated;
  krnl.timestep = timeStepWidth;
  krnl.spaceTimePredictor = stp;

  // The matrix Zinv depends on the timestep
  // If the timestep is not as expected e.g. when approaching a sync point
  // we have to recalculate it

  // beware of float comparison errors with the timestep

  const auto defaultTimestep =
      std::abs((data.template get<LTS::LocalIntegration>().specific.typicalTimeStepWidth -
                timeStepWidth) /
               timeStepWidth) < 1e-7;

  // the members of the Zinv family are stored back to back
  const auto zinvOffset = [](std::size_t i) {
    return yateto::computeFamilySize<tensor::Zinv<Cfg>>(1, i);
  };

  if (!defaultTimestep) {
    auto sourceMatrix = init::ET<Cfg>::view::create(
        data.template get<LTS::LocalIntegration>().specific.sourceMatrix);
    real zinvData[kernels::familySize<tensor::Zinv<Cfg>>()];
    model::ZInvInitializer<Cfg,
                           seissol::model::MaterialOf<Cfg>,
                           0,
                           seissol::model::MaterialOf<Cfg>::NumQuantities,
                           decltype(sourceMatrix)>(zinvData, sourceMatrix, timeStepWidth);
    for (std::size_t i = 0; i < seissol::model::MaterialOf<Cfg>::NumQuantities; i++) {
      krnl.Zinv(i) = zinvData + zinvOffset(i);
    }
    // krnl.execute has to be run here: zinvData is only allocated locally
    krnl.execute();
  } else {
    const real* zinvData = data.template get<LTS::LocalIntegration>().specific.Zinv;
    for (std::size_t i = 0; i < seissol::model::MaterialOf<Cfg>::NumQuantities; i++) {
      krnl.Zinv(i) = zinvData + zinvOffset(i);
    }
    krnl.execute();
  }
}

template <typename Cfg>
void Spacetime<Cfg>::computeAder(const real* coeffs,
                                 double timeStepWidth,
                                 LTS::Ref<Cfg>& data,
                                 LocalTmp<Cfg>& tmp,
                                 real* timeIntegrated,
                                 real* timeDerivatives,
                                 bool updateDisplacement) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.template get<LTS::Dofs>())) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Vectorsize == 0 ||
         timeDerivatives == nullptr);

  alignas(Alignment) real temporaryBuffer[tensor::spaceTimePredictor<Cfg>::size()];
  real* stpBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;
  executeSTP(timeStepWidth, data, timeIntegrated, stpBuffer);
}

template <typename Cfg>
PerformanceEstimate Spacetime<Cfg>::metrics() const {
  auto estimate = PerformanceEstimate::fromKernel<kernel::spaceTimePredictor<Cfg>>();

  estimate.nonzeroFlop += 3 * init::star<Cfg>::size(0);
  estimate.hardwareFlop += 3 * init::star<Cfg>::size(0);

  // legacy memory estimate
  std::uint64_t reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q<Cfg>::size() + 2 * tensor::I<Cfg>::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star<Cfg>>();
  // Zinv
  reals += yateto::computeFamilySize<tensor::Zinv<Cfg>>();
  // G
  reals += 3;

  /// \todo incorporate derivatives

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

template <typename Cfg>
void Spacetime<Cfg>::computeBatchedAder(
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM LocalTmp<Cfg>& tmp,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM bool updateDisplacement,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  kernel::gpu_spaceTimePredictor<Cfg> krnl = deviceKrnlPrototype_;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get<real*>(inner_keys::Wp::Id::Dofs))->getSize();
    krnl.numElements = numElements;

    krnl.I = (entry.get<real*>(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
    krnl.Q =
        const_cast<const real**>((entry.get<real*>(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr());
    krnl.timestep = timeStepWidth;

    krnl.spaceTimePredictor =
        (entry.get<real*>(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, starMatrices);
    for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      krnl.star(i) = const_cast<const real**>(
          (entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      krnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, starMatrices, i);
    }

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.G);
    for (std::size_t i = 0; i < model::MaterialOf<Cfg>::StiffSourceRows.size(); ++i) {
      krnl.Gt(i) = const_cast<const real**>(
          (entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      krnl.extraOffset_Gt(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Cfg>, specific.G, i);
    }

    // checking the first cell should suffice; if we always work on the same cluster.
    // (which we currently always do)
    const auto defaultTimestep =
        std::abs((layer.var<LTS::LocalIntegration>(Cfg())[0].specific.typicalTimeStepWidth -
                  timeStepWidth) /
                 timeStepWidth) < 1e-7;

    if (defaultTimestep) {
      // Zinv is one flat family, so its members are not spaced by the size of a single entry
      SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Cfg>, specific.Zinv);
      for (std::size_t i = 0; i < seissol::model::MaterialOf<Cfg>::NumQuantities; ++i) {
        krnl.Zinv(i) = const_cast<const real**>(
            (entry.get<real*>(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
        krnl.extraOffset_Zinv(i) = SEISSOL_OFFSET(LocalIntegrationData<Cfg>, specific.Zinv) +
                                   yateto::computeFamilySize<tensor::Zinv<Cfg>>(1, i);
      }
    } else {
      auto* layerZinvData = layer.var<LTS::ZinvExtra>(Cfg());
      const auto* layerLocalIntegration = layer.var<LTS::LocalIntegration>(Cfg());
      runtime.enqueueLoop(numElements, [=](std::size_t i) {
        auto* zinvData = layerZinvData + yateto::computeFamilySize<tensor::Zinv<Cfg>>() * i;
        const auto& localIntegration = layerLocalIntegration[i];

        const auto sourceMatrix =
            init::ET<Cfg>::view::create(localIntegration.specific.sourceMatrix);
        model::ZInvInitializer<Cfg,
                               seissol::model::MaterialOf<Cfg>,
                               0,
                               seissol::model::MaterialOf<Cfg>::NumQuantities,
                               decltype(sourceMatrix)>(zinvData, sourceMatrix, timeStepWidth);
      });
      for (std::size_t i = 0; i < seissol::model::MaterialOf<Cfg>::NumQuantities; ++i) {
        krnl.Zinv(i) = const_cast<const real**>(
            (entry.get<real*>(inner_keys::Wp::Id::ZinvExtra))->getDeviceDataPtr());
        krnl.extraOffset_Zinv(i) = yateto::computeFamilySize<tensor::Zinv<Cfg>>(1, i);
      }
    }

    krnl.streamPtr = runtime.stream();

    // TODO: integrate into the following kernel
    device::DeviceInstance::instance().algorithms().setToValue(
        krnl.spaceTimePredictor,
        static_cast<real>(0.0),
        tensor::spaceTimePredictor<Cfg>::size(),
        krnl.numElements,
        krnl.streamPtr);

    krnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
void Time<Cfg>::evaluate(const real* coeffs, const real* timeDerivatives, real* timeEvaluated) {
  kernel::evaluateDOFSAtTimeSTP<Cfg> krnl;
  krnl.spaceTimePredictor = timeDerivatives;
  krnl.QAtTimeSTP = timeEvaluated;
  krnl.timeBasisFunctionsAtPoint = coeffs;
  krnl.execute();
}

template <typename Cfg>
void Time<Cfg>::evaluateBatched(const real* coeffs,
                                const real** timeDerivatives,
                                real** timeIntegratedDofs,
                                std::size_t numElements,
                                seissol::parallel::runtime::StreamRuntime& runtime) {
  // for now, use the Taylor kernel here; since it'll do exactly the same as in the LinearCK case.
  // if there are any errors, check this one again.

#ifdef ACL_DEVICE
  kernel::gpu_derivativeTaylorExpansion<Cfg> krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ<Cfg>>(); ++i) {
    krnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    krnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ<Cfg>>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.streamPtr = runtime.stream();
  krnl.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

template <typename Cfg>
PerformanceEstimate Time<Cfg>::metrics() const {
  return PerformanceEstimate::fromKernel<kernel::evaluateDOFSAtTimeSTP<Cfg>>();
}

template <typename Cfg>
void Time<Cfg>::setGlobalData(const CompoundGlobalData<Cfg>& global) {}

#define SEISSOL_INSTANTIATE(Cfg)                                                                   \
  template class Spacetime<Cfg>;                                                                   \
  template class Time<Cfg>;
SEISSOL_FOR_EACH_CONFIG_STP(SEISSOL_INSTANTIATE)
#undef SEISSOL_INSTANTIATE

} // namespace seissol::kernels::solver::stp
