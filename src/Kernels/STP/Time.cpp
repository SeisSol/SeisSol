// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Time.h"

#include "Common/Marker.h"
#include "Config.h"
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

void Spacetime::setGlobalData(const CompoundGlobalData& global) {
  krnlPrototype_.bindGlobals(*global.onHost);

#ifdef ACL_DEVICE
  deviceKrnlPrototype_.bindGlobals(*global.onDevice);
#endif
}

void Spacetime::executeSTP(double timeStepWidth,
                           LTS::Ref<Config>& data,
                           real* timeIntegrated,
                           real* stp)

{
  assert((reinterpret_cast<uintptr_t>(stp)) % Vectorsize == 0);
  std::fill(stp, stp + tensor::spaceTimePredictor<Config>::size(), 0);
  kernel::spaceTimePredictor<Config> krnl = krnlPrototype_;

  // libxsmm can not generate GEMMs with alpha!=1. As a workaround we multiply the
  // star matrices with dt before we execute the kernel.
  real A_values[init::star<Config>::size(0)];
  real B_values[init::star<Config>::size(1)];
  real C_values[init::star<Config>::size(2)];
  for (std::size_t i = 0; i < init::star<Config>::size(0); i++) {
    A_values[i] = timeStepWidth * data.get<LTS::LocalIntegration>().starMatrices[0][i];
    B_values[i] = timeStepWidth * data.get<LTS::LocalIntegration>().starMatrices[1][i];
    C_values[i] = timeStepWidth * data.get<LTS::LocalIntegration>().starMatrices[2][i];
  }
  krnl.star(0) = A_values;
  krnl.star(1) = B_values;
  krnl.star(2) = C_values;

  for (std::size_t i = 0; i < generated::StiffSourceRowCount; ++i) {
    krnl.G(i) = data.get<LTS::LocalIntegration>().specific.G[i] * timeStepWidth;
  }

  krnl.Q = const_cast<real*>(data.get<LTS::Dofs>());
  krnl.I = timeIntegrated;
  krnl.timestep = timeStepWidth;
  krnl.spaceTimePredictor = stp;

  // The matrix Zinv depends on the timestep
  // If the timestep is not as expected e.g. when approaching a sync point
  // we have to recalculate it

  // beware of float comparison errors with the timestep

  const auto defaultTimestep =
      std::abs((data.get<LTS::LocalIntegration>().specific.typicalTimeStepWidth - timeStepWidth) /
               timeStepWidth) < 1e-7;

  // the members of the Zinv family are stored back to back
  const auto zinvOffset = [](std::size_t i) {
    return yateto::computeFamilySize<tensor::Zinv<Config>>(1, i);
  };

  if (!defaultTimestep) {
    auto sourceMatrix =
        init::ET<Config>::view::create(data.get<LTS::LocalIntegration>().specific.sourceMatrix);
    real zinvData[kernels::familySize<tensor::Zinv<Config>>()];
    model::ZInvInitializer<Config,
                           seissol::model::MaterialT,
                           0,
                           seissol::model::MaterialT::NumQuantities,
                           decltype(sourceMatrix)>(zinvData, sourceMatrix, timeStepWidth);
    for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
      krnl.Zinv(i) = zinvData + zinvOffset(i);
    }
    // krnl.execute has to be run here: zinvData is only allocated locally
    krnl.execute();
  } else {
    const real* zinvData = data.get<LTS::LocalIntegration>().specific.Zinv;
    for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
      krnl.Zinv(i) = zinvData + zinvOffset(i);
    }
    krnl.execute();
  }
}

void Spacetime::computeAder(const real* coeffs,
                            double timeStepWidth,
                            LTS::Ref<Config>& data,
                            LocalTmp& tmp,
                            real* timeIntegrated,
                            real* timeDerivatives,
                            bool updateDisplacement) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>())) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Vectorsize == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Vectorsize == 0 ||
         timeDerivatives == nullptr);

  alignas(Alignment) real temporaryBuffer[tensor::spaceTimePredictor<Config>::size()];
  real* stpBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;
  executeSTP(timeStepWidth, data, timeIntegrated, stpBuffer);
}

PerformanceEstimate Spacetime::metrics() const {
  auto estimate = PerformanceEstimate::fromKernel<kernel::spaceTimePredictor<Config>>();

  estimate.nonzeroFlop += 3 * init::star<Config>::size(0);
  estimate.hardwareFlop += 3 * init::star<Config>::size(0);

  // legacy memory estimate
  std::uint64_t reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q<Config>::size() + 2 * tensor::I<Config>::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star<Config>>();
  // Zinv
  reals += yateto::computeFamilySize<tensor::Zinv<Config>>();
  // G
  reals += 3;

  /// \todo incorporate derivatives

  estimate.bytes = reals * sizeof(real);

  return estimate;
}

void Spacetime::computeBatchedAder(
    SEISSOL_GPU_PARAM const real* coeffs,
    SEISSOL_GPU_PARAM double timeStepWidth,
    SEISSOL_GPU_PARAM LTS::Layer& layer,
    SEISSOL_GPU_PARAM LocalTmp& tmp,
    SEISSOL_GPU_PARAM recording::ConditionalPointersToRealsTable& dataTable,
    SEISSOL_GPU_PARAM bool updateDisplacement,
    SEISSOL_GPU_PARAM seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE

  using namespace seissol::recording;
  kernel::gpu_spaceTimePredictor<Config> krnl = deviceKrnlPrototype_;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    krnl.numElements = numElements;

    krnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
    krnl.Q = const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr());
    krnl.timestep = timeStepWidth;

    krnl.spaceTimePredictor = (entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr();

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Config>, starMatrices);
    for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::star<Config>>(); ++i) {
      krnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      krnl.extraOffset_star(i) =
          SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Config>, starMatrices, i);
    }

    SEISSOL_ARRAY_OFFSET_ASSERT(LocalIntegrationData<Config>, specific.G);
    for (std::size_t i = 0; i < generated::StiffSourceRowCount; ++i) {
      krnl.Gt(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      krnl.extraOffset_Gt(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData<Config>, specific.G, i);
    }

    // checking the first cell should suffice; if we always work on the same cluster.
    // (which we currently always do)
    const auto defaultTimestep =
        std::abs((layer.var<LTS::LocalIntegration>(Config())[0].specific.typicalTimeStepWidth -
                  timeStepWidth) /
                 timeStepWidth) < 1e-7;

    if (defaultTimestep) {
      // Zinv is one flat family, so its members are not spaced by the size of a single entry
      SEISSOL_OFFSET_ASSERT(LocalIntegrationData<Config>, specific.Zinv);
      for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; ++i) {
        krnl.Zinv(i) = const_cast<const real**>(
            (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
        krnl.extraOffset_Zinv(i) = SEISSOL_OFFSET(LocalIntegrationData<Config>, specific.Zinv) +
                                   yateto::computeFamilySize<tensor::Zinv<Config>>(1, i);
      }
    } else {
      auto* layerZinvData = layer.var<LTS::ZinvExtra>(Config());
      const auto* layerLocalIntegration = layer.var<LTS::LocalIntegration>(Config());
      runtime.enqueueLoop(numElements, [=](std::size_t i) {
        auto* zinvData = layerZinvData + yateto::computeFamilySize<tensor::Zinv<Config>>() * i;
        const auto& localIntegration = layerLocalIntegration[i];

        const auto sourceMatrix =
            init::ET<Config>::view::create(localIntegration.specific.sourceMatrix);
        model::ZInvInitializer<Config,
                               seissol::model::MaterialT,
                               0,
                               seissol::model::MaterialT::NumQuantities,
                               decltype(sourceMatrix)>(zinvData, sourceMatrix, timeStepWidth);
      });
      for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; ++i) {
        krnl.Zinv(i) = const_cast<const real**>(
            (entry.get(inner_keys::Wp::Id::ZinvExtra))->getDeviceDataPtr());
        krnl.extraOffset_Zinv(i) = yateto::computeFamilySize<tensor::Zinv<Config>>(1, i);
      }
    }

    krnl.streamPtr = runtime.stream();

    // TODO: integrate into the following kernel
    device.algorithms().setToValue(krnl.spaceTimePredictor,
                                   static_cast<real>(0.0),
                                   tensor::spaceTimePredictor<Config>::size(),
                                   krnl.numElements,
                                   krnl.streamPtr);

    krnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Time::evaluate(const real* coeffs, const real* timeDerivatives, real* timeEvaluated) {
  kernel::evaluateDOFSAtTimeSTP<Config> krnl;
  krnl.spaceTimePredictor = timeDerivatives;
  krnl.QAtTimeSTP = timeEvaluated;
  krnl.timeBasisFunctionsAtPoint = coeffs;
  krnl.execute();
}

void Time::evaluateBatched(const real* coeffs,
                           const real** timeDerivatives,
                           real** timeIntegratedDofs,
                           std::size_t numElements,
                           seissol::parallel::runtime::StreamRuntime& runtime) {
  // for now, use the Taylor kernel here; since it'll do exactly the same as in the LinearCK case.
  // if there are any errors, check this one again.

#ifdef ACL_DEVICE
  kernel::gpu_derivativeTaylorExpansion<Config> krnl;
  krnl.numElements = numElements;
  krnl.I = timeIntegratedDofs;
  for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ<Config>>(); ++i) {
    krnl.dQ(i) = const_cast<const real**>(timeDerivatives);
    krnl.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ<Config>>(1, i);
    krnl.power(i) = coeffs[i];
  }
  krnl.streamPtr = runtime.stream();
  krnl.execute();
#else
  logError() << "No GPU implementation provided";
#endif
}

PerformanceEstimate Time::metrics() const {
  return PerformanceEstimate::fromKernel<kernel::evaluateDOFSAtTimeSTP<Config>>();
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

} // namespace seissol::kernels::solver::stp
