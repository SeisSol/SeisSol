// SPDX-FileCopyrightText: 2022 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "TimeBase.h"

#ifndef NDEBUG
extern long long libxsmm_num_total_flops;
#endif

#include "Kernels/Common.h"
#include "Kernels/MemoryOps.h"

#include <Eigen/Dense>
#include <cassert>
#include <cstring>
#include <stdint.h>

#ifdef ACL_DEVICE
#include "Common/Offset.h"
#endif

#include "Equations/poroelastic/Model/PoroelasticSetup.h"

#include <yateto.h>

GENERATE_HAS_MEMBER(ET)
GENERATE_HAS_MEMBER(sourceMatrix)

namespace seissol::kernels::solver::stp {

void Spacetime::setGlobalData(const CompoundGlobalData& global) {
  for (std::size_t n = 0; n < Cfg::ConvergenceOrder; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        m_krnlPrototype.kDivMTSub(d, n) = init::kDivMTSub<Cfg>::Values[tensor::kDivMTSub<Cfg>::index(d, n)];
      }
    }
    m_krnlPrototype.selectModes(n) = init::selectModes<Cfg>::Values[tensor::selectModes<Cfg>::index(n)];
  }
  for (std::size_t k = 0; k < seissol::model::MaterialT::NumQuantities; k++) {
    m_krnlPrototype.selectQuantity(k) =
        init::selectQuantity<Cfg>::Values[tensor::selectQuantity<Cfg>::index(k)];
    m_krnlPrototype.selectQuantityG(k) =
        init::selectQuantityG<Cfg>::Values[tensor::selectQuantityG<Cfg>::index(k)];
  }
  m_krnlPrototype.timeInt = init::timeInt<Cfg>::Values;
  m_krnlPrototype.wHat = init::wHat<Cfg>::Values;

#ifdef ACL_DEVICE
  // TODO: adjust pointers
  for (std::size_t n = 0; n < Cfg::ConvergenceOrder; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        deviceKrnlPrototype.kDivMTSub(d, n) =
            init::kDivMTSub<Cfg>::Values[tensor::kDivMTSub<Cfg>::index(d, n)];
      }
    }
    deviceKrnlPrototype.selectModes(n) = init::selectModes<Cfg>::Values[tensor::selectModes<Cfg>::index(n)];
  }
  for (std::size_t k = 0; k < seissol::model::MaterialT::NumQuantities; k++) {
    deviceKrnlPrototype.selectQuantity(k) =
        init::selectQuantity<Cfg>::Values[tensor::selectQuantity<Cfg>::index(k)];
    deviceKrnlPrototype.selectQuantityG(k) =
        init::selectQuantityG<Cfg>::Values[tensor::selectQuantityG<Cfg>::index(k)];
  }
  deviceKrnlPrototype.timeInt = init::timeInt<Cfg>::Values;
  deviceKrnlPrototype.wHat = init::wHat<Cfg>::Values;
#endif
}

void Spacetime::executeSTP(double timeStepWidth,
                           LTS::Ref& data,
                           real timeIntegrated[tensor::I<Cfg>::size()],
                           real* stp)

{
  alignas(PagesizeStack) real stpRhs[tensor::spaceTimePredictorRhs<Cfg>::size()];
  assert((reinterpret_cast<uintptr_t>(stp)) % Alignment == 0);
  std::fill(std::begin(stpRhs), std::end(stpRhs), 0);
  std::fill(stp, stp + tensor::spaceTimePredictor<Cfg>::size(), 0);
  kernel::spaceTimePredictor<Cfg> krnl = m_krnlPrototype;

  // libxsmm can not generate GEMMs with alpha!=1. As a workaround we multiply the
  // star matrices with dt before we execute the kernel.
  real A_values[init::star<Cfg>::size(0)];
  real B_values[init::star<Cfg>::size(1)];
  real C_values[init::star<Cfg>::size(2)];
  for (std::size_t i = 0; i < init::star<Cfg>::size(0); i++) {
    A_values[i] = timeStepWidth * data.get<LTS::LocalIntegration>().starMatrices[0][i];
    B_values[i] = timeStepWidth * data.get<LTS::LocalIntegration>().starMatrices[1][i];
    C_values[i] = timeStepWidth * data.get<LTS::LocalIntegration>().starMatrices[2][i];
  }
  krnl.star(0) = A_values;
  krnl.star(1) = B_values;
  krnl.star(2) = C_values;

  krnl.Gk = data.get<LTS::LocalIntegration>().specific.G[10] * timeStepWidth;
  krnl.Gl = data.get<LTS::LocalIntegration>().specific.G[11] * timeStepWidth;
  krnl.Gm = data.get<LTS::LocalIntegration>().specific.G[12] * timeStepWidth;

  krnl.Q = const_cast<real*>(data.get<LTS::Dofs>());
  krnl.I = timeIntegrated;
  krnl.timestep = timeStepWidth;
  krnl.spaceTimePredictor = stp;
  krnl.spaceTimePredictorRhs = stpRhs;

  // The matrix Zinv depends on the timestep
  // If the timestep is not as expected e.g. when approaching a sync point
  // we have to recalculate it
  if (timeStepWidth != data.get<LTS::LocalIntegration>().specific.typicalTimeStepWidth) {
    auto sourceMatrix =
        init::ET<Cfg>::view::create(data.get<LTS::LocalIntegration>().specific.sourceMatrix);
    real ZinvData[seissol::model::MaterialT::NumQuantities][Cfg::ConvergenceOrder * Cfg::ConvergenceOrder];
    model::zInvInitializerForLoop<0,
                                  seissol::model::MaterialT::NumQuantities,
                                  decltype(sourceMatrix)>(ZinvData, sourceMatrix, timeStepWidth);
    for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
      krnl.Zinv(i) = ZinvData[i];
    }
    // krnl.execute has to be run here: ZinvData is only allocated locally
    krnl.execute();
  } else {
    for (std::size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
      krnl.Zinv(i) = data.get<LTS::LocalIntegration>().specific.Zinv[i];
    }
    krnl.execute();
  }
}

void Spacetime::computeAder(const real* coeffs,
                            double timeStepWidth,
                            LTS::Ref& data,
                            LocalTmp<Cfg>& tmp,
                            real timeIntegrated[tensor::I<Cfg>::size()],
                            real* timeDerivatives,
                            bool updateDisplacement) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.get<LTS::Dofs>())) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0 ||
         timeDerivatives == NULL);

  alignas(Alignment) real temporaryBuffer[tensor::spaceTimePredictor<Cfg>::size()];
  real* stpBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;
  executeSTP(timeStepWidth, data, timeIntegrated, stpBuffer);
}

void Spacetime::flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0;
  hardwareFlops = 0;

  nonZeroFlops = kernel::spaceTimePredictor<Cfg>::NonZeroFlops;
  hardwareFlops = kernel::spaceTimePredictor<Cfg>::HardwareFlops;
  // we multiply the star matrices with dt before we execute the kernel
  nonZeroFlops += 3 * init::star<Cfg>::size(0);
  hardwareFlops += 3 * init::star<Cfg>::size(0);
}

std::uint64_t Spacetime::bytesAder() {
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

  return reals * sizeof(real);
}

void Spacetime::computeBatchedAder(const real* coeffs,
                                   double timeStepWidth,
                                   LocalTmp<Cfg>& tmp,
                                   ConditionalPointersToRealsTable& dataTable,
                                   ConditionalMaterialTable& materialTable,
                                   bool updateDisplacement,
                                   seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  kernel::gpu_spaceTimePredictor<Cfg> krnl = deviceKrnlPrototype;

  ConditionalKey timeVolumeKernelKey(KernelNames::Time || KernelNames::Volume);
  if (dataTable.find(timeVolumeKernelKey) != dataTable.end()) {
    auto& entry = dataTable[timeVolumeKernelKey];

    const auto numElements = (entry.get(inner_keys::Wp::Id::Dofs))->getSize();
    krnl.numElements = numElements;

    krnl.I = (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr();
    krnl.Q = const_cast<const real**>((entry.get(inner_keys::Wp::Id::Dofs))->getDeviceDataPtr());
    krnl.timestep = timeStepWidth;

    // TODO: maybe zero init?
    krnl.spaceTimePredictor = (entry.get(inner_keys::Wp::Id::Stp))->getDeviceDataPtr();
    krnl.spaceTimePredictorRhs = (entry.get(inner_keys::Wp::Id::StpRhs))->getDeviceDataPtr();

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star<Cfg>>(); ++i) {
      krnl.star(i) = const_cast<const real**>(
          (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
      krnl.extraOffset_star(i) = SEISSOL_ARRAY_OFFSET(LocalIntegrationData, starMatrices, i);
    }

    krnl.streamPtr = runtime.stream();

    krnl.Gkt = const_cast<const real**>(
        (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
    krnl.Glt = const_cast<const real**>(
        (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
    krnl.Gmt = const_cast<const real**>(
        (entry.get(inner_keys::Wp::Id::LocalIntegrationData))->getDeviceDataPtr());
    krnl.extraOffset_Gkt = SEISSOL_OFFSET(LocalIntegrationData, specific.G[10]);
    krnl.extraOffset_Glt = SEISSOL_OFFSET(LocalIntegrationData, specific.G[11]);
    krnl.extraOffset_Gmt = SEISSOL_OFFSET(LocalIntegrationData, specific.G[12]);

    /*
    // TODO: port

    if (timeStepWidth != data.localIntegration.specific.typicalTimeStepWidth) {
      assert(false && "NYI");
    }

    */
    /*runtime.enqueueOmpFor(numElements, [](std::size_t i) {
      if (timeStepWidth != data.localIntegration.specific.typicalTimeStepWidth) {
        // TODO
      }
    });*/

    std::size_t zinvOffset = SEISSOL_OFFSET(LocalIntegrationData, specific.Zinv);
    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::Zinv<Cfg>>(); i++) {
      krnl.Zinv(i) =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Zinv))->getDeviceDataPtr());
      krnl.extraOffset_Zinv(i) = zinvOffset;
      zinvOffset += tensor::Zinv<Cfg>::size(i);
    }
    krnl.execute();
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

void Time::evaluate(const real* coeffs,
                    const real* timeDerivatives,
                    real timeEvaluated[tensor::I<Cfg>::size()]) {
  kernel::evaluateDOFSAtTimeSTP<Cfg> krnl;
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
  logError() << "No GPU implementation provided";
}

void Time::flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) {
  nonZeroFlops = kernel::evaluateDOFSAtTimeSTP<Cfg>::NonZeroFlops;
  hardwareFlops = kernel::evaluateDOFSAtTimeSTP<Cfg>::HardwareFlops;
}

void Time::setGlobalData(const CompoundGlobalData& global) {}

} // namespace seissol::kernels::solver::stp
