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
#include "Kernels/DenseMatrixOps.h"

#include <Eigen/Dense>
#include <cassert>
#include <cstring>
#include <omp.h>
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
  for (std::size_t n = 0; n < ConvergenceOrder; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        m_krnlPrototype.kDivMTSub(d, n) = init::kDivMTSub::Values[tensor::kDivMTSub::index(d, n)];
      }
    }
    m_krnlPrototype.selectModes(n) = init::selectModes::Values[tensor::selectModes::index(n)];
  }
  for (std::size_t k = 0; k < seissol::model::MaterialT::NumQuantities; k++) {
    m_krnlPrototype.selectQuantity(k) =
        init::selectQuantity::Values[tensor::selectQuantity::index(k)];
    m_krnlPrototype.selectQuantityG(k) =
        init::selectQuantityG::Values[tensor::selectQuantityG::index(k)];
  }
  m_krnlPrototype.timeInt = init::timeInt::Values;
  m_krnlPrototype.wHat = init::wHat::Values;

#ifdef ACL_DEVICE
  // TODO: adjust pointers
  for (std::size_t n = 0; n < ConvergenceOrder; ++n) {
    if (n > 0) {
      for (int d = 0; d < 3; ++d) {
        deviceKrnlPrototype.kDivMTSub(d, n) =
            init::kDivMTSub::Values[tensor::kDivMTSub::index(d, n)];
      }
    }
    deviceKrnlPrototype.selectModes(n) = init::selectModes::Values[tensor::selectModes::index(n)];
  }
  for (std::size_t k = 0; k < seissol::model::MaterialT::NumQuantities; k++) {
    deviceKrnlPrototype.selectQuantity(k) =
        init::selectQuantity::Values[tensor::selectQuantity::index(k)];
    deviceKrnlPrototype.selectQuantityG(k) =
        init::selectQuantityG::Values[tensor::selectQuantityG::index(k)];
  }
  deviceKrnlPrototype.timeInt = init::timeInt::Values;
  deviceKrnlPrototype.wHat = init::wHat::Values;
#endif
}

void Spacetime::executeSTP(double timeStepWidth,
                           LocalData& data,
                           real timeIntegrated[tensor::I::size()],
                           real* stp)

{
  alignas(PagesizeStack) real stpRhs[tensor::spaceTimePredictorRhs::size()];
  assert((reinterpret_cast<uintptr_t>(stp)) % Alignment == 0);
  std::fill(std::begin(stpRhs), std::end(stpRhs), 0);
  std::fill(stp, stp + tensor::spaceTimePredictor::size(), 0);
  kernel::spaceTimePredictor krnl = m_krnlPrototype;

  // libxsmm can not generate GEMMs with alpha!=1. As a workaround we multiply the
  // star matrices with dt before we execute the kernel.
  real A_values[init::star::size(0)];
  real B_values[init::star::size(1)];
  real C_values[init::star::size(2)];
  for (std::size_t i = 0; i < init::star::size(0); i++) {
    A_values[i] = timeStepWidth * data.localIntegration().starMatrices[0][i];
    B_values[i] = timeStepWidth * data.localIntegration().starMatrices[1][i];
    C_values[i] = timeStepWidth * data.localIntegration().starMatrices[2][i];
  }
  krnl.star(0) = A_values;
  krnl.star(1) = B_values;
  krnl.star(2) = C_values;

  krnl.Gk = data.localIntegration().specific.G[10] * timeStepWidth;
  krnl.Gl = data.localIntegration().specific.G[11] * timeStepWidth;
  krnl.Gm = data.localIntegration().specific.G[12] * timeStepWidth;

  krnl.Q = const_cast<real*>(data.dofs());
  krnl.I = timeIntegrated;
  krnl.timestep = timeStepWidth;
  krnl.spaceTimePredictor = stp;
  krnl.spaceTimePredictorRhs = stpRhs;

  // The matrix Zinv depends on the timestep
  // If the timestep is not as expected e.g. when approaching a sync point
  // we have to recalculate it
  if (timeStepWidth != data.localIntegration().specific.typicalTimeStepWidth) {
    auto sourceMatrix = init::ET::view::create(data.localIntegration().specific.sourceMatrix);
    real ZinvData[seissol::model::MaterialT::NumQuantities][ConvergenceOrder * ConvergenceOrder];
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
      krnl.Zinv(i) = data.localIntegration().specific.Zinv[i];
    }
    krnl.execute();
  }
}

void Spacetime::computeAder(double timeStepWidth,
                            LocalData& data,
                            LocalTmp& tmp,
                            real timeIntegrated[tensor::I::size()],
                            real* timeDerivatives,
                            bool updateDisplacement) {
  /*
   * assert alignments.
   */
  assert((reinterpret_cast<uintptr_t>(data.dofs())) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeIntegrated)) % Alignment == 0);
  assert((reinterpret_cast<uintptr_t>(timeDerivatives)) % Alignment == 0 ||
         timeDerivatives == NULL);

  alignas(Alignment) real temporaryBuffer[tensor::spaceTimePredictor::size()];
  real* stpBuffer = (timeDerivatives != nullptr) ? timeDerivatives : temporaryBuffer;
  executeSTP(timeStepWidth, data, timeIntegrated, stpBuffer);
}

void Spacetime::flopsAder(unsigned int& nonZeroFlops, unsigned int& hardwareFlops) {
  // reset flops
  nonZeroFlops = 0;
  hardwareFlops = 0;

  nonZeroFlops = kernel::spaceTimePredictor::NonZeroFlops;
  hardwareFlops = kernel::spaceTimePredictor::HardwareFlops;
  // we multiply the star matrices with dt before we execute the kernel
  nonZeroFlops += 3 * init::star::size(0);
  hardwareFlops += 3 * init::star::size(0);
}

unsigned Spacetime::bytesAder() {
  unsigned reals = 0;

  // DOFs load, tDOFs load, tDOFs write
  reals += tensor::Q::size() + 2 * tensor::I::size();
  // star matrices, source matrix
  reals += yateto::computeFamilySize<tensor::star>();
  // Zinv
  reals += yateto::computeFamilySize<tensor::Zinv>();
  // G
  reals += 3;

  /// \todo incorporate derivatives

  return reals * sizeof(real);
}

void Spacetime::computeBatchedAder(double timeStepWidth,
                                   LocalTmp& tmp,
                                   ConditionalPointersToRealsTable& dataTable,
                                   ConditionalMaterialTable& materialTable,
                                   bool updateDisplacement,
                                   seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  kernel::gpu_spaceTimePredictor krnl = deviceKrnlPrototype;

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

    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::star>(); ++i) {
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
    for (size_t i = 0; i < yateto::numFamilyMembers<tensor::Zinv>(); i++) {
      krnl.Zinv(i) =
          const_cast<const real**>((entry.get(inner_keys::Wp::Id::Zinv))->getDeviceDataPtr());
      krnl.extraOffset_Zinv(i) = zinvOffset;
      zinvOffset += tensor::Zinv::size(i);
    }
    krnl.execute();
  }
#else
  assert(false && "no implementation provided");
#endif
}

} // namespace seissol::kernels::solver::stp
