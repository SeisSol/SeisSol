// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_TIMEBASE_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_TIMEBASE_H_

#include "Common/Constants.h"
#include "generated_code/kernel.h"
#include <Initializer/Typedefs.h>
#include <Kernels/Spacetime.h>
#include <Kernels/Time.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif // ACL_DEVICE

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels::solver::linearck {

class Spacetime : public SpacetimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void computeAder(double timeStepWidth,
                   LocalData& data,
                   LocalTmp& tmp,
                   real timeIntegrated[tensor::I::size()],
                   real* timeDerivativesOrSTP = nullptr,
                   bool updateDisplacement = false) override;
  void computeBatchedAder(double timeStepWidth,
                          LocalTmp& tmp,
                          ConditionalPointersToRealsTable& dataTable,
                          ConditionalMaterialTable& materialTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsAder(unsigned int& nonZeroFlops, unsigned int& hardwareFlops) override;

  unsigned bytesAder() override;

  protected:
  kernel::derivative m_krnlPrototype;
  kernel::projectDerivativeToNodalBoundaryRotated projectDerivativeToNodalBoundaryRotated;

#ifdef ACL_DEVICE
  kernel::gpu_derivative deviceKrnlPrototype;
  kernel::gpu_projectDerivativeToNodalBoundaryRotated deviceDerivativeToNodalBoundaryRotated;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

class Time : public TimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void evaluateAtTime(
      std::shared_ptr<basisFunction::SampledTimeBasisFunctions<real>> evaluatedTimeBasisFunctions,
      const real* timeDerivatives,
      real timeEvaluated[tensor::Q::size()]) override;
  void flopsEvaluateAtTime(long long& nonZeroFlops, long long& hardwareFlops) override;

  void computeIntegral(double expansionPoint,
                       double integrationStart,
                       double integrationEnd,
                       const real* timeDerivatives,
                       real timeIntegrated[tensor::I::size()]) override;

  void computeBatchedIntegral(double expansionPoint,
                              double integrationStart,
                              double integrationEnd,
                              const real** timeDerivatives,
                              real** timeIntegratedDofs,
                              unsigned numElements,
                              seissol::parallel::runtime::StreamRuntime& runtime) override;

  void computeTaylorExpansion(real time,
                              real expansionPoint,
                              const real* timeDerivatives,
                              real timeEvaluated[tensor::Q::size()]) override;

  void computeBatchedTaylorExpansion(real time,
                                     real expansionPoint,
                                     real** timeDerivatives,
                                     real** timeEvaluated,
                                     size_t numElements,
                                     seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops) override;

  protected:
  static void checkGlobalData(const GlobalData* global, size_t alignment);

#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_TIMEBASE_H_
