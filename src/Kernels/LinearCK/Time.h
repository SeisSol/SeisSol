// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_TIME_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_TIME_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Spacetime.h"
#include "Kernels/Time.h"

#ifdef ACL_DEVICE
#include <Device/device.h>
#endif // ACL_DEVICE

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels::solver::linearck {

class Spacetime : public SpacetimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void computeAder(const real* coeffs,
                   double timeStepWidth,
                   LTS::Ref& data,
                   LocalTmp& tmp,
                   real timeIntegrated[tensor::I::size()],
                   real* timeDerivativesOrSTP = nullptr,
                   bool updateDisplacement = false) override;
  void computeBatchedAder(const real* coeffs,
                          double timeStepWidth,
                          LTS::Layer& layer,
                          LocalTmp& tmp,
                          recording::ConditionalPointersToRealsTable& dataTable,
                          recording::ConditionalMaterialTable& materialTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) override;

  std::uint64_t bytesAder() override;

  protected:
  kernel::derivative m_krnlPrototype;
  kernel::projectDerivativeToNodalBoundaryRotated projectDerivativeToNodalBoundaryRotated;

  kernel::fsgKernel fsgKernelPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_derivative deviceKrnlPrototype;
  kernel::gpu_projectDerivativeToNodalBoundaryRotated deviceDerivativeToNodalBoundaryRotated;
  kernel::gpu_fsgKernel deviceFsgKernelPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

class Time : public TimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void evaluate(const real* coeffs,
                const real* timeDerivatives,
                real timeEvaluated[tensor::I::size()]) override;
  void evaluateBatched(const real* coeffs,
                       const real** timeDerivatives,
                       real** timeIntegratedDofs,
                       std::size_t numElements,
                       seissol::parallel::runtime::StreamRuntime& runtime) override;
  void flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) override;
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_KERNELS_LINEARCK_TIME_H_
