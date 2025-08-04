// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_STP_TIMEBASE_H_
#define SEISSOL_SRC_KERNELS_STP_TIMEBASE_H_

#include "Common/Constants.h"
#include "GeneratedCode/kernel.h"
#include <Kernels/Spacetime.h>
#include <Kernels/Time.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif // ACL_DEVICE

namespace seissol::kernels::solver::stp {

template <typename Cfg>
class Spacetime : public SpacetimeKernel<Cfg> {
  public:
  void setGlobalData(const GlobalData& global) override;
  void computeAder(const real* coeffs,
                   double timeStepWidth,
                   LTS::Ref<Cfg>& data,
                   LocalTmp<Cfg>& tmp,
                   real timeIntegrated[tensor::I<Cfg>::size()],
                   real* timeDerivativesOrSTP = nullptr,
                   bool updateDisplacement = false) override;
  void computeBatchedAder(const real* coeffs,
                          double timeStepWidth,
                          LocalTmp<Cfg>& tmp,
                          ConditionalPointersToRealsTable& dataTable,
                          ConditionalMaterialTable& materialTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) override;

  std::uint64_t bytesAder() override;

  private:
  void executeSTP(double timeStepWidth,
                  LTS::Ref<Cfg>& data,
                  real timeIntegrated[tensor::I<Cfg>::size()],
                  real* stp);

  kernel::spaceTimePredictor<Cfg> m_krnlPrototype;
  kernel::projectDerivativeToNodalBoundaryRotated<Cfg> projectDerivativeToNodalBoundaryRotated;

#ifdef ACL_DEVICE
  kernel::gpu_spaceTimePredictor<Cfg> deviceKrnlPrototype;
  kernel::gpu_projectDerivativeToNodalBoundaryRotated<Cfg> deviceDerivativeToNodalBoundaryRotated;
#endif
};

template <typename Cfg>
class Time : public TimeKernel<Cfg> {
  public:
  void setGlobalData(const GlobalData& global) override;
  void evaluate(const real* coeffs,
                const real* timeDerivatives,
                real timeEvaluated[tensor::I<Cfg>::size()]) override;
  void evaluateBatched(const real* coeffs,
                       const real** timeDerivatives,
                       real** timeIntegratedDofs,
                       std::size_t numElements,
                       seissol::parallel::runtime::StreamRuntime& runtime) override;
  void flopsEvaluate(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) override;
};

} // namespace seissol::kernels::solver::stp

#endif // SEISSOL_SRC_KERNELS_STP_TIMEBASE_H_
