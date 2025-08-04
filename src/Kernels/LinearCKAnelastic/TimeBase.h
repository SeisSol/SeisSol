// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIMEBASE_H_
#define SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIMEBASE_H_

#include "GeneratedCode/kernel.h"
#include <Kernels/Spacetime.h>
#include <Kernels/Time.h>

namespace seissol::kernels::solver::linearckanelastic {
class Spacetime : public SpacetimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
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

  protected:
  kernel::derivative<Cfg> m_krnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_derivative<Cfg> deviceKrnlPrototype;
#endif
};

class Time : public TimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
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
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIMEBASE_H_
