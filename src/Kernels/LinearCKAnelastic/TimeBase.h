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

#include "generated_code/kernel.h"
#include <Kernels/Spacetime.h>
#include <Kernels/Time.h>

namespace seissol::kernels::solver::linearckanelastic {
class Spacetime : public SpacetimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void computeAder(double timeStepWidth,
                   LTS::Ref& data,
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

  void flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) override;

  std::uint64_t bytesAder() override;

  protected:
  kernel::derivative m_krnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_derivative deviceKrnlPrototype;
#endif
};

class Time : public TimeKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;
  void evaluateAtTime(
      std::shared_ptr<basisFunction::SampledTimeBasisFunctions<real>> evaluatedTimeBasisFunctions,
      const real* timeDerivatives,
      real timeEvaluated[tensor::Q::size()]) override;
  void flopsEvaluateAtTime(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) override;

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

  void flopsTaylorExpansion(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) override;
};
} // namespace seissol::kernels::solver::linearckanelastic

#endif // SEISSOL_SRC_KERNELS_LINEARCKANELASTIC_TIMEBASE_H_
