// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_TIME_H_
#define SEISSOL_SRC_KERNELS_TIME_H_

#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/TimeBase.h"
#include "generated_code/tensor.h"
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <limits>
#ifdef USE_STP
#include "Numerical/BasisFunction.h"
#include <memory>
#endif

namespace seissol::kernels {

class Time : public TimeBase {
  public:
  void setHostGlobalData(const GlobalData* global);
  void setGlobalData(const CompoundGlobalData& global);

  void computeAder(double timeStepWidth,
                   LocalData& data,
                   LocalTmp& tmp,
                   real timeIntegrated[tensor::I::size()],
                   real* timeDerivatives = nullptr,
                   bool updateDisplacement = false);

#ifdef USE_STP
  void executeSTP(double timeStepWidth,
                  LocalData& data,
                  real timeIntegrated[tensor::I::size()],
                  real* stp);
  void evaluateAtTime(
      std::shared_ptr<basisFunction::SampledTimeBasisFunctions<real>> evaluatedTimeBasisFunctions,
      const real* timeDerivatives,
      real timeEvaluated[tensor::Q::size()]);
  void flopsEvaluateAtTime(long long& nonZeroFlops, long long& hardwareFlops);

#endif
  void computeBatchedAder(double timeStepWidth,
                          LocalTmp& tmp,
                          ConditionalPointersToRealsTable& dataTable,
                          ConditionalMaterialTable& materialTable,
                          bool updateDisplacement,
                          seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsAder(unsigned int& nonZeroFlops, unsigned int& hardwareFlops);

  unsigned bytesAder();

  void computeIntegral(double expansionPoint,
                       double integrationStart,
                       double integrationEnd,
                       const real* timeDerivatives,
                       real timeIntegrated[tensor::I::size()]);

  void computeBatchedIntegral(double expansionPoint,
                              double integrationStart,
                              double integrationEnd,
                              const real** timeDerivatives,
                              real** timeIntegratedDofs,
                              unsigned numElements,
                              seissol::parallel::runtime::StreamRuntime& runtime);

  void computeTaylorExpansion(real time,
                              real expansionPoint,
                              const real* timeDerivatives,
                              real timeEvaluated[tensor::Q::size()]);

  void computeBatchedTaylorExpansion(real time,
                                     real expansionPoint,
                                     real** timeDerivatives,
                                     real** timeEvaluated,
                                     size_t numElements,
                                     seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsTaylorExpansion(long long& nonZeroFlops, long long& hardwareFlops);

  unsigned int* getDerivativesOffsets();
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIME_H_
