// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_SPACETIME_H_
#define SEISSOL_SRC_KERNELS_SPACETIME_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Numerical/BasisFunction.h"
#include <Kernels/Kernel.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <limits>
#include <memory>

namespace seissol::kernels {

template <typename Cfg>
class SpacetimeKernel : public Kernel {
  public:
  using real = Real<Cfg>;

  ~SpacetimeKernel() override = default;
  virtual void computeAder(const real* coeffs,
                           double timeStepWidth,
                           LTS::Ref<Cfg>& data,
                           LocalTmp<Cfg>& tmp,
                           real timeIntegrated[tensor::I<Cfg>::size()],
                           real* timeDerivativesOrSTP = nullptr,
                           bool updateDisplacement = false) = 0;
  virtual void computeBatchedAder(const real* coeffs,
                                  double timeStepWidth,
                                  LocalTmp<Cfg>& tmp,
                                  ConditionalPointersToRealsTable& dataTable,
                                  ConditionalMaterialTable& materialTable,
                                  bool updateDisplacement,
                                  seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) = 0;

  virtual std::uint64_t bytesAder() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_SPACETIME_H_
