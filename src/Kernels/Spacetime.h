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

#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Numerical/BasisFunction.h"
#include "generated_code/tensor.h"
#include <Kernels/Kernel.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <limits>
#include <memory>

namespace seissol::kernels {

class SpacetimeKernel : public Kernel {
  public:
  ~SpacetimeKernel() override = default;
  virtual void computeAder(double timeStepWidth,
                           LocalData& data,
                           LocalTmp& tmp,
                           real timeIntegrated[tensor::I::size()],
                           real* timeDerivativesOrSTP = nullptr,
                           bool updateDisplacement = false) = 0;
  virtual void computeBatchedAder(double timeStepWidth,
                                  LocalTmp& tmp,
                                  ConditionalPointersToRealsTable& dataTable,
                                  ConditionalMaterialTable& materialTable,
                                  bool updateDisplacement,
                                  seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsAder(unsigned int& nonZeroFlops, unsigned int& hardwareFlops) = 0;

  virtual unsigned bytesAder() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_SPACETIME_H_
