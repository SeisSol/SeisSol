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
#include "Kernels/Kernel.h"
#include "Numerical/BasisFunction.h"
#include "Parallel/Runtime/Stream.h"

#include <cassert>
#include <limits>
#include <memory>

namespace seissol::kernels {

class SpacetimeKernel : public Kernel {
  public:
  ~SpacetimeKernel() override = default;

  /**
    @brief Compute the space-time evolution (a.k.a. "derivatives") from the given DOFs.

    Includes an inline time integration; as that's needed for most subsequent operations.

    @param coeffs The time basis coefficients used for the integrated integration.
    @param timeStepWidth The size of the current timestep
    @param timeIntegrated Output: time integration data.
    @param timeDerivativesOrSTP Output: space-time evoluion.
  */
  virtual void computeAder(const real* coeffs,
                           double timeStepWidth,
                           LTS::Ref& data,
                           LocalTmp& tmp,
                           real timeIntegrated[tensor::I::size()],
                           real* timeDerivativesOrSTP = nullptr,
                           bool updateDisplacement = false) = 0;

  virtual void computeBatchedAder(const real* coeffs,
                                  double timeStepWidth,
                                  recording::ConditionalPointersToRealsTable& dataTable,
                                  recording::ConditionalMaterialTable& materialTable,
                                  seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsAder(std::uint64_t& nonZeroFlops, std::uint64_t& hardwareFlops) = 0;

  virtual std::uint64_t bytesAder() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_SPACETIME_H_
