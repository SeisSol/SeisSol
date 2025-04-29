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
#include "generated_code/kernel.h"
#include <Kernels/Spacetime.h>
#include <Kernels/Time.h>

#ifdef ACL_DEVICE
#include <device.h>
#endif // ACL_DEVICE

namespace seissol::kernels::solver::stp {

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

  private:
  void executeSTP(double timeStepWidth,
                  LocalData& data,
                  real timeIntegrated[tensor::I::size()],
                  real* stp);

  kernel::spaceTimePredictor m_krnlPrototype;
  kernel::projectDerivativeToNodalBoundaryRotated projectDerivativeToNodalBoundaryRotated;

#ifdef ACL_DEVICE
  kernel::gpu_spaceTimePredictor deviceKrnlPrototype;
  kernel::gpu_projectDerivativeToNodalBoundaryRotated deviceDerivativeToNodalBoundaryRotated;
#endif
};

} // namespace seissol::kernels::solver::stp

#endif // SEISSOL_SRC_KERNELS_STP_TIMEBASE_H_
