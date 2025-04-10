// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_NEIGHBORBASE_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_NEIGHBORBASE_H_

#include "Common/Constants.h"
#include "generated_code/kernel.h"
#include <Kernels/Neighbor.h>
#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {
struct GlobalData;
} // namespace seissol

namespace seissol::kernels::solver::linearck {

class Neighbor : public NeighborKernel {
  public:
  void setGlobalData(const CompoundGlobalData& global) override;

  void computeNeighborsIntegral(NeighborData& data,
                                const CellDRMapping (&cellDrMapping)[4],
                                real* timeIntegrated[4],
                                real* faceNeighborsPrefetch[4]) override;

  void computeBatchedNeighborsIntegral(ConditionalPointersToRealsTable& table,
                                       seissol::parallel::runtime::StreamRuntime& runtime) override;

  void flopsNeighborsIntegral(const FaceType faceTypes[4],
                              const int neighboringIndices[4][2],
                              const CellDRMapping (&cellDrMapping)[4],
                              unsigned int& nonZeroFlops,
                              unsigned int& hardwareFlops,
                              long long& drNonZeroFlops,
                              long long& drHardwareFlops) override;

  unsigned bytesNeighborsIntegral() override;

  protected:
  kernel::neighboringFlux m_nfKrnlPrototype;
  dynamicRupture::kernel::nodalFlux m_drKrnlPrototype;

#ifdef ACL_DEVICE
  kernel::gpu_neighboringFlux deviceNfKrnlPrototype;
  dynamicRupture::kernel::gpu_nodalFlux deviceDrKrnlPrototype;
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif
};

} // namespace seissol::kernels::solver::linearck

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_NEIGHBORBASE_H_
