// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_KERNELS_LOCAL_H_
#define SEISSOL_SRC_KERNELS_LOCAL_H_

#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/LocalBase.h"
#include "Parallel/Runtime/Stream.h"
#include "generated_code/tensor.h"
#include <cassert>

namespace seissol::kernels {

class Local : public LocalBase {
  public:
  void setHostGlobalData(const GlobalData* global);
  void setGlobalData(const CompoundGlobalData& global);

  void computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I::size()],
                       LocalData& data,
                       LocalTmp& tmp,
                       const CellMaterialData* materialData,
                       const CellBoundaryMapping (*cellBoundaryMapping)[4],
                       double time,
                       double timeStepWidth);

  void computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                              ConditionalMaterialTable& materialTable,
                              ConditionalIndicesTable& indicesTable,
                              kernels::LocalData::Loader& loader,
                              LocalTmp& tmp,
                              double timeStepWidth,
                              seissol::parallel::runtime::StreamRuntime& runtime);

  void evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                      ConditionalIndicesTable& indicesTable,
                                      kernels::LocalData::Loader& loader,
                                      seissol::initializer::Layer& layer,
                                      seissol::initializer::LTS& lts,
                                      double time,
                                      double timeStepWidth,
                                      seissol::parallel::runtime::StreamRuntime& runtime);

  void flopsIntegral(const FaceType faceTypes[4],
                     unsigned int& nonZeroFlops,
                     unsigned int& hardwareFlops);

  unsigned bytesIntegral();
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_LOCAL_H_
