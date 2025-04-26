// SPDX-FileCopyrightText: 2013 SeisSol Group
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
#include "Parallel/Runtime/Stream.h"
#include "generated_code/tensor.h"
#include <Kernels/Kernel.h>
#include <Physics/InitialField.h>
#include <cassert>

namespace seissol::kernels {

class LocalKernel : public Kernel {
  protected:
  double gravitationalAcceleration{9.81};
  const std::vector<std::unique_ptr<physics::InitialField>>* initConds;

  public:
  ~LocalKernel() override = default;
  void setGravitationalAcceleration(double g) { gravitationalAcceleration = g; }
  void setInitConds(decltype(initConds) initConds) { this->initConds = initConds; }

  physics::InitialField* getInitCond(size_t index) {
    const auto& condition = this->initConds->at(index);
    return condition.get();
  }

  virtual void computeIntegral(real timeIntegratedDegreesOfFreedom[tensor::I::size()],
                               LocalData& data,
                               LocalTmp& tmp,
                               const CellMaterialData* materialData,
                               const CellBoundaryMapping (*cellBoundaryMapping)[4],
                               double time,
                               double timeStepWidth) = 0;

  virtual void computeBatchedIntegral(ConditionalPointersToRealsTable& dataTable,
                                      ConditionalMaterialTable& materialTable,
                                      ConditionalIndicesTable& indicesTable,
                                      kernels::LocalData::Loader& loader,
                                      LocalTmp& tmp,
                                      double timeStepWidth,
                                      seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void
      evaluateBatchedTimeDependentBc(ConditionalPointersToRealsTable& dataTable,
                                     ConditionalIndicesTable& indicesTable,
                                     kernels::LocalData::Loader& loader,
                                     seissol::initializer::Layer& layer,
                                     seissol::initializer::LTS& lts,
                                     double time,
                                     double timeStepWidth,
                                     seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void flopsIntegral(const FaceType faceTypes[4],
                             unsigned int& nonZeroFlops,
                             unsigned int& hardwareFlops) = 0;

  virtual unsigned bytesIntegral() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_LOCAL_H_
