// SPDX-FileCopyrightText: 2013 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer

#ifndef SEISSOL_SRC_KERNELS_LOCAL_H_
#define SEISSOL_SRC_KERNELS_LOCAL_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/Interface.h"
#include "Kernels/Kernel.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"
#include "Physics/InitialField.h"

#include <cassert>

namespace seissol::kernels {

class LocalKernel : public Kernel {
  protected:
  double gravitationalAcceleration_{9.81};
  const std::vector<std::unique_ptr<physics::InitialField>>* initConds_{nullptr};

  public:
  ~LocalKernel() override = default;
  void setGravitationalAcceleration(double g) { gravitationalAcceleration_ = g; }
  void setInitConds(decltype(initConds_) initConds) { this->initConds_ = initConds; }

  physics::InitialField* getInitCond(size_t index) {
    const auto& condition = this->initConds_->at(index);
    return condition.get();
  }

  /**
   * @brief Compute the volume integral over the ADER space-time evolution, write the trace of the
   * time integral on every face, and compute the boundary conditions that only need the cell
   * itself.
   *
   * This step equals the first part of the ADER-DG "corrector"; that is, we compute the local
   * space-time volume integral, and write what the corrector of every face reads of the cell: the
   * trace of the time integral on the face (LTS::FaceTraces). The corrector itself runs in the
   * neighbor step, where the trace of the neighbor is there as well. Also, we compute some boundary
   * conditions that can be computed by only using the local cell data.
   *
   * @param timeIntegratedDoFs The time integral over the space-time predictor, as computed by the
   * SpacetimeKernel or the TimeKernel.
   * @param data Cell data reference object (contains references to all stored data arrays for that
   * cell)
   * @param tmp Cell-local temporary data from the SpacetimeKernel
   * @param time The current time in the simulation (needed e.g. for time dependent boundary
   * conditions)
   * @param timeStepWidth The current time step width
   */
  virtual void computeIntegral(real* timeIntegratedDoFs,
                               LTS::Ref& data,
                               LocalTmp& tmp,
                               double time,
                               double timeStepWidth) = 0;

  /**
   * @brief computeIntegral for the batches the recorders noted, without the volume integral: on the
   * device, SpacetimeKernel::computeBatchedAder has applied that already.
   */
  virtual void computeBatchedIntegral(recording::ConditionalPointersToRealsTable& dataTable,
                                      recording::ConditionalIndicesTable& indicesTable,
                                      double timeStepWidth,
                                      seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  virtual void
      evaluateBatchedTimeDependentBc(recording::ConditionalPointersToRealsTable& dataTable,
                                     recording::ConditionalIndicesTable& indicesTable,
                                     LTS::Layer& layer,
                                     double time,
                                     double timeStepWidth,
                                     seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  [[nodiscard]] virtual PerformanceEstimate
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes) const = 0;

  /// What a curved cell adds to the local integral (see CurvedSplit): nothing, where no cell may
  /// be curved.
  [[nodiscard]] virtual PerformanceEstimate
      curvedMetrics(const std::array<FaceType, Cell::NumFaces>& /*faceTypes*/) const {
    return {};
  }
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_LOCAL_H_
