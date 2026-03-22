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
#include "Parallel/Runtime/Stream.h"
#include "Physics/InitialField.h"

#include <cassert>

namespace seissol::kernels {

class LocalKernel : public Kernel {
  protected:
  double gravitationalAcceleration{9.81};
  const std::vector<std::unique_ptr<physics::InitialField>>* initConds{nullptr};

  public:
  ~LocalKernel() override = default;
  void setGravitationalAcceleration(double g) { gravitationalAcceleration = g; }
  void setInitConds(decltype(initConds) initConds) { this->initConds = initConds; }

  physics::InitialField* getInitCond(size_t index) {
    const auto& condition = this->initConds->at(index);
    return condition.get();
  }

  /**
   * @brief Compute the volume integral over the ADER space-time evolution, as well as local
   * boundary conditions as fluxes.
   *
   * This step equals the first part of the ADER-DG "corrector"; that is, we compute the local
   * space-time volume integral, and add the flux contributions that depend on the local element.
   * Also, we compute some boundary conditions that can be computed by only using the local cell
   * data. The dynamic rupture faces are usually ignored in this step.
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

  virtual void computeBatchedIntegral(recording::ConditionalPointersToRealsTable& dataTable,
                                      recording::ConditionalMaterialTable& materialTable,
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

  virtual void flopsIntegral(const std::array<FaceType, Cell::NumFaces>& faceTypes,
                             std::uint64_t& nonZeroFlops,
                             std::uint64_t& hardwareFlops) = 0;

  virtual std::uint64_t bytesIntegral() = 0;
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_LOCAL_H_
