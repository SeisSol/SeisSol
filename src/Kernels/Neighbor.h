// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2014-2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_KERNELS_NEIGHBOR_H_
#define SEISSOL_SRC_KERNELS_NEIGHBOR_H_

#include "Initializer/Typedefs.h"
#include "Kernels/Interface.h"
#include "Kernels/Kernel.h"
#include "Monitoring/Metric.h"
#include "Parallel/Runtime/Stream.h"

namespace seissol::kernels {

class NeighborKernel : public Kernel {
  public:
  ~NeighborKernel() override = default;

  /**
   * @brief Apply the corrector of every face of the cell to its two traces, and the flux the
   * fault supplies on its faces.
   *
   * This step equals the second part of the ADER-DG "corrector": the corrector of a face reads the
   * trace the cell wrote on it after its predictor and the one the neighbor wrote there, applies
   * the flux to the two and lifts the result into the cell; a dynamic rupture face adds the flux of
   * the fault. The neighbor's trace has to be over the time step of the cell: the neighbor wrote it
   * so (same time step, or accumulated over its smaller steps), or it was evaluated from the traces
   * of the neighbor's time derivatives via the TimeKernel (in the case of LTS, if our neighbor is a
   * cell from a larger time cluster).
   *
   * @param data Cell data reference object (contains references to all stored data arrays for that
   * cell)
   * @param neighborTraces The trace of the neighbor on each face, over the time step of the cell;
   * a trace of zeros on a face without neighbor (LTS::ZeroFaceTrace)
   */
  virtual void
      computeNeighborsIntegral(LTS::Ref& data,
                               const std::array<const real*, Cell::NumFaces>& neighborTraces) = 0;

  virtual void
      computeBatchedNeighborsIntegral(recording::ConditionalPointersToRealsTable& table,
                                      seissol::parallel::runtime::StreamRuntime& runtime) = 0;

  [[nodiscard]] virtual std::pair<PerformanceEstimate, PerformanceEstimate>
      metrics(const std::array<FaceType, Cell::NumFaces>& faceTypes,
              const std::array<CellDRMapping, Cell::NumFaces>& cellDrMapping) const = 0;

  /// What a weight-adjusted cell adds to the neighbour integral (see WeightAdjustable): nothing,
  /// where the solver adjusts no cell.
  [[nodiscard]] virtual PerformanceEstimate weightAdjustMetrics() const { return {}; }
};

} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_NEIGHBOR_H_
