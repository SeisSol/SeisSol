// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_TIMECOMMON_H_
#define SEISSOL_SRC_KERNELS_TIMECOMMON_H_

#include "GeneratedCode/tensor.h"
#include "Initializer/LtsSetup.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Solver.h"

namespace seissol::kernels {
struct TimeCommon {
  /**
   * The trace of each neighbor on its face over the time step of the cell, which the corrector of
   * the face reads: what the neighbor wrote there where that is a trace over the time step already
   * (LtsSetup::neighborBuffer is StepIntegrals or AccumulatedIntegrals), or the traces of its time
   * derivatives evaluated over the interval of the cell (Derivatives), which a neighbor writes that
   * takes a larger time step; the zero trace of the layer on a face without neighbor, whose
   * corrector applies the cell's own trace alone.
   *
   * @param ltsSetup the LTS setup of the cell.
   * @param faceTypes face types of the cell.
   * @param subtimeCoeffs the time integration coefficients of the interval of the cell within the
   * time step of a neighbor that takes a larger one
   * @param faceNeighbors what the four neighbors wrote on the faces (LTS::FaceNeighbors).
   * @param zeroTrace a trace of zeros (LTS::ZeroFaceTrace).
   * @param traceBuffers memory where an evaluated trace goes. Ensure thread safety!
   * @param neighborTraces the four traces over the time step of the cell (either in faceNeighbors,
   * in traceBuffers, or zeroTrace).
   **/
  static void computeNeighborTraces(Time& time,
                                    const LtsSetup& ltsSetup,
                                    const std::array<FaceType, Cell::NumFaces>& faceTypes,
                                    const real* subtimeCoeffs,
                                    const std::array<real*, Cell::NumFaces>& faceNeighbors,
                                    const real* zeroTrace,
                                    const std::array<real*, Cell::NumFaces>& traceBuffers,
                                    std::array<const real*, Cell::NumFaces>& neighborTraces);

  /**
   * The traces computeNeighborTraces evaluates, for the batches the recorder noted (see
   * recording::NeighIntegrationRecorder); the others the corrector reads where they are.
   **/
  static void computeBatchedNeighborTraces(Time& time,
                                           const real* subtimeCoeffs,
                                           recording::ConditionalPointersToRealsTable& table,
                                           seissol::parallel::runtime::StreamRuntime& runtime);

  TimeCommon() = delete;
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIMECOMMON_H_
