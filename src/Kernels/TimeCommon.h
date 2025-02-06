// SPDX-FileCopyrightText: 2016-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_KERNELS_TIMECOMMON_H_
#define SEISSOL_SRC_KERNELS_TIMECOMMON_H_

#include "Initializer/Typedefs.h"
#include "Kernels/Time.h"
#include "generated_code/tensor.h"

namespace seissol::kernels {
struct TimeCommon {
  /**
   * Either copies pointers to the DOFs in the time buffer or integrates the DOFs via time
   derivatives.
   *   Evaluation depends on bit 0-3  of the LTS setup.
   *   0 -> copy buffer; 1 -> integrate via time derivatives
   *     Example:

   *     [     4 unused     | copy or int bits  ]
   *     [ -    -    -    - |  0    1    1    0 ]
   *     [ 7    6    5    4 |  3    2    1    0 ]
   *
   *   0 - 0: time integrated DOFs of cell 0 are copied from the buffer.
   *   1 - 1: DOFs of cell 1 are integrated in time via time derivatives.
   *   2 - 1: DOFs of cell 2 are integrated in time via time derivaitves.
   *   3 - 0: time itnegrated DOFs of cell 3 are copied from the buffer.
   *
   * @param ltsSetup bitmask for the LTS setup.
   * @param faceTypes face types of the neighboring cells.
   * @param currentTime current time of the cell [0] and it's four neighbors [1], [2], [3] and [4].
   * @param timeStepWidth time step width of the cell.
   * @param timeDofs pointers to time integrated buffers or time derivatives of the four neighboring
   cells.
   * @param integrationBuffer memory where the time integration goes if derived from derivatives.
   Ensure thread safety!
   * @param timeIntegrated pointers to the time integrated DOFs of the four neighboring cells
   (either local integration buffer or integration buffer of input).
   **/
  static void computeIntegrals(Time& time,
                               unsigned short ltsSetup,
                               const FaceType faceTypes[4],
                               const double currentTime[5],
                               double timeStepWidth,
                               real* const timeDofs[4],
                               real integrationBuffer[4][tensor::I::size()],
                               real* timeIntegrated[4]);

  /**
   * Special case of the computeIntegrals function, which assumes a common "current time" for all
   *face neighbors which provide derivatives.
   *
   * @param ltsSetup bitmask for the LTS setup.
   * @param faceTypes face types of the neighboring cells.
   * @param timeStepStart start time of the current cell with respect to the common point zero: Time
   *of the larger time step width prediction of the face neighbors.
   * @param timeStepWidth time step width of the cell.
   * @param timeDofs pointers to time integrated buffers or time derivatives of the four neighboring
   *cells.
   * @param integrationBuffer memory where the time integration goes if derived from derivatives.
   *Ensure thread safety!
   * @param timeIntegrated pointers to the time integrated DOFs of the four neighboring cells
   *(either local integration buffer or integration buffer of input).
   **/
  static void computeIntegrals(Time& time,
                               unsigned short ltsSetup,
                               const FaceType faceTypes[4],
                               double timeStepStart,
                               double timeStepWidth,
                               real* const timeDofs[4],
                               real integrationBuffer[4][tensor::I::size()],
                               real* timeIntegrated[4]);

  static void computeBatchedIntegrals(Time& time,
                                      double timeStepStart,
                                      double timeStepWidth,
                                      ConditionalPointersToRealsTable& table,
                                      seissol::parallel::runtime::StreamRuntime& runtime);

  TimeCommon() = delete;
};
} // namespace seissol::kernels

#endif // SEISSOL_SRC_KERNELS_TIMECOMMON_H_
