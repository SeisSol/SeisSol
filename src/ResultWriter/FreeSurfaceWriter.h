// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Carsten Uphoff (c.uphoff AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
 * @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITER_H_

#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include "FreeSurfaceWriterExecutor.h"
#include "Geometry/MeshReader.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include <async/Module.h>
#include <utils/logger.h>

namespace seissol {
class SeisSol;
namespace writer {

class FreeSurfaceWriter
    : private async::Module<FreeSurfaceWriterExecutor, FreeSurfaceInitParam, FreeSurfaceParam>,
      public seissol::Module {
  private:
  seissol::SeisSol& seissolInstance;

  /** Is enabled? */
  bool m_enabled{false};

  /** The asynchronous executor */
  FreeSurfaceWriterExecutor m_executor;

  /** Frontend stopwatch */
  Stopwatch m_stopwatch;

  /** free surface integration module. */
  seissol::solver::FreeSurfaceIntegrator* m_freeSurfaceIntegrator{nullptr};

  void constructSurfaceMesh(const seissol::geometry::MeshReader& meshReader,
                            unsigned*& cells,
                            double*& vertices,
                            unsigned& nCells,
                            unsigned& nVertices);

  public:
  FreeSurfaceWriter(seissol::SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

  /**
   * Called by ASYNC on all ranks
   */
  void setUp() override;

  void enable();

  void init(const seissol::geometry::MeshReader& meshReader,
            seissol::solver::FreeSurfaceIntegrator* freeSurfaceIntegrator,
            const char* outputPrefix,
            double interval,
            xdmfwriter::BackendType backend,
            const std::string& backupTimeStamp);

  void write(double time);

  void close() {
    if (m_enabled) {
      wait();
    }

    finalize();

    if (!m_enabled) {
      return;
    }

    m_stopwatch.printTime("Time free surface writer frontend:");
  }

  void tearDown() override { m_executor.finalize(); }

  //
  // Hooks
  //
  void simulationStart() override;

  void syncPoint(double currentTime) override;
};

} // namespace writer

} // namespace seissol


#endif // SEISSOL_SRC_RESULTWRITER_FREESURFACEWRITER_H_

