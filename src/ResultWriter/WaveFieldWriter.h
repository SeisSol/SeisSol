// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de,
 * http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 */

#ifndef SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITER_H_

#include "Parallel/MPI.h"
#include "Parallel/Pin.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <memory>
#include <string>
#include <unordered_set>
#include <vector>

#include "utils/logger.h"

#include "async/Module.h"

#include "Checkpoint/DynStruct.h"
#include "Geometry/Refinement/VariableSubSampler.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"
#include "WaveFieldWriterExecutor.h"

// for OutputBounds
#include "Initializer/Parameters/SeisSolParameters.h"

namespace seissol {
class SeisSol;
namespace refinement {
template <typename T>
class MeshRefiner;
} // namespace refinement

namespace writer {

class WaveFieldWriter
    : private async::Module<WaveFieldWriterExecutor, WaveFieldInitParam, WaveFieldParam>,
      public seissol::Module {
  seissol::SeisSol& seissolInstance;

  /** True if wave field output is enabled */
  bool m_enabled;

  /** The timestep component in the checkpoint header */
  DynStruct::Component<int> m_timestepComp;

  /** False if entire region is to be written */
  bool isExtractRegionEnabled;

  /** The asynchronous executor */
  WaveFieldWriterExecutor m_executor;

  /** Variable buffer ids (high and low order variables) */
  int m_variableBufferIds[2];

  /** The output prefix for the filename */
  std::string m_outputPrefix;

  /** The variable subsampler for the refined mesh */
  std::unique_ptr<refinement::VariableSubsampler<double>> m_variableSubsampler;

  /** The variable subsampler for the refined mesh (plastic strain) */
  std::unique_ptr<refinement::VariableSubsampler<double>> m_variableSubsamplerPStrain;

  /** Number of variables */

  unsigned int m_numVariables;

  /** Number of integrated variables */
  unsigned int m_numIntegratedVariables;

  /** Flag indicated which variables should be written */
  bool* m_outputFlags;

  /** Flag indicated which low variables should be written */
  bool* m_lowOutputFlags;

  /** Refined number of cells */
  unsigned int m_numCells;

  /** Unrefined (low order) number of cells */
  unsigned int m_numLowCells;

  /** Pointer to the degrees of freedom */
  const real* m_dofs;

  /** Pointer to the plastic strain */
  const real* m_pstrain;

  /** Pointer to the integrals */
  const real* m_integrals;

  /** Mapping from the cell order to dofs order */
  unsigned int* m_map;

  /** The stopwatch for the frontend */
  Stopwatch m_stopwatch;

  /** Checks if a vertex given by the vertexCoords lies inside the boxBounds */
  /*   The boxBounds is in the format: xMin, xMax, yMin, yMax, zMin, zMax */
  bool vertexInBox(const double* const boxBounds, const double* const vertexCoords) {
    if (vertexCoords[0] <= boxBounds[1] && vertexCoords[0] >= boxBounds[0] &&
        vertexCoords[1] <= boxBounds[3] && vertexCoords[1] >= boxBounds[2] &&
        vertexCoords[2] <= boxBounds[5] && vertexCoords[2] >= boxBounds[4]) {
      return true;
    } else {
      return false;
    }
  }

  refinement::TetrahedronRefiner<double>* createRefiner(int refinement);

  const unsigned* adjustOffsets(refinement::MeshRefiner<double>* meshRefiner);
  std::vector<unsigned int>
      generateRefinedClusteringData(refinement::MeshRefiner<double>* meshRefiner,
                                    const std::vector<unsigned>& ltsClusteringData,
                                    std::map<int, int>& newToOldCellMap);

  public:
  WaveFieldWriter(seissol::SeisSol& seissolInstance)
      : seissolInstance(seissolInstance), m_enabled(false), isExtractRegionEnabled(false),
        m_numVariables(0), m_outputFlags(0L), m_lowOutputFlags(0L), m_numCells(0), m_numLowCells(0),
        m_dofs(0L), m_pstrain(0L), m_integrals(0L), m_map(0L) {}

  /**
   * Activate the wave field output
   */
  void enable();

  /**
   * @return True if wave field output is enabled, false otherwise
   */
  bool isEnabled() const { return m_enabled; }

  /**
   * Set the output prefix for the filename
   */
  void setFilename(const char* outputPrefix) { m_outputPrefix = outputPrefix; }

  /**
   * Called by ASYNC on all ranks
   */
  void setUp() override;

  void setWaveFieldInterval(double interval) { setSyncInterval(interval); }

  /**
   * Initialize the wave field ouput
   *
   * @param map The mapping from the cell order to dofs order
   * @param timeTolerance The tolerance in the time for ignoring duplicate time steps
   */
  void init(unsigned int numVars,
            int order,
            int numAlignedDOF,
            const seissol::geometry::MeshReader& meshReader,
            const std::vector<unsigned>& ltsClusteringData,
            const real* dofs,
            const real* pstrain,
            const real* integrals,
            unsigned int* map,
            const seissol::initializer::parameters::WaveFieldOutputParameters& parameters,
            xdmfwriter::BackendType backend,
            const std::string& backupTimeStamp);

  /**
   * Write a time step
   */
  void write(double time);

  /**
   * Close wave field writer and free resources
   */
  void close() {
    // Cleanup the executor
    if (m_enabled) {
      wait();
    }

    finalize();

    if (!m_enabled)
      return;

    m_stopwatch.printTime("Time wave field writer frontend:");

    delete[] m_outputFlags;
    m_outputFlags = 0L;
    delete[] m_lowOutputFlags;
    m_lowOutputFlags = 0L;
    if (isExtractRegionEnabled) {
      delete[] m_map;
      m_map = 0L;
    }
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

#endif // SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITER_H_
