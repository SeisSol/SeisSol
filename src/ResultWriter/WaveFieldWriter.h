// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

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
  bool m_enabled{false};

  /** False if entire region is to be written */
  bool isExtractRegionEnabled{false};

  /** The asynchronous executor */
  WaveFieldWriterExecutor m_executor;

  /** Variable buffer ids (high and low order variables) */
  int m_variableBufferIds[2]{};

  /** The output prefix for the filename */
  std::string m_outputPrefix;

  /** The variable subsampler for the refined mesh */
  std::unique_ptr<refinement::VariableSubsampler<double>> m_variableSubsampler;

  /** The variable subsampler for the refined mesh (plastic strain) */
  std::unique_ptr<refinement::VariableSubsampler<double>> m_variableSubsamplerPStrain;

  /** Number of variables */

  unsigned int m_numVariables{0};

  /** Number of integrated variables */
  unsigned int m_numIntegratedVariables{};

  /** Flag indicated which variables should be written */
  bool* m_outputFlags{nullptr};

  /** Flag indicated which low variables should be written */
  bool* m_lowOutputFlags{nullptr};

  /** Refined number of cells */
  unsigned int m_numCells{0};

  /** Unrefined (low order) number of cells */
  unsigned int m_numLowCells{0};

  /** Pointer to the degrees of freedom */
  const real* m_dofs{nullptr};

  /** Pointer to the plastic strain */
  const real* m_pstrain{nullptr};

  /** Pointer to the integrals */
  const real* m_integrals{nullptr};

  /** Mapping from the cell order to dofs order */
  std::vector<unsigned int> m_map;

  /** The stopwatch for the frontend */
  Stopwatch m_stopwatch;

  /** Checks if a vertex given by the vertexCoords lies inside the boxBounds */
  /*   The boxBounds is in the format: xMin, xMax, yMin, yMax, zMin, zMax */
  static bool vertexInBox(const double* const boxBounds, const double* const vertexCoords) {
    return vertexCoords[0] <= boxBounds[1] && vertexCoords[0] >= boxBounds[0] &&
           vertexCoords[1] <= boxBounds[3] && vertexCoords[1] >= boxBounds[2] &&
           vertexCoords[2] <= boxBounds[5] && vertexCoords[2] >= boxBounds[4];
  }

  refinement::TetrahedronRefiner<double>* createRefiner(int refinement);

  const unsigned* adjustOffsets(refinement::MeshRefiner<double>* meshRefiner);
  std::vector<unsigned int>
      generateRefinedClusteringData(refinement::MeshRefiner<double>* meshRefiner,
                                    const std::vector<unsigned>& ltsClusteringData,
                                    std::map<int, int>& newToOldCellMap) const;

  public:
  WaveFieldWriter(seissol::SeisSol& seissolInstance) : seissolInstance(seissolInstance) {}

  /**
   * Activate the wave field output
   */
  void enable();

  /**
   * @return True if wave field output is enabled, false otherwise
   */
  [[nodiscard]] bool isEnabled() const { return m_enabled; }

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
            const unsigned int* map,
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

    if (!m_enabled) {
      return;
    }

    m_stopwatch.printTime("Time wave field writer frontend:");

    delete[] m_outputFlags;
    m_outputFlags = nullptr;
    delete[] m_lowOutputFlags;
    m_lowOutputFlags = nullptr;
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
