// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Sebastian Rettenberger

#ifndef SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITER_H_
#define SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITER_H_

#include "Geometry/Refinement/VariableSubSampler.h"
#include "Modules/Module.h"
#include "Monitoring/Stopwatch.h"
#include "Parallel/MPI.h"
#include "Parallel/Pin.h"
#include "WaveFieldWriterExecutor.h"

#include <algorithm>
#include <array>
#include <async/Module.h>
#include <cassert>
#include <cstddef>
#include <memory>
#include <string>
#include <unordered_set>
#include <utils/logger.h>
#include <vector>

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
  seissol::SeisSol& seissolInstance_;

  /** True if wave field output is enabled */
  bool enabled_{false};

  /** False if entire region is to be written */
  bool isExtractRegionEnabled_{false};

  /** The asynchronous executor */
  WaveFieldWriterExecutor executor_;

  /** Variable buffer ids (high and low order variables) */
  int variableBufferIds_[2]{};

  /** The output prefix for the filename */
  std::string outputPrefix_;

  /** The variable subsampler for the refined mesh */
  std::unique_ptr<refinement::VariableSubsampler<double>> variableSubsampler_;

  /** The variable subsampler for the refined mesh (plastic strain) */
  std::unique_ptr<refinement::VariableSubsampler<double>> variableSubsamplerPStrain_;

  /** Number of variables */

  unsigned int numVariables_{0};

  /** Number of integrated variables */
  unsigned int numIntegratedVariables_{};

  /** Flag indicated which variables should be written */
  bool* outputFlags_{nullptr};

  /** Flag indicated which low variables should be written */
  bool* lowOutputFlags_{nullptr};

  /** Refined number of cells */
  unsigned int numCells_{0};

  /** Unrefined (low order) number of cells */
  unsigned int numLowCells_{0};

  /** Pointer to the degrees of freedom */
  const real* dofs_{nullptr};

  /** Pointer to the plastic strain */
  const real* pstrain_{nullptr};

  /** Pointer to the integrals */
  const real* integrals_{nullptr};

  /** Mapping from the cell order to dofs order */
  std::vector<unsigned int> map_;

  /** The stopwatch for the frontend */
  Stopwatch stopwatch_;

  /** Checks if a vertex given by the vertexCoords lies inside the boxBounds */
  /*   The boxBounds is in the format: xMin, xMax, yMin, yMax, zMin, zMax */
  static bool vertexInBox(const double* const boxBounds, const double* const vertexCoords) {
    return vertexCoords[0] <= boxBounds[1] && vertexCoords[0] >= boxBounds[0] &&
           vertexCoords[1] <= boxBounds[3] && vertexCoords[1] >= boxBounds[2] &&
           vertexCoords[2] <= boxBounds[5] && vertexCoords[2] >= boxBounds[4];
  }

  const refinement::TetrahedronRefiner<double>* createRefiner(int refinement);

  const unsigned* adjustOffsets(refinement::MeshRefiner<double>* meshRefiner);
  std::vector<unsigned int>
      generateRefinedClusteringData(refinement::MeshRefiner<double>* meshRefiner,
                                    const std::vector<unsigned>& ltsClusteringData,
                                    std::map<int, int>& newToOldCellMap) const;

  public:
  explicit WaveFieldWriter(seissol::SeisSol& seissolInstance) : seissolInstance_(seissolInstance) {}

  /**
   * Activate the wave field output
   */
  void enable();

  /**
   * @return True if wave field output is enabled, false otherwise
   */
  [[nodiscard]] bool isEnabled() const { return enabled_; }

  /**
   * Set the output prefix for the filename
   */
  void setFilename(const char* outputPrefix) { outputPrefix_ = outputPrefix; }

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
            const std::vector<unsigned>& ltsIdData,
            const real* dofs,
            const real* pstrain,
            const real* integrals,
            const std::size_t* map,
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
    if (enabled_) {
      wait();
    }

    finalize();

    if (!enabled_) {
      return;
    }

    stopwatch_.printTime("Time wave field writer frontend:");

    delete[] outputFlags_;
    outputFlags_ = nullptr;
    delete[] lowOutputFlags_;
    lowOutputFlags_ = nullptr;
  }

  void tearDown() override { executor_.finalize(); }

  //
  // Hooks
  //
  void simulationStart(std::optional<double> checkpointTime) override;

  void syncPoint(double currentTime) override;
};

} // namespace writer

} // namespace seissol

#endif // SEISSOL_SRC_RESULTWRITER_WAVEFIELDWRITER_H_
