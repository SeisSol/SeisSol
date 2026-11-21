// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SIMULATION_PARAMETERS_OUTPUTPARAMETERS_H_
#define SEISSOL_SRC_SIMULATION_PARAMETERS_OUTPUTPARAMETERS_H_

#include "Core/Common/ConfigRegistry.h"
#include "Core/Common/InputAux.h"
#include "Core/IO/WriterParameters.h"
#include "Features/DynamicRupture/Output/FaultOutputParameters.h"
#include "Features/Plasticity/Physics/PlasticityQuantities.h"
#include "ParameterReader.h"
#include "WavePropagation/Output/ReceiverParameters.h"

#include <array>
#include <cstdint>
#include <optional>
#include <string>
#include <unordered_set>
#include <utils/logger.h>
#include <variant>
#include <vector>

namespace seissol::initializer::parameters {

enum class OutputFormat : int { None = 10, Xdmf = 6 };

enum class VolumeRefinement : int { NoRefine = 0, Refine4 = 1, Refine8 = 2, Refine32 = 3 };

/**
 * How volumetric/surface data is transferred onto the output points: either by evaluating the
 * solution there, or by an L2 projection onto the output (Lagrange) space. The latter is
 * conservative on each output subcell, the former is cheaper and reproduces the solution exactly
 * at the sampled points.
 *
 * The default differs per output, so that each keeps the behavior it had before the output
 * modules were unified: the wavefield output sampled the solution at the subcell barycenters,
 * while the free-surface output averaged over each subcell.
 */
enum class ProjectionMethod : int { Pointwise, L2 };

struct CheckpointParameters {
  double interval{0};
};

struct EnergyOutputParameters {
  int computeVolumeEnergiesEveryOutput{0};
  double interval{0};
  bool terminalOutput{false};
  int terminalPrecision{0};
  double terminatorMaxTimePostRupture{0};
  double terminatorMomentRateThreshold{0};
};

struct FreeSurfaceOutputParameters {
  unsigned refinement{0};
  double interval{0};
  int vtkorder{-1};
  // per configuration, by its id (empty for those the run does not compute in), per quantity of
  // its material
  std::vector<std::vector<bool>> outputMask;
  ProjectionMethod projection{ProjectionMethod::L2};
  // none if the parameters choose none (see TimeSeriesMode)
  std::optional<TimeSeriesMode> timeSeries;
};

struct OutputInterval {
  double lower{0};
  double upper{0};

  [[nodiscard]] bool contains(double value) const { return value >= lower && value <= upper; }
};

struct OutputBounds {
  bool enabled{false};
  OutputInterval boundsX, boundsY, boundsZ;

  OutputBounds() = default;
  OutputBounds(bool enabled,
               OutputInterval intervalX,
               OutputInterval intervalY,
               OutputInterval intervalZ)
      : enabled(enabled), boundsX(intervalX), boundsY(intervalY), boundsZ(intervalZ) {};

  [[nodiscard]] bool contains(double x, double y, double z) const {
    if (enabled) {
      return boundsX.contains(x) && boundsY.contains(y) && boundsZ.contains(z);
    } else {
      return true;
    }
  }
};

struct WaveFieldOutputParameters {
  int vtkorder{-1};
  double interval{0};
  VolumeRefinement refinement{VolumeRefinement::NoRefine};
  OutputBounds bounds;
  // per configuration, by its id (empty for those the run does not compute in), per quantity of
  // its material
  std::vector<std::vector<bool>> outputMask;
  std::array<bool, seissol::model::PlasticityQuantityCount> plasticityMask{};
  // as outputMask
  std::vector<std::vector<bool>> integrationMask;
  std::unordered_set<int> groups;
  bool computeRotation{false};
  bool computeStrain{false};
  ProjectionMethod projection{ProjectionMethod::Pointwise};
  // none if the parameters choose none (see TimeSeriesMode)
  std::optional<TimeSeriesMode> timeSeries;
};

/// One output of a run, and what it writes.
struct OutputEntry {
  //! unique among the outputs of the run; the outputs of the wave field and of the surface name
  //! their writer and files after it
  std::string name;
  std::variant<WaveFieldOutputParameters,
               FreeSurfaceOutputParameters,
               ReceiverOutputParameters,
               ElementwiseFaultParameters,
               PickpointParameters>
      settings;
};

struct OutputParameters {
  XdmfBackend xdmfWriterBackend{XdmfBackend::Posix};
  uint32_t hdfcompress{0};
  //! the directory of the outputs; empty for the working directory
  std::string outputPath;
  //! what the names of the output files start with, in their directory
  std::string outputName;
  //! whether each run writes into a directory of its own in outputPath, named by outputName and
  //! the ID of the run (see initprocedure::createOutputDirectory)
  bool outputDirPerSim{false};
  //! none if the run writes no checkpoints
  std::optional<CheckpointParameters> checkpoint;
  //! none if the run writes no energies
  std::optional<EnergyOutputParameters> energy;
  //! the outputs of the run; empty if it writes none
  std::vector<OutputEntry> outputs;
  //! the simulated time after which a simulation epoch ends at the latest; zero for no limit
  double epochInterval{0};
  //! whether the run measures what its cells cost, and writes it to <prefix>-vertexweights.h5 at
  //! the end, for a later run to partition by (discretization.vertexweightfile)
  bool vertexWeightOutput{false};

  /// The settings of the outputs of the kind `T`, in the order of `outputs`.
  template <typename T>
  [[nodiscard]] std::vector<const T*> outputsOf() const {
    std::vector<const T*> found;
    for (const auto& entry : outputs) {
      if (const auto* settings = std::get_if<T>(&entry.settings); settings != nullptr) {
        found.push_back(settings);
      }
    }
    return found;
  }

  /// The settings of the output of the kind `T`, or none. The time clusters sample the outputs at
  /// points (the receivers, those of the fault) into one place each, so a run can have at most one
  /// of each of these kinds for now; `kind` names it in the error otherwise.
  template <typename T>
  [[nodiscard]] const T* singleOutputOf(const std::string& kind) const {
    const auto found = outputsOf<T>();
    if (found.size() > 1) {
      logError() << "A run can have only one" << kind << "output for now, but this one has"
                 << found.size() << ".";
    }
    return found.empty() ? nullptr : found.front();
  }
};

void warnIntervalAndDisable(bool& enabled,
                            double interval,
                            const std::string& valName,
                            const std::string& intName);
// Each of the following reads all fields of its output, and gives none if they switch it off.
std::optional<CheckpointParameters> readCheckpointParameters(ParameterReader* baseReader);
std::optional<EnergyOutputParameters> readEnergyParameters(ParameterReader* baseReader);
std::optional<FreeSurfaceOutputParameters>
    readFreeSurfaceParameters(ParameterReader* baseReader,
                              std::optional<TimeSeriesMode> defaultTimeSeries,
                              const std::vector<ConfigId>& configs);
std::optional<ReceiverOutputParameters> readReceiverParameters(ParameterReader* baseReader);
std::optional<WaveFieldOutputParameters>
    readWaveFieldParameters(ParameterReader* baseReader,
                            std::optional<TimeSeriesMode> defaultTimeSeries,
                            const std::vector<ConfigId>& configs);
// The section of the rupture switches the outputs of the fault on (see readOutputParameters).
ElementwiseFaultParameters
    readElementwiseParameters(ParameterReader* baseReader,
                              std::optional<TimeSeriesMode> defaultTimeSeries);
PickpointParameters readPickpointParameters(ParameterReader* baseReader);
/// The output parameters of a run whose cells compute in the configurations `configs`; with the
/// outputs of its fault only if it has one (`withFault`).
OutputParameters readOutputParameters(ParameterReader* baseReader,
                                      const std::vector<ConfigId>& configs,
                                      bool withFault);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_OUTPUTPARAMETERS_H_
