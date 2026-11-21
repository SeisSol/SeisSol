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
#include <list>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
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
  bool enabled{false};
  double interval{0};
};

struct EnergyOutputParameters {
  bool enabled{false};
  int computeVolumeEnergiesEveryOutput{0};
  double interval{0};
  bool terminalOutput{false};
  int terminalPrecision{0};
  double terminatorMaxTimePostRupture{0};
  double terminatorMomentRateThreshold{0};
};

struct FreeSurfaceOutputParameters {
  bool enabled{false};
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
  bool enabled{false};
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

struct OutputParameters {
  OutputFormat format{OutputFormat::None};
  XdmfBackend xdmfWriterBackend{};
  uint32_t hdfcompress{0};
  //! the directory of the outputs; empty for the working directory
  std::string outputPath;
  //! what the names of the output files start with, in their directory
  std::string outputName;
  //! whether each run writes into a directory of its own in outputPath, named by outputName and
  //! the ID of the run (see initprocedure::createOutputDirectory)
  bool outputDirPerSim{false};
  CheckpointParameters checkpointParameters;
  ElementwiseFaultParameters elementwiseParameters;
  EnergyOutputParameters energyParameters;
  FreeSurfaceOutputParameters freeSurfaceParameters;
  PickpointParameters pickpointParameters;
  ReceiverOutputParameters receiverParameters;
  WaveFieldOutputParameters waveFieldParameters;
  //! the simulated time after which a simulation epoch ends at the latest; zero for no limit
  double epochInterval{0};
  //! whether the run measures what its cells cost, and writes it to <prefix>-vertexweights.h5 at
  //! the end, for a later run to partition by (discretization.vertexweightfile)
  bool vertexWeightOutput{false};

  OutputParameters() = default;
  OutputParameters(OutputFormat format,
                   XdmfBackend xdmfWriterBackend,
                   uint32_t hdfcompress,
                   std::string outputPath,
                   std::string outputName,
                   bool outputDirPerSim,
                   const CheckpointParameters& checkpointParameters,
                   const ElementwiseFaultParameters& elementwiseParameters,
                   const EnergyOutputParameters& energyParameters,
                   FreeSurfaceOutputParameters freeSurfaceParameters,
                   PickpointParameters pickpointParameters,
                   ReceiverOutputParameters receiverParameters,
                   WaveFieldOutputParameters waveFieldParameters)
      : format(format), xdmfWriterBackend(xdmfWriterBackend), hdfcompress(hdfcompress),
        outputPath(std::move(outputPath)), outputName(std::move(outputName)),
        outputDirPerSim(outputDirPerSim), checkpointParameters(checkpointParameters),
        elementwiseParameters(elementwiseParameters), energyParameters(energyParameters),
        freeSurfaceParameters(std::move(freeSurfaceParameters)),
        pickpointParameters(std::move(pickpointParameters)),
        receiverParameters(std::move(receiverParameters)),
        waveFieldParameters(std::move(waveFieldParameters)) {}
};

void warnIntervalAndDisable(bool& enabled,
                            double interval,
                            const std::string& valName,
                            const std::string& intName);
CheckpointParameters readCheckpointParameters(ParameterReader* baseReader);
ElementwiseFaultParameters
    readElementwiseParameters(ParameterReader* baseReader,
                              std::optional<TimeSeriesMode> defaultTimeSeries);
EnergyOutputParameters readEnergyParameters(ParameterReader* baseReader);
FreeSurfaceOutputParameters
    readFreeSurfaceParameters(ParameterReader* baseReader,
                              std::optional<TimeSeriesMode> defaultTimeSeries,
                              const std::vector<ConfigId>& configs);
PickpointParameters readPickpointParameters(ParameterReader* baseReader);
ReceiverOutputParameters readReceiverParameters(ParameterReader* baseReader);
WaveFieldOutputParameters readWaveFieldParameters(ParameterReader* baseReader,
                                                  std::optional<TimeSeriesMode> defaultTimeSeries,
                                                  const std::vector<ConfigId>& configs);
/// The output parameters of a run whose cells compute in the configurations `configs`.
OutputParameters readOutputParameters(ParameterReader* baseReader,
                                      const std::vector<ConfigId>& configs);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_SIMULATION_PARAMETERS_OUTPUTPARAMETERS_H_
