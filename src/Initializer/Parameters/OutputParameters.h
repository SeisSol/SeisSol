// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_OUTPUTPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_OUTPUTPARAMETERS_H_

#include "Common/ConfigRegistry.h"
#include "Initializer/InputAux.h"
#include "Model/PlasticityQuantities.h"
#include "ParameterReader.h"

#include <array>
#include <list>
#include <optional>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace seissol::initializer::parameters {

constexpr double VeryLongTime = 1.0e100;

enum class FaultRefinement { Triple = 1, Quad = 2, None = 3 };

enum class OutputFormat : int { None = 10, Xdmf = 6 };

enum class VolumeRefinement : int { NoRefine = 0, Refine4 = 1, Refine8 = 2, Refine32 = 3 };

enum class XdmfBackend : int { Posix, Hdf5 };

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

/**
 * @brief How the steps of a mesh output are distributed over files.
 *
 * `Snapshot` writes one file per step, each complete in itself. `Incremental` keeps the geometry
 * -- which does not move -- in a file of its own and lets every snapshot link to it, which is the
 * cheaper choice whenever the mesh is large against the fields written on it. `Monolith` puts
 * every step into a single file as a VTKHDF time series, so a run leaves one file behind instead
 * of one per output step.
 *
 * Only the VTKHDF format can do more than `Snapshot`; an Xdmf output (that is, one written
 * without a `vtkorder`) references its unchanging parts from every step anyway. An output whose
 * parameters choose none is written as `Monolith` in VTKHDF and as a snapshot per step in Xdmf.
 */
enum class TimeSeriesMode : int { Snapshot, Incremental, Monolith };

struct CheckpointParameters {
  bool enabled{false};
  double interval{0};
};

struct ElementwiseFaultParameters {
  double printTimeIntervalSec{1.0};
  std::array<bool, 12> outputMask{true, true, true, true};
  FaultRefinement refinementStrategy{FaultRefinement::Quad};
  int refinement{2};
  int vtkorder{-1};
  // none if the parameters choose none (see TimeSeriesMode)
  std::optional<TimeSeriesMode> timeSeries;
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
  // a derived-output program (sderiv or Lua) whose outputs are written along
  std::string script;
};

enum class ReceiverOutputFormat { Csv, Hdf5 };

struct PickpointParameters {
  int printTimeInterval{1};
  double writeInterval{VeryLongTime};
  std::array<bool, 12> outputMask{true, true, true};
  std::optional<std::string> pickpointFileName;
  bool collectiveio{false};
  ReceiverOutputFormat format{ReceiverOutputFormat::Csv};
  //! @brief As ReceiverOutputParameters::samplechunk, for the on-fault receivers.
  std::size_t samplechunk{0};
  /**
   * @brief The time between two samples of the on-fault receivers (`pickdt`); zero samples every
   * `printTimeInterval` time steps of the finest time cluster with fault faces instead.
   */
  double samplingInterval{0};
};

struct ReceiverOutputParameters {
  bool enabled{false};
  ReceiverOutputFormat format{ReceiverOutputFormat::Csv};
  bool computeRotation{false};
  bool computeStrain{false};
  double writeInterval{VeryLongTime};
  double samplingInterval{0};
  std::string fileName;
  bool collectiveio{false};
  /**
   * @brief How far a storage chunk of the HDF5 receiver output reaches along the sample axis.
   *
   * Zero lets the writer decide. The chunking is settled when the file is created, while the
   * number of samples a write carries varies as soon as the sampling interval does not divide
   * the write interval, so which value pays off is a property of how the file is read afterwards.
   */
  std::size_t samplechunk{0};
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
  // a derived-output program (sderiv or Lua) whose outputs are written along
  std::string script;
};

struct OutputParameters {
  OutputFormat format{OutputFormat::None};
  XdmfBackend xdmfWriterBackend{};
  uint32_t hdfcompress{0};
  std::string prefix;
  CheckpointParameters checkpointParameters;
  ElementwiseFaultParameters elementwiseParameters;
  EnergyOutputParameters energyParameters;
  FreeSurfaceOutputParameters freeSurfaceParameters;
  PickpointParameters pickpointParameters;
  ReceiverOutputParameters receiverParameters;
  WaveFieldOutputParameters waveFieldParameters;
  /// The model (an easi file, or a script) whose `wavefield` and `surface` restrict the cells and
  /// the faces of the free surface written (cf. OutputRegions); empty for none.
  std::string regionFileName;
  //! the simulated time after which a simulation epoch ends at the latest; zero for no limit
  double epochInterval{0};
  //! whether the run measures what its cells cost, and writes it to <prefix>-vertexweights.h5 at
  //! the end, for a later run to partition by (discretization.vertexweightfile)
  bool vertexWeightOutput{false};

  OutputParameters() = default;
  OutputParameters(OutputFormat format,
                   XdmfBackend xdmfWriterBackend,
                   uint32_t hdfcompress,
                   std::string prefix,
                   const CheckpointParameters& checkpointParameters,
                   const ElementwiseFaultParameters& elementwiseParameters,
                   const EnergyOutputParameters& energyParameters,
                   FreeSurfaceOutputParameters freeSurfaceParameters,
                   PickpointParameters pickpointParameters,
                   ReceiverOutputParameters receiverParameters,
                   WaveFieldOutputParameters waveFieldParameters,
                   std::string regionFileName)
      : format(format), xdmfWriterBackend(xdmfWriterBackend), hdfcompress(hdfcompress),
        prefix(std::move(prefix)), checkpointParameters(checkpointParameters),
        elementwiseParameters(elementwiseParameters), energyParameters(energyParameters),
        freeSurfaceParameters(std::move(freeSurfaceParameters)),
        pickpointParameters(std::move(pickpointParameters)),
        receiverParameters(std::move(receiverParameters)),
        waveFieldParameters(std::move(waveFieldParameters)),
        regionFileName(std::move(regionFileName)) {}
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

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_OUTPUTPARAMETERS_H_
