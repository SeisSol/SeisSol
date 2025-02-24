// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_OUTPUTPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_OUTPUTPARAMETERS_H_

#include <list>
#include <string>
#include <unordered_set>

#include <xdmfwriter/backends/Backend.h>

#include "Equations/Datastructures.h"
#include "Initializer/InputAux.h"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

constexpr double VeryLongTime = 1.0e100;

enum class FaultRefinement { Triple = 1, Quad = 2, None = 3 };

enum class OutputFormat : int { None = 10, Xdmf = 6 };

enum class VolumeRefinement : int { NoRefine = 0, Refine4 = 1, Refine8 = 2, Refine32 = 3 };

struct CheckpointParameters {
  bool enabled;
  double interval;
};

struct ElementwiseFaultParameters {
  double printTimeIntervalSec{1.0};
  std::array<bool, 12> outputMask{true, true, true, true};
  FaultRefinement refinementStrategy{FaultRefinement::Quad};
  int refinement{2};
  int vtkorder{-1};
};

struct EnergyOutputParameters {
  bool enabled;
  int computeVolumeEnergiesEveryOutput;
  double interval;
  bool terminalOutput;
  int terminalPrecision;
  double terminatorMaxTimePostRupture;
  double terminatorMomentRateThreshold;
};

struct FreeSurfaceOutputParameters {
  bool enabled;
  unsigned refinement;
  double interval;
  int vtkorder{-1};
};

struct PickpointParameters {
  int printTimeInterval{1};
  int maxPickStore{50};
  std::array<bool, 12> outputMask{true, true, true};
  std::string pickpointFileName;
  bool collectiveio{false};
};

struct ReceiverOutputParameters {
  bool enabled;
  bool computeRotation;
  bool computeStrain;
  double interval;
  double samplingInterval;
  std::string fileName;
  bool collectiveio{false};
};

struct OutputInterval {
  double lower;
  double upper;

  [[nodiscard]] bool contains(double value) const { return value >= lower && value <= upper; }
};

struct OutputBounds {
  bool enabled;
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
  bool enabled;
  int vtkorder;
  double interval;
  VolumeRefinement refinement;
  OutputBounds bounds;
  std::array<bool, seissol::model::MaterialT::NumQuantities> outputMask;
  std::array<bool, 7> plasticityMask;
  std::array<bool, 9> integrationMask;
  std::unordered_set<int> groups;
};

struct OutputParameters {
  bool loopStatisticsNetcdfOutput;
  OutputFormat format;
  xdmfwriter::BackendType xdmfWriterBackend;
  std::string prefix;
  CheckpointParameters checkpointParameters;
  ElementwiseFaultParameters elementwiseParameters;
  EnergyOutputParameters energyParameters;
  FreeSurfaceOutputParameters freeSurfaceParameters;
  PickpointParameters pickpointParameters;
  ReceiverOutputParameters receiverParameters;
  WaveFieldOutputParameters waveFieldParameters;

  OutputParameters() = default;
  OutputParameters(bool loopStatisticsNetcdfOutput,
                   OutputFormat format,
                   xdmfwriter::BackendType xdmfWriterBackend,
                   const std::string& prefix,
                   const CheckpointParameters& checkpointParameters,
                   const ElementwiseFaultParameters& elementwiseParameters,
                   const EnergyOutputParameters& energyParameters,
                   const FreeSurfaceOutputParameters& freeSurfaceParameters,
                   const PickpointParameters& pickpointParameters,
                   const ReceiverOutputParameters& receiverParameters,
                   const WaveFieldOutputParameters& waveFieldParameters)
      : loopStatisticsNetcdfOutput(loopStatisticsNetcdfOutput), format(format),
        xdmfWriterBackend(xdmfWriterBackend), prefix(prefix),
        checkpointParameters(checkpointParameters), elementwiseParameters(elementwiseParameters),
        energyParameters(energyParameters), freeSurfaceParameters(freeSurfaceParameters),
        pickpointParameters(pickpointParameters), receiverParameters(receiverParameters),
        waveFieldParameters(waveFieldParameters) {}
};

void warnIntervalAndDisable(bool& enabled,
                            double interval,
                            const std::string& valName,
                            const std::string& intName);
CheckpointParameters readCheckpointParameters(ParameterReader* baseReader);
ElementwiseFaultParameters readElementwiseParameters(ParameterReader* baseReader);
EnergyOutputParameters readEnergyParameters(ParameterReader* baseReader);
FreeSurfaceOutputParameters readFreeSurfaceParameters(ParameterReader* baseReader);
PickpointParameters readPickpointParameters(ParameterReader* baseReader);
ReceiverOutputParameters readReceiverParameters(ParameterReader* baseReader);
WaveFieldOutputParameters readWaveFieldParameters(ParameterReader* baseReader);
OutputParameters readOutputParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_OUTPUTPARAMETERS_H_
