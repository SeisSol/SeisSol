
#ifndef INPUT_PARAMETERS_HPP_
#define INPUT_PARAMETERS_HPP_

#include <cstdint>
#include <string>
#include <array>
#include <yaml-cpp/yaml.h>

#include <xdmfwriter/XdmfWriter.h>

#include "Geometry/MeshReader.h"
#include "SourceTerm/Manager.h"
#include "Checkpoint/Backend.h"
#include "time_stepping/LtsWeights/WeightsFactory.h"

namespace seissol::initializer::parameters {

constexpr bool isModelAnelastic() { return NUMBER_OF_RELAXATION_MECHANISMS > 0; }

constexpr bool isModelElastic() {
#ifdef USE_ELASTIC
  return true;
#else
  return false;
#endif
}

constexpr bool isModelViscoelastic() {
#if defined(USE_VISCOELASTIC) || defined(USE_VISCOELASTIC2)
  return true;
#else
  return false;
#endif
}

constexpr bool isModelPoroelastic() {
#ifdef USE_POROELASTIC
  return true;
#else
  return false;
#endif
}

constexpr bool isModelAnisotropic() {
#ifdef USE_ANISOTROPIC
  return true;
#else
  return false;
#endif
}

struct ModelParameters {
  double gravitationalAcceleration = 9.81;
  double tv = 0.1;
  bool plasticity = false;
  bool useCellHomogenizedMaterial = false;
  double freqCentral = 0.0;
  double freqRatio = 0.0;
  std::string materialFileName = "";
  std::string boundaryFileName = "";
  bool hasBoundaryFile = false;
};

enum class InitializationType : int {
  Zero,
  Planarwave,
  SuperimposedPlanarwave,
  Travelling,
  Scholte,
  Snell,
  Ocean0,
  Ocean1,
  Ocean2
};

struct InitializationParameters {
  InitializationType type;
  std::array<double, 3> origin;
  std::array<double, 3> kVec;
  std::array<double, NUMBER_OF_QUANTITIES> ampField;
};

struct DynamicRuptureParameters {
  bool hasFault;
  // TODO(David): port rest of the DR parameters here?
};

enum class OutputFormat : int { None = 10, Xdmf = 6 };

enum class OutputRefinement : int { NoRefine = 0, Refine4 = 1, Refine8 = 2, Refine32 = 3 };

struct VertexWeightParameters {
  int weightElement;
  int weightDynamicRupture;
  int weightFreeSurfaceWithGravity;
};

struct MeshParameters {
  std::string meshFileName;
  MeshFormat meshFormat;
  std::array<double, 3> displacement;
  std::array<std::array<double, 3>, 3> scaling;
};

struct OutputInterval {
  double lower;
  double upper;

  bool contains(double value) const { return value >= lower && value <= upper; }
};

struct OutputBounds {
  bool enabled = false;
  OutputInterval boundsX, boundsY, boundsZ;

  bool contains(double x, double y, double z) const {
    if (enabled) {
      return boundsX.contains(x) && boundsY.contains(y) && boundsZ.contains(z);
    } else {
      return true;
    }
  }
};

struct ReceiverOutputParameters {
  bool enabled = false;
  double interval = 1.0e100;
  bool computeRotation = false;
  std::string fileName = "";
  double samplingInterval = 0.0;
};

struct FreeSurfaceOutputParameters {
  bool enabled = false;
  double interval = 1.0e100;
  unsigned refinement;
};

struct EnergyOutputParameters {
  bool enabled = false;
  double interval = 1.0e100;
  bool terminalOutput = false;
  int computeVolumeEnergiesEveryOutput = 1;
};

struct CheckpointParameters {
  bool enabled = false;
  double interval = 1.0e100;
  std::string fileName;
  seissol::checkpoint::Backend backend;
};

struct WaveFieldOutputParameters {
  bool enabled = false;
  double interval = 1.0e100;
  OutputRefinement refinement = OutputRefinement::NoRefine;
  OutputBounds bounds;
  std::array<bool, NUMBER_OF_QUANTITIES> outputMask;
  std::array<bool, 7> plasticityMask;
  std::array<bool, 9> integrationMask;
  std::unordered_set<int> groups;
};

struct OutputParameters {
  std::string prefix = "data";
  OutputFormat format = OutputFormat::None;
  xdmfwriter::BackendType xdmfWriterBackend;
  CheckpointParameters checkpointParameters;
  WaveFieldOutputParameters waveFieldParameters;
  ReceiverOutputParameters receiverParameters;
  FreeSurfaceOutputParameters freeSurfaceParameters;
  EnergyOutputParameters energyParameters;
  bool faultOutput = false;
  bool loopStatisticsNetcdfOutput = false;
};

struct LtsParameters {
  unsigned rate = 2;
  seissol::initializers::time_stepping::LtsWeightsTypes weighttype =
      seissol::initializers::time_stepping::LtsWeightsTypes::ExponentialWeights;
};

struct TimesteppingParameters {
  double cfl = 0.5;
  double maxTimestepWidth = 5000;
  LtsParameters lts;
  VertexWeightParameters vertexWeight;
};

struct SourceParameters {
  seissol::sourceterm::SourceType type = seissol::sourceterm::SourceType::None;
  std::string fileName = "";
};

struct EndParameters {
  double endTime = 15.0;
};

struct SeisSolParameters {
  ModelParameters model;
  DynamicRuptureParameters dynamicRupture;
  MeshParameters mesh;
  InitializationParameters initialization;
  OutputParameters output;
  TimesteppingParameters timestepping;
  SourceParameters source;
  EndParameters end;

  void readPar(const YAML::Node& baseNode);
};
} // namespace seissol::initializer::parameters

#endif
