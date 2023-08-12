
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
#include "Model/datastructures.hpp"

namespace seissol::initializer::parameters {

struct ModelParameters {
  double gravitationalAcceleration;
  double tv;
  bool plasticity;
  bool useCellHomogenizedMaterial;
  double freqCentral;
  double freqRatio;
  std::string materialFileName;
  std::string configFileName;
  std::string boundaryFileName;
  bool hasBoundaryFile;
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
  std::array<double, seissol::model::Material_t::NumberOfQuantities> ampField;
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
  std::string partitioningLib;
  seissol::geometry::MeshFormat meshFormat;
  std::array<double, 3> displacement;
  std::array<std::array<double, 3>, 3> scaling;
};

struct OutputInterval {
  double lower;
  double upper;

  bool contains(double value) const { return value >= lower && value <= upper; }
};

struct OutputBounds {
  bool enabled;
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
  bool enabled;
  double interval;
  bool computeRotation;
  std::string fileName;
  double samplingInterval;
};

struct FreeSurfaceOutputParameters {
  bool enabled;
  double interval;
  unsigned refinement;
};

struct EnergyOutputParameters {
  bool enabled;
  double interval;
  bool terminalOutput;
  int computeVolumeEnergiesEveryOutput;
};

struct CheckpointParameters {
  bool enabled;
  double interval;
  std::string fileName;
  seissol::checkpoint::Backend backend;
};

struct WaveFieldOutputParameters {
  bool enabled;
  double interval;
  OutputRefinement refinement;
  OutputBounds bounds;
  std::array<bool, NUMBER_OF_QUANTITIES> outputMask;
  std::array<bool, 7> plasticityMask;
  std::array<bool, 9> integrationMask;
  std::unordered_set<int> groups;
};

struct OutputParameters {
  std::string prefix;
  OutputFormat format;
  xdmfwriter::BackendType xdmfWriterBackend;
  CheckpointParameters checkpointParameters;
  WaveFieldOutputParameters waveFieldParameters;
  ReceiverOutputParameters receiverParameters;
  FreeSurfaceOutputParameters freeSurfaceParameters;
  EnergyOutputParameters energyParameters;
  bool faultOutput;
  bool loopStatisticsNetcdfOutput;
};

struct LtsParameters {
  unsigned rate;
  seissol::initializers::time_stepping::LtsWeightsTypes weighttype;
};

struct TimeSteppingParameters {
  double cfl;
  double maxTimestepWidth;
  LtsParameters lts;
  VertexWeightParameters vertexWeight;
};

struct SourceParameters {
  seissol::sourceterm::SourceType type;
  std::string fileName;
};

struct EndParameters {
  double endTime;
};

struct SeisSolParameters {
  ModelParameters model;
  DynamicRuptureParameters dynamicRupture;
  MeshParameters mesh;
  InitializationParameters initialization;
  OutputParameters output;
  TimeSteppingParameters timeStepping;
  SourceParameters source;
  EndParameters end;

  void readParameters(const YAML::Node& baseNode);
};
} // namespace seissol::initializer::parameters

#endif
