
#ifndef INPUT_PARAMETERS_HPP_
#define INPUT_PARAMETERS_HPP_

#include <cstdint>
#include <string>
#include <array>
#include <unordered_set>
#include <yaml-cpp/yaml.h>

#include <xdmfwriter/XdmfWriter.h>

#include "Geometry/MeshReader.h"
#include "Geometry/CubeGenerator.h"
#include "SourceTerm/typedefs.hpp"
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
  double gravitationalAcceleration;
  double tv;
  bool plasticity;
  bool useCellHomogenizedMaterial;
  double freqCentral;
  double freqRatio;
  std::string materialFileName;
  std::string boundaryFileName;
  bool hasBoundaryFile;
};

enum class ReflectionType : int { bothwaves = 1, bothwaves_velocity = 2, pwave = 3, swave = 4 };

struct ITMParameters {
  double ITMTime;
  double ITMVelocityScalingFactor;
  double ITMStartingTime;
  bool ITMToggle;
  ReflectionType reflectionType;
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
  Ocean2,
  PressureInjection,
  AcousticTravellingwithITM,
};

struct InitializationParameters {
  InitializationType type;
  std::array<double, 3> origin;
  std::array<double, 3> kVec;
  std::array<double, NUMBER_OF_QUANTITIES> ampField;
  double magnitude;
  double width;
  double k;
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
  bool showEdgeCutStatistics;
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
  MeshParameters mesh;
  seissol::geometry::CubeGeneratorParameters cubeGenerator;
  InitializationParameters initialization;
  OutputParameters output;
  TimeSteppingParameters timeStepping;
  SourceParameters source;
  EndParameters end;
  ITMParameters itmParameters;

  void readParameters(const YAML::Node& baseNode);
};
} // namespace seissol::initializer::parameters

#endif
