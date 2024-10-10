#ifndef SEISSOL_MODEL_PARAMETERS_H
#define SEISSOL_MODEL_PARAMETERS_H

#include <string>

#include "ParameterReader.h"

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

enum class ReflectionType { BothWaves = 1, BothWavesVelocity, Pwave, Swave };

struct ITMParameters {
  bool itmEnabled;
  double itmStartingTime;
  double itmDuration;
  double itmVelocityScalingFactor;
  ReflectionType itmReflectionType;
};

enum class NumericalFlux { Godunov, Rusanov };

struct ModelParameters {
  bool hasBoundaryFile;
  bool plasticity;
  bool useCellHomogenizedMaterial;
  double freqCentral;
  double freqRatio;
  double gravitationalAcceleration;
  double tv;
  std::string boundaryFileName;
  std::string materialFileName;
  ITMParameters itmParameters;
  NumericalFlux flux;
  NumericalFlux fluxNearFault;
};

ModelParameters readModelParameters(ParameterReader* baseReader);
ITMParameters readITMParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif
