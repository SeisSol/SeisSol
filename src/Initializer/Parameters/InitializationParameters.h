#ifndef SEISSOL_INITIALIZATION_PARAMETERS_H
#define SEISSOL_INITIALIZATION_PARAMETERS_H

#include <Eigen/Dense>

#include "Equations/Datastructures.h"
#include "Initializer/InputAux.h"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

enum class InitializationType : int {
  Zero,
  Planarwave,
  SuperimposedPlanarwave,
  Travelling,
  AcousticTravellingWithITM,
  Scholte,
  Snell,
  Ocean0,
  Ocean1,
  Ocean2,
  PressureInjection,
  Easi
};

struct InitializationParameters {
  InitializationType type;
  Eigen::Vector3d origin;
  Eigen::Vector3d kVec;
  Eigen::Vector<double, seissol::model::MaterialT::NumQuantities> ampField;
  double magnitude;
  double width;
  double k;
  std::string filename;
  bool hasTime;
};

InitializationParameters readInitializationParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif
