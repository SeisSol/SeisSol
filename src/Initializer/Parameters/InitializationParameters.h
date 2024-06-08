#ifndef SEISSOL_INITIALIZATION_PARAMETERS_H
#define SEISSOL_INITIALIZATION_PARAMETERS_H

#include <Eigen/Dense>

#include "Equations/datastructures.hpp"
#include "Initializer/InputAux.hpp"
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
  PressureInjection
};

struct InitializationParameters {
  InitializationType type;
  Eigen::Vector3d origin;
  Eigen::Vector3d kVec;
  Eigen::Vector<double, seissol::model::Material_t::NumberOfQuantities> ampField;
  double magnitude;
  double width;
  double k;
};

InitializationParameters readInitializationParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif
