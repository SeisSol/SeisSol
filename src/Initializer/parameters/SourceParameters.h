#ifndef SEISSOL_SOURCE_PARAMETERS_H
#define SEISSOL_SOURCE_PARAMETERS_H

#include <string>

#include "ParameterReader.h"

namespace seissol::initializers::parameters {

enum class PointSourceType : int { None = 0, NrfSource = 42, FsrmSource = 50 };

struct SourceParameters {
  PointSourceType type;
  std::string fileName;
};

SourceParameters readSourceParameters(ParameterReader& baseReader);
} // namespace seissol::initializers::parameters

#endif
