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

inline SourceParameters readSourceParameters(ParameterReader& baseReader) {
  auto reader = baseReader.readSubNode("sourcetype");

  const auto type = reader.readWithDefaultEnum(
      "type",
      PointSourceType::None,
      {PointSourceType::None, PointSourceType::FsrmSource, PointSourceType::NrfSource});
  auto readFilename = [&reader](bool enabled) {
    std::string fileName;
    if (enabled) {
      fileName = reader.readOrFail<std::string>("fileName", "No source file specified.");
    } else {
      reader.markUnused("fileName");
    }

    return fileName;
  };
  const auto fileName = readFilename(type != PointSourceType::None);
  reader.warnDeprecated({"rtype", "ndirac", "npulsesource", "nricker"});
  reader.warnUnknown();

  return SourceParameters{type, fileName};
}

} // namespace seissol::initializers::parameters

#endif
