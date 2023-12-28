#include "SourceParameters.h"

namespace seissol::initializers::parameters {

SourceParameters readSourceParameters(ParameterReader& baseReader) {
  auto& reader = baseReader.readSubNode("sourcetype");

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

  return SourceParameters{type, fileName};
}

} // namespace seissol::initializers::parameters
