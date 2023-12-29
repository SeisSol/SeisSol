#include "SourceParameters.h"

namespace seissol::initializers::parameters {

SourceParameters readSourceParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("sourcetype");

  const auto type = reader->readWithDefaultEnum(
      "type",
      PointSourceType::None,
      {PointSourceType::None, PointSourceType::FsrmSource, PointSourceType::NrfSource});
  const auto fileName = reader->readIfRequired<std::string>("filename", type != PointSourceType::None);
  reader->warnDeprecated({"rtype", "ndirac", "npulsesource", "nricker"});

  return SourceParameters{type, fileName};
}

} // namespace seissol::initializers::parameters
