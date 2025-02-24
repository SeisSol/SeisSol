// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SourceParameters.h"
#include <Initializer/Parameters/ParameterReader.h>

namespace seissol::initializer::parameters {

SourceParameters readSourceParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("sourcetype");

  const auto type = reader->readWithDefaultEnum(
      "type",
      PointSourceType::None,
      {PointSourceType::None, PointSourceType::FsrmSource, PointSourceType::NrfSource});
  const auto fileName = [&]() -> std::string {
    if (type == PointSourceType::None) {
      reader->markUnused({"filename"});
      return "";
    } else {
      return reader->readPathOrFail("filename",
                                    "Point sources were enabled but no file was given.");
    }
  }();
  reader->warnDeprecated({"rtype", "ndirac", "npulsesource", "nricker"});

  return SourceParameters{type, fileName};
}

} // namespace seissol::initializer::parameters
