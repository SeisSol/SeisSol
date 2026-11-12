// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "SourceParameters.h"

#include "Initializer/Parameters/ParameterReader.h"

#include <string>
#include <utils/logger.h>

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
#if !defined(USE_NETCDF) || defined(NETCDF_PASSIVE)
  if (type == PointSourceType::NrfSource) {
    logError() << "NRF sources (type 42) need SeisSol to be linked with an (active) Netcdf "
                  "library. However, this is not the case for this build.";
  }
#endif
  reader->warnDeprecated({"rtype", "ndirac", "npulsesource", "nricker"});

  return SourceParameters{.type = type, .fileName = fileName};
}

} // namespace seissol::initializer::parameters
