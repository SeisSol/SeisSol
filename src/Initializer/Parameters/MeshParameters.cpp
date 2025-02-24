// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MeshParameters.h"
#include <Initializer/InputAux.h>
#include <Initializer/Parameters/ParameterReader.h>

namespace seissol::initializer::parameters {

MeshParameters readMeshParameters(ParameterReader* baseReader) {
  auto* reader = baseReader->readSubNode("meshnml");

  const auto meshFormat =
      reader->readWithDefaultStringEnum<MeshFormat>("meshgenerator",
                                                    "puml",
                                                    {{"netcdf", MeshFormat::Netcdf},
                                                     {"puml", MeshFormat::PUML},
                                                     {"cubegenerator", MeshFormat::CubeGenerator}});
  const std::string meshFileName = reader->readPathOrFail("meshfile", "No mesh file given.");
  const std::string partitioningLib =
      reader->readWithDefault("partitioninglib", std::string("Default"));
  const auto pumlBoundaryFormat =
      reader->readWithDefaultStringEnum<BoundaryFormat>("pumlboundaryformat",
                                                        "auto",
                                                        {
                                                            {"auto", BoundaryFormat::Auto},
                                                            {"i32", BoundaryFormat::I32},
                                                            {"i64", BoundaryFormat::I64},
                                                            {"i32x4", BoundaryFormat::I32x4},
                                                        });

  const auto displacementRaw = seissol::initializer::convertStringToArray<double, 3>(
      reader->readWithDefault("displacement", std::string("0.0 0.0 0.0")));
  const Eigen::Vector3d displacement(displacementRaw.data());
  const auto scalingXRaw = seissol::initializer::convertStringToArray<double, 3>(
      reader->readWithDefault("scalingmatrixx", std::string("1.0 0.0 0.0")));
  const auto scalingYRaw = seissol::initializer::convertStringToArray<double, 3>(
      reader->readWithDefault("scalingmatrixy", std::string("0.0 1.0 0.0")));
  const auto scalingZRaw = seissol::initializer::convertStringToArray<double, 3>(
      reader->readWithDefault("scalingmatrixz", std::string("0.0 0.0 1.0")));
  Eigen::Matrix3d scaling;
  scaling << scalingXRaw[0], scalingXRaw[1], scalingXRaw[2], scalingYRaw[0], scalingYRaw[1],
      scalingYRaw[2], scalingZRaw[0], scalingZRaw[1], scalingZRaw[2];

  const bool showEdgeCutStatistics = reader->readWithDefault("showedgecutstatistics", false);

  reader->warnDeprecated({"periodic", "periodic_direction"});

  return MeshParameters{showEdgeCutStatistics,
                        pumlBoundaryFormat,
                        meshFormat,
                        meshFileName,
                        partitioningLib,
                        displacement,
                        scaling};
}
} // namespace seissol::initializer::parameters
