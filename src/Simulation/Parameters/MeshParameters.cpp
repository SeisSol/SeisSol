// SPDX-FileCopyrightText: 2023 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "MeshParameters.h"

#include "Core/Common/InputAux.h"
#include "Core/Geometry/MeshParameters.h"
#include "Simulation/Parameters/ParameterReader.h"

#include <Eigen/Core>
#include <optional>

namespace seissol::initializer::parameters {

MeshParameters readMeshParameters(ParameterReader* baseReader) {
  const MeshParameters defaults{};
  auto* reader = baseReader->readSubNode("meshnml");

  const auto meshFormat = reader->readWithDefaultStringEnum<MeshFormat>(
      "meshgenerator",
      defaults.meshFormat,
      {{"puml", MeshFormat::PUML}, {"cubegenerator", MeshFormat::CubeGenerator}});
  const std::string meshFileName = reader->readPathOrFail("meshfile", "No mesh file given.");
  const std::string partitioningLib =
      reader->readWithDefault("partitioninglib", defaults.partitioningLib);
  const auto pumlBoundaryFormat =
      reader->readWithDefaultStringEnum<BoundaryFormat>("pumlboundaryformat",
                                                        defaults.pumlBoundaryFormat,
                                                        {
                                                            {"auto", BoundaryFormat::Auto},
                                                            {"i32", BoundaryFormat::I32},
                                                            {"i64", BoundaryFormat::I64},
                                                            {"i32x4", BoundaryFormat::I32x4},
                                                        });
  const auto pumlTopologyFormat = reader->readWithDefaultStringEnum<TopologyFormat>(
      "pumltopologyformat",
      defaults.pumlTopologyFormat,
      {
          {"auto", TopologyFormat::Auto},
          {"geometric", TopologyFormat::Geometric},
          {"identify-face", TopologyFormat::IdentifyFace},
          {"identify-vertex", TopologyFormat::IdentifyVertex},
      });

  const auto toVector = [](const std::string& text) {
    const auto raw = seissol::initializer::convertStringToArray<double, 3>(text);
    return Eigen::Vector3d(raw.data());
  };
  const Eigen::Vector3d displacement =
      reader->readWithDefaultConverted("displacement", defaults.displacement, toVector);
  // the scaling matrix, row by row
  Eigen::Matrix3d scaling;
  Eigen::Index row = 0;
  for (const auto* field : {"scalingmatrixx", "scalingmatrixy", "scalingmatrixz"}) {
    const Eigen::Vector3d defaultRow = defaults.scaling.row(row).transpose();
    scaling.row(row) = reader->readWithDefaultConverted(field, defaultRow, toVector).transpose();
    ++row;
  }

  const bool showEdgeCutStatistics =
      reader->readWithDefault("showedgecutstatistics", defaults.showEdgeCutStatistics);

  const auto faceMapPre = reader->readPath("facemapfile");
  const auto faceMap = faceMapPre == "" ? std::optional<std::string>{} : faceMapPre;

  reader->warnDeprecated({"periodic", "periodic_direction"});

  return MeshParameters{.showEdgeCutStatistics = showEdgeCutStatistics,
                        .pumlBoundaryFormat = pumlBoundaryFormat,
                        .pumlTopologyFormat = pumlTopologyFormat,
                        .meshFormat = meshFormat,
                        .meshFileName = meshFileName,
                        .partitioningLib = partitioningLib,
                        .faceMapFile = faceMap,
                        .displacement = displacement,
                        .scaling = scaling};
}
} // namespace seissol::initializer::parameters
