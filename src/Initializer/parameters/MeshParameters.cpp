#include "MeshParameters.h"

namespace seissol::initializers::parameters {

MeshParameters readMeshParameters(ParameterReader& baseReader) {
  auto& reader = baseReader.readSubNode("meshnml");

  const MeshFormat meshFormat =
      reader.readWithDefaultStringEnum<MeshFormat>("meshgenerator",
                                                   "puml",
                                                   {{"netcdf", MeshFormat::Netcdf},
                                                    {"puml", MeshFormat::PUML},
                                                    {"cubegenerator", MeshFormat::CubeGenerator}});
  const std::string meshFileName =
      reader.readOrFail<std::string>("meshfile", "No mesh file given.");
  const std::string partitioningLib =
      reader.readWithDefault("partitioninglib", std::string("Default"));

  const auto displacementRaw = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("displacement", std::string("0.0 0.0 0.0")));
  Eigen::Vector3d displacement(displacementRaw.data());
  const auto scalingXRaw = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("scalingmatrixx", std::string("1.0 0.0 0.0")));
  const auto scalingYRaw = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("scalingmatrixy", std::string("0.0 1.0 0.0")));
  const auto scalingZRaw = seissol::initializers::convertStringToArray<double, 3>(
      reader.readWithDefault("scalingmatrixz", std::string("0.0 0.0 1.0")));
  Eigen::Matrix3d scaling;
  scaling << scalingXRaw[0], scalingXRaw[1], scalingXRaw[2], scalingYRaw[0], scalingYRaw[1],
      scalingYRaw[2], scalingZRaw[0], scalingZRaw[1], scalingZRaw[2];

  const bool showEdgeCutStatistics = reader.readWithDefault("showedgecutstatistics", false);

  reader.warnDeprecated({"periodic", "periodic_direction"});

  return MeshParameters{
      showEdgeCutStatistics, meshFormat, meshFileName, partitioningLib, displacement, scaling};
}
} // namespace seissol::initializers::parameters
