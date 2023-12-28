#ifndef SEISSOL_MESH_PARAMETERS_H
#define SEISSOL_MESH_PARAMETERS_H

#include <string>
#include <Eigen/Dense>

#include "Initializer/InputAux.hpp"
#include "ParameterReader.h"

namespace seissol::initializers::parameters {

enum class MeshFormat : int { Netcdf, PUML, CubeGenerator };

struct MeshParameters {
  bool showEdgeCutStatistics;
  MeshFormat meshFormat;
  std::string meshFileName;
  std::string partitioningLib;
  Eigen::Vector3d displacement;
  Eigen::Matrix3d scaling;
};

MeshParameters readMeshParameters(ParameterReader& baseReader);
} // namespace seissol::initializers::parameters

#endif
