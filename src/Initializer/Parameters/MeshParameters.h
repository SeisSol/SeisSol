#ifndef SEISSOL_MESH_PARAMETERS_H
#define SEISSOL_MESH_PARAMETERS_H

#include <string>
#include <Eigen/Dense>

#include "Initializer/InputAux.hpp"
#include "ParameterReader.h"

namespace seissol::initializer::parameters {

enum class MeshFormat : int { Netcdf, PUML, CubeGenerator };

enum class BoundaryFormat : int { Auto, I32, I64, I32x4 };

struct MeshParameters {
  bool showEdgeCutStatistics;
  BoundaryFormat pumlBoundaryFormat;
  MeshFormat meshFormat;
  std::string meshFileName;
  std::string partitioningLib;
  Eigen::Vector3d displacement;
  Eigen::Matrix3d scaling;
};

MeshParameters readMeshParameters(ParameterReader* baseReader);
} // namespace seissol::initializer::parameters

#endif
