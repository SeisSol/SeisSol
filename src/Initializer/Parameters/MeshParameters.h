// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_MESHPARAMETERS_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_MESHPARAMETERS_H_

#include <Eigen/Dense>
#include <string>

#include "Initializer/InputAux.h"
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

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_MESHPARAMETERS_H_
