// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOCOMMUNICATION_H_
#define SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOCOMMUNICATION_H_

#include "Common/Real.h"
#include "Initializer/TimeStepping/ClusterLayout.h"
#include "Initializer/Typedefs.h"

#include <cstddef>
#include <limits>
#include <vector>
namespace seissol::solver {

struct RemoteCluster {
  void* data{nullptr};
  std::size_t size{0};
  RealType datatype{RealType::F64};
  int rank{-1};
  std::size_t tag{std::numeric_limits<std::size_t>::max()};

  RemoteCluster() = default;

  RemoteCluster(void* data, std::size_t size, RealType datatype, int rank, std::size_t tag)
      : data(data), size(size), datatype(datatype), rank(rank), tag(tag) {}
};

struct RemoteClusterPair {
  std::vector<RemoteCluster> copy;
  std::vector<RemoteCluster> ghost;
};

using HaloCommunication = std::vector<std::vector<RemoteClusterPair>>;
} // namespace seissol::solver

#endif // SEISSOL_SRC_SOLVER_TIMESTEPPING_HALOCOMMUNICATION_H_
