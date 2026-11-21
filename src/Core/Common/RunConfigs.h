// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_CORE_COMMON_RUNCONFIGS_H_
#define SEISSOL_SRC_CORE_COMMON_RUNCONFIGS_H_

#include "Core/Common/ConfigRegistry.h"

#include <cstddef>
#include <unordered_map>
#include <vector>

namespace seissol::initializer::parameters {

/// The configurations the cells of a run compute in, and the number of its simulations.
struct RunConfigs {
  // the configuration the cells of the run compute in, unless their mesh group has its own
  ConfigId config{};
  // the mesh groups whose cells compute in a configuration of their own
  std::unordered_map<int, ConfigId> groupConfigs;
  // the number of simulations of the run; the fused simulations of the cells hold them, and as
  // many copies of the cells as needed
  std::size_t simulationCount{1};

  /// The configuration the cells of a mesh group compute in.
  [[nodiscard]] ConfigId configOfGroup(int group) const;
  /// Every configuration cells of the run may compute in: `config` first, then those of the mesh
  /// groups in the order of their ids, each once.
  [[nodiscard]] std::vector<ConfigId> configs() const;
};

} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_CORE_COMMON_RUNCONFIGS_H_
