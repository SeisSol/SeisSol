// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_LATTICEDPSEARCH_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_LATTICEDPSEARCH_H_

#include "Initializer/Clustering/ClusteringEvaluator.h"
#include "Initializer/Clustering/LadderSearch.h"

namespace seissol::initializer {

/// Chooses the ladder as well as the wiggle factor.
///
/// For every wiggle factor on the same grid the legacy search uses, the cell timesteps are
/// binned by update factor and the cheapest divisibility chain is found by dynamic
/// programming. Candidate evaluation therefore costs one histogram reduction plus local
/// work, instead of a smoothing fixed point per candidate; only the winner is realized.
///
/// With the default cost model this degenerates to "as many clusters as possible", because
/// pure update counting never charges for a cluster. Setting a launch cost or a fill
/// threshold is what makes coarser ladders competitive.
class LatticeDpSearch : public LadderSearch {
  public:
  SearchResult run(ClusteringEvaluator& evaluator, const SearchConstraints& constraints) override;

  [[nodiscard]] std::size_t reductions() const { return reductions_; }

  private:
  std::size_t reductions_{0};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_LATTICEDPSEARCH_H_
