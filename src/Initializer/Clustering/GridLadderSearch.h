// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_CLUSTERING_GRIDLADDERSEARCH_H_
#define SEISSOL_SRC_INITIALIZER_CLUSTERING_GRIDLADDERSEARCH_H_

#include "Initializer/Clustering/ClusteringEvaluator.h"
#include "Initializer/Clustering/LadderSearch.h"

#include <optional>

namespace seissol::initializer {

/// The wiggle factor grid sweep with a top-merge scan per candidate, as SeisSol has always
/// done it. Kept as the default so that existing parameter files keep their behavior.
class GridLadderSearch : public LadderSearch {
  public:
  SearchResult run(ClusteringEvaluator& evaluator, const SearchConstraints& constraints) override;

  /// Total number of demotions performed while searching, for reporting.
  [[nodiscard]] int reductions() const { return reductions_; }

  private:
  SearchResult sweep(ClusteringEvaluator& evaluator,
                     const SearchConstraints& constraints,
                     std::optional<double> baselineCost,
                     bool autoMerge);

  int reductions_{0};
};

} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_CLUSTERING_GRIDLADDERSEARCH_H_
