// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#ifndef LTSSTRUCTURE_20231019_HPP
#define LTSSTRUCTURE_20231019_HPP

#include <vector>
#include <iosfwd>

namespace seissol {

struct Region {
  long neighborRank;
  long neighborTimeClusterId;
  long neighborTimeStepRate;
  double neighborCflTimeStepWidth;
  long numberOfGhostRegionCells;
  long numberOfGhostRegionDerivatives;
  long ghostRegionSize;
  long numberOfCopyRegionCells;
  long numberOfCommunicatedCopyRegionDerivatives;
  long copyRegionSize;
  long sendIdentifier;
  long receiveIdentifier;
};

struct TimeCluster {
  long timeClusterId;
  long timeStepRate;
  double cflTimeStepWidth;
  long numberOfInteriorCells;
  long interiorLayerSize;
  std::vector<Region> regions;
};

} // namespace seissol

namespace std {
std::ostream& operator<<(std::ostream& os, seissol::Region const& region);
std::ostream& operator<<(std::ostream& os, seissol::TimeCluster const& tc);
}; // namespace std

#endif // LTSSTRUCTURE_20231019_HPP
