// Copyright (C) 2023 Intel Corporation
// SPDX-License-Identifier: BSD-3-Clause

#include "LtsStructure.hpp"

#include <ostream>

namespace std {
std::ostream& operator<<(std::ostream& os, seissol::Region const& region) {
  os << "Region(" << std::endl;
  os << "  neighborRank=" << region.neighborRank << std::endl;
  os << "  neighborTimeClusterId=" << region.neighborTimeClusterId << std::endl;
  os << "  numberOfGhostRegionCells=" << region.numberOfGhostRegionCells << std::endl;
  os << "  numberOfGhostRegionDerivatives=" << region.numberOfGhostRegionDerivatives << std::endl;
  os << "  ghostRegionSize=" << region.ghostRegionSize << std::endl;
  os << "  numberOfCopyRegionCells=" << region.numberOfCopyRegionCells << std::endl;
  os << "  numberOfCommunicatedCopyRegionDerivatives="
     << region.numberOfCommunicatedCopyRegionDerivatives << std::endl;
  os << "  copyRegionSize=" << region.copyRegionSize << std::endl;
  os << "  sendIdentifier=" << region.sendIdentifier << std::endl;
  os << "  receiveIdentifier=" << region.receiveIdentifier << std::endl;
  os << ")";
  return os;
}
std::ostream& operator<<(std::ostream& os, seissol::TimeCluster const& tc) {
  os << "TimeCluster(" << std::endl;
  os << "  timeClusterId=" << tc.timeClusterId << std::endl;
  os << "  timeStepRate=" << tc.timeStepRate << std::endl;
  os << "  cflTimeStepWidth=" << tc.cflTimeStepWidth << std::endl;
  os << "  numberOfInteriorCells=" << tc.numberOfInteriorCells << std::endl;
  os << "[" << std::endl;
  auto it = tc.regions.begin();
  os << *it++;
  for (; it != tc.regions.end(); ++it) {
    os << "," << std::endl << *it;
  }
  os << "])";
  return os;
}
} // namespace std
