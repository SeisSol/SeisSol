// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

#include "ODEVector.h"

namespace seissol::ode {
template <typename RealT>
std::pair<std::size_t, std::size_t> ODEVector<RealT>::index(std::size_t idx) const {
  for (std::size_t i = 0; i < storages.size(); ++i) {
    const auto begin = offsets[i];
    const auto end = begin + sizes[i];
    assert(begin <= idx);
    if (idx < end) {
      return {i, idx - begin};
    }
  }
  std::abort(); // Unreachable!
}

template <typename RealT>
ODEVector<RealT>::ODEVector(std::vector<RealT*> storages, std::vector<std::size_t> sizes)
    : storages(std::move(storages)), sizes(std::move(sizes)) {
  std::size_t curOffset = 0;
  for (const unsigned long size : this->sizes) {
    offsets.push_back(curOffset);
    curOffset += size;
  }
}

template <typename RealT>
void ODEVector<RealT>::updateStoragesAndSizes(std::vector<RealT*> newStorages,
                                              std::vector<std::size_t> newSizes) {
  storages = std::move(newStorages);
  sizes = std::move(newSizes);
  offsets.clear();
  std::size_t curOffset = 0;
  for (const unsigned long size : this->sizes) {
    offsets.push_back(curOffset);
    curOffset += size;
  }
}

template <typename RealT>
std::pair<RealT*, size_t> ODEVector<RealT>::getSubvector(size_t storageIdx) {
  return {storages[storageIdx], sizes[storageIdx]};
}

template <typename RealT>
RealT& ODEVector<RealT>::operator[](std::size_t idx) {
  const auto idxPair = index(idx);
  return storages[idxPair.first][idxPair.second];
}

template <typename RealT>
const RealT& ODEVector<RealT>::operator[](std::size_t idx) const {
  const auto idxPair = index(idx);
  return storages[idxPair.first][idxPair.second];
}

template <typename RealT>
ODEVector<RealT>& ODEVector<RealT>::operator+=(ODEVector<RealT>& rhs) {
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == rhs.sizes[i]);
#pragma omp simd
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      storages[i][j] += rhs.storages[i][j];
    }
  }
  return *this;
}

template <typename RealT>
ODEVector<RealT>& ODEVector<RealT>::operator*=(RealT scalar) {
  for (std::size_t i = 0; i < storages.size(); ++i) {
#pragma omp simd
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      storages[i][j] *= scalar;
    }
  }
  return *this;
}

template <typename RealT>
ODEVector<RealT>& ODEVector<RealT>::copyFrom(const ODEVector<RealT>& other) {
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == other.sizes[i]);
    std::copy_n(other.storages[i], sizes[i], storages[i]);
  }
  return *this;
}

template <typename RealT>
void ODEVector<RealT>::weightedAddInplace(RealT weight, const ODEVector<RealT>& rhs) {
  if (weight == 0.0) {
    return;
  }
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == rhs.sizes[i]);
#pragma omp simd
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      storages[i][j] += weight * rhs.storages[i][j];
    }
  }
}

template <typename RealT>
RealT ODEVector<RealT>::normDifferenceTo(ODEVector<RealT>& other, bool useLInfNorm) {
  // Computes the L2 or LInf norm of the difference between two vectors.
  RealT error = 0.0;
  RealT maxError = -1;
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == other.sizes[i]);
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      const RealT curDiff = storages[i][j] - other.storages[i][j];
      error += curDiff * curDiff;
      maxError = std::max(std::abs(curDiff), maxError);
    }
  }
  if (useLInfNorm) {
    return maxError;
  }
  return std::sqrt(error);
}

template <typename RealT>
RealT ODEVector<RealT>::l2Norm() {
  // Computes the L2 norm.
  RealT norm = 0.0;
  for (std::size_t i = 0; i < storages.size(); ++i) {
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      norm += storages[i][j] * storages[i][j];
    }
  }
  return std::sqrt(norm);
}

template <typename RealT>
void ODEVector<RealT>::print() {
  const auto* const delim = "----------- print() -----------";
  for (std::size_t i = 0; i < storages.size(); ++i) {
    std::cout << delim << std::endl;
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      std::cout << storages[i][j] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << delim << std::endl;
}

template class ODEVector<float>;
template class ODEVector<double>;

} // namespace seissol::ode
