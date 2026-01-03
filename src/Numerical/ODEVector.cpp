// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ODEVector.h"

#include "Kernels/Precision.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstddef>
#include <cstdlib>
#include <iostream>
#include <utility>
#include <vector>

namespace seissol::ode {
std::pair<std::size_t, std::size_t> ODEVector::index(std::size_t idx) const {
  for (std::size_t i = 0; i < storages_.size(); ++i) {
    const auto begin = offsets_[i];
    const auto end = begin + sizes_[i];
    assert(begin <= idx);
    if (idx < end) {
      return {i, idx - begin};
    }
  }
  std::abort(); // Unreachable!
}

ODEVector::ODEVector(std::vector<real*> storages, std::vector<std::size_t> sizes)
    : storages_(std::move(storages)), sizes_(std::move(sizes)) {
  std::size_t curOffset = 0;
  for (const unsigned long size : this->sizes_) {
    offsets_.push_back(curOffset);
    curOffset += size;
  }
}

void ODEVector::updateStoragesAndSizes(std::vector<real*> newStorages,
                                       std::vector<std::size_t> newSizes) {
  storages_ = std::move(newStorages);
  sizes_ = std::move(newSizes);
  offsets_.clear();
  std::size_t curOffset = 0;
  for (const unsigned long size : this->sizes_) {
    offsets_.push_back(curOffset);
    curOffset += size;
  }
}

std::pair<real*, size_t> ODEVector::getSubvector(size_t storageIdx) {
  return {storages_[storageIdx], sizes_[storageIdx]};
}

real& ODEVector::operator[](std::size_t idx) {
  const auto idxPair = index(idx);
  return storages_[idxPair.first][idxPair.second];
}

const real& ODEVector::operator[](std::size_t idx) const {
  const auto idxPair = index(idx);
  return storages_[idxPair.first][idxPair.second];
}

ODEVector& ODEVector::operator+=(ODEVector& rhs) {
  for (std::size_t i = 0; i < storages_.size(); ++i) {
    assert(sizes[i] == rhs.sizes[i]);
#pragma omp simd
    for (std::size_t j = 0; j < sizes_[i]; ++j) {
      storages_[i][j] += rhs.storages_[i][j];
    }
  }
  return *this;
}

ODEVector& ODEVector::operator*=(real scalar) {
  for (std::size_t i = 0; i < storages_.size(); ++i) {
#pragma omp simd
    for (std::size_t j = 0; j < sizes_[i]; ++j) {
      storages_[i][j] *= scalar;
    }
  }
  return *this;
}

ODEVector& ODEVector::copyFrom(const ODEVector& other) {
  for (std::size_t i = 0; i < storages_.size(); ++i) {
    assert(sizes[i] == other.sizes[i]);
    std::copy_n(other.storages_[i], sizes_[i], storages_[i]);
  }
  return *this;
}

void ODEVector::weightedAddInplace(real weight, const ODEVector& rhs) {
  if (weight == 0.0) {
    return;
  }
  for (std::size_t i = 0; i < storages_.size(); ++i) {
    assert(sizes[i] == rhs.sizes[i]);
#pragma omp simd
    for (std::size_t j = 0; j < sizes_[i]; ++j) {
      storages_[i][j] += weight * rhs.storages_[i][j];
    }
  }
}

real ODEVector::normDifferenceTo(ODEVector& other, bool useLInfNorm) {
  // Computes the L2 or LInf norm of the difference between two vectors.
  real error = 0.0;
  real maxError = -1;
  for (std::size_t i = 0; i < storages_.size(); ++i) {
    assert(sizes[i] == other.sizes[i]);
    for (std::size_t j = 0; j < sizes_[i]; ++j) {
      const real curDiff = storages_[i][j] - other.storages_[i][j];
      error += curDiff * curDiff;
      maxError = std::max(std::abs(curDiff), maxError);
    }
  }
  if (useLInfNorm) {
    return maxError;
  }
  return std::sqrt(error);
}

real ODEVector::l2Norm() {
  // Computes the L2 norm.
  real norm = 0.0;
  for (std::size_t i = 0; i < storages_.size(); ++i) {
    for (std::size_t j = 0; j < sizes_[i]; ++j) {
      norm += storages_[i][j] * storages_[i][j];
    }
  }
  return std::sqrt(norm);
}

void ODEVector::print() {
  const auto* const delim = "----------- print() -----------";
  for (std::size_t i = 0; i < storages_.size(); ++i) {
    std::cout << delim << std::endl;
    for (std::size_t j = 0; j < sizes_[i]; ++j) {
      std::cout << storages_[i][j] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << delim << std::endl;
}

} // namespace seissol::ode
