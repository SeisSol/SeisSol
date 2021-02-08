#include <cassert>
#include <cmath>
#include <iostream>
#include <utility>
#include <tuple>
#include <vector>

#include "ODEVector.h"

namespace seissol::ode {
std::pair<std::size_t, std::size_t> ODEVector::index(std::size_t idx) const {
  for (std::size_t i = 0; i < storages.size(); ++i) {
    const auto begin = offsets[i];
    const auto end = begin + sizes[i];
    assert(begin <= idx);
    // TODO(Lukas) Verify this logic!
    if (idx < end) {
      return {i, idx - begin};
    }
  }
  assert(false); // Unreachable!
}

ODEVector::ODEVector(std::vector<real*> storages,
                     std::vector<std::size_t> sizes)
    : storages(std::move(storages)), sizes(std::move(sizes)) {
  std::size_t curOffset = 0;
  for (unsigned long size : this->sizes) {
    offsets.push_back(curOffset);
    curOffset += size;
  }
}

void ODEVector::updateStoragesAndSizes(std::vector<real*> newStorages, std::vector<std::size_t> newSizes) {
  storages = std::move(newStorages);
  sizes = std::move(newSizes);
  offsets.clear();
  std::size_t curOffset = 0;
  for (unsigned long size : this->sizes) {
    offsets.push_back(curOffset);
    curOffset += size;
  }
}

std::pair<real*, size_t> ODEVector::getSubvector(size_t storageIdx) {
  return {storages[storageIdx], sizes[storageIdx]};
}

real& ODEVector::operator[](std::size_t idx) {
  const auto idxPair = index(idx);
  return storages[idxPair.first][idxPair.second];
}

const real& ODEVector::operator[](std::size_t idx) const {
  const auto idxPair = index(idx);
  return storages[idxPair.first][idxPair.second];
}

ODEVector& ODEVector::operator+=(ODEVector& rhs) {
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == rhs.sizes[i]);
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      storages[i][j] += rhs.storages[i][j];
    }
  }
  return *this;
}

ODEVector& ODEVector::operator*=(real rhs) {
  for (std::size_t i = 0; i < storages.size(); ++i) {
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      storages[i][j] *= rhs;
    }
  }
  return *this;
}

void ODEVector::weightedAddInplace(real weight, const ODEVector& rhs) {
  if (weight == 0.0) return;
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == rhs.sizes[i]);
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      storages[i][j] += weight * rhs.storages[i][j];
    }
  }
}

ODEVector& ODEVector::operator=(const ODEVector& other) {
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == other.sizes[i]);
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      storages[i][j] = other.storages[i][j];
    }
  }
  return *this;
}

real ODEVector::normDifferenceTo(ODEVector& other, bool useLInfNorm) {
  // Computes the L2 or LInf norm of the difference between two vectors.
  real error = 0.0;
  real maxError = -1;
  for (std::size_t i = 0; i < storages.size(); ++i) {
    assert(sizes[i] == other.sizes[i]);
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      const real curDiff = storages[i][j] - other.storages[i][j];
      error += curDiff * curDiff;
      maxError = std::max(std::abs(curDiff), maxError);
    }
  }
  if (useLInfNorm) return maxError;
  return std::sqrt(error);
}

real ODEVector::norm() {
  // Computes the L2 norm.
  real norm = 0.0;
  for (std::size_t i = 0; i < storages.size(); ++i) {
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      norm += storages[i][j] * storages[i][j];
    }
  }
  return std::sqrt(norm);
}

void ODEVector::print() {
  const auto delim = "----------- print() -----------";
  for (std::size_t i = 0; i < storages.size(); ++i) {
    std::cout << delim << std::endl;
    for (std::size_t j = 0; j < sizes[i]; ++j) {
      std::cout << storages[i][j] << ", ";
    }
    std::cout << std::endl;
  }
  std::cout << delim << std::endl;
}

} // namespace seissol::ode
