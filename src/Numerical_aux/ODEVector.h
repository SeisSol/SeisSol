#ifndef SEISSOL_ODEVECTOR_H
#define SEISSOL_ODEVECTOR_H
#include <utility>
#include <vector>
#include <cassert>
#include <tuple>

class ODEVector {
  using real = double;
  // A simple vector, patched together from multiple smaller vectors.
  std::vector<real*> storages;
  std::vector<std::size_t> sizes;
  std::vector<std::size_t> offsets;

  [[nodiscard]] std::pair<std::size_t, std::size_t> index(std::size_t idx) const {
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
public:
  ODEVector(std::vector<real*> storages,
            std::vector<std::size_t> sizes)
      : storages(std::move(storages)), sizes(std::move(sizes)) {
    std::size_t curOffset = 0;
    for (unsigned long size : this->sizes) {
      offsets.push_back(curOffset);
      curOffset += size;
    }
  }

  std::pair<real*, size_t> getSubvector(size_t storageIdx) {
    return {storages[storageIdx], sizes[storageIdx]};
  }

  real& operator[](std::size_t idx) {
    const auto idxPair = index(idx);
    return storages[idxPair.first][idxPair.second];
  }
  const real& operator[](std::size_t idx) const {
    const auto idxPair = index(idx);
    return storages[idxPair.first][idxPair.second];
  }

  ODEVector& operator+=(ODEVector& rhs) {
    for (std::size_t i = 0; i < storages.size(); ++i) {
      assert(sizes[i] == rhs.sizes[i]);
      for (std::size_t j = 0; j < sizes[i]; ++j) {
        storages[i][j] += rhs.storages[i][j];
      }
    }
    return *this;
  }

  ODEVector& operator*=(real rhs) {
    for (std::size_t i = 0; i < storages.size(); ++i) {
      for (std::size_t j = 0; j < sizes[i]; ++j) {
        storages[i][j] *= rhs;
      }
    }
    return *this;
  }

  ODEVector& operator=(const ODEVector& other) {
    for (std::size_t i = 0; i < storages.size(); ++i) {
      assert(sizes[i] == other.sizes[i]);
      for (std::size_t j = 0; j < sizes[i]; ++j) {
        storages[i][j] = other.storages[i][j];
      }
    }
    return *this;
  }


  real normDifferenceTo(ODEVector& other, bool useLInfNorm = true) {
    // Computes the L2 or LInf norm of the difference between two vectors.
    real error = 0.0;
    real maxError = -1;
    for (std::size_t i = 0; i < storages.size(); ++i) {
      assert(sizes[i] == other.sizes[i]);
      for (std::size_t j = 0; j < sizes[i]; ++j) {
        const double curDiff = storages[i][j] - other.storages[i][j];
        error += curDiff * curDiff;
        maxError = std::max(std::abs(curDiff), maxError);
      }
    }
    if (useLInfNorm) return maxError;
    return std::sqrt(error);
  }

  real norm() {
    // Computes the L2 norm.
    real norm = 0.0;
    for (std::size_t i = 0; i < storages.size(); ++i) {
      for (std::size_t j = 0; j < sizes[i]; ++j) {
        norm += storages[i][j] * storages[i][j];
      }
    }
    return std::sqrt(norm);
  }

  void print() {
    for (std::size_t i = 0; i < storages.size(); ++i) {
      std::cout << "-------------------------------------" << std::endl;
      for (std::size_t j = 0; j < sizes[i]; ++j) {
        std::cout << storages[i][j] << ", ";
      }
      std::cout << std::endl;
    }
  }
};


#endif //SEISSOL_ODEVECTOR_H
