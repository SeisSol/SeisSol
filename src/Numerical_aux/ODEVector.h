#ifndef SEISSOL_ODEVECTOR_H
#define SEISSOL_ODEVECTOR_H
#include <vector>
#include <cassert>
#include <tuple>

class ODEVector {
  using real = double;
  // A simple vector, patched together from multiple smaller vectors.
  std::vector<real*> storages;
  std::vector<std::size_t> sizes;
  std::vector<std::size_t> offsets;

  std::pair<std::size_t, std::size_t> index(std::size_t idx) const {
    for (std::size_t i = 0; i < storages.size(); ++i) {
      const auto begin = offsets[i];
      const auto end = begin + sizes[i];
      assert(begin < idx);
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
      : storages(storages){
    std::size_t curOffset = 0;
    for (std::size_t i = 0; i < sizes.size(); ++i) {
      offsets.push_back(curOffset);
      curOffset += sizes[i];
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
      for (std::size_t j = 0; j < sizes[i]; ++i) {
        storages[i][j] += rhs.storages[i][j];
      }
    }
    return *this;
  }

  ODEVector& operator*=(real rhs) {
    for (std::size_t i = 0; i < storages.size(); ++i) {
      for (std::size_t j = 0; j < sizes[i]; ++i) {
        storages[i][j] *= rhs;
      }
    }
    return *this;
  }

  ODEVector& operator=(const ODEVector& other) {
    for (std::size_t i = 0; i < storages.size(); ++i) {
      for (std::size_t j = 0; j < sizes[i]; ++i) {
        storages[i][j] = other.storages[i][j];
      }
    }
    return *this;
  }


  real normDifferenceTo(ODEVector& other) {
    // Computes the L2 norm difference.
    real error = 0.0;
    for (std::size_t i = 0; i < storages.size(); ++i) {
      assert(sizes[i] == other.sizes[i]);
      for (std::size_t j = 0; j < sizes[i]; ++i) {
        const double curDiff = storages[i][j] - other.storages[i][j];
        error += curDiff * curDiff;
      }
    }
    return error; // TODO(Lukas) Maybe use sqrt?

  }
};


#endif //SEISSOL_ODEVECTOR_H
