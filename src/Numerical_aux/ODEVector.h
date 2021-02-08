#ifndef SEISSOL_ODEVECTOR_H
#define SEISSOL_ODEVECTOR_H

#include "Kernels/precision.hpp"

namespace seissol::ode {

class ODEVector {
  // A simple vector, patched together from multiple smaller vectors.
  std::vector<real*> storages{};
  std::vector<std::size_t> sizes{};
  std::vector<std::size_t> offsets{};

  [[nodiscard]] std::pair<std::size_t, std::size_t> index(std::size_t idx) const;

public:
  ODEVector() = default;

  ODEVector(std::vector<real*> storages,
            std::vector<std::size_t> sizes);


  void updateStoragesAndSizes(std::vector<real*> newStorages, std::vector<std::size_t> newSizes);

  std::pair<real*, size_t> getSubvector(size_t storageIdx);

  real& operator[](std::size_t idx);

  const real& operator[](std::size_t idx) const;

  ODEVector& operator+=(ODEVector& rhs);

  ODEVector& operator*=(real rhs);

  // this += weight * rhs
  void weightedAddInplace(real weight, const ODEVector& rhs);

  ODEVector& operator=(const ODEVector& other);

  real normDifferenceTo(ODEVector& other, bool useLInfNorm = true);

  real norm();

  void print();
};

} // namespace seissol::ode

#endif //SEISSOL_ODEVECTOR_H
