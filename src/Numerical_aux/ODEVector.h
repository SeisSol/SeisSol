#ifndef SEISSOL_ODEVECTOR_H
#define SEISSOL_ODEVECTOR_H

#include "Kernels/precision.hpp"

namespace seissol::ode {

/**
 * A simple vector, patched together from multiple smaller vectors.
 * Can be used as data type for different ODE solvers.
 * It does not allocate but rather uses pre-allocated storages.
 * For examples, see unit test ODEInt.t.h.
 */
class ODEVector {
  std::vector<real*> storages{};
  std::vector<std::size_t> sizes{};
  std::vector<std::size_t> offsets{};

  [[nodiscard]] std::pair<std::size_t, std::size_t> index(std::size_t idx) const;

public:
  ODEVector() = default;

  /**
   *
   * @param storages a vector of pointers to arrays of real which are used as storage
   * @param sizes a vector is size_t that gives the size of the storages
   */
  ODEVector(std::vector<real*> storages,
            std::vector<std::size_t> sizes);


  /**
   * Updates storages and sizes
   * @param newStorages
   * @param newSizes
   */
  void updateStoragesAndSizes(std::vector<real*> newStorages, std::vector<std::size_t> newSizes);

  /**
   * @param storageIdx an index of the subarray
   * @return pair of storage pointer and size
   */
  std::pair<real*, size_t> getSubvector(size_t storageIdx);

  /**
   * @param idx an array index (of combined vector)
   * @return reference to entry at index idx (considering all storages)
   */
  real& operator[](std::size_t idx);

  /**
   * @param idx an array index (of combined vector)
   * @return const reference to entry at index idx (considering all storages)
   */
  const real& operator[](std::size_t idx) const;

  /**
   * Computes this += rhs in-place.
   * @param rhs
   * @return reference to new ODEVector
   */
  ODEVector& operator+=(ODEVector& rhs);

  /**
   * Computes this *= scalar in-place.
   * @param scalar a real scaling factor
   * @return reference to new ODEVector
   */
  ODEVector& operator*=(real scalar);

  //
  /**
   * Computes this += weight * rhs inplace
   * @param weight a scalar factor
   * @param rhs a reference to another ODEVector
   */
  void weightedAddInplace(real weight, const ODEVector& rhs);

  /**
   * Sets ODEVector to values from other vector.
   * Warning: Only shallow copy, points to same storages
   * @param other
   * @return a reference to updated ODEVEctor
   */
  ODEVector& operator=(const ODEVector& other);

  /**
   *
   * @param other
   * @param useLInfNorm
   * @return returns either the L2 difference between this and other
   * or the LInf difference if useLInfNorm is true
   */
  real normDifferenceTo(ODEVector& other, bool useLInfNorm = true);

  /**
   * @return the L2 norm of this
   */
  real l2Norm();

  /**
   * Prints out entries of this. For debugging.
   */
  void print();
};

} // namespace seissol::ode

#endif //SEISSOL_ODEVECTOR_H
