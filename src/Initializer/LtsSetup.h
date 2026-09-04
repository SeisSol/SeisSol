// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_LTSSETUP_H_

#include "Common/Constants.h"
#include "Common/Literals.h"

#include <cassert>
#include <cstdint>

namespace seissol {

/**
  All possible buffer types that are stored and MPI-transferred
  between ADER-DG prediction and correction.
  */
enum class BufferType : int8_t {
  /// the full space-time evolution (for e.g. DR, clusterLocal > clusterNeighbor)
  Derivatives = 0,

  /// time-integrated dofs, over one step only (GTS-like, e.g. for clusterLocal == clusterNeighbor)
  StepIntegrals = 1,

  /// time-integrated dofs, buffered over multiple timesteps (for clusterLocal < clusterNeighbor)
  AccumulatedIntegrals = 2,
};

constexpr static std::uint32_t BufferCount = 3;
constexpr static std::uint32_t BufferCountBits = 2;

static_assert(1_U32 << BufferCountBits >= BufferCount, "BufferCountBits is too small.");

/**
  Encapsules LTS-relevant data, i.e. if pre-time-integrated buffers
  or a full space-time evolution should be used.
 */
class LtsSetup {
  private:
  constexpr static std::uint32_t IndexNeighborBuffer = 0;
  constexpr static std::uint32_t IndexNeighborGTSRelation = BufferCountBits * Cell::NumFaces;
  constexpr static std::uint32_t IndexHasBuffer = IndexNeighborGTSRelation + Cell::NumFaces;
  constexpr static std::uint32_t CountIndex = IndexHasBuffer + BufferCount;

  public:
  using BitmapType = std::uint16_t;
  static_assert(CountIndex <= sizeof(BitmapType) * 8, "Capacity of the LtsSetup exceeded");

  LtsSetup() = default;

  explicit LtsSetup(BitmapType data) : data_(data) {}

  constexpr auto setNeighborBuffer(std::uint32_t face, BufferType type) -> LtsSetup& {
    assert(face < Cell::NumFaces);
    return setBits(face * BufferCountBits + IndexNeighborBuffer,
                   static_cast<std::uint32_t>(type),
                   BufferCountBits);
  }

  /**
    Return the BufferType that is supplied from this face neighbor.
   */
  [[nodiscard]] constexpr auto neighborBuffer(std::uint32_t face) const -> BufferType {
    assert(face < Cell::NumFaces);
    return static_cast<BufferType>(
        getBits(face * BufferCountBits + IndexNeighborBuffer, BufferCountBits));
  }

  constexpr auto setNeighborGTSRelation(std::uint32_t face, bool gts) -> LtsSetup& {
    assert(face < Cell::NumFaces);
    return set(face + IndexNeighborGTSRelation, gts);
  }

  /**
    Iff true, we have the same timestep as our neighbor over `face`.
   */
  [[nodiscard]] constexpr auto neighborGTSRelation(std::uint32_t face) const -> bool {
    assert(face < Cell::NumFaces);
    return test(face + IndexNeighborGTSRelation);
  }

  constexpr auto setHasBuffer(bool val, BufferType type) -> LtsSetup& {
    return set(IndexHasBuffer + static_cast<std::uint32_t>(type), val);
  }

  /**
    True, iff this cell pre-computes and stores the given buffer type.
   */
  [[nodiscard]] constexpr auto hasBuffer(BufferType type) const -> bool {
    return test(IndexHasBuffer + static_cast<std::uint32_t>(type));
  }

  [[nodiscard]] constexpr auto test(std::uint32_t index) const -> bool {
    return (data_ & (1 << index)) != 0;
  }

  constexpr auto set(std::uint32_t index, bool value) -> LtsSetup& {
    if (value) {
      data_ |= 1 << index;
    } else {
      data_ &= ~(1 << index);
    }
    return *this;
  }

  [[nodiscard]] constexpr auto getBits(std::uint32_t index, std::uint32_t size) const
      -> std::uint32_t {
    const auto mask = (1_U32 << size) - 1;
    return (data_ >> index) & mask;
  }

  constexpr auto setBits(std::uint32_t index, std::uint32_t bits, std::uint32_t size) -> LtsSetup& {
    const auto mask = (1_U32 << size) - 1;
    data_ &= ~(mask << index);
    data_ |= (bits & mask) << index;
    return *this;
  }

  [[nodiscard]] constexpr auto unwrap() const -> BitmapType { return data_; }

  private:
  BitmapType data_{0};
};

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
