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
#include <optional>

namespace seissol {

/**
  All possible buffer types a cell writes on a face between ADER-DG prediction and correction, for
  the neighbour across that face, and that are MPI-transferred where the neighbour is on another
  rank. Each is a face trace: what the predictor gave the cell, restricted to the face, in the face
  basis.
  */
enum class BufferType : int8_t {
  /// the trace of the time integral, over one step only (GTS-like, e.g. for clusterLocal ==
  /// clusterNeighbor)
  StepIntegrals = 0,

  /// the traces of the time derivatives, the full time evolution on the face (for e.g. DR,
  /// clusterLocal > clusterNeighbor)
  Derivatives = 1,

  /// the trace of the time integral, accumulated over multiple timesteps (for clusterLocal <
  /// clusterNeighbor)
  AccumulatedIntegrals = 2,
};

constexpr static std::uint32_t BufferCount = 3;
constexpr static std::uint32_t BufferCountBits = 2;

static_assert(1_U32 << BufferCountBits > BufferCount,
              "BufferCountBits leaves no value for a face without a buffer.");

/**
  Encapsules LTS-relevant data per face: what the neighbour across a face writes there for this
  cell, whether it runs at the same time step, and what this cell writes there for the neighbour.
 */
class LtsSetup {
  private:
  constexpr static std::uint32_t IndexNeighborBuffer = 0;
  constexpr static std::uint32_t IndexNeighborGTSRelation = BufferCountBits * Cell::NumFaces;
  constexpr static std::uint32_t IndexFaceBuffer = IndexNeighborGTSRelation + Cell::NumFaces;
  constexpr static std::uint32_t CountIndex = IndexFaceBuffer + BufferCountBits * Cell::NumFaces;
  /// the value of the bits of a face buffer where the cell writes none
  constexpr static std::uint32_t NoFaceBuffer = (1_U32 << BufferCountBits) - 1;

  public:
  using BitmapType = std::uint32_t;
  static_assert(CountIndex <= sizeof(BitmapType) * 8, "Capacity of the LtsSetup exceeded");

  constexpr LtsSetup() {
    for (std::uint32_t face = 0; face < Cell::NumFaces; ++face) {
      setFaceBuffer(face, std::nullopt);
    }
  }

  explicit constexpr LtsSetup(BitmapType data) : data_(data) {}

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

  constexpr auto setFaceBuffer(std::uint32_t face, std::optional<BufferType> type) -> LtsSetup& {
    assert(face < Cell::NumFaces);
    return setBits(face * BufferCountBits + IndexFaceBuffer,
                   type.has_value() ? static_cast<std::uint32_t>(type.value()) : NoFaceBuffer,
                   BufferCountBits);
  }

  /**
    The BufferType this cell writes on `face` for the neighbor across it, or none where nothing
    reads it there.
   */
  [[nodiscard]] constexpr auto faceBuffer(std::uint32_t face) const -> std::optional<BufferType> {
    assert(face < Cell::NumFaces);
    const auto bits = getBits(face * BufferCountBits + IndexFaceBuffer, BufferCountBits);
    if (bits == NoFaceBuffer) {
      return std::nullopt;
    }
    return static_cast<BufferType>(bits);
  }

  /**
    True, iff this cell writes the given buffer type on at least one of its faces.
   */
  [[nodiscard]] constexpr auto hasBuffer(BufferType type) const -> bool {
    for (std::uint32_t face = 0; face < Cell::NumFaces; ++face) {
      if (faceBuffer(face) == type) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] constexpr auto hasAnyBuffer() const -> bool {
    for (std::uint32_t face = 0; face < Cell::NumFaces; ++face) {
      if (faceBuffer(face).has_value()) {
        return true;
      }
    }
    return false;
  }

  [[nodiscard]] constexpr auto test(std::uint32_t index) const -> bool {
    return (data_ & (1_U32 << index)) != 0;
  }

  constexpr auto set(std::uint32_t index, bool value) -> LtsSetup& {
    if (value) {
      data_ |= 1_U32 << index;
    } else {
      data_ &= ~(1_U32 << index);
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
