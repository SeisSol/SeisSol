// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_LTSSETUP_H_

#include "Common/Constants.h"
#include <cassert>
#include <cstdint>

namespace seissol {

/**
  Encapsules LTS-relevant data, i.e. if pre-time-integrated buffers
  or a full space-time evolution should be used.
 */
class LtsSetup {
  private:
  constexpr static std::uint32_t IndexNeighborHasDerivatives = 0;
  constexpr static std::uint32_t IndexNeighborGTSRelation = Cell::NumFaces;
  constexpr static std::uint32_t IndexBuffers = Cell::NumFaces * 2;
  constexpr static std::uint32_t IndexDerivatives = IndexBuffers + 1;
  constexpr static std::uint32_t IndexCache = IndexDerivatives + 1;
  constexpr static std::uint32_t CountIndex = IndexCache + 1;

  public:
  using BitmapType = std::uint16_t;
  static_assert(CountIndex <= sizeof(BitmapType) * 8, "Capacity of the LtsSetup exceeded");

  LtsSetup() = default;

  LtsSetup(BitmapType data) : data(data) {}

  constexpr auto setNeighborHasDerivatives(std::uint32_t face, bool derivatives) -> LtsSetup& {
    assert(face < Cell::NumFaces);
    return set(face + IndexNeighborHasDerivatives, derivatives);
  }

  [[nodiscard]] constexpr auto neighborHasDerivatives(std::uint32_t face) const -> bool {
    assert(face < Cell::NumFaces);
    return test(face + IndexNeighborHasDerivatives);
  }

  constexpr auto setNeighborGTSRelation(std::uint32_t face, bool gts) -> LtsSetup& {
    assert(face < Cell::NumFaces);
    return set(face + IndexNeighborGTSRelation, gts);
  }

  [[nodiscard]] constexpr auto neighborGTSRelation(std::uint32_t face) const -> bool {
    assert(face < Cell::NumFaces);
    return test(face + IndexNeighborGTSRelation);
  }

  constexpr auto setHasBuffers(bool val) -> LtsSetup& { return set(IndexBuffers, val); }

  [[nodiscard]] constexpr auto hasBuffers() const -> bool { return test(IndexBuffers); }

  constexpr auto setHasDerivatives(bool val) -> LtsSetup& { return set(IndexDerivatives, val); }

  [[nodiscard]] constexpr auto hasDerivatives() const -> bool { return test(IndexDerivatives); }

  constexpr auto setAccumulateBuffers(bool val) -> LtsSetup& { return set(IndexCache, val); }

  [[nodiscard]] constexpr auto accumulateBuffers() const -> bool { return test(IndexCache); }

  [[nodiscard]] constexpr auto test(std::uint32_t index) const -> bool {
    return (data & (1 << index)) != 0;
  }

  constexpr auto set(std::uint32_t index, bool value) -> LtsSetup& {
    if (value) {
      data |= 1 << index;
    } else {
      data &= ~(1 << index);
    }
    return *this;
  }

  [[nodiscard]] constexpr auto unwrap() const -> BitmapType { return data; }

  private:
  BitmapType data{0};
};

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
