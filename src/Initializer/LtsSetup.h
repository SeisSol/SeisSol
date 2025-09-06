// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_LTSSETUP_H_

#include <Common/Constants.h>
#include <cassert>
#include <cstdint>

namespace seissol {

class LtsSetup {
  private:
  constexpr static int IndexNeighborHasDerivatives = 0;
  constexpr static int IndexNeighborHasGTS = Cell::NumFaces;
  constexpr static int IndexBuffers = Cell::NumFaces * 2;
  constexpr static int IndexDerivatives = IndexBuffers + 1;
  constexpr static int IndexCache = IndexDerivatives + 1;

  public:
  LtsSetup() = default;

  LtsSetup(uint16_t data) : data(data) {}

  constexpr auto setNeighborHasDerivatives(int face, bool derivatives) -> LtsSetup& {
    assert(face < Cell::NumFaces);
    return set(face + IndexNeighborHasDerivatives, derivatives);
  }

  [[nodiscard]] constexpr auto neighborHasDerivatives(int face) const -> bool {
    assert(face < Cell::NumFaces);
    return test(face + IndexNeighborHasDerivatives);
  }

  constexpr auto setNeighborGTS(int face, bool gts) -> LtsSetup& {
    assert(face < Cell::NumFaces);
    return set(face + IndexNeighborHasGTS, gts);
  }

  [[nodiscard]] constexpr auto neighborGTS(int face) const -> bool {
    assert(face < Cell::NumFaces);
    return test(face + IndexNeighborHasGTS);
  }

  constexpr auto setHasBuffers(bool val) -> LtsSetup& { return set(IndexBuffers, val); }

  [[nodiscard]] constexpr auto hasBuffers() const -> bool { return test(IndexBuffers); }

  constexpr auto setHasDerivatives(bool val) -> LtsSetup& { return set(IndexDerivatives, val); }

  [[nodiscard]] constexpr auto hasDerivatives() const -> bool { return test(IndexDerivatives); }

  constexpr auto setCacheBuffers(bool val) -> LtsSetup& { return set(IndexCache, val); }

  [[nodiscard]] constexpr auto cacheBuffers() const -> bool { return test(IndexCache); }

  [[nodiscard]] constexpr auto test(int index) const -> bool { return (data & (1 << index)) != 0; }

  constexpr auto set(int index, bool value) -> LtsSetup& {
    if (value) {
      data |= 1 << index;
    } else {
      data &= ~(1 << index);
    }
    return *this;
  }

  [[nodiscard]] constexpr auto unwrap() const -> uint16_t { return data; }

  private:
  uint16_t data{0};
};

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
