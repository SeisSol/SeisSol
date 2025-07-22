// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_LTSSETUP_H_

#include <cstdint>

namespace seissol {

class LtsSetup {
public:
    LtsSetup() = default;

    LtsSetup(uint16_t data) : data(data) {}

    constexpr auto setNeighborHasDerivatives(int face, bool derivatives) -> LtsSetup& {
        return set(face, derivatives);
    }

    [[nodiscard]] constexpr auto neighborHasDerivatives(int face) const -> bool {
        return test(face);
    }

    constexpr auto setNeighborGTS(int face, bool gts) -> LtsSetup& {
        return set(face + 4, gts);
    }

    [[nodiscard]] constexpr auto neighborGTS(int face) const -> bool {
        return test(face + 4);
    }

    constexpr auto setHasBuffers(bool val) -> LtsSetup& {
        return set(8, val);
    }

    [[nodiscard]] constexpr auto hasBuffers() const -> bool {
        return test(8);
    }

    constexpr auto setHasDerivatives(bool val) -> LtsSetup& {
        return set(9, val);
    }

    [[nodiscard]] constexpr auto hasDerivatives() const -> bool {
        return test(9);
    }

    constexpr auto setCacheBuffers(bool val) -> LtsSetup& {
        return set(10, val);
    }

    [[nodiscard]] constexpr auto cacheBuffers() const -> bool {
        return test(10);
    }

    [[nodiscard]] constexpr auto test(int index) const -> bool {
        return (data & (1 << index)) != 0;
    }

    constexpr auto set(int index, bool value) -> LtsSetup& {
        if (value) {
            data |= 1 << index;
        }
        else {
            data &= ~(1 << index);
        }
        return *this;
    }

    [[nodiscard]] constexpr auto unwrap() const -> uint16_t {
        return data;
    }
private:
    uint16_t data;
};

} // namespace seissol
#endif // SEISSOL_SRC_INITIALIZER_LTSSETUP_H_
