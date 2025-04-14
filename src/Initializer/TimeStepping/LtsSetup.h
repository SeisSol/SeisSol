// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
#ifndef SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSSETUP_H_
#define SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSSETUP_H_

#include <bitset>
#include <cstdint>
namespace seissol::initializer {

class LtsSetup {
public:
    LtsSetup() = default;

    LtsSetup(uint16_t data) : data(data) {}

    auto setNeighborHasBuffers(int face, bool buffers) -> LtsSetup& {
        return set(face, buffers);
    }

    auto neighborHasBuffers(int face) const -> bool {
        return test(face);
    }

    auto setNeighborGTS(int face, bool gts) -> LtsSetup& {
        return set(face + 4, gts);
    }

    auto neighborGTS(int face) const -> bool {
        return test(face + 4);
    }

    auto setHasBuffers(bool val) -> LtsSetup& {
        return set(8, val);
    }

    auto hasBuffers() const -> bool {
        return test(8);
    }

    auto setHasDerivatives(bool val) -> LtsSetup& {
        return set(9, val);
    }

    auto hasDerivatives() const -> bool {
        return test(9);
    }

    auto setCacheBuffers(bool val) -> LtsSetup& {
        return set(10, val);
    }

    auto cacheBuffers() const -> bool {
        return test(10);
    }

    auto test(int index) const -> bool {
        return (data & (1 << index)) != 0;
    }

    auto set(int index, bool value) -> LtsSetup& {
        if (value) {
            data |= 1 << index;
        }
        else {
            data &= ~(1 << index);
        }
        return *this;
    }

    auto unwrap() const -> uint16_t {
        return data;
    }
private:
    uint16_t data;
};

} // namespace seissol::initializer
#endif // SEISSOL_SRC_INITIALIZER_TIMESTEPPING_LTSSETUP_H_
