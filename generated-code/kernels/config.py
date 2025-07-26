# SPDX-FileCopyrightText: 2025 SeisSol Group
#
# SPDX-License-Identifier: BSD-3-Clause
# SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
#
# SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff


def make_config(index, config):
    return """
struct Config{id} {{
    static constexpr std::size_t ConvergenceOrder = {order};
    static constexpr std::size_t RelaxationMechanisms = {mechanisms};
    static constexpr model::MaterialType MaterialType = model::MaterialType::{equation};
    static constexpr RealType Precision = RealType::{precision};
    static constexpr ViscoImplementation ViscoMode = ViscoMode::{viscomode};
    static constexpr DRQuadRuleType DRQuadRule = DRQuadRuleType::{drquadrule};
    static constexpr std::size_t NumSimulations = {numsims};
}};
    """.format(
        **config, id=index
    )


def make_configfile(configs):
    configcode = "".join(make_config(i, config) for i, config in enumerate(configs))
    confignames = ", ".join(f"Config{i}" for i, _ in enumerate(configs))

    return f"""
#ifndef SEISSOL_CONFIGS_H_
#define SEISSOL_CONFIGS_H_

#include <Common/Real.h>
#include <Common/Typedefs.h>
#include <Model/CommonDatastructures.h>
#include <cstddef>
#include <variant>

namespace seissol {{

{configcode}

using ConfigVariant = std::variant<{confignames}>;

// for now
using Config = Config0;

}} // namespace seissol

#endif
"""
