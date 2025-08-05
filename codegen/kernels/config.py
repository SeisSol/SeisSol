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
    static constexpr ViscoImplementation ViscoMode = ViscoImplementation::{viscomode};
    static constexpr DRQuadRuleType DRQuadRule = DRQuadRuleType::{drquadrule};
    static constexpr std::size_t NumSimulations = {numsims};
}};
    """.format(
        **config, id=index
    )


def make_configfile(configs):
    configcode = "".join(make_config(i, config) for i, config in enumerate(configs))
    confignames = ", ".join(config["name"] for config in configs)

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

}} // namespace seissol

#endif
"""


def make_configincludefile(configs):
    configincludes = "\n".join(f"_H_({config['name']})" for config in configs)

    return f"""
// (include guard knowlingly omitted)

// Config instantiation file.
// Needed due to C++ not allowing to define templates in source files.
// I know, I'm not a fan of it as well.

// Usage: define _H_(x) in our source file with the wanted value,
// then include this file

{configincludes}

#undef _H_
"""
