#pragma once

#include "cellconfig.hpp"
#include "constants.hpp"
#include "Kernels/precision.hpp"
#include "Equations/datastructures.hpp"
#include <variant>

namespace seissol {

// This is it where it all comes together.

// TODO(David): auto-generate this very type.
/*using SupportedConfigs =
    std::variant<CellConfig<seissol::model::Material_t, real, ConvergenceOrder, true>,
                 CellConfig<seissol::model::Material_t, real, ConvergenceOrder, false>>;*/

using SupportedConfigs =
    std::variant<CellConfig<seissol::model::ElasticMaterial, real, ConvergenceOrder, true>,
                 CellConfig<seissol::model::ElasticMaterial, real, ConvergenceOrder, false>,
                 CellConfig<seissol::model::ViscoElasticMaterial<3>, real, ConvergenceOrder, true>,
                 CellConfig<seissol::model::ViscoElasticMaterial<3>, real, ConvergenceOrder, false>,
                 CellConfig<seissol::model::AnisotropicMaterial, real, ConvergenceOrder, true>,
                 CellConfig<seissol::model::AnisotropicMaterial, real, ConvergenceOrder, false>>;

} // namespace seissol
