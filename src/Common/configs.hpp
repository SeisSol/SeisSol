#pragma once

#include "cellconfig.hpp"
#include "constants.hpp"
#include "Kernels/precision.hpp"
#include <variant>

namespace seissol {

// TODO(David): forward-declare all models
namespace model {
class ElasticMaterial;
template <std::size_t>
class ViscoElasticMaterial;
class PoroElasticMaterial;
class AnisotropicMaterial;
} // namespace model

// This is it where it all comes together.

// TODO(David): auto-generate this very type.
/*using SupportedConfigs =
    std::variant<CellConfig<seissol::model::Material_t, real, ConvergenceOrder, true>,
                 CellConfig<seissol::model::Material_t, real, ConvergenceOrder, false>>;*/

using SupportedConfigs = std::variant<
    CellConfig<seissol::model::ElasticMaterial, real, GivenConvergenceOrder, true>,
    CellConfig<seissol::model::ElasticMaterial, real, GivenConvergenceOrder, false>,
    CellConfig<seissol::model::ViscoElasticMaterial<3>, real, GivenConvergenceOrder, true>,
    CellConfig<seissol::model::ViscoElasticMaterial<3>, real, GivenConvergenceOrder, false>,
    CellConfig<seissol::model::AnisotropicMaterial, real, GivenConvergenceOrder, true>,
    CellConfig<seissol::model::AnisotropicMaterial, real, GivenConvergenceOrder, false>>;

} // namespace seissol
