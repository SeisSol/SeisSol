// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigDispatch.h"
#include "Common/ConfigLayout.h"
#include "Common/ConfigRegistry.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/runtime.h"
#include "GeneratedCode/tensor.h"
#include "GeneratedCode/variant.h"
#include "Solver/MultipleSimulations.h"

#include <cstddef>
#include <utility>
#include <variant>
#include <vector>

namespace seissol {

namespace {

template <std::size_t... Ids>
constexpr bool idsAreVariants(std::index_sequence<Ids...> /*ids*/) {
  return ((runtime::variantOf<std::variant_alternative_t<Ids, ConfigVariant>>() == Ids) && ...);
}

// The id of a configuration is also the variant of its kernels in runtime.h, so that whatever
// knows the configuration of its data by id calls the kernels of that configuration.
static_assert(std::variant_size_v<ConfigVariant> == runtime::VariantCount,
              "Every variant of the kernels in runtime.h belongs to one built configuration.");
static_assert(idsAreVariants(std::make_index_sequence<std::variant_size_v<ConfigVariant>>()),
              "The configurations are listed in the order of the variants of their kernels.");

/// The layout of the configuration the generated code of this translation unit belongs to.
ConfigLayout generatedLayout() {
  ConfigLayout layout;
  layout.config = Config::Value;

  layout.numQuantities = model::MaterialT::NumQuantities;
  layout.numBasisFunctions = tensor::Q<Config>::Shape[multisim::BasisFunctionDimension];
  layout.basisFunctionDimension = multisim::BasisFunctionDimension;
  layout.dofsSize = tensor::Q<Config>::size();

  layout.drNumPoints = dr::misc::NumBoundaryGaussPoints<Config>;
  layout.drNumPaddedPoints = dr::misc::NumPaddedPoints<Config>;
  layout.drNumQuantities = dr::misc::NumQuantities<Config>;
  layout.drNumTimePoints = dr::misc::TimeSteps<Config>;
  return layout;
}

} // namespace

const std::vector<ConfigLayout>& builtConfigLayouts() {
  static const std::vector<ConfigLayout> Layouts{generatedLayout()};
  return Layouts;
}

} // namespace seissol
