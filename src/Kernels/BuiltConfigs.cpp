// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/ConfigLayout.h"
#include "Common/ConfigRegistry.h"
#include "Config.h"
#include "DynamicRupture/Misc.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "Solver/MultipleSimulations.h"

#include <vector>

namespace seissol {

namespace {

/// The layout of the configuration the generated code of this translation unit belongs to.
ConfigLayout generatedLayout() {
  ConfigLayout layout;
  layout.config = Config::Value;

  layout.numQuantities = model::MaterialT::NumQuantities;
  layout.numBasisFunctions = tensor::Q::Shape[multisim::BasisFunctionDimension];
  layout.basisFunctionDimension = multisim::BasisFunctionDimension;
  layout.dofsSize = tensor::Q::size();

  layout.drNumPoints = dr::misc::NumBoundaryGaussPoints;
  layout.drNumPaddedPoints = dr::misc::NumPaddedPoints;
  layout.drNumQuantities = dr::misc::NumQuantities;
  layout.drNumTimePoints = dr::misc::TimeSteps;
  return layout;
}

} // namespace

const std::vector<ConfigLayout>& builtConfigLayouts() {
  static const std::vector<ConfigLayout> Layouts{generatedLayout()};
  return Layouts;
}

} // namespace seissol
