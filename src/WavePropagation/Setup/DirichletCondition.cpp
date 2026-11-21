// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#include "DirichletCondition.h"

#include "Config.h"
#include "Core/Common/MultipleSimulations.h"
#include "Core/Common/Real.h"
#include "Core/Data/DataReader.h"
#include "Core/Data/DataTable.h"
#include "GeneratedCode/init.h"
#include "WavePropagation/Physics/Datastructures.h"

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_set>
#include <utility>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {

DirichletCondition::DirichletCondition(std::unique_ptr<data::DataReader> model)
    : model_(std::move(model)) {}

template <typename Cfg>
BoundaryFrame DirichletCondition::query(const double* barycenter,
                                        Real<Cfg>* mapTermsData,
                                        Real<Cfg>* constantTermsData,
                                        const std::vector<std::size_t>& simulations) const {
  using real = Real<Cfg>;
  if (model_ == nullptr) {
    logError() << "Model for easi-provided boundary is not initialized.";
  }
  assert(mapTermsData != nullptr);
  assert(constantTermsData != nullptr);
  assert(!simulations.empty());

  // The boundary condition is constant over the face, so it is sampled at the
  // face barycenter.
  data::DataTable table(1);
  table.bindViewConst("x", data::Direction::In, barycenter, 3, 0);
  table.bindViewConst("y", data::Direction::In, barycenter, 3, 1);
  table.bindViewConst("z", data::Direction::In, barycenter, 3, 2);
  table.bindConstant("group", std::int32_t{1});

  const auto supplied = data::suppliedParameters(*model_);

  // The ghost cell state is an affine function of the interior state, given in
  // global coordinates: q_ghost = A q_inside + b. The entries of A are named
  // map_{to}_{from}, those of b const_{to}, where the quantity names are the
  // ones of the material at hand. Mirroring the x velocity at the ghost cell is
  // therefore map_v1_v1: -1.
  const auto& varNames = model::MaterialOf<Cfg>::Quantities;

  auto mapTerms = init::dirichletMapGlobal<Cfg>::view::create(mapTermsData);
  auto constantTerms = init::dirichletOffsetGlobal<Cfg>::view::create(constantTermsData);

  std::unordered_set<std::string> known;

  // a model supplies numbers, so the frame is stated as one: 0 for global, 1 for face-aligned.
  real frame = 0.0;
  known.insert("frame");
  if (supplied.contains("frame")) {
    table.bindView("frame", data::Direction::Out, &frame);
  }

  for (size_t i = 0; i < varNames.size(); ++i) {
    const auto termName = std::string{"const_"} + varNames[i];
    known.insert(termName);
    auto& term = multisim::multisimWrap<Cfg>(constantTerms, simulations.front(), i);
    if (supplied.count(termName) > 0) {
      table.bindView(termName, data::Direction::Out, &term);
    } else {
      term = 0.0;
    }
  }
  for (size_t i = 0; i < varNames.size(); ++i) {
    for (size_t j = 0; j < varNames.size(); ++j) {
      auto termName = std::string{"map_"};
      termName += varNames[i];
      termName += '_';
      termName += varNames[j];
      known.insert(termName);
      if (supplied.contains(termName)) {
        table.bindView(termName, data::Direction::Out, &mapTerms(i, j));
      } else {
        // Default: Extrapolate
        mapTerms(i, j) = (i == j) ? 1.0 : 0.0;
      }
    }
  }

  for (const auto& termName : supplied) {
    if (!known.contains(termName)) {
      std::ostringstream valid;
      for (size_t i = 0; i < varNames.size(); ++i) {
        valid << (i == 0 ? "" : ", ") << varNames[i];
      }
      logError() << "The boundary condition file supplies" << termName
                 << "which is not a term of the boundary condition. Terms are named"
                 << "map_{to}_{from} and const_{to}, where both quantity names are one of:"
                 << valid.str() << ".";
    }
  }

  data::evaluate(*model_, table, "Dirichlet BC data");

  if (frame != 0.0 && frame != 1.0) {
    logError() << "The boundary condition file supplies a frame of" << frame
               << "-- it has to be 0 for a condition stated in global coordinates, or 1 for one "
                  "stated in the face-aligned basis.";
  }

  // every fused simulation that takes the condition gets its constant terms
  for (const auto simulation : simulations) {
    for (size_t i = 0; i < varNames.size(); ++i) {
      multisim::multisimWrap<Cfg>(constantTerms, simulation, i) =
          multisim::multisimWrap<Cfg>(constantTerms, simulations.front(), i);
    }
  }

  return frame == 0.0 ? BoundaryFrame::Global : BoundaryFrame::FaceAligned;
}

#define SEISSOL_CONFIG_INSTANTIATE(Cfg)                                                            \
  template BoundaryFrame DirichletCondition::query<Cfg>(                                           \
      const double*, Real<Cfg>*, Real<Cfg>*, const std::vector<std::size_t>&) const;
SEISSOL_FOR_EACH_CONFIG(SEISSOL_CONFIG_INSTANTIATE)
#undef SEISSOL_CONFIG_INSTANTIATE

} // namespace seissol::initializer
