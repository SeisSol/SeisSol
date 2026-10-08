// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/quantities.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/SolverSelector.h"
#include "Model/Quantities.h"
#include "TestConfigs.h"

#include <cstddef>

namespace seissol::unit_test {

namespace materialinvariants {

/// The mechanism-dependent assertions live in a template so that `if
/// constexpr` actually discards the branch that does not apply; in a plain
/// function both branches are still evaluated and every static_assert fires.
template <typename Cfg>
void checkRelaxation() {
  using MaterialT = model::MaterialOf<Cfg>;
  static_assert(MaterialT::Mechanisms == Cfg::RelaxationMechanisms,
                "the material and the configuration disagree on the mechanism count");
  if constexpr (MaterialT::Mechanisms > 0) {
    static_assert(MaterialT::NumberPerMechanism > 0,
                  "a mechanism that occupies nothing cannot carry a memory variable");
    static_assert(MaterialT::NumQuantities ==
                      MaterialT::NumElasticQuantities +
                          MaterialT::Mechanisms * MaterialT::NumberPerMechanism,
                  "the quantity count has to account for every mechanism");
  } else {
    static_assert(MaterialT::NumQuantities == MaterialT::NumElasticQuantities);
  }
}

/// Only a material with mechanisms has relaxation frequencies to check; the
/// loop cannot even be written for one without them.
template <typename MaterialT>
void checkRelaxationIsInert(const MaterialT& material) {
  if constexpr (MaterialT::Mechanisms > 0) {
    for (std::size_t mech = 0; mech < MaterialT::Mechanisms; ++mech) {
      CHECK(material.omega[mech] == 0.0);
    }
  }
}

} // namespace materialinvariants

/// Whatever a configuration is, these have to hold for its material. Checking
/// them for the material of every configuration built into the executable
/// rather than a fixed one means every configuration the CI builds gets them.
TEST_CASE_TEMPLATE("Configured material is self consistent" * doctest::test_suite("equations"),
                   Cfg,
                   SEISSOL_CONFIG_TYPES) {
  using namespace materialinvariants;

  using MaterialT = model::MaterialOf<Cfg>;
  using SolverT = kernels::SolverOf<Cfg>;

  SUBCASE("the quantity groups account for the whole layout") {
    static_assert(model::totalExtent(MaterialT::PrimaryGroups) <= MaterialT::NumQuantities,
                  "the primary groups cannot cover more than the material has");
    static_assert(model::quantitiesWellFormed(MaterialT::template RotationGroups<SolverT>,
                                              tensor::T<Cfg>::Shape[0]));
    static_assert(model::quantitiesWellFormed(MaterialT::template InverseRotationGroups<SolverT>,
                                              tensor::Tinv<Cfg>::Shape[0]));
  }

  SUBCASE("the face roles are where the rest of the code expects them") {
    // Everything that reaches for a velocity component goes through this
    // offset, and everything that counts stress components through the other.
    static_assert(MaterialT::VelocityOffset ==
                  model::roleOffset(MaterialT::PrimaryGroups, model::FaceRole::Velocity));
    static_assert(MaterialT::TractionComponents ==
                  model::roleExtent(MaterialT::PrimaryGroups, model::FaceRole::Traction));
    static_assert(MaterialT::TractionComponents > 0, "a material needs a traction");
    static_assert(MaterialT::VelocityOffset + 3 <= MaterialT::NumQuantities,
                  "the three velocity components have to fit");
  }

  SUBCASE("the traction comes before the velocities") {
    // Both readings of the old single constant -- the velocity offset and the
    // traction component count -- were the same number only because of this.
    // Several call sites still rely on it.
    static_assert(MaterialT::TractionComponents == MaterialT::VelocityOffset);
  }

  SUBCASE("relaxation is configured consistently") { checkRelaxation<Cfg>(); }

  SUBCASE("the stiff source rows are inside the material") {
    static_assert(MaterialT::StiffSourceRows.size() ==
                      generated::Quantities<Cfg>::StiffSourceRowCount,
                  "the material and the generated kernels disagree on the stiff rows");
    for (const auto& row : MaterialT::StiffSourceRows) {
      CHECK(row.quantity < MaterialT::NumQuantities);
      CHECK(row.target < MaterialT::NumQuantities);
      // The predictor works down from the last quantity and substitutes into
      // rows it has not reached yet, so a row may only feed a lower one.
      CHECK(row.target < row.quantity);
    }
  }

  SUBCASE("a default material is inert") {
    const MaterialT material{};
    CHECK(material.rho == 0.0);
    checkRelaxationIsInert(material);
  }
}

} // namespace seissol::unit_test
