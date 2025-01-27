// SPDX-FileCopyrightText: 2015-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff
// SPDX-FileContributor: Sebastian Wolf

#ifndef SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_DATASTRUCTURES_H_
#define SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_DATASTRUCTURES_H_

#include "Common/Constants.h"
#include "Equations/elastic/Model/Datastructures.h"
#include "Initializer/PreProcessorMacros.h"
#include "Model/CommonDatastructures.h"
#include "generated_code/tensor.h"
#include <array>
#include <cstddef>
#include <string>

namespace seissol::model {
class ViscoElasticLocalData;
class ViscoElasticNeighborData;

template <std::size_t MechanismsP>
struct ViscoElasticMaterialParametrized : public ElasticMaterial {
  static constexpr std::size_t NumberPerMechanism = 6;
  static constexpr std::size_t NumElasticQuantities = 9;
  static constexpr std::size_t NumQuantities =
      NumElasticQuantities + MechanismsP * NumberPerMechanism;
  static constexpr std::size_t Mechanisms = MechanismsP;
  static constexpr MaterialType Type = MaterialType::Viscoelastic;
  static constexpr LocalSolver Solver = LocalSolver::CauchyKovalevskiAnelastic;
  static inline const std::string Text = "viscoelastic-" + std::to_string(MechanismsP);
  static inline const std::array<std::string, NumElasticQuantities> Quantities{
      "s_xx", "s_yy", "s_zz", "s_xy", "s_yz", "s_xz", "v1", "v2", "v3"};

  using LocalSpecificData = ViscoElasticLocalData;
  using NeighborSpecificData = ViscoElasticNeighborData;

  //! Relaxation frequencies
  double omega[zeroLengthArrayHandler(Mechanisms)];
  /** Entries of the source matrix (E)
   * theta[0] = -(lambda * Y_lambda + 2.0 * mu * Y_mu)
   * theta[1] = -lambda * Y_lambda
   * theta[2] = -2.0 * mu * Y_mu
   **/
  double theta[zeroLengthArrayHandler(Mechanisms)][3];
  double Qp;
  double Qs;

  ViscoElasticMaterialParametrized() = default;
  ViscoElasticMaterialParametrized(const std::vector<double>& materialValues)
      : ElasticMaterial(materialValues) {
    for (int mech = 0; mech < Mechanisms; ++mech) {
      this->omega[mech] = materialValues.at(3 + 4 * mech);
      for (unsigned i = 1; i < 4; ++i) {
        this->theta[mech][i - 1] = materialValues.at(3 + 4 * mech + i);
      }
    }
    // This constructor is used to initialize a ViscoElasticMaterial
    // from the values in Fortran. Qp and Qs are not part of the
    // material in Fortran, so we set these to NaN.
    Qp = std::numeric_limits<double>::signaling_NaN();
    Qs = std::numeric_limits<double>::signaling_NaN();
  }

  ~ViscoElasticMaterialParametrized() override = default;

  [[nodiscard]] MaterialType getMaterialType() const override { return Type; }
};

using ViscoElasticMaterial = ViscoElasticMaterialParametrized<NUMBER_OF_RELAXATION_MECHANISMS>;
} // namespace seissol::model

#endif // SEISSOL_SRC_EQUATIONS_VISCOELASTIC2_MODEL_DATASTRUCTURES_H_
