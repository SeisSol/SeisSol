// SPDX-FileCopyrightText: 2015 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Common/Typedefs.h"
#include "Config.h"
#include "Equations/Datastructures.h"
#include "GeneratedCode/tensor.h"
#include "Kernels/Precision.h"
#include "Model/CommonDatastructures.h"
#include "SourceTerm/PointSource.h"
#include "TestHelper.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <limits>
#include <memory>

namespace seissol::unit_test {

TEST_CASE("Transform moment tensor" * doctest::test_suite("sourceterm")) {
  constexpr double Epsilon = 100 * std::numeric_limits<real>::epsilon();

  // the acoustic and the viscoacoustic material carry a single isotropic stress, the pressure,
  // so only the first diagonal entry of the moment tensor makes it into the source
  constexpr bool ScalarStress = model::MaterialT::TractionComponents == 1;

  // strike = dip = rake = pi / 3
  double strike = M_PI / 3.0;
  double dip = M_PI / 3.0;
  double rake = M_PI / 3.0;

  // M_xy = M_yx = 1, others are zero
  const double localMomentTensorXY[3][3] = {
      {0.0, 1.0, 0.0},
      {1.0, 0.0, 0.0},
      {0.0, 0.0, 0.0},
  };
  const double localSolidVelocityComponent[3] = {0.0, 0.0, 0.0};
  const double localPressureComponent = 0.0;
  const double localFluidVelocityComponent[3] = {0.0, 0.0, 0.0};

  auto momentTensor = seissol::memory::AlignedArray<real, tensor::update::Size>{};

  seissol::sourceterm::transformMomentTensor(localMomentTensorXY,
                                             localSolidVelocityComponent,
                                             localPressureComponent,
                                             localFluidVelocityComponent,
                                             strike,
                                             dip,
                                             rake,
                                             momentTensor.data());

  // Compare to hand-computed reference solution
  CHECK(momentTensor[0] == AbsApprox(-5.0 * std::sqrt(3.0) / 32.0).epsilon(Epsilon));
  if (!ScalarStress) {
    CHECK(momentTensor[1] == AbsApprox(-7.0 * std::sqrt(3.0) / 32.0).epsilon(Epsilon));
    CHECK(momentTensor[2] == AbsApprox(3.0 * std::sqrt(3.0) / 8.0).epsilon(Epsilon));
    CHECK(momentTensor[3] == AbsApprox(19.0 / 32.0).epsilon(Epsilon));
    CHECK(momentTensor[4] == AbsApprox(-9.0 / 16.0).epsilon(Epsilon));
    CHECK(momentTensor[5] == AbsApprox(-std::sqrt(3.0) / 16.0).epsilon(Epsilon));
    CHECK(momentTensor[6] == 0);
    CHECK(momentTensor[7] == 0);
    CHECK(momentTensor[8] == 0);
  } else {
    CHECK(momentTensor[1] == 0);
    CHECK(momentTensor[2] == 0);
    CHECK(momentTensor[3] == 0);
  }

  // strike = dip = rake = pi / 3
  strike = -1.349886940156521;
  dip = 3.034923466331855;
  rake = 0.725404224946106;

  // Random M
  const double localMomentTensorXZ[3][3] = {
      {1.833885014595086, -0.970040810572334, 0.602398893453385},
      {-0.970040810572334, -1.307688296305273, 1.572402458710038},
      {0.602398893453385, 1.572402458710038, 2.769437029884877},
  };

  seissol::sourceterm::transformMomentTensor(localMomentTensorXZ,
                                             localSolidVelocityComponent,
                                             localPressureComponent,
                                             localFluidVelocityComponent,
                                             strike,
                                             dip,
                                             rake,
                                             momentTensor.data());

  // Compare to hand-computed reference solution
  CHECK(momentTensor[0] == AbsApprox(-0.415053502680640).epsilon(Epsilon));
  if (!ScalarStress) {
    CHECK(momentTensor[1] == AbsApprox(0.648994284092410).epsilon(Epsilon));
    CHECK(momentTensor[2] == AbsApprox(3.061692966762920).epsilon(Epsilon));
    CHECK(momentTensor[3] == AbsApprox(1.909053142737053).epsilon(Epsilon));
    CHECK(momentTensor[4] == AbsApprox(0.677535767462651).epsilon(Epsilon));
    CHECK(momentTensor[5] == AbsApprox(-1.029826812214912).epsilon(Epsilon));
    CHECK(momentTensor[6] == 0.0);
    CHECK(momentTensor[7] == 0.0);
    CHECK(momentTensor[8] == 0.0);
  } else {
    CHECK(momentTensor[1] == 0);
    CHECK(momentTensor[2] == 0);
    CHECK(momentTensor[3] == 0);
  }
}

TEST_CASE("Transform moment tensor into the quantity layout" * doctest::test_suite("sourceterm")) {
  // the memory variables share the quantity axis with the other quantities in the fused layout only
  constexpr bool Split = Config::Solver == SolverType::LinearCKAnelastic;
  static_assert(tensor::update::Size == (Split ? model::MaterialT::NumElasticQuantities
                                               : model::MaterialT::NumQuantities),
                "The point source update has to span the quantity axis of Q.");

  constexpr auto Family = model::MaterialT::Type;
  constexpr bool AcousticFamily =
      Family == model::MaterialType::Acoustic || Family == model::MaterialType::Viscoacoustic;
  constexpr bool Poroelastic = Family == model::MaterialType::Poroelastic;

  // the sources lie back to back in one buffer, so nothing behind the update may be touched
  constexpr std::size_t Guard = 16;
  constexpr real Sentinel = 12345;
  auto buffer = std::array<real, tensor::update::Size + Guard>{};
  buffer.fill(Sentinel);

  const double localSolidVelocityComponent[3] = {7.0, 8.0, 9.0};
  const double localPressureComponent = 10.0;
  const double localFluidVelocityComponent[3] = {11.0, 12.0, 13.0};

  SUBCASE("every entry lands on its quantity") {
    const double localMomentTensor[3][3] = {
        {1.0, 4.0, 6.0},
        {4.0, 2.0, 5.0},
        {6.0, 5.0, 3.0},
    };

    // strike = dip = rake = 0 is the identity rotation, so every input arrives exactly
    seissol::sourceterm::transformMomentTensor(localMomentTensor,
                                               localSolidVelocityComponent,
                                               localPressureComponent,
                                               localFluidVelocityComponent,
                                               0.0,
                                               0.0,
                                               0.0,
                                               buffer.data());

    // Spelled out per material family on purpose, instead of being derived from the quantity
    // groups the implementation reads. The array has room for every index named below in any
    // build; the memory variables of a fused layout stay zero.
    auto expected = std::array<real, tensor::update::Size + Guard>{};
    std::fill(expected.begin() + tensor::update::Size, expected.end(), Sentinel);
    if constexpr (AcousticFamily) {
      // (pprime, v1, v2, v3); the pressure takes M_xx only
      expected[0] = 1.0;
      expected[1] = 7.0;
      expected[2] = 8.0;
      expected[3] = 9.0;
    } else {
      // (s_xx, s_yy, s_zz, s_xy, s_yz, s_xz, v1, v2, v3)
      expected[0] = 1.0;
      expected[1] = 2.0;
      expected[2] = 3.0;
      expected[3] = 4.0;
      expected[4] = 5.0;
      expected[5] = 6.0;
      expected[6] = 7.0;
      expected[7] = 8.0;
      expected[8] = 9.0;
      if constexpr (Poroelastic) {
        // (p, v1_f, v2_f, v3_f)
        expected[9] = 10.0;
        expected[10] = 11.0;
        expected[11] = 12.0;
        expected[12] = 13.0;
      }
    }

    for (std::size_t i = 0; i < buffer.size(); ++i) {
      CAPTURE(i);
      CHECK(buffer[i] == expected[i]);
    }
  }

  SUBCASE("the forces are rotated") {
    constexpr double Epsilon = 100 * std::numeric_limits<real>::epsilon();
    const double localMomentTensor[3][3] = {
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
        {0.0, 0.0, 0.0},
    };

    // strike = pi / 2 maps (x, y, z) to (y, -x, z)
    seissol::sourceterm::transformMomentTensor(localMomentTensor,
                                               localSolidVelocityComponent,
                                               localPressureComponent,
                                               localFluidVelocityComponent,
                                               M_PI / 2.0,
                                               0.0,
                                               0.0,
                                               buffer.data());

    constexpr std::size_t VelocityOffset = AcousticFamily ? 1 : 6;
    CHECK(buffer[VelocityOffset + 0] == AbsApprox(8.0).epsilon(Epsilon));
    CHECK(buffer[VelocityOffset + 1] == AbsApprox(-7.0).epsilon(Epsilon));
    CHECK(buffer[VelocityOffset + 2] == AbsApprox(9.0).epsilon(Epsilon));
    if constexpr (Poroelastic) {
      CHECK(buffer[10] == AbsApprox(12.0).epsilon(Epsilon));
      CHECK(buffer[11] == AbsApprox(-11.0).epsilon(Epsilon));
      CHECK(buffer[12] == AbsApprox(13.0).epsilon(Epsilon));
    }
    for (std::size_t i = tensor::update::Size; i < buffer.size(); ++i) {
      CAPTURE(i);
      CHECK(buffer[i] == Sentinel);
    }
  }
}

} // namespace seissol::unit_test
