// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Equations/poroelastic/Model/Datastructures.h"
#include "Equations/poroelastic/Model/Helper.h"
#include "Equations/poroelastic/Model/Setup.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/pool.h"
#include "GeneratedCode/tensor.h"
#include "Geometry/CellTransform.h"
#include "Initializer/Typedefs.h"
#include "Kernels/Common.h"
#include "Kernels/STP/Setup.h"
#include "Kernels/StarOperands.h"
#include "Model/Common.h"
#include "Model/OperatorLayout.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <random>
#include <type_traits>

namespace seissol::unit_test {

class SpaceTimePredictorTestFixture {
  protected:
  constexpr static const double Epsilon = std::numeric_limits<real>::epsilon();
  constexpr static const double Dt = 1.05109e-06;
  real starMatrices0[tensor::star::size(0)];
  real starMatrices1[tensor::star::size(1)];
  real starMatrices2[tensor::star::size(2)];
  real sourceMatrix[tensor::ET::size()];
  real zMatrix[seissol::model::MaterialT::NumQuantities][tensor::Zinv::size(0)];
  LocalIntegrationData localIntegration;

  void setStarMatrix(const real* at,
                     const real* bt,
                     const real* ct,
                     const std::array<double, Cell::Dim>& grad,
                     real* starMatrix) {
    for (unsigned idx = 0; idx < seissol::tensor::star::size(0); ++idx) {
      starMatrix[idx] = grad[0] * at[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star::size(1); ++idx) {
      starMatrix[idx] += grad[1] * bt[idx];
    }

    for (unsigned idx = 0; idx < seissol::tensor::star::size(2); ++idx) {
      starMatrix[idx] += grad[2] * ct[idx];
    }
  }

  void prepareModel() {
    // prepare Material
    const auto materialVals =
        std::vector<double>{{40.0e9, 2500, 12.0e9, 10.0e9, 0.2, 600.0e-15, 3, 2.5e9, 1040, 0.001}};
    const model::PoroElasticMaterial material(materialVals);

    // NOLINTNEXTLINE (-cert-dcl59-cpp)
    std::mt19937 generator(20210109); // Standard mersenne_twister_engine seeded with today's date
    std::uniform_real_distribution<real> distribution(0, 1);
    std::array<CoordinateT, Cell::NumVertices> vertices{};
    for (auto& vertex : vertices) {
      vertex =
          CoordinateT{distribution(generator), distribution(generator), distribution(generator)};
    }

    const auto grad = seissol::geometry::AffineTransform(vertices).refToSpaceJacobianInverse(
        seissol::geometry::CellTransform::VectorEigenT(Cell::ReferenceBarycenter.data()));

    std::array<double, Cell::Dim> gradXi{};
    std::array<double, Cell::Dim> gradEta{};
    std::array<double, Cell::Dim> gradZeta{};
    for (std::size_t i = 0; i < Cell::Dim; ++i) {
      gradXi[i] = grad(0, i);
      gradEta[i] = grad(1, i);
      gradZeta[i] = grad(2, i);
    }

    // prepare starmatrices
    real atData[tensor::star::size(0)];
    real btData[tensor::star::size(1)];
    real ctData[tensor::star::size(2)];
    auto at = init::star::view<0>::create(atData);
    auto bt = init::star::view<0>::create(btData);
    auto ct = init::star::view<0>::create(ctData);
    model::getTransposedCoefficientMatrix(material, 0, at);
    model::getTransposedCoefficientMatrix(material, 1, bt);
    model::getTransposedCoefficientMatrix(material, 2, ct);
    setStarMatrix(atData, btData, ctData, gradXi, starMatrices0);
    setStarMatrix(atData, btData, ctData, gradEta, starMatrices1);
    setStarMatrix(atData, btData, ctData, gradZeta, starMatrices2);

    // The three matrices above state the system the predictor has to solve;
    // what a cell hands the kernel is whatever that build has it carry, so
    // fill both and let the binding pick.
    if constexpr (FactoredStar) {
      for (std::size_t dim = 0; dim < Cell::Dim; ++dim) {
        for (std::size_t component = 0; component < Cell::Dim; ++component) {
          localIntegration.referenceGradients[dim][component] = grad(dim, component);
        }
      }
      const auto coefficients = seissol::model::getStarCoefficients(material);
      for (std::size_t i = 0; i < coefficients.size(); ++i) {
        // a material that does not vary carries the same value at every sample
        for (std::size_t point = 0; point < MaterialSampleCount; ++point) {
          localIntegration.materialCoefficients[i][point] = coefficients[i];
        }
      }
    } else {
      const real* const matrices[3] = {starMatrices0, starMatrices1, starMatrices2};
      for (std::size_t dim = 0; dim < 3; ++dim) {
        std::copy_n(matrices[dim], tensor::star::size(dim), localIntegration.starMatrices[dim]);
      }
    }

    // prepare sourceterm
    auto et = init::ET::view::create(sourceMatrix);
    model::getTransposedSourceCoefficientTensor(material, et);

    // prepare Zinv
    auto zinv0 = init::Zinv::view<0>::create(zMatrix[0]);
    model::calcZinv(zinv0, et, 0, model::isStiffRow<model::PoroElasticMaterial>(0), Dt);
    auto zinv1 = init::Zinv::view<1>::create(zMatrix[1]);
    model::calcZinv(zinv1, et, 1, model::isStiffRow<model::PoroElasticMaterial>(1), Dt);
    auto zinv2 = init::Zinv::view<2>::create(zMatrix[2]);
    model::calcZinv(zinv2, et, 2, model::isStiffRow<model::PoroElasticMaterial>(2), Dt);
    auto zinv3 = init::Zinv::view<3>::create(zMatrix[3]);
    model::calcZinv(zinv3, et, 3, model::isStiffRow<model::PoroElasticMaterial>(3), Dt);
    auto zinv4 = init::Zinv::view<4>::create(zMatrix[4]);
    model::calcZinv(zinv4, et, 4, model::isStiffRow<model::PoroElasticMaterial>(4), Dt);
    auto zinv5 = init::Zinv::view<5>::create(zMatrix[5]);
    model::calcZinv(zinv5, et, 5, model::isStiffRow<model::PoroElasticMaterial>(5), Dt);
    auto zinv6 = init::Zinv::view<6>::create(zMatrix[6]);
    model::calcZinv(zinv6, et, 6, model::isStiffRow<model::PoroElasticMaterial>(6), Dt);
    auto zinv7 = init::Zinv::view<7>::create(zMatrix[7]);
    model::calcZinv(zinv7, et, 7, model::isStiffRow<model::PoroElasticMaterial>(7), Dt);
    auto zinv8 = init::Zinv::view<8>::create(zMatrix[8]);
    model::calcZinv(zinv8, et, 8, model::isStiffRow<model::PoroElasticMaterial>(8), Dt);
    auto zinv9 = init::Zinv::view<9>::create(zMatrix[9]);
    model::calcZinv(zinv9, et, 9, model::isStiffRow<model::PoroElasticMaterial>(9), Dt);
    auto zinv10 = init::Zinv::view<10>::create(zMatrix[10]);
    model::calcZinv(zinv10, et, 10, model::isStiffRow<model::PoroElasticMaterial>(10), Dt);
    auto zinv11 = init::Zinv::view<11>::create(zMatrix[11]);
    model::calcZinv(zinv11, et, 11, model::isStiffRow<model::PoroElasticMaterial>(11), Dt);
    auto zinv12 = init::Zinv::view<12>::create(zMatrix[12]);
    model::calcZinv(zinv12, et, 12, model::isStiffRow<model::PoroElasticMaterial>(12), Dt);
  }

  // Which constant matrices a kernel reads is the generator's business, so let
  // it hand them over: the operator a cell carries decides how many there are.
  template <typename KernelT>
  void prepareKernel(KernelT& krnlPrototype) {
    krnlPrototype.bindGlobals(seissol::Pool::host());
  }

  void prepareQ(real* qData) {
    // scale quantities to make it more realistic
    std::array<real, 13> factor = {{1e9, 1e9, 1e9, 1e9, 1e9, 1e9, 1, 1, 1, 1e9, 1, 1, 1}};
    auto q = init::Q::view::create(qData);

    // NOLINTNEXTLINE (-cert-dcl59-cpp)
    std::mt19937 rnggen(1234);
    std::uniform_real_distribution<> rngdist(0.0, 1.0);

    for (std::size_t qi = 0; qi < q.shape(1); ++qi) {
      for (std::size_t bf = 0; bf < q.shape(0); ++bf) {
        q(bf, qi) = rngdist(rnggen) * factor.at(qi);
      }
    }
  }

  void solveWithKernel(real stp[], const real* qData) {
    real timeIntegrated[seissol::tensor::I::size()];

    seissol::kernel::spaceTimePredictor krnl;
    prepareKernel(krnl);

    // the predictor reads the operator the way a cell carries it, and scales
    // it by the timestep itself. A material that does not vary asks nothing of
    // the source beyond what the cell carries, so that operand stays zero.
    kernels::bindStarOperands(krnl, localIntegration);
    kernels::bindSourceDeviationOperands(krnl, localIntegration);

    for (size_t i = 0; i < seissol::model::MaterialT::NumQuantities; i++) {
      krnl.Zinv(i) = zMatrix[i];
    }

    auto sourceView = init::ET::view::create(sourceMatrix);
    for (std::size_t i = 0; i < model::PoroElasticMaterial::StiffSourceRows.size(); ++i) {
      const auto& row = model::PoroElasticMaterial::StiffSourceRows[i];
      krnl.G(i) = sourceView(row.quantity, row.target) * Dt;
    }

    krnl.Q = qData;
    krnl.I = timeIntegrated;
    krnl.timestep = Dt;
    krnl.spaceTimePredictor = stp;
    krnl.execute();
  }

  void computeLhs(const real* stp, real* lhs) {
    kernel::stpTestLhs testLhsKrnl;
    prepareKernel(testLhsKrnl);
    testLhsKrnl.ET = sourceMatrix;
    testLhsKrnl.spaceTimePredictor = stp;
    testLhsKrnl.testLhs = lhs;
    testLhsKrnl.minus = -Dt;
    testLhsKrnl.execute();
  };

  void computeRhs(const real* stp, const real* qData, real* rhs) {
    kernel::stpTestRhs testRhsKrnl;
    prepareKernel(testRhsKrnl);
    testRhsKrnl.Q = qData;
    testRhsKrnl.star(0) = starMatrices0;
    testRhsKrnl.star(1) = starMatrices1;
    testRhsKrnl.star(2) = starMatrices2;
    testRhsKrnl.spaceTimePredictor = stp;
    testRhsKrnl.timestep = Dt;
    testRhsKrnl.testRhs = rhs;
    testRhsKrnl.execute();
  };

  public:
  SpaceTimePredictorTestFixture() { prepareModel(); };
};

TEST_CASE_FIXTURE(SpaceTimePredictorTestFixture,
                  "Solve Space Time Predictor" * doctest::test_suite("kernel")) {
  alignas(PagesizeStack) real stp[seissol::tensor::spaceTimePredictor::size()];
  alignas(PagesizeStack) real rhs[seissol::tensor::testLhs::size()];
  alignas(PagesizeStack) real lhs[seissol::tensor::testRhs::size()];
  alignas(PagesizeStack) real qData[seissol::tensor::Q::size()];
  std::fill(std::begin(stp), std::end(stp), 0);
  std::fill(std::begin(rhs), std::end(rhs), 0);
  std::fill(std::begin(lhs), std::end(lhs), 0);
  std::fill(std::begin(qData), std::end(qData), 0);

  prepareQ(qData);

  solveWithKernel(stp, qData);

  computeLhs(stp, lhs);
  computeRhs(stp, qData, rhs);

  double diffNorm = 0;
  double refNorm = 0;

  auto lhsView = init::testLhs::view::create(lhs);
  auto rhsView = init::testRhs::view::create(rhs);

  for (size_t b = 0; b < tensor::spaceTimePredictor::Shape[0]; b++) {
    for (size_t q = 0; q < tensor::spaceTimePredictor::Shape[1]; q++) {
      for (size_t o = 0; o < tensor::spaceTimePredictor::Shape[2]; o++) {
        const double d = std::abs(lhsView(b, q, o) - rhsView(b, q, o));
        const double a = std::abs(lhsView(b, q, o));
        diffNorm += d * d;
        refNorm += a * a;
      }
    }
  }

  CHECK(diffNorm / refNorm < Epsilon);
}

} // namespace seissol::unit_test
