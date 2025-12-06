// SPDX-FileCopyrightText: 2020 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_KERNELS_LINEARCK_GRAVITATIONALFREESURFACEBC_H_
#define SEISSOL_SRC_KERNELS_LINEARCK_GRAVITATIONALFREESURFACEBC_H_

#include "Equations/Datastructures.h"
#include "GeneratedCode/init.h"
#include "GeneratedCode/kernel.h"
#include "GeneratedCode/tensor.h"
#include "Initializer/Typedefs.h"
#include "Numerical/ODEInt.h"
#include "Numerical/Quadrature.h"
#include "Parallel/Runtime/Stream.h"
#include "Solver/MultipleSimulations.h"

#include <utility>

#ifdef ACL_DEVICE
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "Kernels/LinearCK/DeviceAux/KernelsAux.h"

#include <Device/device.h>
#include <tuple>
#endif

namespace seissol {

class GravitationalFreeSurfaceBc {
  private:
  double gravitationalAcceleration;

  public:
  explicit GravitationalFreeSurfaceBc(double gravitationalAcceleration)
      : gravitationalAcceleration(gravitationalAcceleration) {};

  static std::pair<std::uint64_t, std::uint64_t>
      getFlopsDisplacementFace(unsigned face, [[maybe_unused]] FaceType faceType);

  template <typename TimeKrnl, typename MappingKrnl>
  void evaluate(unsigned faceIdx,
                MappingKrnl&& projectKernelPrototype,
                const CellBoundaryMapping& boundaryMapping,
                real* displacementNodalData,
                real* integratedDisplacementNodalData,
                TimeKrnl& /*timeKernel*/,
                const real* derivatives,
                double timeStepWidth,
                CellMaterialData& materialData,
                FaceType /*faceType*/) {
    // This function does two things:
    // 1: Compute eta (for all three dimensions) at the end of the timestep
    // 2: Compute the integral of eta in normal direction over the timestep
    // We do this by building up the Taylor series of eta.
    // Eta is defined by the ODE eta_t = u^R - 1/Z * (rho g eta - p^R)
    // Compute coefficients by differentiating ODE recursively, e.g.:
    // eta_tt = u^R_t - 1/Z * (rho eta_t g - p^R_t)
    // and substituting the previous coefficient eta_t
    // This implementation sums up the Taylor series directly without storing
    // all coefficients.
    if constexpr (multisim::MultisimEnabled) {
      logError()
          << "The Free Surface Gravity BC kernel does not work with multiple simulations yet.";
    } else {
      constexpr int PIdx = 0;
      constexpr int UIdx = model::MaterialT::TractionQuantities;

      // Prepare kernel that projects volume data to face and rotates it to face-nodal basis.
      assert(boundaryMapping.nodes != nullptr);
      assert(boundaryMapping.dataTinv != nullptr);
      assert(boundaryMapping.dataT != nullptr);
      auto tinv = init::Tinv::view::create(boundaryMapping.dataTinv);
      auto t = init::Tinv::view::create(boundaryMapping.dataT);
      auto projectKernel = std::forward<MappingKrnl>(projectKernelPrototype);
      projectKernel.Tinv = tinv.data();

      // Prepare projection of displacement/velocity to face-nodal basis.
      alignas(Alignment)
          real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];
      auto rotateDisplacementToFaceNormal =
          init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
      alignas(Alignment)
          real rotateDisplacementToGlobalData[init::displacementRotationMatrix::Size];
      auto rotateDisplacementToGlobal =
          init::displacementRotationMatrix::view::create(rotateDisplacementToGlobalData);
      for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
          // Extract part that rotates velocity from T
          rotateDisplacementToFaceNormal(i, j) = tinv(i + UIdx, j + UIdx);
          rotateDisplacementToGlobal(i, j) = t(i + UIdx, j + UIdx);
        }
      }
      static_assert(init::rotatedFaceDisplacement::Size == init::faceDisplacement::Size);
      alignas(Alignment) real rotatedFaceDisplacementData[init::rotatedFaceDisplacement::Size];

      auto integratedDisplacementNodal =
          init::averageNormalDisplacement::view::create(integratedDisplacementNodalData);
      auto rotatedFaceDisplacement =
          init::faceDisplacement::view::create(rotatedFaceDisplacementData);

      // Rotate face displacement to face-normal coordinate system in which the computation is
      // more convenient.
      auto rotateFaceDisplacementKrnl = kernel::rotateFaceDisplacement();
      rotateFaceDisplacementKrnl.faceDisplacement = displacementNodalData;
      rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
      rotateFaceDisplacementKrnl.rotatedFaceDisplacement = rotatedFaceDisplacementData;
      rotateFaceDisplacementKrnl.execute();

      // Temporary buffer to store nodal face dofs at some time t
      alignas(Alignment) real dofsFaceNodalStorage[tensor::INodal::size()];
      auto dofsFaceNodal = init::INodal::view::create(dofsFaceNodalStorage);

      // Temporary buffer to store nodal face coefficients at some time t
      alignas(Alignment) std::array<real, nodal::tensor::nodes2D::Shape[0]> prevCoefficients{};

      const double deltaT = timeStepWidth;
      const double deltaTInt = timeStepWidth;

      // Initialize first component of Taylor series
      for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
        const auto localCoeff = rotatedFaceDisplacement(i, 0);
        prevCoefficients[i] = localCoeff;
        // This is clearly a zeroth order approximation of the integral!
        integratedDisplacementNodal(i) = deltaTInt * localCoeff; // 1 FLOP
      }

      // Coefficients for Taylor series
      double factorEvaluated = 1;
      double factorInt = deltaTInt;

      projectKernel.INodal = dofsFaceNodal.data();
      for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        projectKernel.dQ(i) = derivatives + yateto::computeFamilySize<tensor::dQ>(1, i);
      }

      const double rho = materialData.local->getDensity();
      const double g = gravitationalAcceleration; // [m/s^2]
      const double z = std::sqrt(materialData.local->getLambdaBar() * rho);

      // Note: Probably need to increase ConvergenceOrderby 1 here!
      for (std::size_t order = 1; order < ConvergenceOrder + 1; ++order) {
        dofsFaceNodal.setZero();

        projectKernel.execute(order - 1, faceIdx);

        factorEvaluated *= deltaT / (1.0 * order);
        factorInt *= deltaTInt / (order + 1.0);

#pragma omp simd
        for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
          // Derivatives of interior variables
          const auto uInside = dofsFaceNodal(i, UIdx + 0);
          const auto vInside = dofsFaceNodal(i, UIdx + 1);
          const auto wInside = dofsFaceNodal(i, UIdx + 2);
          const auto pressureInside = dofsFaceNodal(i, PIdx);

          const double curCoeff =
              uInside - (1.0 / z) * (rho * g * prevCoefficients[i] + pressureInside);
          // Basically uInside - C_1 * (c_2 * prevCoeff[i] + pressureInside)
          // 2 add, 2 mul = 4 flops

          prevCoefficients[i] = curCoeff;

          // 2 * 3 = 6 flops for updating displacement
          rotatedFaceDisplacement(i, 0) += factorEvaluated * curCoeff;
          rotatedFaceDisplacement(i, 1) += factorEvaluated * vInside;
          rotatedFaceDisplacement(i, 2) += factorEvaluated * wInside;

          // 2 flops for updating integral of displacement
          integratedDisplacementNodal(i) += factorInt * curCoeff;
        }
      }

      // Rotate face displacement back to global coordinate system which we use as storage
      // coordinate system
      rotateFaceDisplacementKrnl.faceDisplacement = rotatedFaceDisplacementData;
      rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToGlobalData;
      rotateFaceDisplacementKrnl.rotatedFaceDisplacement = displacementNodalData;
      rotateFaceDisplacementKrnl.execute();
    }
  }

#ifdef ACL_DEVICE
  template <typename TimeKrnl, typename MappingKrnl>
  void evaluateOnDevice(unsigned faceIdx,
                        MappingKrnl&& projectKernelPrototype,
                        TimeKrnl& timeKernel,
                        recording::ConditionalPointersToRealsTable& dataTable,
                        recording::ConditionalMaterialTable& materialTable,
                        double timeStepWidth,
                        device::DeviceInstance& device,
                        seissol::parallel::runtime::StreamRuntime& runtime) {

    using namespace seissol::recording;
    auto* deviceStream = runtime.stream();
    ConditionalKey key(
        *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, faceIdx);
    if (dataTable.find(key) != dataTable.end()) {
      const size_t numElements{dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getSize()};
      // const double g = gravitationalAcceleration;

      auto* invImpedances =
          materialTable[key].get(inner_keys::Material::Id::InvImpedances)->getDeviceDataPtr();

      auto* rhos = materialTable[key].get(inner_keys::Material::Id::Rho)->getDeviceDataPtr();
      auto* lambdas = materialTable[key].get(inner_keys::Material::Id::Lambda)->getDeviceDataPtr();
      kernels::time::aux::computeInvAcousticImpedance(
          invImpedances, rhos, lambdas, numElements, deviceStream);

      auto** TinvDataPtrs = dataTable[key].get(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();
      auto** TDataPtrs = dataTable[key].get(inner_keys::Wp::Id::T)->getDeviceDataPtr();
      auto** derivativesPtrs =
          dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getDeviceDataPtr();

      auto** displacementsPtrs =
          dataTable[key].get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
      auto** integratedDisplacementNodalPtrs =
          dataTable[key].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();

      kernel::gpu_fsgKernel kernel;
      for (std::size_t i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        kernel.dQ(i) = const_cast<const real**>(derivativesPtrs);
        kernel.extraOffset_dQ(i) = yateto::computeFamilySize<tensor::dQ>(1, i);
      }
      kernel.Tinv = const_cast<const real**>(TinvDataPtrs);
      kernel.T = const_cast<const real**>(TDataPtrs);
      kernel.faceDisplacement = displacementsPtrs;
      kernel.Iint = integratedDisplacementNodalPtrs;

      double coeffInt = timeStepWidth;
      double coeff = 1;

      kernel.power(0) = timeStepWidth;
      for (std::size_t i = 1; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
        coeffInt *= timeStepWidth / (i + 1);
        coeff *= timeStepWidth / i;

        kernel.coeff(i) = coeff;
        kernel.power(i) = coeffInt;
      }

      kernel.streamPtr = deviceStream;
      kernel.execute(faceIdx);
    }
  }
#endif // ACL_DEVICE
};

} // namespace seissol

#endif // SEISSOL_SRC_KERNELS_LINEARCK_GRAVITATIONALFREESURFACEBC_H_
