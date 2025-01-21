// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_GRAVITATIONALFREESURFACEBC_H_
#define SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_GRAVITATIONALFREESURFACEBC_H_

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/Typedefs.h"

#include "Numerical/ODEInt.h"
#include "Numerical/Quadrature.h"
#include <Parallel/Runtime/Stream.h>

#include <utility>

#ifdef ACL_DEVICE
#include "Equations/elastic/Kernels/DeviceAux/KernelsAux.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.h"
#include "device.h"
#include <tuple>
#endif

namespace seissol {

class GravitationalFreeSurfaceBc {
  private:
  const double gravitationalAcceleration;

  public:
  GravitationalFreeSurfaceBc(double gravitationalAcceleration)
      : gravitationalAcceleration(gravitationalAcceleration) {};

  static std::pair<long long, long long>
      getFlopsDisplacementFace(unsigned face, [[maybe_unused]] FaceType faceType);

  template <typename TimeKrnl, typename MappingKrnl>
  void evaluate(unsigned faceIdx,
                MappingKrnl&& projectKernelPrototype,
                const CellBoundaryMapping& boundaryMapping,
                real* displacementNodalData,
                real* integratedDisplacementNodalData,
                TimeKrnl& timeKernel,
                const real* derivatives,
                double timeStepWidth,
                CellMaterialData& materialData,
                FaceType faceType) {
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

    // Prepare kernel that projects volume data to face and rotates it to face-nodal basis.
    assert(boundaryMapping.nodes != nullptr);
    assert(boundaryMapping.TinvData != nullptr);
    assert(boundaryMapping.TData != nullptr);
    auto tinv = init::Tinv::view::create(boundaryMapping.TinvData);
    auto t = init::Tinv::view::create(boundaryMapping.TData);
    auto projectKernel = std::forward<MappingKrnl>(projectKernelPrototype);
    projectKernel.Tinv = tinv.data();

    // Prepare projection of displacement/velocity to face-nodal basis.
    alignas(Alignment)
        real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];
    auto rotateDisplacementToFaceNormal =
        init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
    alignas(Alignment) real rotateDisplacementToGlobalData[init::displacementRotationMatrix::Size];
    auto rotateDisplacementToGlobal =
        init::displacementRotationMatrix::view::create(rotateDisplacementToGlobalData);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        // Extract part that rotates velocity from T
        rotateDisplacementToFaceNormal(i, j) = tinv(i + 6, j + 6);
        rotateDisplacementToGlobal(i, j) = t(i + 6, j + 6);
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
    alignas(Alignment) std::array<real, nodal::tensor::nodes2D::Shape[0]> prevCoefficients;

    const double deltaT = timeStepWidth;
    const double deltaTInt = timeStepWidth;

    // Initialize first component of Taylor series
    for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
      prevCoefficients[i] = rotatedFaceDisplacement(i, 0);
      // This is clearly a zeroth order approximation of the integral!
      integratedDisplacementNodal(i) = deltaTInt * rotatedFaceDisplacement(i, 0); // 1 FLOP
    }

    // Coefficients for Taylor series
    double factorEvaluated = 1;
    double factorInt = deltaTInt;

    auto* derivativesOffsets = timeKernel.getDerivativesOffsets();
    projectKernel.INodal = dofsFaceNodal.data();
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      projectKernel.dQ(i) = derivatives + derivativesOffsets[i];
    }

    const double rho = materialData.local.rho;
    const double g = gravitationalAcceleration; // [m/s^2]
    const double z = std::sqrt(materialData.local.getLambdaBar() * rho);

    // Note: Probably need to increase ConvergenceOrderby 1 here!
    for (std::size_t order = 1; order < ConvergenceOrder + 1; ++order) {
      dofsFaceNodal.setZero();

      projectKernel.execute(order - 1, faceIdx);

      factorEvaluated *= deltaT / (1.0 * order);
      factorInt *= deltaTInt / (order + 1.0);

#pragma omp simd
      for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
        // Derivatives of interior variables
        constexpr int PIdx = 0;
        constexpr int UIdx = 6;

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

#ifdef ACL_DEVICE
  template <typename TimeKrnl, typename MappingKrnl>
  void evaluateOnDevice(unsigned faceIdx,
                        MappingKrnl&& projectKernelPrototype,
                        TimeKrnl& timeKernel,
                        ConditionalPointersToRealsTable& dataTable,
                        ConditionalMaterialTable& materialTable,
                        double timeStepWidth,
                        device::DeviceInstance& device,
                        seissol::parallel::runtime::StreamRuntime& runtime) {

    auto* deviceStream = runtime.stream();
    ConditionalKey key(
        *KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, faceIdx);
    if (dataTable.find(key) != dataTable.end()) {

      real** rotateDisplacementToFaceNormalPtrs{nullptr};
      real** rotateDisplacementToGlobalPtrs{nullptr};
      real** rotatedFaceDisplacementPtrs{nullptr};
      real** dofsFaceNodalPtrs{nullptr};
      real** prevCoefficientsPtrs{nullptr};

      real* rotateDisplacementToFaceNormalMemory{nullptr};
      real* rotateDisplacementToGlobalMemory{nullptr};
      real* rotatedFaceDisplacementMemory{nullptr};
      real* dofsFaceNodalMemory{nullptr};
      real* prevCoefficientsMemory{nullptr};

      std::array<std::tuple<real***, real**, unsigned long>, 5> dataPack{
          std::make_tuple(&rotateDisplacementToFaceNormalPtrs,
                          &rotateDisplacementToFaceNormalMemory,
                          init::displacementRotationMatrix::Size),
          std::make_tuple(&rotateDisplacementToGlobalPtrs,
                          &rotateDisplacementToGlobalMemory,
                          init::displacementRotationMatrix::Size),
          std::make_tuple(&rotatedFaceDisplacementPtrs,
                          &rotatedFaceDisplacementMemory,
                          init::rotatedFaceDisplacement::Size),
          std::make_tuple(&dofsFaceNodalPtrs, &dofsFaceNodalMemory, tensor::INodal::size()),
          std::make_tuple(
              &prevCoefficientsPtrs, &prevCoefficientsMemory, nodal::tensor::nodes2D::Shape[0])};

      const size_t numElements{dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getSize()};
      for (auto& item : dataPack) {
        auto& ptrs = *std::get<0>(item);
        auto& memory = *std::get<1>(item);
        auto elementSize = std::get<2>(item);
        memory = reinterpret_cast<real*>(
            device.api->getStackMemory(elementSize * numElements * sizeof(real)));
        ptrs = reinterpret_cast<real**>(device.api->getStackMemory(numElements * sizeof(real*)));
        device.algorithms.incrementalAdd(ptrs, memory, elementSize, numElements, deviceStream);
      }
      size_t memCounter{2 * dataPack.size()};

      auto** TinvDataPtrs = dataTable[key].get(inner_keys::Wp::Id::Tinv)->getDeviceDataPtr();
      auto** TDataPtrs = dataTable[key].get(inner_keys::Wp::Id::T)->getDeviceDataPtr();
      kernels::time::aux::extractRotationMatrices(rotateDisplacementToFaceNormalPtrs,
                                                  rotateDisplacementToGlobalPtrs,
                                                  TDataPtrs,
                                                  TinvDataPtrs,
                                                  numElements,
                                                  deviceStream);

      auto rotateFaceDisplacementKrnl = kernel::gpu_rotateFaceDisplacement();
      const auto auxTmpMemSize =
          yateto::getMaxTmpMemRequired(rotateFaceDisplacementKrnl, projectKernelPrototype);
      auto* auxTmpMem =
          reinterpret_cast<real*>(device.api->getStackMemory(auxTmpMemSize * numElements));
      ++memCounter;

      auto** displacementsPtrs =
          dataTable[key].get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
      rotateFaceDisplacementKrnl.numElements = numElements;
      rotateFaceDisplacementKrnl.faceDisplacement = const_cast<const real**>(displacementsPtrs);
      rotateFaceDisplacementKrnl.displacementRotationMatrix =
          const_cast<const real**>(rotateDisplacementToFaceNormalPtrs);
      rotateFaceDisplacementKrnl.rotatedFaceDisplacement = rotatedFaceDisplacementPtrs;
      rotateFaceDisplacementKrnl.linearAllocator.initialize(auxTmpMem);
      rotateFaceDisplacementKrnl.streamPtr = deviceStream;
      rotateFaceDisplacementKrnl.execute();

      const double deltaT = timeStepWidth;
      const double deltaTInt = timeStepWidth;

      auto** integratedDisplacementNodalPtrs =
          dataTable[key].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();
      kernels::time::aux::initializeTaylorSeriesForGravitationalBoundary(
          prevCoefficientsPtrs,
          integratedDisplacementNodalPtrs,
          rotatedFaceDisplacementPtrs,
          deltaTInt,
          numElements,
          deviceStream);

      auto* invImpedances =
          reinterpret_cast<double*>(device.api->getStackMemory(sizeof(double) * numElements));
      ++memCounter;

      auto* rhos = materialTable[key].get(inner_keys::Material::Id::Rho)->getDeviceDataPtr();
      auto* lambdas = materialTable[key].get(inner_keys::Material::Id::Lambda)->getDeviceDataPtr();
      kernels::time::aux::computeInvAcousticImpedance(
          invImpedances, rhos, lambdas, numElements, deviceStream);

      double factorEvaluated = 1;
      double factorInt = deltaTInt;
      const double g = gravitationalAcceleration;

      auto** derivativesPtrs =
          dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getDeviceDataPtr();
      for (std::size_t order = 1; order < ConvergenceOrder + 1; ++order) {

        factorEvaluated *= deltaT / (1.0 * order);
        factorInt *= deltaTInt / (order + 1.0);

        device.algorithms.setToValue(
            dofsFaceNodalPtrs, 0.0, tensor::INodal::size(), numElements, deviceStream);

        auto projectKernel = projectKernelPrototype;
        projectKernel.numElements = numElements;
        projectKernel.Tinv = const_cast<const real**>(TinvDataPtrs);
        projectKernel.INodal = dofsFaceNodalPtrs;
        projectKernel.linearAllocator.initialize(auxTmpMem);
        projectKernel.streamPtr = deviceStream;

        auto* derivativesOffsets = timeKernel.getDerivativesOffsets();
        for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
          projectKernel.dQ(i) = const_cast<const real**>(derivativesPtrs);
          projectKernel.extraOffset_dQ(i) = derivativesOffsets[i];
        }

        projectKernel.execute(order - 1, faceIdx);

        kernels::time::aux::updateRotatedFaceDisplacement(rotatedFaceDisplacementPtrs,
                                                          prevCoefficientsPtrs,
                                                          integratedDisplacementNodalPtrs,
                                                          dofsFaceNodalPtrs,
                                                          invImpedances,
                                                          rhos,
                                                          g,
                                                          factorEvaluated,
                                                          factorInt,
                                                          numElements,
                                                          deviceStream);
      }

      rotateFaceDisplacementKrnl.numElements = numElements;
      rotateFaceDisplacementKrnl.faceDisplacement =
          const_cast<const real**>(rotatedFaceDisplacementPtrs);
      rotateFaceDisplacementKrnl.displacementRotationMatrix =
          const_cast<const real**>(rotateDisplacementToGlobalPtrs);
      rotateFaceDisplacementKrnl.rotatedFaceDisplacement = displacementsPtrs;
      rotateFaceDisplacementKrnl.linearAllocator.initialize(auxTmpMem);
      rotateFaceDisplacementKrnl.streamPtr = deviceStream;
      rotateFaceDisplacementKrnl.execute();

      for (size_t i = 0; i < memCounter; ++i) {
        device.api->popStackMemory();
      }
    }
  }
#endif // ACL_DEVICE
};

} // namespace seissol

#endif // SEISSOL_SRC_EQUATIONS_ELASTIC_KERNELS_GRAVITATIONALFREESURFACEBC_H_
