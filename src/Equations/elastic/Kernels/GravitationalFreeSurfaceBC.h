#ifndef SEISSOL_GRAVIATIONALFREESURFACEBC_H
#define SEISSOL_GRAVIATIONALFREESURFACEBC_H

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/typedefs.hpp"

#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/ODEInt.h"

#ifdef ACL_DEVICE
#include "device.h"
#include "Initializer/BatchRecorders/DataTypes/ConditionalTable.hpp"
#include "Equations/elastic/Kernels/DeviceAux/KernelsAux.h"
#include <tuple>
#endif

namespace seissol {

// Used to avoid including SeisSo.h here as this leads to all sorts of issues
double getGravitationalAcceleration();

template<typename Config>
class GravitationalFreeSurfaceBc {
public:
  GravitationalFreeSurfaceBc() = default;

  using RealT = typename Config::RealT;
  using MaterialT = typename Config::MaterialT;

  static std::pair<long long, long long> getFlopsDisplacementFace(unsigned face,
                                                                  [[maybe_unused]] FaceType faceType){
    long long hardwareFlops = 0;
    long long nonZeroFlops = 0;

    constexpr auto numberOfNodes = nodal::tensor::nodes2D::Shape[0];

    // initialize integral of displacement
    hardwareFlops += 1 * numberOfNodes;
    nonZeroFlops += 1 * numberOfNodes;

    // Before adjusting the range of the loop, check range of loop in computation!
    for (int order = 1; order < ConvergenceOrder+ 1; ++order) {
  #ifdef USE_ELASTIC
      constexpr auto flopsPerQuadpoint =
          4 + // Computing coefficient
          6 + // Updating displacement
          2; // Updating integral of displacement
  #else
      constexpr auto flopsPerQuadpoint = 0;
  #endif
      constexpr auto flopsUpdates = flopsPerQuadpoint * numberOfNodes;

      nonZeroFlops += kernel::projectDerivativeToNodalBoundaryRotated::nonZeroFlops(order - 1, face) + flopsUpdates;
      hardwareFlops += kernel::projectDerivativeToNodalBoundaryRotated::hardwareFlops(order - 1, face) + flopsUpdates;
    }

    // Two rotations: One to face-aligned, one to global
    hardwareFlops += 2 * kernel::rotateFaceDisplacement::HardwareFlops;
    nonZeroFlops += 2 * kernel::rotateFaceDisplacement::NonZeroFlops;

    return {nonZeroFlops, hardwareFlops};
  }

  template<typename TimeKrnl, typename MappingKrnl, typename MaterialT>
  void evaluate(unsigned faceIdx,
                MappingKrnl&& projectKernelPrototype,
                const CellBoundaryMapping& boundaryMapping,
                RealT* displacementNodalData,
                RealT* integratedDisplacementNodalData,
                TimeKrnl& timeKernel,
                RealT* derivatives,
                double timeStepWidth,
                const MaterialT& materialData,
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
    auto Tinv = init::Tinv::view::create(boundaryMapping.TinvData);
    auto T = init::Tinv::view::create(boundaryMapping.TData);
    auto projectKernel = projectKernelPrototype;
    projectKernel.Tinv = Tinv.data();

    // Prepare projection of displacement/velocity to face-nodal basis.
    alignas(Alignment) RealT rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];
    auto rotateDisplacementToFaceNormal = init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
    alignas(Alignment) RealT rotateDisplacementToGlobalData[init::displacementRotationMatrix::Size];
    auto rotateDisplacementToGlobal = init::displacementRotationMatrix::view::create(rotateDisplacementToGlobalData);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        // Extract part that rotates velocity from T
        rotateDisplacementToFaceNormal(i, j) = Tinv(i + 6, j + 6);
        rotateDisplacementToGlobal(i,j) = T(i+6, j+6);
      }
    }
    static_assert(init::rotatedFaceDisplacement::Size == init::faceDisplacement::Size);
    alignas(Alignment) RealT rotatedFaceDisplacementData[init::rotatedFaceDisplacement::Size];

    auto integratedDisplacementNodal = init::averageNormalDisplacement::view::create(integratedDisplacementNodalData);
    auto rotatedFaceDisplacement = init::faceDisplacement::view::create(rotatedFaceDisplacementData);

    // Rotate face displacement to face-normal coordinate system in which the computation is
    // more convenient.
    auto rotateFaceDisplacementKrnl = kernel::rotateFaceDisplacement();
    rotateFaceDisplacementKrnl.faceDisplacement = displacementNodalData;
    rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
    rotateFaceDisplacementKrnl.rotatedFaceDisplacement = rotatedFaceDisplacementData;
    rotateFaceDisplacementKrnl.execute();

    // Temporary buffer to store nodal face dofs at some time t
    alignas(Alignment) RealT dofsFaceNodalStorage[tensor::INodal::size()];
    auto dofsFaceNodal = init::INodal::view::create(dofsFaceNodalStorage);

    // Temporary buffer to store nodal face coefficients at some time t
    alignas(Alignment) std::array<RealT, nodal::tensor::nodes2D::Shape[0]> prevCoefficients;

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

    const double rho = materialData.rho;
    const double Z = MaterialT::Type == seissol::model::MaterialType::elastic ? std::sqrt(materialData.lambda * rho) : 1; // TODO non-elastic
    const double g = getGravitationalAcceleration(); // [m/s^2]

    // Note: Probably need to increase ConvergenceOrdery 1 here!
    for (int order = 1; order < ConvergenceOrder; ++order) {
      dofsFaceNodal.setZero();

      projectKernel.execute(order - 1, faceIdx);

      factorEvaluated *= deltaT / (1.0 * order);
      factorInt *= deltaTInt / (order + 1.0);

#pragma omp simd
      for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
        // Derivatives of interior variables
        constexpr int pIdx = 0;
        constexpr int uIdx = 6;

        const auto uInside = dofsFaceNodal(i, uIdx + 0);
        const auto vInside = dofsFaceNodal(i, uIdx + 1);
        const auto wInside = dofsFaceNodal(i, uIdx + 2);
        const auto pressureInside = dofsFaceNodal(i, pIdx);

        const double curCoeff = MaterialT::Type == seissol::model::MaterialType::elastic ? (uInside - (1.0/Z) * (rho * g * prevCoefficients[i] + pressureInside)) : uInside;
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
  template<typename TimeKrnl, typename MappingKrnl>
  void evaluateOnDevice(unsigned faceIdx,
                        MappingKrnl&& projectKernelPrototype,
                        TimeKrnl& timeKernel,
                        ConditionalPointersToRealsTable &dataTable,
                        ConditionalMaterialTable &materialTable,
                        double timeStepWidth,
                        device::DeviceInstance& device) {

    auto* deviceStream = device.api->getDefaultStream();
    ConditionalKey key(*KernelNames::BoundaryConditions, *ComputationKind::FreeSurfaceGravity, faceIdx);
    if(dataTable.find(key) != dataTable.end()) {

      RealT** rotateDisplacementToFaceNormalPtrs{nullptr};
      RealT** rotateDisplacementToGlobalPtrs{nullptr};
      RealT** rotatedFaceDisplacementPtrs{nullptr};
      RealT** dofsFaceNodalPtrs{nullptr};
      RealT** prevCoefficientsPtrs{nullptr};

      RealT* rotateDisplacementToFaceNormalMemory{nullptr};
      RealT* rotateDisplacementToGlobalMemory{nullptr};
      RealT* rotatedFaceDisplacementMemory{nullptr};
      RealT* dofsFaceNodalMemory{nullptr};
      RealT* prevCoefficientsMemory{nullptr};

      std::array<std::tuple<RealT***, RealT**, unsigned long>, 5> dataPack{
          std::make_tuple(&rotateDisplacementToFaceNormalPtrs, &rotateDisplacementToFaceNormalMemory, init::displacementRotationMatrix::Size),
          std::make_tuple(&rotateDisplacementToGlobalPtrs, &rotateDisplacementToGlobalMemory, init::displacementRotationMatrix::Size),
          std::make_tuple(&rotatedFaceDisplacementPtrs, &rotatedFaceDisplacementMemory, init::rotatedFaceDisplacement::Size),
          std::make_tuple(&dofsFaceNodalPtrs, &dofsFaceNodalMemory, tensor::INodal::size()),
          std::make_tuple(&prevCoefficientsPtrs, &prevCoefficientsMemory, nodal::tensor::nodes2D::Shape[0])
      };

      const size_t numElements{dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getSize()};
      for (auto& item: dataPack) {
        auto& ptrs = *std::get<0>(item);
        auto& memory = *std::get<1>(item);
        auto elementSize = std::get<2>(item);
        memory = reinterpret_cast<RealT*>(device.api->getStackMemory(elementSize * numElements * sizeof(RealT)));
        ptrs  = reinterpret_cast<RealT**>(device.api->getStackMemory(numElements * sizeof(RealT*)));
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
      const auto auxTmpMemSize = yateto::getMaxTmpMemRequired(rotateFaceDisplacementKrnl, projectKernelPrototype);
      auto* auxTmpMem = reinterpret_cast<RealT*>(device.api->getStackMemory(auxTmpMemSize * numElements));
      ++memCounter;

      auto** displacementsPtrs = dataTable[key].get(inner_keys::Wp::Id::FaceDisplacement)->getDeviceDataPtr();
      rotateFaceDisplacementKrnl.numElements = numElements;
      rotateFaceDisplacementKrnl.faceDisplacement = const_cast<const RealT**>(displacementsPtrs);
      rotateFaceDisplacementKrnl.displacementRotationMatrix = const_cast<const RealT**>(rotateDisplacementToFaceNormalPtrs);
      rotateFaceDisplacementKrnl.rotatedFaceDisplacement = rotatedFaceDisplacementPtrs;
      rotateFaceDisplacementKrnl.linearAllocator.initialize(auxTmpMem);
      rotateFaceDisplacementKrnl.streamPtr = deviceStream;
      rotateFaceDisplacementKrnl.execute();

      const double deltaT = timeStepWidth;
      const double deltaTInt = timeStepWidth;

      auto** integratedDisplacementNodalPtrs = dataTable[key].get(inner_keys::Wp::Id::NodalAvgDisplacements)->getDeviceDataPtr();
      kernels::time::aux::initializeTaylorSeriesForGravitationalBoundary(
        prevCoefficientsPtrs,
        integratedDisplacementNodalPtrs,
        rotatedFaceDisplacementPtrs,
        deltaTInt,
        numElements,
        deviceStream);

      auto* invImpedances = reinterpret_cast<double*>(device.api->getStackMemory(sizeof(double) * numElements));
      ++memCounter;

      auto* rhos = materialTable[key].get(inner_keys::Material::Id::Rho)->getDeviceDataPtr();
      auto* lambdas = materialTable[key].get(inner_keys::Material::Id::Lambda)->getDeviceDataPtr();
      kernels::time::aux::computeInvAcousticImpedance(invImpedances,
                                                      rhos,
                                                      lambdas,
                                                      numElements,
                                                      deviceStream);

      double factorEvaluated = 1;
      double factorInt = deltaTInt;
      const double g = getGravitationalAcceleration();

      auto** derivativesPtrs = dataTable[key].get(inner_keys::Wp::Id::Derivatives)->getDeviceDataPtr();
      for (int order = 1; order < ConvergenceOrder+1; ++order) {

        factorEvaluated *= deltaT / (1.0 * order);
        factorInt *= deltaTInt / (order + 1.0);

        device.algorithms.setToValue(dofsFaceNodalPtrs,
                                     0.0,
                                     tensor::INodal::size(),
                                     numElements,
                                     deviceStream);

        auto projectKernel = projectKernelPrototype;
        projectKernel.numElements = numElements;
        projectKernel.Tinv = const_cast<const RealT**>(TinvDataPtrs);
        projectKernel.INodal = dofsFaceNodalPtrs;
        projectKernel.linearAllocator.initialize(auxTmpMem);
        projectKernel.streamPtr = deviceStream;

        auto* derivativesOffsets = timeKernel.getDerivativesOffsets();
        for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
          projectKernel.dQ(i) = const_cast<const RealT **>(derivativesPtrs);
          projectKernel.extraOffset_dQ(i) = derivativesOffsets[i];
        }

        projectKernel.execute(order - 1, faceIdx);

        kernels::time::aux::updateRotatedFaceDisplacement(
            rotatedFaceDisplacementPtrs,
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
      rotateFaceDisplacementKrnl.faceDisplacement = const_cast<const RealT**>(rotatedFaceDisplacementPtrs);
      rotateFaceDisplacementKrnl.displacementRotationMatrix = const_cast<const RealT**>(rotateDisplacementToGlobalPtrs);
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

#endif //SEISSOL_GRAVIATIONALFREESURFACEBC_H
