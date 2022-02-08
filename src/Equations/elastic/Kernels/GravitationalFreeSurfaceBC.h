#ifndef SEISSOL_GRAVIATIONALFREESURFACEBC_H
#define SEISSOL_GRAVIATIONALFREESURFACEBC_H

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/typedefs.hpp"

#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/ODEInt.h"

namespace seissol {

// Used to avoid including SeisSo.h here as this leads to all sorts of issues
double getGravitationalAcceleration();

class GravitationalFreeSurfaceBc {
public:
  GravitationalFreeSurfaceBc() { }

  template <typename TimeKrnl>
  static std::pair<long long, long long> getFlopsDisplacementFace(unsigned face,
                                                                  [[maybe_unused]] FaceType faceType,
                                                                  TimeKrnl& timeKrnl) {
    const auto numDofs = init::faceDisplacement::size() + init::averageNormalDisplacement::size();
    long long hardwareFlops = 0;
    long long nonZeroFlops = 0;

    // Note: This neglects roughly 10 * CONVERGENCE_ORDER * numNodes2D flops
    // Before adjusting the range of the loop, check range of loop in computation!
    for (int order = 1; order < CONVERGENCE_ORDER + 1; ++order) {
      nonZeroFlops += kernel::projectDerivativeToNodalBoundaryRotated::nonZeroFlops(order - 1, face);
      hardwareFlops += kernel::projectDerivativeToNodalBoundaryRotated::hardwareFlops(order - 1, face);
    }

    return {nonZeroFlops, hardwareFlops};
  }

  template<typename TimeKrnl, typename MappingKrnl>
  void evaluate(unsigned faceIdx,
                MappingKrnl&& projectKernelPrototype,
                const CellBoundaryMapping& boundaryMapping,
                real* displacementNodalData,
                real* integratedDisplacementNodalData,
                TimeKrnl& timeKernel,
                real* derivatives,
                double startTime,
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
    auto Tinv = init::Tinv::view::create(boundaryMapping.TinvData);
    auto T = init::Tinv::view::create(boundaryMapping.TData);
    auto projectKernel = projectKernelPrototype;
    projectKernel.Tinv = Tinv.data();

    // Prepare projection of displacement/velocity to face-nodal basis.
    alignas(ALIGNMENT) real rotateDisplacementToFaceNormalData[init::displacementRotationMatrix::Size];
    auto rotateDisplacementToFaceNormal = init::displacementRotationMatrix::view::create(rotateDisplacementToFaceNormalData);
    alignas(ALIGNMENT) real rotateDisplacementToGlobalData[init::displacementRotationMatrix::Size];
    auto rotateDisplacementToGlobal = init::displacementRotationMatrix::view::create(rotateDisplacementToGlobalData);
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        // Extract part that rotates velocity from T
        rotateDisplacementToFaceNormal(i, j) = Tinv(i + 6, j + 6);
        rotateDisplacementToGlobal(i,j) = T(i+6, j+6);
      }
    }
    static_assert(init::rotatedFaceDisplacement::Size == init::faceDisplacement::Size);
    alignas(ALIGNMENT) real rotatedFaceDisplacementData[init::rotatedFaceDisplacement::Size];

    auto integratedDisplacementNodal = init::averageNormalDisplacement::view::create(integratedDisplacementNodalData);
    auto rotatedFaceDisplacement = init::faceDisplacement::view::create(rotatedFaceDisplacementData);

    // Rotate face displacement to face-normal coordinate system in which the computation is
    // more convenient.
    auto rotateFaceDisplacementKrnl = kernel::rotateFaceDisplacement();
    rotateFaceDisplacementKrnl.faceDisplacement = displacementNodalData;
    rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
    rotateFaceDisplacementKrnl.rotatedFaceDisplacement = rotatedFaceDisplacementData;
    rotateFaceDisplacementKrnl.execute();

    // Temporary buffer to store dofs at some time t
    alignas(ALIGNMENT) real dofsVolumeInteriorModalStorage[tensor::I::size()];
    auto dofsVolumeInteriorModal = init::I::view::create(dofsVolumeInteriorModalStorage);

    // Temporary buffer to store nodal face dofs at some time t
    alignas(ALIGNMENT) real dofsFaceNodalStorage[tensor::INodal::size()];
    auto dofsFaceNodal = init::INodal::view::create(dofsFaceNodalStorage);

    // Temporary buffer to store nodal face coefficients at some time t
    alignas(ALIGNMENT) std::array<real, NUMBER_OF_ALIGNED_BASIS_FUNCTIONS> prevCoefficients;

    const double deltaT = timeStepWidth;
    const double deltaTInt = timeStepWidth;

    // Initialize first component of Taylor series
    for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
      prevCoefficients[i] = rotatedFaceDisplacement(i, 0);
      // This is clearly a zeroth order approximation of the integral!
      integratedDisplacementNodal(i) = deltaTInt * rotatedFaceDisplacement(i, 0);
    }

    // Coefficients for Taylor series
    double factorEvaluated = 1;
    double factorInt = deltaTInt;

    auto* derivativesOffsets = timeKernel.getDerivativesOffsets();
    projectKernel.INodal = dofsFaceNodal.data();
    for (unsigned i = 0; i < yateto::numFamilyMembers<tensor::dQ>(); ++i) {
      projectKernel.dQ(i) = derivatives + derivativesOffsets[i];
    }

#ifdef USE_ELASTIC
    const double rho = materialData.local.rho;
    const double g = getGravitationalAcceleration(); // [m/s^2]
    const double Z = std::sqrt(materialData.local.lambda * rho) ;
#endif

    // Note: Probably need to increase CONVERGENCE_ORDER by 1 here!
    for (int order = 1; order < CONVERGENCE_ORDER+1; ++order) {
      dofsVolumeInteriorModal.setZero();
      dofsFaceNodal.setZero();

      projectKernel.execute(order - 1, faceIdx);

      factorEvaluated *= deltaT / (1.0 * order);
      factorInt *= deltaTInt / (order + 1.0);

      for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
        // Derivatives of interior variables
        constexpr int pIdx = 0;
        constexpr int uIdx = 6;

        const auto uInside = dofsFaceNodal(i, uIdx + 0);
        const auto vInside = dofsFaceNodal(i, uIdx + 1);
        const auto wInside = dofsFaceNodal(i, uIdx + 2);
        const auto pressureInside = dofsFaceNodal(i, pIdx);

#ifdef USE_ELASTIC
        const double curCoeff = uInside - (1.0/Z) * (rho * g * prevCoefficients[i] + pressureInside);
#else
        const double curCoeff = uInside;
#endif
        prevCoefficients[i] = curCoeff;

        rotatedFaceDisplacement(i, 0) += factorEvaluated * curCoeff;
        rotatedFaceDisplacement(i, 1) += factorEvaluated * vInside;
        rotatedFaceDisplacement(i, 2) += factorEvaluated * wInside;

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
};

} // namespace seissol

#endif //SEISSOL_GRAVIATIONALFREESURFACEBC_H
