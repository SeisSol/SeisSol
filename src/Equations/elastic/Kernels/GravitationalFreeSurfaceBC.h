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
    const auto config = ode::ODESolverConfig(1.0); // Use default config
    int numStages = ode::getNumberOfStages(config.solver);
    const auto numDofs = init::faceDisplacement::size() + init::averageNormalDisplacement::size();

    Eigen::MatrixXd a;
    Eigen::VectorXd b;
    Eigen::VectorXd c;
    initializeRungeKuttaScheme(config.solver, numStages, a, b, c);

    long long nonZeroFlopsTaylor, hardwareFlopsTaylor;
    timeKrnl.flopsTaylorExpansion(nonZeroFlopsTaylor, hardwareFlopsTaylor);

    const auto nonZeroFlopsFunctionEvaluation =
        nonZeroFlopsTaylor +
        kernel::projectToNodalBoundaryRotated::NonZeroFlops[face];
    const auto hardwareFlopsFunctionEvaluation =
        hardwareFlopsTaylor +
        kernel::projectToNodalBoundaryRotated::HardwareFlops[face];

    const auto nonZeroFlopsFunctionEvaluations = numStages * nonZeroFlopsFunctionEvaluation;
    const auto hardwareFlopsFunctionEvaluations = numStages * hardwareFlopsFunctionEvaluation;

    const auto intermediateStages = a.count();
    const auto flopsRKStages = intermediateStages * numDofs * 2; // One mul to scale with a_{ij} \Delta t, one add
    const auto flopsRKFinalValue = numStages * 2 * numDofs; // One mul to scale with b \Delta t, one add

    const auto hardwareFlopsRK = flopsRKStages + flopsRKFinalValue;
    const auto nonZeroFlopsRK = hardwareFlopsRK;

    const auto nonZeroFlops = nonZeroFlopsFunctionEvaluations + nonZeroFlopsRK;
    const auto hardwareFlops = hardwareFlopsFunctionEvaluations + hardwareFlopsRK;

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
    // 2: Compute the integral of eta in normal direction (called H) over the timestep

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
    alignas(ALIGNMENT) real rotatedFaceDisplacement[init::rotatedFaceDisplacement::Size];

    auto integratedDisplacementNodal = init::averageNormalDisplacement::view::create(integratedDisplacementNodalData);
    auto displacementNodal = init::faceDisplacement::view::create(rotatedFaceDisplacement);

    // Rotate face displacement to face-normal coordinate system in which the computation is
    // more convenient.
    auto rotateFaceDisplacementKrnl = kernel::rotateFaceDisplacement();
    rotateFaceDisplacementKrnl.faceDisplacement = displacementNodalData;
    rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
    rotateFaceDisplacementKrnl.rotatedFaceDisplacement = rotatedFaceDisplacement;
    rotateFaceDisplacementKrnl.execute();

    // Temporary buffer to store dofs at some time t
    alignas(ALIGNMENT) real dofsVolumeInteriorModalStorage[tensor::I::size()];
    auto dofsVolumeInteriorModal = init::I::view::create(dofsVolumeInteriorModalStorage);

    // Temporary buffer to store nodal face dofs at some time t
    alignas(ALIGNMENT) real dofsFaceNodalStorage[tensor::INodal::size()];
    auto dofsFaceNodal = init::INodal::view::create(dofsFaceNodalStorage);

    // Temporary buffer to store nodal face coefficients at some time t
    // TODO(Lukas) We actually only need one coeff per node! -> std::array
    alignas(ALIGNMENT) real prevCoefficientsStorage[tensor::INodal::size()];
    auto prevCoefficients = init::INodal::view::create(prevCoefficientsStorage);

    const double deltaT = timeStepWidth;
    const double deltaTInt = timeStepWidth;
    // Initialize first component of series
    //evaluated = eta_ic
    //evaluated_int =  delta_t_int * eta_ic

    for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
      prevCoefficients(i,0) = displacementNodal(i, 0);
      integratedDisplacementNodal(i) = deltaTInt * displacementNodal(i,0);
    }

    // Coefficients for Taylor series
    double factorEvaluated = 1;
    double factorInt = deltaTInt;

    // Note: Probably need to increase CONVERGENCE_ORDER by 1 here!
    for (int order = 1; order < CONVERGENCE_ORDER+1; ++order) {
      // TODO(Lukas): Actually just picks out the nth coefficient...
      dofsVolumeInteriorModal.setZero();
      dofsFaceNodal.setZero();
      timeKernel.computeDerivativeTaylorExpansion(startTime,
                                        startTime,
                                        derivatives,
                                        dofsVolumeInteriorModal.data(),
                                        order-1);

      projectKernel.I = dofsVolumeInteriorModal.data();
      projectKernel.INodal = dofsFaceNodal.data();
      projectKernel.execute(faceIdx);

      factorEvaluated *= deltaT / (1.0 * order);
      factorInt *= deltaTInt / (order + 1.0);
      for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
        const double rho = materialData.local.rho;
        const double g = getGravitationalAcceleration(); // [m/s^2]
        const double Z = std::sqrt(materialData.local.lambda * rho) ;

        // Derivatives of interior variables
        constexpr int pIdx = 0;
        constexpr int uIdx = 6;

        const auto uInside = dofsFaceNodal(i, uIdx + 0);
        const auto vInside = dofsFaceNodal(i, uIdx + 1);
        const auto wInside = dofsFaceNodal(i, uIdx + 2);
        const auto pressureInside = dofsFaceNodal(i, pIdx);

        const double prevCoeff = prevCoefficients(i, 0);
        // TODO(Lukas) Check sign
        const double curCoeff = uInside - (1.0/Z) * (rho * g * prevCoeff + pressureInside);
        prevCoefficients(i,0) = curCoeff;

        displacementNodal(i,0) += factorEvaluated * curCoeff;
        //TODO(Lukas) Check if Taylor series for displ is correct
        displacementNodal(i,1) += factorEvaluated * vInside;
        displacementNodal(i,2) += factorEvaluated * wInside;

        integratedDisplacementNodal(i) += factorInt * curCoeff;
      }

    }

    // Rotate face displacement back to global coordinate system which we use as storage
    // coordinate system
    rotateFaceDisplacementKrnl.faceDisplacement = rotatedFaceDisplacement;
    rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToGlobalData;
    rotateFaceDisplacementKrnl.rotatedFaceDisplacement = displacementNodalData;
    rotateFaceDisplacementKrnl.execute();
  }
};

} // namespace seissol

#endif //SEISSOL_GRAVIATIONALFREESURFACEBC_H
