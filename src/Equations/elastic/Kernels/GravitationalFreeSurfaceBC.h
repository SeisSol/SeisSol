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
  GravitationalFreeSurfaceBc()
      : odeSolver(ode::RungeKuttaODESolver(
      {init::averageNormalDisplacement::size(),
       init::faceDisplacement::size()}, ode::ODESolverConfig(1.0)
  )) {}

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
                real* displacementNodal,
                real* integratedDisplacementNodal,
                TimeKrnl& timeKernel,
                real* derivatives,
                double startTime,
                double timeStepWidth,
                CellMaterialData& materialData,
                FaceType faceType) {
    // This function does two things:
    // 1: Compute eta (for all three dimensions) at the end of the timestep
    // 2: Compute the integral of eta in normal direction (called H) over the timestep
    // TODO(Lukas) Values below only correct for acoustic!
    // We do this by solving the following set of ODEs (described only for 1 point here!)
    // pressureAtBnd = 0 (for free surface)
    // pressureAtBnd = -\rho \eta g (for gravitational free surface)
    // dH/dt = etaN,
    // dEtaN/dt =  u_r + 1/Z ( pressureAtBnd - pressureInside ),
    // dEtaS/dt = v_r
    // dEtaT/dt = w_r
    // where the vales u_r/p_r are boundary extrapolated from the interior
    // The initial conditions are H(t^n) = 0, eta(t^n) = eta_0

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
    // Rotate face displacement to face-normal coordinate system in which the computation is
    // more convenient.
    auto rotateFaceDisplacementKrnl = kernel::rotateFaceDisplacement();
    rotateFaceDisplacementKrnl.faceDisplacement = displacementNodal;
    rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToFaceNormalData;
    rotateFaceDisplacementKrnl.rotatedFaceDisplacement = rotatedFaceDisplacement;
    rotateFaceDisplacementKrnl.execute();

    // Temporary buffer to store dofs at some time t
    alignas(ALIGNMENT) real dofsVolumeInteriorModalStorage[tensor::I::size()];
    auto dofsVolumeInteriorModal = init::I::view::create(dofsVolumeInteriorModalStorage);

    // Temporary buffer to store nodal face dofs at some time t
    alignas(ALIGNMENT) real dofsFaceNodalStorage[tensor::INodal::size()];
    auto dofsFaceNodal = init::INodal::view::create(dofsFaceNodalStorage);

    auto f = [&](
        ode::ODEVector& du,
        ode::ODEVector& u,
        double time) {
      // Evaluate Taylor series at time
      dofsVolumeInteriorModal.setZero();
      dofsFaceNodal.setZero();

      timeKernel.computeTaylorExpansion(time,
                                        startTime,
                                        derivatives,
                                        dofsVolumeInteriorModal.data());

      // Project to face and rotate into face-normal aligned coordinate system.
      projectKernel.I = dofsVolumeInteriorModal.data();
      projectKernel.INodal = dofsFaceNodal.data();
      projectKernel.execute(faceIdx);

      // Unpack du
      auto[dEtaIntegratedStorage, dEtaIntegratedSize] = du.getSubvector(0);
      assert(dEtaIntegratedSize == init::averageNormalDisplacement::size());
      auto dEtaIntegrated = init::averageNormalDisplacement::view::create(dEtaIntegratedStorage);
      auto[dEtaStorage, dEtaSize] = du.getSubvector(1);
      assert(dEtaSize == init::faceDisplacement::size());
      auto dEta = init::faceDisplacement::view::create(dEtaStorage);

      // Unpack u
      auto[etaStorage, etaSize] = u.getSubvector(1);
      assert(etaSize == init::faceDisplacement::size());
      auto eta = init::faceDisplacement::view::create(etaStorage);

      constexpr int pIdx = 0;
      constexpr int uIdx = 6;

      dEta.setZero();
      dEtaIntegrated.setZero();
      for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
        const double rho = materialData.local.rho;
        const double g = getGravitationalAcceleration(); // [m/s^2]
        double pressureAtBnd = 0;
        if (faceType == FaceType::freeSurfaceGravity) {
          pressureAtBnd = -1 * rho * g * eta(i, 0);
        }
#ifdef USE_ELASTIC
        const double Z = std::sqrt(materialData.local.lambda * rho) ;
        bool isAcoustic = std::abs(materialData.local.mu) < 1e-15;
#else
        const real Z = 0.0; // Sets penalty term effectively to zero
        bool isAcoustic = false;
#endif

        // dH/dt = etaN,
        // dEtaN/dt =  u_r - 1/Z ( p_r - pressureAtBnd ),
        const auto uInside = dofsFaceNodal(i, uIdx + 0);
        const auto vInside = dofsFaceNodal(i, uIdx + 1);
        const auto wInside = dofsFaceNodal(i, uIdx + 2);
        const auto pressureInside = dofsFaceNodal(i, pIdx);
        if (isAcoustic
            && (faceType == FaceType::freeSurface || faceType == FaceType::freeSurfaceGravity)) {
          dEta(i, 0) = uInside + 1 / Z * (pressureAtBnd - pressureInside);
        } else {
          // If not on acoustic boundary, just use values inside
          dEta(i, 0) = uInside;
        }
        dEta(i, 1) = vInside;
        dEta(i, 2) = wInside;
        dEtaIntegrated(i) = eta(i, 0);
      }
    };

    constexpr auto integratedEtaSize = init::averageNormalDisplacement::size();
    constexpr auto etaSize = init::faceDisplacement::size();

    auto curValue = ode::ODEVector{{integratedDisplacementNodal, rotatedFaceDisplacement},
                                   {integratedEtaSize,           etaSize}};

    // Apply initial condition to integrated displacement (start from 0 each PDE timestep)
    std::fill_n(integratedDisplacementNodal, integratedEtaSize, 0.0);

    // Setup ODE solver
    ode::TimeSpan timeSpan = {startTime, startTime + timeStepWidth};
    auto odeSolverConfig = ode::ODESolverConfig(timeStepWidth);
    odeSolverConfig.initialDt = timeStepWidth;

    odeSolver.setConfig(odeSolverConfig);
    odeSolver.solve(f, curValue, timeSpan);

    // Rotate face displacement back to global coordinate system which we use as storage
    // coordinate system
    rotateFaceDisplacementKrnl.faceDisplacement = rotatedFaceDisplacement;
    rotateFaceDisplacementKrnl.displacementRotationMatrix = rotateDisplacementToGlobalData;
    rotateFaceDisplacementKrnl.rotatedFaceDisplacement = displacementNodal;
    rotateFaceDisplacementKrnl.execute();
  }

private:
  ode::RungeKuttaODESolver odeSolver;
};

} // namespace seissol

#endif //SEISSOL_GRAVIATIONALFREESURFACEBC_H
