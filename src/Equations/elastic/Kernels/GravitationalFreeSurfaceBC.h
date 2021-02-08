#ifndef SEISSOL_GRAVIATIONALFREESURFACEBC_H
#define SEISSOL_GRAVIATIONALFREESURFACEBC_H

#include "generated_code/init.h"
#include "generated_code/kernel.h"
#include "generated_code/tensor.h"

#include "Initializer/typedefs.hpp"

#include "Numerical_aux/Quadrature.h"
#include "Numerical_aux/ODEInt.h"
#include "Kernels/Time.h"

namespace seissol {
class GravitationalFreeSurfaceBc {
public:
  GravitationalFreeSurfaceBc() = default;


  template<typename MappingKrnl>
  void evaluate(unsigned faceIdx,
                MappingKrnl&& projectKernelPrototype,
                const CellBoundaryMapping& boundaryMapping,
                real* displacementNodal,
                real* integratedDisplacementNodal,
                seissol::kernels::Time& timeKernel,
                real* derivatives,
                double startTime,
                double timeStepWidth,
                CellMaterialData& materialData,
                FaceType faceType) const {
    // This function does two things:
    // 1: Compute eta (for all three dimensions) at the end of the timestep
    // 2: Compute the integral of eta in normal direction (called H) over the timestep
    // TODO(Lukas) Values below only correct for acoustic!
    // We do this by solving the following set of ODEs (described only for 1 point here!)
    // pressureAtBnd = 0 (for free surface)
    // pressureAtBnd = \rho \eta g (for gravitational free surface)
    // dH/dt = etaN,
    // dEtaN/dt =  u_r - 1/Z ( p_r - pressureAtBnd ),
    // dEtaS/dt = v_r
    // dEtaT/dt = w_r
    // where the vales u_r/p_r are boundary extrapolated from the interior
    // The initial conditions are H(t^n) = 0, eta(t^n) = eta_0

    // Prepare kernel that projects volume data to face and rotates it to face-nodal basis.
    assert(boundaryMapping.nodes != nullptr);
    assert(boundaryMapping.TinvData != nullptr);
    auto projectKernel = projectKernelPrototype;
    projectKernel.Tinv = boundaryMapping.TinvData;

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
        const double g = 9.81; // [m/s^2]
        double pressureAtBnd = 0;
        if (faceType == FaceType::freeSurfaceGravity) {
          pressureAtBnd = -1 * rho * g * eta(i, 0);
        }
#ifdef USE_ELASTIC
        const double Z = std::sqrt(materialData.local.lambda / rho);
        bool isAcoustic = std::abs(materialData.local.mu) < 1e-15;
#else
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
          // TODO(Lukas) Check sign
          dEta(i, 0) = uInside + 1 / Z * (pressureInside - pressureAtBnd);
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

    auto curValue = ode::ODEVector{{integratedDisplacementNodal, displacementNodal},
                                   {integratedEtaSize,           etaSize}};

    // Apply boundary condition to integrated displacement (start from 0 each PDE timestep)
    std::fill_n(integratedDisplacementNodal, integratedEtaSize, 0.0);

    // Setup ODE solver
    ode::TimeSpan timeSpan = {startTime, startTime + timeStepWidth};
    auto odeSolverConfig = ode::ODESolverConfig(timeStepWidth);
    odeSolverConfig.initialDt = timeStepWidth;

    auto solver = ode::RungeKuttaODESolver({integratedEtaSize, etaSize}, odeSolverConfig);
    solver.solve(f, curValue, timeSpan);
  }

private:
  //ode::ODESolver odeSolver;
};

} // namespace seissol

#endif //SEISSOL_GRAVIATIONALFREESURFACEBC_H
