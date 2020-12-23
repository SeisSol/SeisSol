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
  void evaluate(int faceIdx,
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
    // Prepare kernel that projects volume data to face and rotates it to face-nodal basis.
    static_assert(nodal::tensor::nodes2D::Shape[0] == tensor::INodal::Shape[0],
                  "Need evaluation at all nodes!");
    assert(boundaryMapping.nodes != nullptr);
    assert(boundaryMapping.TinvData != nullptr);
    auto projectKernel = projectKernelPrototype;
    projectKernel.Tinv = boundaryMapping.TinvData;
    alignas(ALIGNMENT) real dofsFaceNodal[tensor::INodal::size()];
    auto boundaryDofs = init::INodal::view::create(dofsFaceNodal);

    // This function does two things:
    // 1: Compute eta (for all three dimensions) at the end of the timestep
    // 2: Compute the integral of eta in normal direction (called H) over the timestep
    //alignas(ALIGNMENT) real faceDisplacementData[tensor::faceDisplacement::size()];
    //auto displacement = init::faceDisplacement::view::create(faceDisplacementData);

    //alignas(ALIGNMENT) real averageFaceNormalDisplacementData[tensor::averageNormalDisplacement::size()];
    //auto averageNormalFaceDisplacement = init::averageNormalDisplacement::view::create(averageFaceNormalDisplacementData);
    //

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

    // TODO(Lukas) Optimize
    // Right now this is hacked together:
    // Our ODE solver interface only allows to update one vector at the same time
    // This means that we have to glue both parts together.
    // We can avoid many of the copies by being smarter.
    // TODO(Lukas) Wrong
    //alignas(ALIGNMENT) real odeSolverOutputData[tensor::INodal::size()];
    //auto odeSolverOutputView = init::INodal::view::create(dofsFaceNodal);

    alignas(ALIGNMENT) real dofsVolumeInteriorModal[tensor::I::size()];

    // TODO: Improve capturing

    auto f = [&](
        ODEVector& du,
        ODEVector& u,
        double time) {
        // Evaluate Taylor series at time t
        timeKernel.computeTaylorExpansion(time,
                                          startTime,
                                          derivatives,
                                          dofsVolumeInteriorModal);

        projectKernel.I = dofsVolumeInteriorModal;
        projectKernel.INodal = dofsFaceNodal;
        projectKernel.execute(faceIdx);
        auto dofsFaceNodalView = init::INodal::view::create(dofsFaceNodal);

        // Unpack du
        auto [dEtaIntegratedStorage, dEtaIntegratedSize] = du.getSubvector(0);
        assert(dEtaIntegratedSize == init::averageNormalDisplacement::size());
        auto dEtaIntegrated = init::averageNormalDisplacement::view::create(dEtaIntegratedStorage);
        auto [dEtaStorage, dEtaSize] = du.getSubvector(1);
        assert(dEtaSize == init::faceDisplacement::size());
        auto dEta = init::faceDisplacement::view::create(dEtaStorage);

        // Unpack u
        auto [etaIntegratedStorage, etaIntegratedSize] = u.getSubvector(0);
        assert(etaIntegratedSize == init::averageNormalDisplacement::size());
        auto etaIntegrated = init::averageNormalDisplacement::view::create(etaIntegratedStorage);
        auto [etaStorage, etaSize] = u.getSubvector(1);
        assert(etaSize == init::faceDisplacement::size());
        auto eta = init::faceDisplacement::view::create(dEtaStorage);

        // Evaluate function evaluation for eta
        // TODO(Lukas) Check indices
        constexpr int uIdx = 6;
        constexpr int pIdx = 0;

        for (unsigned int i = 0; i < nodal::tensor::nodes2D::Shape[0]; ++i) {
          const double rho = materialData.local.rho;
          const double g = 9.81; // [m/s^2]
          double pressureAtBnd = 0;
          if (faceType == FaceType::freeSurfaceGravity) {
            pressureAtBnd = -1 * rho * g * eta(i, 0);
          }
          const double Z = std::sqrt(materialData.local.lambda / rho);

          // dH/dt = etaN,
          // dEtaN/dt =  u_r - 1/Z ( p_r - pressureAtBnd ),
          const auto uInside = dofsFaceNodalView(i, uIdx+0);
          const auto vInside = dofsFaceNodalView(i, uIdx+1);
          const auto wInside = dofsFaceNodalView(i, uIdx+2);
          const auto pressureInside = dofsFaceNodalView(i, pIdx);
          if (faceType == FaceType::freeSurface || faceType == FaceType::freeSurfaceGravity) {
            dEta(i,0) = uInside - 1/Z * (pressureInside - pressureAtBnd);
          } else {
            // Elastic-acoustic interface => no penalty term
            dEta(i,0) = uInside;
          }
          dEta(i,1) = vInside;
          dEta(i,2) = wInside;
          dEtaIntegrated(i) = eta(i,0);
      }


    };

    // Allocate storage
    // TODO(Lukas) Do this in ODEInt!
    constexpr auto integratedEtaSize = init::averageNormalDisplacement::size();
    constexpr auto etaSize = init::faceDisplacement::size();

    alignas(ALIGNMENT) real fEvalIntegrated[integratedEtaSize];
    alignas(ALIGNMENT) real stage1Integrated[integratedEtaSize];
    alignas(ALIGNMENT) real fEvalHeunIntegrated[integratedEtaSize];
    alignas(ALIGNMENT) real updateHeunIntegrated[integratedEtaSize];

    alignas(ALIGNMENT) real fEvalEta[etaSize];
    alignas(ALIGNMENT) real stage1Eta[etaSize];
    alignas(ALIGNMENT) real fEvalHeunEta[etaSize];
    alignas(ALIGNMENT) real updateHeunEta[etaSize];

    auto fEval = ODEVector{{fEvalIntegrated, fEvalEta}, {integratedEtaSize, etaSize}};
    auto stage1 = ODEVector{{stage1Integrated, stage1Eta}, {integratedEtaSize, etaSize}};
    auto fEvalHeun = ODEVector{{fEvalHeunIntegrated, fEvalHeunEta}, {integratedEtaSize, etaSize}};
    auto updateHeun = ODEVector{{updateHeunIntegrated, updateHeunEta}, {integratedEtaSize, etaSize}};

    // TODO(Lukas) Don't allocate this, we should actually update the value...
    //alignas(ALIGNMENT) real curValueIntegrated[integratedEtaSize];
    //alignas(ALIGNMENT) real curValueEta[etaSize];
    //auto curValue = ODEVector{{curValueIntegrated, curValueEta}, {integratedEtaSize, etaSize}};
    auto curValue = ODEVector{{ integratedDisplacementNodal, displacementNodal }, {integratedEtaSize, etaSize}};
    auto integratedDisplacementNodalView = init::averageNormalDisplacement::view::create(integratedDisplacementNodal);
    integratedDisplacementNodalView.setZero();

    // Setup ODE solver
    ode::TimeSpan timeSpan = {startTime, startTime + timeStepWidth};
    auto odeSolverConfig = ode::ODESolverConfig(timeStepWidth);
    odeSolverConfig.initialDt = timeStepWidth / 2; // TODO(Lukas) Be smarter!
    odeSolverConfig.acceptableError = 1e-5;
    odeSolverConfig.minimumDt = 1e-8;

    auto solver = ode::ODESolver(fEval, stage1, fEvalHeun, updateHeun, odeSolverConfig);
    solver.solve(f, curValue, timeSpan);


  }

private:
  //ode::ODESolver odeSolver;
};

} // namespace seissol

#endif //SEISSOL_GRAVIATIONALFREESURFACEBC_H
