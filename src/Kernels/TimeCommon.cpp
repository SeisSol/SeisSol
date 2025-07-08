// SPDX-FileCopyrightText: 2016 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Carsten Uphoff

#include "TimeCommon.h"
#include <Common/Constants.h>
#include <DataTypes/ConditionalTable.h>
#include <Initializer/BasicTypedefs.h>
#include <Kernels/Precision.h>
#include <Kernels/Solver.h>
#include <Parallel/Runtime/Stream.h>
#include <cassert>
#include <cstddef>
#include <stdint.h>
#include <tensor.h>

#include "utils/logger.h"

#ifdef ACL_DEVICE
#include <DataTypes/ConditionalKey.h>
#include <DataTypes/EncodedConstants.h>
#endif

#ifndef NDEBUG
#include "Alignment.h"
#include <cstdint>
#endif

namespace seissol::kernels {

#ifdef USE_DAMAGE
void TimeCommon::computeNonIntegrals(Time& time,
                                  LocalData& data,
                                  unsigned short ltsSetup,
                                  const FaceType faceTypes[4],
                                  const double currentTime[5],
                                  double timeStepWidth,
                                  real* const timeDofs[4],
                                  real integrationBuffer[4][seissol::kernels::Solver::BufferSize],
                                  real* timeIntegrated[4]) {
  /*
   * assert valid input.
   */
  // only lower 10 bits are used for lts encoding
  assert(ltsSetup < 2048);

#ifndef NDEBUG
  // alignment of the time derivatives/integrated dofs and the buffer
  for (int dofneighbor = 0; dofneighbor < 4; dofneighbor++) {
    assert(reinterpret_cast<uintptr_t>(timeDofs[dofneighbor]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(integrationBuffer[dofneighbor]) % Alignment == 0);
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for (unsigned dofneighbor = 0; dofneighbor < 4; ++dofneighbor) {
    // collect information only in the case that neighboring element contributions are required
    if (faceTypes[dofneighbor] != FaceType::Outflow &&
        faceTypes[dofneighbor] != FaceType::DynamicRupture) {
      // check if the time integration is already done (-> copy pointer)
      if ((ltsSetup >> dofneighbor) % 2 == 0) {
        timeIntegrated[dofneighbor] = timeDofs[dofneighbor];
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      // Zihua: In nonlinear case, integrate Q, Fx, Fy and Fz again from the 
      // derivatives, which is timeDofs[dofneighbor];
      // Store it in integrationBuffer[dofneighbor] first,
      // Then pass it to timeIntegrated[dofneighbor]
      else {
        // Step 1: Convert from Modal to Nodal for each temporal quadrature point;:
        // Meanwhile, compute the info for Rusanov fluxes for neighboring cells,
        // including the integrated Fx, Fy, Fz
        // Note: need to consider start and end of temporal integration, as in 
        // the original time.computeIntegral() for I here
        
        // Step 1.1: Integrate I as in the linear case.
        time.computeIntegral(currentTime[dofneighbor + 1],
                             currentTime[0],
                             currentTime[0] + timeStepWidth,
                             timeDofs[dofneighbor],
                             integrationBuffer[dofneighbor]);

        // Step 1.2: Integrate Fx, Fy, 
        // The integration needs to be sub-timestep
        double timePoints[ConvergenceOrder];
        double timeWeights[ConvergenceOrder];
        seissol::quadrature::GaussLegendre(timePoints, timeWeights, ConvergenceOrder);
        // time difference between substep begins and expansion point:
        // currentTime[0] is subTimeStart, which is already the dt from lastSubTIme
        double dtStart = currentTime[0] - currentTime[dofneighbor + 1];
        for (unsigned int point = 0; point < ConvergenceOrder; ++point) {
          timePoints[point] = dtStart + 0.5 * (timeStepWidth * timePoints[point] + timeStepWidth);
          timeWeights[point] = 0.5 * timeStepWidth * timeWeights[point];
        }

        alignas(PagesizeStack) real QInterpolatedBody[ConvergenceOrder][tensor::Q::size()] = {{0.0}};
        real* QInterpolatedBodyi;
        alignas(PagesizeStack) real QInterpolatedBodyNodal[ConvergenceOrder][tensor::QNodal::size()];
        real* QInterpolatedBodyNodali;

        kernel::damageConvertToNodal d_converToKrnl;
        for (unsigned int timeInterval = 0; timeInterval < ConvergenceOrder; ++timeInterval){
          QInterpolatedBodyi = QInterpolatedBody[timeInterval];
          QInterpolatedBodyNodali = QInterpolatedBodyNodal[timeInterval];
          time.computeTaylorExpansion(timePoints[timeInterval], 0.0, timeDofs[dofneighbor], QInterpolatedBodyi);
          /// Convert Q_{lp}(tau_z) in modal basis to QN_{ip}(tau_z) in nodal basis
          d_converToKrnl.v = init::v::Values;
          d_converToKrnl.QNodal = QInterpolatedBodyNodali;
          d_converToKrnl.Q = QInterpolatedBodyi;
          d_converToKrnl.execute();
        }

        // Step 2: Do time integration for the Rusanov flux
        real* integratedFx = integrationBuffer[dofneighbor] + 1*tensor::I::size();
        real* integratedFy = integrationBuffer[dofneighbor] + 2*tensor::I::size();
        real* integratedFz = integrationBuffer[dofneighbor] + 3*tensor::I::size();

        alignas(Alignment) real integratedFxNodal[tensor::QNodal::size()] = {0.0};
        alignas(Alignment) real integratedFyNodal[tensor::QNodal::size()] = {0.0};
        alignas(Alignment) real integratedFzNodal[tensor::QNodal::size()] = {0.0};

        alignas(Alignment) real sxxNodal[tensor::QNodal::Shape[0]] = {0};
        alignas(Alignment) real syyNodal[tensor::QNodal::Shape[0]] = {0};
        alignas(Alignment) real szzNodal[tensor::QNodal::Shape[0]] = {0};
        alignas(Alignment) real sxyNodal[tensor::QNodal::Shape[0]] = {0};
        alignas(Alignment) real syzNodal[tensor::QNodal::Shape[0]] = {0};
        alignas(Alignment) real szxNodal[tensor::QNodal::Shape[0]] = {0};

        real mat_mu0 = data.material().local.mu0;
        real mat_lambda0 = data.material().local.lambda0;
        real mat_Cd = data.material().local.Cd;
        real mat_gammaR = data.material().local.gammaR;

        real epsInitxx = data.material().local.epsInit_xx;
        real epsInityy = data.material().local.epsInit_yy;
        real epsInitzz = data.material().local.epsInit_zz;
        real epsInitxy = data.material().local.epsInit_xy;
        real epsInityz = data.material().local.epsInit_yz;
        real epsInitzx = data.material().local.epsInit_xz;

        for (unsigned timeInterval = 0; timeInterval < ConvergenceOrder; ++timeInterval) {
          // Compute rusanovFluxMinusNodal at each time step and integrate in time
          // with temporal quadrature
          auto weight = timeWeights[timeInterval];
          real* exxNodal = (QInterpolatedBodyNodal[timeInterval] + 0*tensor::Q::Shape[0]);
          real* eyyNodal = (QInterpolatedBodyNodal[timeInterval] + 1*tensor::Q::Shape[0]);
          real* ezzNodal = (QInterpolatedBodyNodal[timeInterval] + 2*tensor::Q::Shape[0]);
          real* alphaNodal = (QInterpolatedBodyNodal[timeInterval] + 9*tensor::Q::Shape[0]);
          real* exyNodal = (QInterpolatedBodyNodal[timeInterval] + 3*tensor::Q::Shape[0]);
          real* eyzNodal = (QInterpolatedBodyNodal[timeInterval] + 4*tensor::Q::Shape[0]);
          real* ezxNodal = (QInterpolatedBodyNodal[timeInterval] + 5*tensor::Q::Shape[0]);
          real* vxNodal = (QInterpolatedBodyNodal[timeInterval] + 6*tensor::Q::Shape[0]);
          real* vyNodal = (QInterpolatedBodyNodal[timeInterval] + 7*tensor::Q::Shape[0]);
          real* vzNodal = (QInterpolatedBodyNodal[timeInterval] + 8*tensor::Q::Shape[0]);
          // time quadrature for each nodal values
          for (unsigned q = 0; q < tensor::Q::Shape[0]; ++q){
            real EspI = (exxNodal[q]+epsInitxx) + (eyyNodal[q]+epsInityy) + (ezzNodal[q]+epsInitzz);
            real EspII = (exxNodal[q]+epsInitxx)*(exxNodal[q]+epsInitxx)
              + (eyyNodal[q]+epsInityy)*(eyyNodal[q]+epsInityy)
              + (ezzNodal[q]+epsInitzz)*(ezzNodal[q]+epsInitzz)
              + 2*(exyNodal[q]+epsInitxy)*(exyNodal[q]+epsInitxy)
              + 2*(eyzNodal[q]+epsInityz)*(eyzNodal[q]+epsInityz)
              + 2*(ezxNodal[q]+epsInitzx)*(ezxNodal[q]+epsInitzx);
            
            // Compute nonlinear flux term
            real mu_eff = mat_mu0 - alphaNodal[q]*mat_mu0;
            sxxNodal[q] = (1-alphaNodal[q])*mat_lambda0*EspI
                          + 2*mu_eff*(exxNodal[q]+epsInitxx);
            syyNodal[q] = (1-alphaNodal[q])*mat_lambda0*EspI
                          + 2*mu_eff*(eyyNodal[q]+epsInityy);
            szzNodal[q] = (1-alphaNodal[q])*mat_lambda0*EspI
                          + 2*mu_eff*(ezzNodal[q]+epsInitzz);
            sxyNodal[q] = 2*mu_eff*(exyNodal[q]+epsInitxy);
            syzNodal[q] = 2*mu_eff*(eyzNodal[q]+epsInityz);
            szxNodal[q] = 2*mu_eff*(ezxNodal[q]+epsInitzx);

            integratedFxNodal[0*tensor::Q::Shape[0] + q] += -vxNodal[q]*weight;
            integratedFxNodal[1*tensor::Q::Shape[0] + q] += 0;
            integratedFxNodal[2*tensor::Q::Shape[0] + q] += 0;
            integratedFxNodal[3*tensor::Q::Shape[0] + q] += -0.5*vyNodal[q]*weight;
            integratedFxNodal[4*tensor::Q::Shape[0] + q] += 0;
            integratedFxNodal[5*tensor::Q::Shape[0] + q] += -0.5*vzNodal[q]*weight;
            integratedFxNodal[6*tensor::Q::Shape[0] + q] += -sxxNodal[q]/data.material().local.rho*weight;
            integratedFxNodal[7*tensor::Q::Shape[0] + q] += -sxyNodal[q]/data.material().local.rho*weight;
            integratedFxNodal[8*tensor::Q::Shape[0] + q] += -szxNodal[q]/data.material().local.rho*weight;
            integratedFxNodal[9*tensor::Q::Shape[0] + q] += 0;

            integratedFyNodal[0*tensor::Q::Shape[0] + q] += 0;
            integratedFyNodal[1*tensor::Q::Shape[0] + q] += -vyNodal[q]*weight;
            integratedFyNodal[2*tensor::Q::Shape[0] + q] += 0;
            integratedFyNodal[3*tensor::Q::Shape[0] + q] += -0.5*vxNodal[q]*weight;
            integratedFyNodal[4*tensor::Q::Shape[0] + q] += -0.5*vzNodal[q]*weight;
            integratedFyNodal[5*tensor::Q::Shape[0] + q] += 0;
            integratedFyNodal[6*tensor::Q::Shape[0] + q] += -sxyNodal[q]/data.material().local.rho*weight;
            integratedFyNodal[7*tensor::Q::Shape[0] + q] += -syyNodal[q]/data.material().local.rho*weight;
            integratedFyNodal[8*tensor::Q::Shape[0] + q] += -syzNodal[q]/data.material().local.rho*weight;
            integratedFyNodal[9*tensor::Q::Shape[0] + q] += 0;

            integratedFzNodal[0*tensor::Q::Shape[0] + q] += 0;
            integratedFzNodal[1*tensor::Q::Shape[0] + q] += 0;
            integratedFzNodal[2*tensor::Q::Shape[0] + q] += -vzNodal[q]*weight;
            integratedFzNodal[3*tensor::Q::Shape[0] + q] += 0;
            integratedFzNodal[4*tensor::Q::Shape[0] + q] += -0.5*vyNodal[q]*weight;
            integratedFzNodal[5*tensor::Q::Shape[0] + q] += -0.5*vxNodal[q]*weight;
            integratedFzNodal[6*tensor::Q::Shape[0] + q] += -szxNodal[q]/data.material().local.rho*weight;
            integratedFzNodal[7*tensor::Q::Shape[0] + q] += -syzNodal[q]/data.material().local.rho*weight;
            integratedFzNodal[8*tensor::Q::Shape[0] + q] += -szzNodal[q]/data.material().local.rho*weight;
            integratedFzNodal[9*tensor::Q::Shape[0] + q] += 0;
          }
        } // end of time integration

        // // Step 3: Convert the INTEGRATED Fx, Fy and Fz flux from Nodal to Modal space
        // Store them in integratedFx, integratedFy, and integratedFz
        kernel::damageAssignFToDQ d_convertBackKrnl;
        d_convertBackKrnl.vInv = init::vInv::Values;
        d_convertBackKrnl.FNodal = integratedFxNodal;
        d_convertBackKrnl.dQModal = integratedFx;
        d_convertBackKrnl.execute();

        d_convertBackKrnl.FNodal = integratedFyNodal;
        d_convertBackKrnl.dQModal = integratedFy;
        d_convertBackKrnl.execute();

        d_convertBackKrnl.FNodal = integratedFzNodal;
        d_convertBackKrnl.dQModal = integratedFz;
        d_convertBackKrnl.execute();

        timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
      }
    }
  }
}

void TimeCommon::computeNonIntegrals(Time& time,
                                  LocalData& data,
                                  unsigned short ltsSetup,
                                  const FaceType faceTypes[4],
                                  const double timeStepStart,
                                  const double timeStepWidth,
                                  real* const timeDofs[4],
                                  real integrationBuffer[4][seissol::kernels::Solver::BufferSize],
                                  real* timeIntegrated[4]) {
  double startTimes[5];
  startTimes[0] = timeStepStart;
  startTimes[1] = startTimes[2] = startTimes[3] = startTimes[4] = 0;

  // adjust start times for GTS on derivatives
  for (unsigned int face = 0; face < 4; face++) {
    if (((ltsSetup >> (face + 4)) % 2) != 0) {
      startTimes[face + 1] = timeStepStart;
    }
  }

  // call the more general assembly
  computeNonIntegrals(time,
                   data,
                   ltsSetup,
                   faceTypes,
                   startTimes,
                   timeStepWidth,
                   timeDofs,
                   integrationBuffer,
                   timeIntegrated);
}
#endif

void TimeCommon::computeIntegrals(Time& time,
                                  unsigned short ltsSetup,
                                  const FaceType faceTypes[4],
                                  const double currentTime[5],
                                  double timeStepWidth,
                                  real* const timeDofs[4],
                                  real integrationBuffer[4][seissol::kernels::Solver::BufferSize],
                                  real* timeIntegrated[4]) {
  /*
   * assert valid input.
   */
  // only lower 10 bits are used for lts encoding
  assert(ltsSetup < 2048);

#ifndef NDEBUG
  // alignment of the time derivatives/integrated dofs and the buffer
  for (std::size_t dofneighbor = 0; dofneighbor < Cell::NumFaces; dofneighbor++) {
    assert(reinterpret_cast<uintptr_t>(timeDofs[dofneighbor]) % Alignment == 0);
    assert(reinterpret_cast<uintptr_t>(integrationBuffer[dofneighbor]) % Alignment == 0);
  }
#endif

  /*
   * set/compute time integrated DOFs.
   */
  for (std::size_t dofneighbor = 0; dofneighbor < Cell::NumFaces; ++dofneighbor) {
    // collect information only in the case that neighboring element contributions are required
    if (faceTypes[dofneighbor] != FaceType::Outflow &&
        faceTypes[dofneighbor] != FaceType::DynamicRupture) {
      // check if the time integration is already done (-> copy pointer)
      if ((ltsSetup >> dofneighbor) % 2 == 0) {
        timeIntegrated[dofneighbor] = timeDofs[dofneighbor];
      }
      // integrate the DOFs in time via the derivatives and set pointer to local buffer
      else {
        time.computeIntegral(currentTime[dofneighbor + 1],
                             currentTime[0],
                             currentTime[0] + timeStepWidth,
                             timeDofs[dofneighbor],
                             integrationBuffer[dofneighbor]);

        timeIntegrated[dofneighbor] = integrationBuffer[dofneighbor];
      }
    }
  }
}

void TimeCommon::computeIntegrals(Time& time,
                                  unsigned short ltsSetup,
                                  const FaceType faceTypes[4],
                                  const double timeStepStart,
                                  const double timeStepWidth,
                                  real* const timeDofs[4],
                                  real integrationBuffer[4][seissol::kernels::Solver::BufferSize],
                                  real* timeIntegrated[4]) {
  double startTimes[5];
  startTimes[0] = timeStepStart;
  startTimes[1] = startTimes[2] = startTimes[3] = startTimes[4] = 0;

  // adjust start times for GTS on derivatives
  for (std::size_t face = 0; face < Cell::NumFaces; face++) {
    if (((ltsSetup >> (face + 4)) % 2) != 0) {
      startTimes[face + 1] = timeStepStart;
    }
  }

  // call the more general assembly
  computeIntegrals(time,
                   ltsSetup,
                   faceTypes,
                   startTimes,
                   timeStepWidth,
                   timeDofs,
                   integrationBuffer,
                   timeIntegrated);
}

void TimeCommon::computeBatchedIntegrals(Time& time,
                                         const double timeStepStart,
                                         const double timeStepWidth,
                                         ConditionalPointersToRealsTable& table,
                                         seissol::parallel::runtime::StreamRuntime& runtime) {
#ifdef ACL_DEVICE
  // Compute time integrated dofs using neighbors derivatives using the GTS relation,
  // i.e. the expansion point is around 'timeStepStart'
  ConditionalKey key(*KernelNames::NeighborFlux, *ComputationKind::WithGtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.computeBatchedIntegral(
        timeStepStart,
        timeStepStart,
        timeStepStart + timeStepWidth,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }

  // Compute time integrated dofs using neighbors derivatives using the LTS relation,
  // i.e. the expansion point is around '0'
  key = ConditionalKey(*KernelNames::NeighborFlux, *ComputationKind::WithLtsDerivatives);
  if (table.find(key) != table.end()) {
    auto& entry = table[key];
    time.computeBatchedIntegral(
        0.0,
        timeStepStart,
        timeStepStart + timeStepWidth,
        const_cast<const real**>((entry.get(inner_keys::Wp::Id::Derivatives))->getDeviceDataPtr()),
        (entry.get(inner_keys::Wp::Id::Idofs))->getDeviceDataPtr(),
        (entry.get(inner_keys::Wp::Id::Idofs))->getSize(),
        runtime);
  }
#else
  logError() << "No GPU implementation provided";
#endif
}

} // namespace seissol::kernels
