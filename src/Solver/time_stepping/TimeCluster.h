/******************************************************************************
** Copyright (c) 2015, Intel Corporation                                     **
** All rights reserved.                                                      **
**                                                                           **
** Redistribution and use in source and binary forms, with or without        **
** modification, are permitted provided that the following conditions        **
** are met:                                                                  **
** 1. Redistributions of source code must retain the above copyright         **
**    notice, this list of conditions and the following disclaimer.          **
** 2. Redistributions in binary form must reproduce the above copyright      **
**    notice, this list of conditions and the following disclaimer in the    **
**    documentation and/or other materials provided with the distribution.   **
** 3. Neither the name of the copyright holder nor the names of its          **
**    contributors may be used to endorse or promote products derived        **
**    from this software without specific prior written permission.          **
**                                                                           **
** THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS       **
** "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT         **
** LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR     **
** A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT      **
** HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,    **
** SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED  **
** TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR    **
** PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF    **
** LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING      **
** NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS        **
** SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.              **
******************************************************************************/
/* Alexander Heinecke (Intel Corp.)
******************************************************************************/
/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2013-2015, SeisSol Group
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 *
 * 2. Redistributions in binary form must reproduce the above copyright notice,
 *    this list of conditions and the following disclaimer in the documentation
 *    and/or other materials provided with the distribution.
 *
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from this
 *    software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * LTS cluster in SeisSol.
 **/

#ifndef TIMECLUSTER_H_
#define TIMECLUSTER_H_

#ifdef USE_MPI
#include <mpi.h>
#include <list>
#endif

#include <sstream>
#include "DynamicRupture/Misc.h"

#include <Initializer/typedefs.hpp>
#include <SourceTerm/typedefs.hpp>
#include <utils/logger.h>
#include <Initializer/LTS.h>
#include <Initializer/tree/LTSTree.hpp>

#include <Kernels/Time.h>
#include <Kernels/Local.h>
#include <Kernels/Neighbor.h>
#include <Kernels/DynamicRupture.h>
#include <Kernels/Plasticity.h>
#include <Kernels/TimeCommon.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include <Monitoring/LoopStatistics.h>
#include <Monitoring/ActorStateStatistics.h>

#include "AbstractTimeCluster.h"

#include <generated_code/kernel.h>

#include <algorithm>

#ifdef ACL_DEVICE
#include <device.h>
#include <Solver/Pipeline/DrPipeline.h>
#endif

namespace seissol {
  namespace time_stepping {
    class TimeCluster;
  }

  namespace kernels {
    class ReceiverCluster;
  }
}

/**
 * Time cluster, which represents a collection of elements having the same time step width.
 **/
class seissol::time_stepping::TimeCluster : public seissol::time_stepping::AbstractTimeCluster
{
private:
    // Last correction time of the neighboring cluster with higher dt
    double lastSubTime;

    void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
    void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
    void start() override {}
    void predict() override;
    void correct() override;
    bool usePlasticity;

    //! number of time steps
    unsigned long m_numberOfTimeSteps;

    /*
     * integrators
     */
    //! time kernel
    kernels::Time m_timeKernel;

    //! local kernel
    kernels::Local m_localKernel;

    //! neighbor kernel
    kernels::Neighbor m_neighborKernel;

    // //! cell average kernel
    // kernel::cellAve m_cellAverageKernel;

    kernels::DynamicRupture m_dynamicRuptureKernel;

    kernel::nonlEvaluateAndRotateQAtInterpolationPoints m_nonlinearInterpolation;

    kernel::nonlinearSurfaceIntegral m_nonlSurfIntPrototype;

    kernel::nonlinearVolumeIntegration m_krnlNonlVolPrototype;

  /*
   * global data
   */
     //! global data structures
    GlobalData *m_globalDataOnHost{nullptr};
    GlobalData *m_globalDataOnDevice{nullptr};
#ifdef ACL_DEVICE
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    dr::pipeline::DrPipeline drPipeline;
#endif

    /*
     * element data
     */
    seissol::initializers::Layer* m_clusterData;
    seissol::initializers::Layer* dynRupInteriorData;
    seissol::initializers::Layer* dynRupCopyData;
    seissol::initializers::LTS*         m_lts;
    seissol::initializers::DynamicRupture* m_dynRup;
    dr::friction_law::FrictionSolver* frictionSolver;
    dr::output::OutputManager* faultOutputManager;

    //! Mapping of cells to point sources
    sourceterm::CellToPointSourcesMapping const* m_cellToPointSources;

    //! Number of mapping of cells to point sources
    unsigned m_numberOfCellToPointSourcesMappings;

    //! Point sources
    sourceterm::PointSources const* m_pointSources;

    enum class ComputePart {
      Local = 0,
      Neighbor,
      DRNeighbor,
      DRFrictionLawInterior,
      DRFrictionLawCopy,
      PlasticityCheck,
      PlasticityYield,
      NUM_COMPUTE_PARTS
    };

    long long m_flops_nonZero[static_cast<int>(ComputePart::NUM_COMPUTE_PARTS)];
    long long m_flops_hardware[static_cast<int>(ComputePart::NUM_COMPUTE_PARTS)];

    //! Tv parameter for plasticity
    double m_tv;

    //! Relax time for plasticity
    double m_oneMinusIntegratingFactor;

    //! Stopwatch of TimeManager
    LoopStatistics* m_loopStatistics;
    ActorStateStatistics* actorStateStatistics;
    unsigned        m_regionComputeLocalIntegration;
    unsigned        m_regionComputeNeighboringIntegration;
    unsigned        m_regionComputeDynamicRupture;

    kernels::ReceiverCluster* m_receiverCluster;

    /**
     * Writes the receiver output if applicable (receivers present, receivers have to be written).
     **/
    void writeReceivers();

    /**
     * Computes the source terms if applicable.
     **/
    void computeSources();

    /**
     * Computes dynamic rupture.
     **/
    void computeDynamicRupture( seissol::initializers::Layer&  layerData );


    /**
     * Update all cell local material properties and the resulted matrices.
     *
     * This are:
     *  * time integration
     *  * volume integration
     *  * local boundary integration
     *
     * Remark: After this step the DOFs are only updated half with the boundary contribution
     *         of the neighborings cells missing.
     *
     * @param i_layerData number of cells.
     **/
    void updateMaterialLocal( seissol::initializers::Layer&  i_layerData);

    /**
     * Computes all cell local integration.
     *
     * This are:
     *  * time integration
     *  * volume integration
     *  * local boundary integration
     *
     * Remark: After this step the DOFs are only updated half with the boundary contribution
     *         of the neighborings cells missing.
     *
     * @param i_numberOfCells number of cells.
     * @param i_cellInformation cell local information.
     * @param i_cellData cell data.
     * @param io_buffers time integration buffers.
     * @param io_derivatives time derivatives.
     * @param io_dofs degrees of freedom.
     **/
    void computeLocalIntegration( seissol::initializers::Layer&  i_layerData, bool resetBuffers);

    /**
     * Computes the contribution of the neighboring cells to the boundary integral.
     *
     * Remark: After this step (in combination with the local integration) the DOFs are at the next time step.
     * TODO: This excludes dynamic rupture contribution.
     *
     * @param i_numberOfCells number of cells.
     * @param i_cellInformation cell local information.
     * @param i_cellData cell data.
     * @param i_faceNeighbors pointers to neighboring time buffers or derivatives.
     * @param io_dofs degrees of freedom.
     **/
    void computeNeighboringIntegration( seissol::initializers::Layer&  i_layerData, double subTimeStart );

    void computeLocalIntegrationFlops(seissol::initializers::Layer& layerData);
#ifndef ACL_DEVICE
    template<bool usePlasticity>
    std::pair<long, long> computeNeighboringIntegrationImplementation(seissol::initializers::Layer& i_layerData,
                                                                      double subTimeStart) {
      SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

      m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

      real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
      CellDRMapping (*drMapping)[4] = i_layerData.var(m_lts->drMapping);
      CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);
      PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
      real (*pstrain)[7 * NUMBER_OF_ALIGNED_BASIS_FUNCTIONS] = i_layerData.var(m_lts->pstrain);
      unsigned numberOTetsWithPlasticYielding = 0;

      kernels::NeighborData::Loader loader;
      loader.load(*m_lts, i_layerData);

      real *l_timeIntegrated[4];
      real *l_faceNeighbors_prefetch[4];

      real** derivatives = i_layerData.var(m_lts->derivatives);
      LocalIntegrationData* localIntegration = i_layerData.var(m_lts->localIntegration);
      CellMaterialData* materialData = i_layerData.var(m_lts->material);

      #if defined USE_DAMAGEDELASTIC
      double timePoints[CONVERGENCE_ORDER];
      double timeWeights[CONVERGENCE_ORDER];

      seissol::quadrature::GaussLegendre(timePoints, timeWeights, CONVERGENCE_ORDER);
      for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
        timePoints[point] = 0.5 * (timeStepSize() * timePoints[point] + timeStepSize());
        timeWeights[point] = 0.5 * timeStepSize() * timeWeights[point];
      }
      #endif

      // for (int i_out = 0; i_out<50; ++i_out){
      //   std::cout << faceNeighbors[i_out][0][20*6+0] << " ";
      // }
      // std::cout << faceNeighbors[50][0][20*9] << " "<< std::endl;

      // // Print nodes corresponding to exx and vx:
      // auto faceNeighborsView = init::Q::view::create(faceNeighbors[100][0]);
      // auto data1 = loader.entry(100);
      // // real* derivative = derivatives[100];
      // for (int i_out = 0; i_out<2; ++i_out){
      //   std::cout
      //             << faceNeighbors[100][0][20*i_out*6+0]/timeStepSize() << " "
      //             << faceNeighborsView(0,i_out*6)/timeStepSize() << " "
      //             << data1.dofs[20*i_out*6+0] << " "
      //             // << derivative[0]
      //             // << derivatives[l_cell][20*i_out+0]
      //             << "s ";
      // }
      // std::cout << faceNeighbors[100][0][20*6+0]/faceNeighbors[100][0][20*0+0]
      //           << " " << data1.dofs[20*6+0]/data1.dofs[20*0+0] << std::endl;

      #ifdef _OPENMP
      #pragma omp parallel for collapse(1) schedule(static) default(none) private(l_timeIntegrated, l_faceNeighbors_prefetch) shared(std::cout, materialData, localIntegration, cellInformation, loader, faceNeighbors, derivatives, pstrain, i_layerData, plasticity, drMapping, subTimeStart, timePoints, timeWeights) reduction(+:numberOTetsWithPlasticYielding)
      #endif
      for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
        auto data = loader.entry(l_cell);

        #if defined USE_DAMAGEDELASTIC
          // Nonlinear surface flux integration -- for regular/periodic flux
          /// TODO: Check if it works for periodic BCs.
          // Here, plus side is actually minus (or local solution side),
          // minus side is neighbor solution side.
        for (unsigned int side = 0; side < 4; side++ ){
          if (cellInformation[l_cell].faceTypes[side] == FaceType::regular
          || cellInformation[l_cell].faceTypes[side] == FaceType::periodic){
            // Compute local integrals with derivatives and Rusanov flux
            /// S1: compute the space-time interpolated Q on both side of 4 faces
            /// S2: at the same time rotate the field to face-aligned coord.
            alignas(PAGESIZE_STACK) real QInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};
            alignas(PAGESIZE_STACK) real QInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};

            for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
              // alignas(PAGESIZE_STACK)
              real degreesOfFreedomPlus[tensor::Q::size()];
              // alignas(PAGESIZE_STACK)
              real degreesOfFreedomMinus[tensor::Q::size()];

              for (unsigned i_f = 0; i_f < tensor::Q::size(); i_f++){
                degreesOfFreedomPlus[i_f] = static_cast<real>(0.0);
                degreesOfFreedomMinus[i_f] = static_cast<real>(0.0);
              }

              // !!! Make sure every time after entering this function, the last input should be reinitialized to zero
              m_timeKernel.computeTaylorExpansion(timePoints[timeInterval], 0.0, derivatives[l_cell], degreesOfFreedomPlus);
              m_timeKernel.computeTaylorExpansion(timePoints[timeInterval], 0.0, faceNeighbors[l_cell][side], degreesOfFreedomMinus);

              /// Prototype is necessary for openmp
              kernel::nonlEvaluateAndRotateQAtInterpolationPoints m_nonLinInter
                = m_nonlinearInterpolation;

              m_nonLinInter.QInterpolated = &QInterpolatedPlus[timeInterval][0];
              m_nonLinInter.Q = degreesOfFreedomPlus;
              // m_nonlinearInterpolation.TinvT = localIntegration[l_cell].TinvT[side];
              // m_nonlinearInterpolation._prefetch.QInterpolated = plusPrefetch;
              m_nonLinInter.execute(side, 0);

              m_nonLinInter.QInterpolated = &QInterpolatedMinus[timeInterval][0];
              m_nonLinInter.Q = degreesOfFreedomMinus;
              // m_nonlinearInterpolation.TinvT = localIntegration[l_cell].TinvT[side];
              // m_nonlinearInterpolation._prefetch.QInterpolated = minusPrefetch;
              m_nonLinInter.execute(cellInformation[l_cell].faceRelations[side][0]
                            , cellInformation[l_cell].faceRelations[side][1]+1);
            }

            /// S3: Construct matrices to store Rusanov flux on surface quadrature nodes.
            //// Reshape the interpolated results
            using QInterpolatedShapeT = const real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];
            // std::cout << seissol::dr::misc::numQuantities << std::endl;
            auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedPlus));
            auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedMinus));

            //// The arrays to store time integrated flux
            alignas(PAGESIZE_STACK) real rusanovFluxPlus[tensor::QInterpolated::size()] = {0.0};
            // alignas(PAGESIZE_STACK) real rusanovFluxMinus[tensor::QInterpolated::size()] = {0.0};

            for (unsigned i_f = 0; i_f < tensor::QInterpolated::size(); i_f++){
              rusanovFluxPlus[i_f] = static_cast<real>(0.0);
            }

            using rusanovFluxShape = real(*)[seissol::dr::misc::numPaddedPoints];
            auto* rusanovFluxP = reinterpret_cast<rusanovFluxShape>(rusanovFluxPlus);
            // auto* rusanovFluxM = reinterpret_cast<rusanovFluxShape>(rusanovFluxMinus);

            /// Checked that, after reshaping, it still uses the same memory address
            /// S4: Integration in time the Rusanov flux on surface quadrature nodes.
            using namespace seissol::dr::misc::quantity_indices;
            unsigned DAM = 9;

            real epsInitxx = -1e-2; // eps_xx0
            real epsInityy = -0e-1; // eps_yy0
            real epsInitzz = -0e-1; // eps_zz0
            real lambda0 = materialData[l_cell].local.lambda0;
            real mu0 = materialData[l_cell].local.mu0;
            real rho0 = materialData[l_cell].local.rho;

            real lambda_max = 1.0*std::sqrt( (lambda0+2*mu0)/rho0 ) ;
            real sxxP, syyP, szzP, sxyP, syzP, szxP
            ,sxxM, syyM, szzM, sxyM, syzM, szxM;

            /// In this time loop, use "qIPlus" and "qIMinus" to interpolate "rusanovFluxP"
            for (unsigned o = 0; o < CONVERGENCE_ORDER; ++o) {
              auto weight = timeWeights[o];

              for (unsigned i = 0; i < seissol::dr::misc::numPaddedPoints;
                  i ++) {
                lambda_max = std::max(
                  std::sqrt( (1- qIPlus[o][DAM][i]) * (lambda0+2*mu0)/rho0 ),
                  std::sqrt( (1-qIMinus[o][DAM][i]) * (lambda0+2*mu0)/rho0 )
                );

                real EspIp = (qIPlus[o][XX][i]+epsInitxx) + (qIPlus[o][YY][i]+epsInityy) + (qIPlus[o][ZZ][i]+epsInitzz);
                real EspIIp = (qIPlus[o][XX][i]+epsInitxx)*(qIPlus[o][XX][i]+epsInitxx)
                  + (qIPlus[o][YY][i]+epsInityy)*(qIPlus[o][YY][i]+epsInityy)
                  + (qIPlus[o][ZZ][i]+epsInitzz)*(qIPlus[o][ZZ][i]+epsInitzz)
                  + 2*qIPlus[o][XY][i]*qIPlus[o][XY][i]
                  + 2*qIPlus[o][YZ][i]*qIPlus[o][YZ][i]
                  + 2*qIPlus[o][XZ][i]*qIPlus[o][XZ][i];
                real alphap = qIPlus[o][DAM][i];
                real xip;
                if (EspIIp > 1e-30){
                  xip = EspIp / std::sqrt(EspIIp);
                } else{
                  xip = 0.0;
                }

                sxxP = (lambda0*EspIp - alphap*materialData[l_cell].local.gammaR*std::sqrt(EspIIp))
                      + (2*(mu0 - alphap*materialData[l_cell].local.gammaR*materialData[l_cell].local.xi0)
                          - alphap*materialData[l_cell].local.gammaR*xip)
                        *(qIPlus[o][XX][i]+epsInitxx);

                syyP = (lambda0*EspIp - alphap*materialData[l_cell].local.gammaR*std::sqrt(EspIIp))
                      + (2*(mu0 - alphap*materialData[l_cell].local.gammaR*materialData[l_cell].local.xi0)
                          - alphap*materialData[l_cell].local.gammaR*xip)
                        *(qIPlus[o][YY][i]+epsInityy);

                szzP = (lambda0*EspIp - alphap*materialData[l_cell].local.gammaR*std::sqrt(EspIIp))
                      + (2*(mu0 - alphap*materialData[l_cell].local.gammaR*materialData[l_cell].local.xi0)
                          - alphap*materialData[l_cell].local.gammaR*xip)
                        *(qIPlus[o][ZZ][i]+epsInitzz);

                sxyP = 0
                      + (2*(mu0 - alphap*materialData[l_cell].local.gammaR*materialData[l_cell].local.xi0)
                          - alphap*materialData[l_cell].local.gammaR*xip)
                        *qIPlus[o][XY][i];

                syzP = 0
                      + (2*(mu0 - alphap*materialData[l_cell].local.gammaR*materialData[l_cell].local.xi0)
                          - alphap*materialData[l_cell].local.gammaR*xip)
                        *qIPlus[o][YZ][i];

                szxP = 0
                      + (2*(mu0 - alphap*materialData[l_cell].local.gammaR*materialData[l_cell].local.xi0)
                          - alphap*materialData[l_cell].local.gammaR*xip)
                        *qIPlus[o][XZ][i];

                real EspIm = (qIMinus[o][XX][i]+epsInitxx) + (qIMinus[o][YY][i]+epsInityy) + (qIMinus[o][ZZ][i]+epsInitzz);
                real EspIIm = (qIMinus[o][XX][i]+epsInitxx)*(qIMinus[o][XX][i]+epsInitxx)
                  + (qIMinus[o][YY][i]+epsInityy)*(qIMinus[o][YY][i]+epsInityy)
                  + (qIMinus[o][ZZ][i]+epsInitzz)*(qIMinus[o][ZZ][i]+epsInitzz)
                  + 2*qIMinus[o][XY][i]*qIMinus[o][XY][i]
                  + 2*qIMinus[o][YZ][i]*qIMinus[o][YZ][i]
                  + 2*qIMinus[o][XZ][i]*qIMinus[o][XZ][i];
                real alpham = qIMinus[o][DAM][i];
                real xim;
                if (EspIIm > 1e-30){
                  xim = EspIm / std::sqrt(EspIIm);
                } else{
                  xim = 0.0;
                }

                sxxM = (lambda0*EspIm - alpham*materialData[l_cell].neighbor[side].gammaR*std::sqrt(EspIIm))
                      + (2*(mu0 - alpham*materialData[l_cell].neighbor[side].gammaR*materialData[l_cell].neighbor[side].xi0)
                          - alpham*materialData[l_cell].neighbor[side].gammaR*xim)
                        *(qIMinus[o][XX][i]+epsInitxx);

                syyM = (lambda0*EspIm - alpham*materialData[l_cell].neighbor[side].gammaR*std::sqrt(EspIIm))
                      + (2*(mu0 - alpham*materialData[l_cell].neighbor[side].gammaR*materialData[l_cell].neighbor[side].xi0)
                          - alpham*materialData[l_cell].neighbor[side].gammaR*xim)
                        *(qIMinus[o][YY][i]+epsInityy);

                szzM = (lambda0*EspIm - alpham*materialData[l_cell].neighbor[side].gammaR*std::sqrt(EspIIm))
                      + (2*(mu0 - alpham*materialData[l_cell].neighbor[side].gammaR*materialData[l_cell].neighbor[side].xi0)
                          - alpham*materialData[l_cell].neighbor[side].gammaR*xim)
                        *(qIMinus[o][ZZ][i]+epsInitzz);

                sxyM = 0
                      + (2*(mu0 - alpham*materialData[l_cell].neighbor[side].gammaR*materialData[l_cell].neighbor[side].xi0)
                          - alpham*materialData[l_cell].neighbor[side].gammaR*xim)
                        *qIMinus[o][XY][i];

                syzM = 0
                      + (2*(mu0 - alpham*materialData[l_cell].neighbor[side].gammaR*materialData[l_cell].neighbor[side].xi0)
                          - alpham*materialData[l_cell].neighbor[side].gammaR*xim)
                        *qIMinus[o][YZ][i];

                szxM = 0
                      + (2*(mu0 - alpham*materialData[l_cell].neighbor[side].gammaR*materialData[l_cell].neighbor[side].xi0)
                          - alpham*materialData[l_cell].neighbor[side].gammaR*xim)
                        *qIMinus[o][XZ][i];

                rusanovFluxP[XX][i] += weight * (
                  (
                    0.5*(-qIPlus[o][U][i]) + 0.5*(-qIMinus[o][U][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][XX][i]+epsInitxx) - 0.5*lambda_max*(qIMinus[o][XX][i]+epsInitxx)
                );

                rusanovFluxP[YY][i] += weight * (
                  (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-qIPlus[o][V][i]) + 0.5*(-qIMinus[o][V][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][YY][i]+epsInityy) - 0.5*lambda_max*(qIMinus[o][YY][i]+epsInityy)
                );

                rusanovFluxP[ZZ][i] += weight * (
                  (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-qIPlus[o][W][i]) + 0.5*(-qIMinus[o][W][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][ZZ][i]+epsInitzz) - 0.5*lambda_max*(qIMinus[o][ZZ][i]+epsInitzz)
                );

                rusanovFluxP[XY][i] += weight * (
                  (
                    0.5*(-0.5*qIPlus[o][V][i]) + 0.5*(-0.5*qIMinus[o][V][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-0.5*qIPlus[o][U][i]) + 0.5*(-0.5*qIMinus[o][U][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][XY][i]) - 0.5*lambda_max*(qIMinus[o][XY][i])
                );

                rusanovFluxP[YZ][i] += weight * (
                  (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-0.5*qIPlus[o][W][i]) + 0.5*(-0.5*qIMinus[o][W][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-0.5*qIPlus[o][V][i]) + 0.5*(-0.5*qIMinus[o][V][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][YZ][i]) - 0.5*lambda_max*(qIMinus[o][YZ][i])
                );

                rusanovFluxP[XZ][i] += weight * (
                  (
                    0.5*(-0.5*qIPlus[o][W][i]) + 0.5*(-0.5*qIMinus[o][W][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-0.5*qIPlus[o][U][i]) + 0.5*(-0.5*qIMinus[o][U][i])
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][XZ][i]) - 0.5*lambda_max*(qIMinus[o][XZ][i])
                );

                rusanovFluxP[U][i] += weight * (
                  (
                    0.5*(-sxxP/rho0) + 0.5*(-sxxM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-sxyP/rho0) + 0.5*(-sxyM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-szxP/rho0) + 0.5*(-szxM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][U][i]) - 0.5*lambda_max*(qIMinus[o][U][i])
                );

                rusanovFluxP[V][i] += weight * (
                  (
                    0.5*(-sxyP/rho0) + 0.5*(-sxyM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-syyP/rho0) + 0.5*(-syyM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-syzP/rho0) + 0.5*(-syzM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][V][i]) - 0.5*lambda_max*(qIMinus[o][V][i])
                );

                rusanovFluxP[W][i] += weight * (
                  (
                    0.5*(-szxP/rho0) + 0.5*(-szxM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-syzP/rho0) + 0.5*(-syzM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-szzP/rho0) + 0.5*(-szzM/rho0)
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][W][i]) - 0.5*lambda_max*(qIMinus[o][W][i])
                );

                rusanovFluxP[DAM][i] += weight * (
                  (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][0]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][1]
                  + (
                    0.5*(-0) + 0.5*(-0)
                  ) * localIntegration[l_cell].surfaceNormal[side][2]
                  + 0.5*lambda_max*(qIPlus[o][DAM][i]) - 0.5*lambda_max*(qIMinus[o][DAM][i])
                );
              }
            } // time integration loop
            // Tested that after this step, rusanovFluxPlus is also changed

            /// S5: Integrate in space using quadrature.
            kernel::nonlinearSurfaceIntegral m_surfIntegral = m_nonlSurfIntPrototype;
            m_surfIntegral.Q = data.dofs;
            m_surfIntegral.Flux = rusanovFluxPlus;
            // m_surfIntegral.TT = localIntegration[l_cell].TT[side];
            m_surfIntegral.fluxScale = localIntegration[l_cell].fluxScales[side];
            // m_surfIntegral._prefetch.I = &QInterpolatedPlus[0][0];
            m_surfIntegral.execute(side, 0);
          } // if (faceTypes)
        } // for (side)

        #else

        // for (int i_out = 0; i_out<50; ++i_out){
        //   std::cout << faceNeighbors[l_cell][0][20*i_out+0] << " ";
        // }
        // std::cout << faceNeighbors[l_cell][0][50] << " "<< std::endl;
        seissol::kernels::TimeCommon::computeIntegrals(m_timeKernel,
                                                       data.cellInformation.ltsSetup,
                                                       data.cellInformation.faceTypes,
                                                       subTimeStart,
                                                       timeStepSize(),
                                                       faceNeighbors[l_cell],
#ifdef _OPENMP
                                                       *reinterpret_cast<real (*)[4][tensor::I::size()]>(&(m_globalDataOnHost->integrationBufferLTS[omp_get_thread_num()*4*tensor::I::size()])),
#else
            *reinterpret_cast<real (*)[4][tensor::I::size()]>(m_globalData->integrationBufferLTS),
#endif
                                                       l_timeIntegrated);

#ifdef ENABLE_MATRIX_PREFETCH
        l_faceNeighbors_prefetch[0] = (cellInformation[l_cell].faceTypes[1] != FaceType::dynamicRupture) ?
                                      faceNeighbors[l_cell][1] :
                                      drMapping[l_cell][1].godunov;
        l_faceNeighbors_prefetch[1] = (cellInformation[l_cell].faceTypes[2] != FaceType::dynamicRupture) ?
                                      faceNeighbors[l_cell][2] :
                                      drMapping[l_cell][2].godunov;
        l_faceNeighbors_prefetch[2] = (cellInformation[l_cell].faceTypes[3] != FaceType::dynamicRupture) ?
                                      faceNeighbors[l_cell][3] :
                                      drMapping[l_cell][3].godunov;

        // fourth face's prefetches
        if (l_cell < (i_layerData.getNumberOfCells()-1) ) {
          l_faceNeighbors_prefetch[3] = (cellInformation[l_cell+1].faceTypes[0] != FaceType::dynamicRupture) ?
                                        faceNeighbors[l_cell+1][0] :
                                        drMapping[l_cell+1][0].godunov;
        } else {
          l_faceNeighbors_prefetch[3] = faceNeighbors[l_cell][3];
        }
#endif

        m_neighborKernel.computeNeighborsIntegral( data,
                                                   drMapping[l_cell],
#ifdef ENABLE_MATRIX_PREFETCH
                                                   l_timeIntegrated, l_faceNeighbors_prefetch
#else
            l_timeIntegrated
#endif
        );
        #endif // end of DAMAGE conditions

        if constexpr (usePlasticity) {
          updateRelaxTime();
          numberOTetsWithPlasticYielding += seissol::kernels::Plasticity::computePlasticity( m_oneMinusIntegratingFactor,
                                                                                             timeStepSize(),
                                                                                             m_tv,
                                                                                             m_globalDataOnHost,
                                                                                             &plasticity[l_cell],
                                                                                             data.dofs,
                                                                                             pstrain[l_cell] );
        }
#ifdef INTEGRATE_QUANTITIES
        seissol::SeisSol::main.postProcessor().integrateQuantities( m_timeStepWidth,
                                                              i_layerData,
                                                              l_cell,
                                                              dofs[l_cell] );
#endif // INTEGRATE_QUANTITIES
      }

      const long long nonZeroFlopsPlasticity =
          i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)] +
          numberOTetsWithPlasticYielding * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)];
      const long long hardwareFlopsPlasticity =
          i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)] +
          numberOTetsWithPlasticYielding * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)];

      m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_globalClusterId);

      return {nonZeroFlopsPlasticity, hardwareFlopsPlasticity};
    }
#endif // ACL_DEVICE

    void computeLocalIntegrationFlops(unsigned numberOfCells,
                                      CellLocalInformation const* cellInformation,
                                      long long& nonZeroFlops,
                                      long long& hardwareFlops);

    void computeNeighborIntegrationFlops(seissol::initializers::Layer &layerData);

    void computeDynamicRuptureFlops(seissol::initializers::Layer &layerData,
                                    long long& nonZeroFlops,
                                    long long& hardwareFlops);

    void computeFlops();

    //! Update relax time for plasticity
    void updateRelaxTime() {
      m_oneMinusIntegratingFactor = (m_tv > 0.0) ? 1.0 - exp(-timeStepSize() / m_tv) : 1.0;
    }

  const LayerType layerType;
  //! time of the next receiver output
  double m_receiverTime;

  //! print status every 100th timestep
  bool printProgress;
  //! cluster id on this rank
  const unsigned int m_clusterId;

  //! global cluster cluster id
  const unsigned int m_globalClusterId;

  DynamicRuptureScheduler* dynamicRuptureScheduler;

  void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;

public:
  ActResult act() override;

  /**
   * Constructs a new LTS cluster.
   *
   * @param i_clusterId id of this cluster with respect to the current rank.
   * @param i_globalClusterId global id of this cluster.
   * @param usePlasticity true if using plasticity
   **/
  TimeCluster(unsigned int i_clusterId, unsigned int i_globalClusterId, bool usePlasticity,
              LayerType layerType, double maxTimeStepSize,
              long timeStepRate, bool printProgress,
              DynamicRuptureScheduler* dynamicRuptureScheduler, CompoundGlobalData i_globalData,
              seissol::initializers::Layer *i_clusterData, seissol::initializers::Layer* dynRupInteriorData,
              seissol::initializers::Layer* dynRupCopyData, seissol::initializers::LTS* i_lts,
              seissol::initializers::DynamicRupture* i_dynRup,
              seissol::dr::friction_law::FrictionSolver* i_FrictionSolver,
              dr::output::OutputManager* i_faultOutputManager, LoopStatistics* i_loopStatistics,
              ActorStateStatistics* actorStateStatistics);

  /**
   * Destructor of a LTS cluster.
   * TODO: Currently prints only statistics in debug mode.
   **/
  ~TimeCluster() override;

  /**
   * Sets the pointer to the cluster's point sources
   *
   * @param i_cellToPointSources Contains mappings of 1 cell offset to m point sources
   * @param i_numberOfCellToPointSourcesMappings Size of i_cellToPointSources
   * @param i_pointSources pointer to all point sources used on this cluster
   */
  void setPointSources( sourceterm::CellToPointSourcesMapping const* i_cellToPointSources,
                        unsigned i_numberOfCellToPointSourcesMappings,
                        sourceterm::PointSources const* i_pointSources );

  void setReceiverCluster( kernels::ReceiverCluster* receiverCluster) {
    m_receiverCluster = receiverCluster;
  }

  void setFaultOutputManager(dr::output::OutputManager* outputManager) {
    faultOutputManager = outputManager;
  }

  /**
   * Set Tv constant for plasticity.
   */
  void setTv(double tv) {
    m_tv = tv;
    updateRelaxTime();
  }

  /**
   * Initialize the derivatives with data.dofs.
   *
   * Remark: This is called before entering time loop.
   *
   **/
  void updateDerivatives();


  void reset() override;

  [[nodiscard]] unsigned int getClusterId() const;
  [[nodiscard]] unsigned int getGlobalClusterId() const;
  [[nodiscard]] LayerType getLayerType() const;
  void setReceiverTime(double receiverTime);
};

#endif
