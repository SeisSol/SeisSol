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

#ifndef TIMECLUSTER_H_ // okay
#define TIMECLUSTER_H_

#ifdef USE_MPI // okay
#include <mpi.h>
#include <list>
#endif // USE_MPI

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
#include <Kernels/PointSourceCluster.h>
#include <Kernels/TimeCommon.h>
#include <Solver/FreeSurfaceIntegrator.h>
#include <Monitoring/LoopStatistics.h>
#include <Monitoring/ActorStateStatistics.h>
#include "Initializer/DynamicRupture.h"
#include "DynamicRupture/FrictionLaws/FrictionSolver.h"
#include "DynamicRupture/Output/OutputManager.hpp"

#include "AbstractTimeCluster.h"

#ifdef ACL_DEVICE // okay
#include <device.h>
#include <Solver/Pipeline/DrPipeline.h>
#endif // ACL_DEVICE

namespace seissol {
  class SeisSol;
  namespace time_stepping {
    class TimeCluster;
  } // namespace time_stepping

  namespace kernels {
    class ReceiverCluster;
  } // namespace kernels
} // namespace seissol

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

    seissol::SeisSol& seissolInstance;
    /*
     * integrators
     */
    //! time kernel
    kernels::Time m_timeKernel;

    //! local kernel
    kernels::Local m_localKernel;

    //! neighbor kernel
    kernels::Neighbor m_neighborKernel;
    
    kernels::DynamicRupture m_dynamicRuptureKernel;

#ifdef USE_DAMAGEDELASTIC
    kernel::nonlEvaluateAndRotateQAtInterpolationPoints m_nonlinearInterpolation;
    kernel::nonlinearSurfaceIntegral m_nonlSurfIntPrototype;
    kernel::nonlinearVolumeIntegration m_krnlNonlVolPrototype;
#endif

  /*
   * global data
   */
     //! global data structures
    GlobalData *m_globalDataOnHost{nullptr};
    GlobalData *m_globalDataOnDevice{nullptr};
#ifdef ACL_DEVICE // okay
    device::DeviceInstance& device = device::DeviceInstance::getInstance();
    dr::pipeline::DrPipeline drPipeline;
#endif // ACL_DEVICE

    /*
     * element data
     */     
    seissol::initializer::Layer* m_clusterData;
    seissol::initializer::Layer* dynRupInteriorData;
    seissol::initializer::Layer* dynRupCopyData;
    seissol::initializer::LTS*         m_lts;
    seissol::initializer::DynamicRupture* m_dynRup;
    dr::friction_law::FrictionSolver* frictionSolver;
    dr::output::OutputManager* faultOutputManager;

    std::unique_ptr<kernels::PointSourceCluster> m_sourceCluster;

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
    unsigned        m_regionComputePointSources;

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
    void computeDynamicRupture( seissol::initializer::Layer&  layerData );

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
    void updateMaterialLocal( seissol::initializer::Layer&  i_layerData);

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
    void computeLocalIntegration( seissol::initializer::Layer&  i_layerData, bool resetBuffers);

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
    void computeNeighboringIntegration( seissol::initializer::Layer&  i_layerData, double subTimeStart );

    void computeLocalIntegrationFlops(seissol::initializer::Layer& layerData);
#ifndef ACL_DEVICE
#ifdef USE_DAMAGEDELASTIC
    template<bool usePlasticity>
    std::pair<long, long> computeNeighboringIntegrationImplementation(seissol::initializer::Layer& i_layerData,
                                                                      double subTimeStart) {
      if (i_layerData.getNumberOfCells() == 0) return {0,0};
      SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

      m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

      real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
      CellDRMapping (*drMapping)[4] = i_layerData.var(m_lts->drMapping);
      CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);
      PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
      auto* pstrain = i_layerData.var(m_lts->pstrain);
      unsigned numberOTetsWithPlasticYielding = 0;

      kernels::NeighborData::Loader loader;
      loader.load(*m_lts, i_layerData);

      real *l_timeIntegrated[4];
      real *l_faceNeighbors_prefetch[4];

      real** derivatives = i_layerData.var(m_lts->derivatives);
      LocalIntegrationData* localIntegration = i_layerData.var(m_lts->localIntegration);
      CellMaterialData* materialData = i_layerData.var(m_lts->material);

      double timePoints[CONVERGENCE_ORDER];
      double timeWeights[CONVERGENCE_ORDER];

      seissol::quadrature::GaussLegendre(timePoints, timeWeights, CONVERGENCE_ORDER);
      for (unsigned point = 0; point < CONVERGENCE_ORDER; ++point) {
        timePoints[point] = 0.5 * (timeStepSize() * timePoints[point] + timeStepSize());
        timeWeights[point] = 0.5 * timeStepSize() * timeWeights[point];
      }

      
      #ifdef _OPENMP // okay
      #pragma omp parallel for collapse(1) schedule(static) default(none) private(l_timeIntegrated, l_faceNeighbors_prefetch) shared( materialData, localIntegration, cellInformation, loader, faceNeighbors, derivatives, pstrain, i_layerData, plasticity, drMapping, subTimeStart, timePoints, timeWeights) reduction(+:numberOTetsWithPlasticYielding)
      #endif // _OPENMP
      for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
        auto data = loader.entry(l_cell);

        // Nonlinear surface flux integration -- for regular/periodic flux
        // TODO: Check if it works for periodic BCs.
        // Here, plus side is actually minus (or local solution side),
        // minus side is neighbor solution side.
        for (unsigned int side = 0; side < 4; side++ ){
          if (cellInformation[l_cell].faceTypes[side] == FaceType::regular
          || cellInformation[l_cell].faceTypes[side] == FaceType::periodic){
            // Compute local integrals with derivatives and Rusanov flux
            /// S1: compute the space-time interpolated Q on both side of 4 faces
            /// S2: at the same time rotate the field to face-aligned coord.
            alignas(ALIGNMENT) real QInterpolatedPlus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};
            alignas(ALIGNMENT) real QInterpolatedMinus[CONVERGENCE_ORDER][seissol::tensor::QInterpolated::size()] = {{0.0}};

            for (unsigned timeInterval = 0; timeInterval < CONVERGENCE_ORDER; ++timeInterval) {
              alignas(ALIGNMENT) real degreesOfFreedomPlus[tensor::Q::size()];
              alignas(ALIGNMENT) real degreesOfFreedomMinus[tensor::Q::size()];

              for (unsigned i_f = 0; i_f < tensor::Q::size(); i_f++){
                degreesOfFreedomPlus[i_f] = static_cast<real>(0.0);
                degreesOfFreedomMinus[i_f] = static_cast<real>(0.0);
              }

              // !!! Make sure every time after entering this function, the last input should be reinitialized to zero
              m_timeKernel.computeTaylorExpansion(timePoints[timeInterval], 0.0, derivatives[l_cell], degreesOfFreedomPlus);
              m_timeKernel.computeTaylorExpansion(timePoints[timeInterval], 0.0, faceNeighbors[l_cell][side], degreesOfFreedomMinus);

              // Prototype is necessary for openmp
              kernel::nonlEvaluateAndRotateQAtInterpolationPoints m_nonLinInter
                = m_nonlinearInterpolation;

              m_nonLinInter.QInterpolated = &QInterpolatedPlus[timeInterval][0];
              m_nonLinInter.Q = degreesOfFreedomPlus;
              m_nonLinInter.execute(side, 0);

              m_nonLinInter.QInterpolated = &QInterpolatedMinus[timeInterval][0];
              m_nonLinInter.Q = degreesOfFreedomMinus;
              m_nonLinInter.execute(cellInformation[l_cell].faceRelations[side][0]
                            , cellInformation[l_cell].faceRelations[side][1]+1);
            }

            // S3: Construct matrices to store Rusanov flux on surface quadrature nodes.
            // Reshape the interpolated results
            using QInterpolatedShapeT = const real(*)[seissol::dr::misc::numQuantities][seissol::dr::misc::numPaddedPoints];

            auto* qIPlus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedPlus));
            auto* qIMinus = (reinterpret_cast<QInterpolatedShapeT>(QInterpolatedMinus));

            // The arrays to store time integrated flux
            alignas(ALIGNMENT) real rusanovFluxPlus[tensor::QInterpolated::size()] = {0.0};
            for (unsigned i_f = 0; i_f < tensor::QInterpolated::size(); i_f++){
              rusanovFluxPlus[i_f] = static_cast<real>(0.0);
            }
            using rusanovFluxShape = real(*)[seissol::dr::misc::numPaddedPoints];
            auto* rusanovFluxP = reinterpret_cast<rusanovFluxShape>(rusanovFluxPlus);
            m_timeKernel.computeNonLinearRusanovFlux(materialData, l_cell, side, timeWeights, *qIPlus[0], 
                *qIMinus[0], *rusanovFluxP, localIntegration);
            /// S5: Integrate in space using quadrature.
            kernel::nonlinearSurfaceIntegral m_surfIntegral = m_nonlSurfIntPrototype;
            m_surfIntegral.Q = data.dofs;
            m_surfIntegral.Flux = rusanovFluxPlus;
            m_surfIntegral.fluxScale = localIntegration[l_cell].fluxScales[side];
            m_surfIntegral.execute(side, 0);
          }
          else if (cellInformation[l_cell].faceTypes[side] == FaceType::dynamicRupture) {
            // No neighboring cell contribution, interior bc.
            assert(reinterpret_cast<uintptr_t>(drMapping[l_cell][side].godunov) % ALIGNMENT == 0);

            kernel::nonlinearSurfaceIntegral m_drIntegral = m_nonlSurfIntPrototype;
            m_drIntegral.Q = data.dofs;
            m_drIntegral.Flux = drMapping[l_cell][side].godunov;
            m_drIntegral.fluxScale = localIntegration[l_cell].fluxScales[side];
            m_drIntegral.execute(side, drMapping[l_cell][side].faceRelation);

          } // if (faceTypes)
        } // for (side)
      }

      const long long nonZeroFlopsPlasticity =
          i_layerData.getNumberOfCells() * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityCheck)] +
          numberOTetsWithPlasticYielding * m_flops_nonZero[static_cast<int>(ComputePart::PlasticityYield)];
      const long long hardwareFlopsPlasticity =
          i_layerData.getNumberOfCells() * m_flops_hardware[static_cast<int>(ComputePart::PlasticityCheck)] +
          numberOTetsWithPlasticYielding * m_flops_hardware[static_cast<int>(ComputePart::PlasticityYield)];

      m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_profilingId);

      return {nonZeroFlopsPlasticity, hardwareFlopsPlasticity};
    }
#else // USE_DAMAGEDELASTIC
    template<bool usePlasticity>
    std::pair<long, long> computeNeighboringIntegrationImplementation(seissol::initializer::Layer& i_layerData,
                                                                      double subTimeStart) {
      if (i_layerData.getNumberOfCells() == 0) return {0,0};
      SCOREP_USER_REGION( "computeNeighboringIntegration", SCOREP_USER_REGION_TYPE_FUNCTION )

      m_loopStatistics->begin(m_regionComputeNeighboringIntegration);

      real* (*faceNeighbors)[4] = i_layerData.var(m_lts->faceNeighbors);
      CellDRMapping (*drMapping)[4] = i_layerData.var(m_lts->drMapping);
      CellLocalInformation* cellInformation = i_layerData.var(m_lts->cellInformation);
      PlasticityData* plasticity = i_layerData.var(m_lts->plasticity);
      auto* pstrain = i_layerData.var(m_lts->pstrain);
      unsigned numberOTetsWithPlasticYielding = 0;

      kernels::NeighborData::Loader loader;
      loader.load(*m_lts, i_layerData);

      real *l_timeIntegrated[4];
      real *l_faceNeighbors_prefetch[4];

      real** derivatives = i_layerData.var(m_lts->derivatives);
      LocalIntegrationData* localIntegration = i_layerData.var(m_lts->localIntegration);
      CellMaterialData* materialData = i_layerData.var(m_lts->material);
      
      #ifdef _OPENMP // okay
      #pragma omp parallel for schedule(static) default(none) private(l_timeIntegrated, l_faceNeighbors_prefetch) shared(cellInformation, loader, faceNeighbors, pstrain, i_layerData, plasticity, drMapping, subTimeStart) reduction(+:numberOTetsWithPlasticYielding)
      #endif // _OPENMP
      for( unsigned int l_cell = 0; l_cell < i_layerData.getNumberOfCells(); l_cell++ ) {
        auto data = loader.entry(l_cell);

        seissol::kernels::TimeCommon::computeIntegrals(m_timeKernel,
                                                       data.cellInformation.ltsSetup,
                                                       data.cellInformation.faceTypes,
                                                       subTimeStart,
                                                       timeStepSize(),
                                                       faceNeighbors[l_cell],
#ifdef _OPENMP //okay
                                                       *reinterpret_cast<real (*)[4][tensor::I::size()]>(&(m_globalDataOnHost->integrationBufferLTS[omp_get_thread_num()*4*tensor::I::size()])),
#else // _OPENMP
            *reinterpret_cast<real (*)[4][tensor::I::size()]>(m_globalData->integrationBufferLTS),
#endif // _OPENMP
                                                       l_timeIntegrated);

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

        m_neighborKernel.computeNeighborsIntegral( data,
                                                   drMapping[l_cell],
                                                   l_timeIntegrated, l_faceNeighbors_prefetch
        );

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
#ifdef INTEGRATE_QUANTITIES // okay
        seissolInstance.postProcessor().integrateQuantities( m_timeStepWidth,
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

      m_loopStatistics->end(m_regionComputeNeighboringIntegration, i_layerData.getNumberOfCells(), m_profilingId);

      return {nonZeroFlopsPlasticity, hardwareFlopsPlasticity};
    }
#endif // USE_DAMAGEDELASTIC
#endif // ACL_DEVICE

    void computeLocalIntegrationFlops(unsigned numberOfCells,
                                      CellLocalInformation const* cellInformation,
                                      long long& nonZeroFlops,
                                      long long& hardwareFlops);

    void computeNeighborIntegrationFlops(seissol::initializer::Layer &layerData);

    void computeDynamicRuptureFlops(seissol::initializer::Layer &layerData,
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

  //! id used to identify this cluster (including layer type) when profiling
  const unsigned int m_profilingId;

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
  TimeCluster(unsigned int i_clusterId,
      unsigned int i_globalClusterId,
      unsigned int profilingId,
      bool usePlasticity,
      LayerType layerType,
      double maxTimeStepSize,
      long timeStepRate,
      bool printProgress,
      DynamicRuptureScheduler* dynamicRuptureScheduler,
      CompoundGlobalData i_globalData,
      seissol::initializer::Layer *i_clusterData,
      seissol::initializer::Layer* dynRupInteriorData,
      seissol::initializer::Layer* dynRupCopyData,
      seissol::initializer::LTS* i_lts,
      seissol::initializer::DynamicRupture* i_dynRup,
      seissol::dr::friction_law::FrictionSolver* i_FrictionSolver,
      dr::output::OutputManager* i_faultOutputManager,
      seissol::SeisSol& seissolInstance,
      LoopStatistics* i_loopStatistics,
      ActorStateStatistics* actorStateStatistics);

  /**
   * Destructor of a LTS cluster.
   * TODO: Currently prints only statistics in debug mode.
   **/
  ~TimeCluster() override;

  /**
   * Sets the the cluster's point sources
   *
   * @param sourceCluster Contains point sources for cluster
   */
  void setPointSources(std::unique_ptr<kernels::PointSourceCluster> sourceCluster);
  void freePointSources() { m_sourceCluster.reset(nullptr); }

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

  std::vector<NeighborCluster>* getNeighborClusters();
};

#endif // TIMECLUSTER_H_
