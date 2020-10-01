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
#include <Solver/FreeSurfaceIntegrator.h>
#include <Monitoring/LoopStatistics.h>
#include "AbstractTimeCluster.h"

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
    double lastSubTime;

    void handleAdvancedPredictionTimeMessage(const NeighborCluster& neighborCluster) override;
    void handleAdvancedCorrectionTimeMessage(const NeighborCluster& neighborCluster) override;
    void start() override {}
    void predict() override;
    void correct() override;

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

  /*
   * global data
   */
     //! global data structures
    struct GlobalData *m_globalData;

    /*
     * element data
     */     
    seissol::initializers::Layer* m_clusterData;
    seissol::initializers::Layer* m_dynRupClusterData;
    seissol::initializers::LTS*         m_lts;
    seissol::initializers::DynamicRupture* m_dynRup;

    //! Mapping of cells to point sources
    sourceterm::CellToPointSourcesMapping const* m_cellToPointSources;

    //! Number of mapping of cells to point sources
    unsigned m_numberOfCellToPointSourcesMappings;

    //! Point sources
    sourceterm::PointSources const* m_pointSources;

    //! true if dynamic rupture faces are present
    bool m_dynamicRuptureFaces;

    enum class ComputePart {
      Local = 0,
      Neighbor,
      DRNeighbor,
      DRFrictionLaw,
      PlasticityCheck,
      PlasticityYield,
      NUM_COMPUTE_PARTS
    };

    long long m_flops_nonZero[static_cast<int>(ComputePart::NUM_COMPUTE_PARTS)];
    long long m_flops_hardware[static_cast<int>(ComputePart::NUM_COMPUTE_PARTS)];
    
    //! Tv parameter for plasticity
    double m_tv;
    
    //! Relax time for plasticity
    double m_relaxTime;
    
    //! Stopwatch of TimeManager
    LoopStatistics* m_loopStatistics;
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

    void computeNeighborIntegrationFlops(seissol::initializers::Layer &layerData);

    void computeDynamicRuptureFlops(seissol::initializers::Layer& layerData);
                                          
    void computeFlops();
    
    //! Update relax time for plasticity
    void updateRelaxTime() {
      m_relaxTime = (m_tv > 0.0) ? 1.0 - exp(-timeStepSize() / m_tv) : 1.0;
    }

public:
    const LayerType layerType;

    //! time of the next receiver output
    double m_receiverTime;

    //! print status every 100th timestep
    bool printProgress;

    /**
     * Constructs a new LTS cluster.
     *
     * @param i_clusterId id of this cluster with respect to the current rank.
     * @param i_globalClusterId global id of this cluster.
     * @param i_timeKernel time integration kernel.
     * @param i_volumeKernel volume integration kernel.
     * @param i_boundaryKernel boundary integration kernel.
     * @param i_copyCellInformation cell information in the copy layer.
     * @param i_interiorCellInformation cell information in the interior.
     * @param i_globalData global data.
     * @param i_copyCellData cell data in the copy layer.
     * @param i_interiorCellData cell data in the interior.
     * @param i_cells degrees of freedom, time buffers, time derivatives.
     **/
    TimeCluster(unsigned int i_clusterId,
                unsigned int i_globalClusterId,
                LayerType layerType,
                double maxTimeStepSize,
                long timeStepRate,
                double timeTolerance,
                bool printProgress,
                struct GlobalData *i_globalData,
                seissol::initializers::Layer *i_clusterData,
                seissol::initializers::Layer *i_dynRupClusterData,
                seissol::initializers::LTS *i_lts,
                seissol::initializers::DynamicRupture *i_dynRup,
                LoopStatistics *i_loopStatistics);

    /**
     * Destructor of a LTS cluster.
     * TODO: Currently prints only statistics in debug mode.
     **/
    ~TimeCluster() override;

    /**
     * Adds a source to the cluster.
     *
     * @param i_meshId mesh id of the point of interest.
     **/
    void addSource( unsigned int i_meshId );
    
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

    /**
     * Set Tv constant for plasticity.
     */
    void setTv(double tv) {
      m_tv = tv;
      updateRelaxTime();
    }

    //! cluster id on this rank
    const unsigned int m_clusterId;

    //! global cluster cluster id
    const unsigned int m_globalClusterId;

    void printTimeoutMessage(std::chrono::seconds timeSinceLastUpdate) override;
    void reset() override;
};

#endif
