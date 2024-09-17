// SPDX-FileCopyrightText: 2013-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

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
 * @author Alex Breuer (breuer AT mytum.de,
 *http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
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

#include <Initializer/Tree/Layer.h>
#include <SeisSol.h>
#ifdef USE_MPI
#include <list>
#include <mpi.h>
#endif

#include "Initializer/LTS.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Typedefs.h"
#include "SourceTerm/Typedefs.h"
#include <utils/logger.h>

#include "Initializer/DynamicRupture.h"
#include "Kernels/Local.h"
#include "Kernels/Neighbor.h"
#include "Kernels/Plasticity.h"
#include "Kernels/PointSourceCluster.h"
#include "Kernels/Time.h"
#include "Kernels/TimeCommon.h"
#include "Monitoring/ActorStateStatistics.h"
#include "Monitoring/LoopStatistics.h"
#include "Solver/FreeSurfaceIntegrator.h"
#include <Common/Executor.h>

#include "Solver/Clustering/AbstractTimeCluster.h"

#ifdef ACL_DEVICE
#include <device.h>
#endif

namespace seissol {

namespace kernels {
class ReceiverCluster;
} // namespace kernels

namespace solver::clustering::computation {

/**
 * Time cluster, which represents a collection of elements having the same time step width.
 **/
class TimeCluster : public CellCluster {

  protected:
  void runCompute(ComputeStep step) override;

  private:
  // Last correction time of the neighboring cluster with higher dt
  double lastSubTime;

  void handleAdvancedComputeTimeMessage(ComputeStep,
                                        const NeighborCluster& neighborCluster) override;
  void start() override {}
  void predict();
  void correct();
  bool usePlasticity;
  //! number of time steps
  unsigned long numberOfTimeSteps;

  seissol::SeisSol& seissolInstance;
  /*
   * integrators
   */
  //! time kernel
  kernels::Time timeKernel;

  //! local kernel
  kernels::Local localKernel;

  //! neighbor kernel
  kernels::Neighbor neighborKernel;

  /*
   * global data
   */
  //! global data structures
  GlobalData* globalDataOnHost{nullptr};
  GlobalData* globalDataOnDevice{nullptr};
#ifdef ACL_DEVICE
  device::DeviceInstance& device = device::DeviceInstance::getInstance();
#endif

  /*
   * element data
   */
  seissol::initializer::Layer* layer;
  seissol::initializer::LTS* lts;

  seissol::kernels::PointSourceClusterPair sourceCluster;

  enum class ComputePart {
    Local = 0,
    Neighbor,
    DRNeighbor,
    DRFrictionLawInterior,
    DRFrictionLawCopy,
    PlasticityCheck,
    PlasticityYield,
    NumComputeParts
  };

  long long flops_nonZero[static_cast<int>(ComputePart::NumComputeParts)];
  long long flops_hardware[static_cast<int>(ComputePart::NumComputeParts)];

  //! Tv parameter for plasticity
  double tv;

  //! Stopwatch of TimeManager
  LoopStatistics* loopStatistics;
  ActorStateStatistics* actorStateStatistics;
  unsigned regionComputeLocalIntegration;
  unsigned regionComputeNeighboringIntegration;
  unsigned regionComputePointSources;

  kernels::ReceiverCluster* receiverCluster;

  /**
   * Writes the receiver output if applicable (receivers present, receivers have to be written).
   **/
  void writeReceivers();

  /**
   * Computes the source terms if applicable.
   **/
  void computeSources();

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
  void computeLocalIntegration(bool resetBuffers);

  /**
   * Computes the contribution of the neighboring cells to the boundary integral.
   *
   * Remark: After this step (in combination with the local integration) the DOFs are at the next
   *time step.
   * TODO: This excludes dynamic rupture contribution.
   *
   * @param i_numberOfCells number of cells.
   * @param i_cellInformation cell local information.
   * @param i_cellData cell data.
   * @param i_faceNeighbors pointers to neighboring time buffers or derivatives.
   * @param io_dofs degrees of freedom.
   **/
  void computeNeighboringIntegration(double subTimeStart);

#ifdef ACL_DEVICE
  void computeLocalIntegrationDevice(bool resetBuffers);
  void computeNeighboringIntegrationDevice(double subTimeStart);
#endif

  void computeLocalIntegrationFlops();

  template <bool UsePlasticity>
  std::pair<long, long> computeNeighboringIntegrationImplementation(double subTimeStart);

  void computeLocalIntegrationFlops(unsigned numberOfCells,
                                    const CellLocalInformation* cellInformation,
                                    long long& nonZeroFlops,
                                    long long& hardwareFlops);

  void computeNeighborIntegrationFlops();

  void computeFlops();

  //! Update relax time for plasticity
  double getRelaxTime() { return (tv > 0.0) ? 1.0 - exp(-timeStepSize() / tv) : 1.0; }

  const LayerType layerType;
  //! time of the next receiver output
  double receiverTime;

  //! print status every 100th timestep
  bool printProgress;
  //! cluster id on this rank
  const unsigned int clusterId;

  //! global cluster cluster id
  const unsigned int globalClusterId;

  //! id used to identify this cluster (including layer type) when profiling
  const unsigned int profilingId;

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
  TimeCluster(unsigned int clusterId,
              unsigned int globalClusterId,
              unsigned int profilingId,
              bool usePlasticity,
              double maxTimeStepSize,
              long timeStepRate,
              bool printProgress,
              CompoundGlobalData globalData,
              seissol::initializer::Layer* clusterData,
              seissol::initializer::LTS* lts,
              seissol::SeisSol& seissolInstance,
              LoopStatistics* loopStatistics,
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
  void setPointSources(seissol::kernels::PointSourceClusterPair sourceCluster);
  void freePointSources() {
    sourceCluster.host.reset(nullptr);
    sourceCluster.device.reset(nullptr);
  }

  void setReceiverCluster(kernels::ReceiverCluster* receiverCluster) {
    this->receiverCluster = receiverCluster;
  }

  /**
   * Set Tv constant for plasticity.
   */
  void setTv(double tv) { tv = tv; }

  void finalize() override;

  [[nodiscard]] unsigned int getClusterId() const;
  [[nodiscard]] unsigned int getGlobalClusterId() const;
  LayerType getLayerType() const override { return Interior; }
  void setReceiverTime(double receiverTime);

  std::vector<NeighborCluster>* getNeighborClusters();

  void synchronizeTo(seissol::initializer::AllocationPlace place, void* stream);

  std::string description() const override { return "cell-cluster"; }
};

} // namespace solver::clustering::computation
} // namespace seissol

#endif
