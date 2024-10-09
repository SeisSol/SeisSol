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
 * Memory management of SeisSol.
 **/

#ifndef MEMORYMANAGER_H_
#define MEMORYMANAGER_H_

#include "Parameters/DRParameters.h"
#include "Tree/Layer.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include <DynamicRupture/FrictionLaws/FrictionSolver.h>
#include <array>
#ifdef USE_MPI
#include <mpi.h>
#endif

#include <utils/logger.h>

#include "Initializer/Typedefs.h"
#include "MemoryAllocator.h"

#include "Initializer/LTS.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Lut.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/InputAux.h"
#include "Initializer/Boundary.h"
#include "Initializer/ParameterDB.h"

#include "Physics/InitialField.h"

#include <vector>
#include <memory>

#include "DynamicRupture/Factory.h"
#include <yaml-cpp/yaml.h>

namespace seissol {
  class SeisSol;
  namespace initializer {

/**
 * Memory manager of SeisSol.
 **/
class MemoryManager {
  private: // explicit private for unit tests
    seissol::SeisSol& seissolInstance;

    //! memory allocator
    seissol::memory::ManagedAllocator m_memoryAllocator;

    //! LTS mesh structure
    struct MeshStructure *m_meshStructure;

    unsigned int* ltsToFace;

    /*
     * Interior
     */
    //! number of buffers in the interior per cluster
    unsigned int *m_numberOfInteriorBuffers;

    //! number of derivatives in the interior per cluster
    unsigned int *m_numberOfInteriorDerivatives;

#ifdef USE_MPI
    /*
     * Ghost layer
     */
    //! number of buffers in the ghost layer per cluster
    unsigned int  *m_numberOfGhostBuffers;

    //! number of buffers in the ghost regions per cluster
    unsigned int **m_numberOfGhostRegionBuffers;

    //! number of derivatives in the ghost layer per cluster
    unsigned int  *m_numberOfGhostDerivatives;

    //! number of derivatives in the ghost regions per cluster
    unsigned int **m_numberOfGhostRegionDerivatives;

    /*
     * Copy Layer
     */
    //! number of buffers in the copy layer per cluster
    unsigned int  *m_numberOfCopyBuffers;

    //! number of buffers in the copy regions per cluster
    unsigned int **m_numberOfCopyRegionBuffers;

    //! number of derivatives in the copy layer per cluster
    unsigned int  *m_numberOfCopyDerivatives;

    //! number of derivatives in the copy regionsper cluster
    unsigned int **m_numberOfCopyRegionDerivatives;
#endif

    /*
     * Cross-cluster
     */
    //! global data
    GlobalData            m_globalDataOnHost;
    GlobalData            m_globalDataOnDevice;

    //! Memory organization tree
    LTSTree               m_ltsTree;
    LTS                   m_lts;
    Lut                   m_ltsLut;

    std::vector<std::unique_ptr<physics::InitialField>> m_iniConds;

    std::array<LTSTree*, MULTIPLE_SIMULATIONS> m_dynRupTree = {nullptr};
    // LTSTree m_dynRupTree;
    std::array<std::shared_ptr<DynamicRupture>, MULTIPLE_SIMULATIONS> m_dynRup = {nullptr};
    // std::unique_ptr<DynamicRupture> m_dynRup = nullptr;
    std::array<std::shared_ptr<dr::initializer::BaseDRInitializer>, MULTIPLE_SIMULATIONS> m_DRInitializer = {nullptr};
    // std::unique_ptr<dr::initializer::BaseDRInitializer> m_DRInitializer = nullptr;
    std::array<std::shared_ptr<dr::friction_law::FrictionSolver>, MULTIPLE_SIMULATIONS> m_FrictionLaw = {nullptr};
    std::array<std::shared_ptr<dr::friction_law::FrictionSolver>, MULTIPLE_SIMULATIONS> m_FrictionLawDevice = {nullptr};
    // std::unique_ptr<dr::friction_law::FrictionSolver> m_FrictionLaw = nullptr;
    std::array<std::shared_ptr<dr::output::OutputManager>, MULTIPLE_SIMULATIONS> m_faultOutputManager = {nullptr};
    // std::unique_ptr<dr::output::OutputManager> m_faultOutputManager = nullptr;
    std::shared_ptr<seissol::initializer::parameters::SeisSolParameters> m_seissolParams = nullptr;

    LTSTree m_boundaryTree;
    Boundary m_boundary;

    EasiBoundary m_easiBoundary;

    /**
     * Corrects the LTS Setups (buffer or derivatives, never both) in the ghost region
     **/
    void correctGhostRegionSetups(); 

    /**
     * Derives the layouts -- number of buffers and derivatives -- of the layers.
     **/
    void deriveLayerLayouts();

    /**
     * Initializes the face neighbor pointers of the internal state.
     **/
    void initializeFaceNeighbors( unsigned    cluster,
                                  Layer& layer);

    /**
     * Initializes the pointers of the internal state.
     **/
    void initializeBuffersDerivatives();

    /**
    * Derives the size of the displacement accumulation buffer.
    */
    void deriveFaceDisplacementsBucket();

    /**
     * Derives the size of the displacement accumulation buffer.
     */
    void deriveDisplacementsBucket();

    /**
     * Initializes the displacement accumulation buffer.
     */
    void initializeDisplacements();

    /**
     * Initializes the displacement accumulation buffer.
     */
  void initializeFaceDisplacements();

#ifdef USE_MPI
    /**
     * Initializes the communication structure.
     **/
    void initializeCommunicationStructure();
#endif

  public:
    /**
     * Constructor
     **/
    MemoryManager(seissol::SeisSol& instance) : seissolInstance(instance) {};

    /**
     * Destructor, memory is freed by managed allocator
     **/
    ~MemoryManager() = default;
    
    /**
     * Initialization function, which allocates memory for the global matrices and initializes them.
     **/
    void initialize();
    
    /**
     * Sets the number of cells in each leaf of the lts tree, fixates the variables, and allocates memory.
     * Afterwards the tree cannot be changed anymore.
     *
     * @param i_meshStructrue mesh structure.
     **/
    void fixateLtsTree(struct TimeStepping& i_timeStepping,
                       struct MeshStructure*i_meshStructure,
                       unsigned* numberOfDRCopyFaces,
                       unsigned* numberOfDRInteriorFaces,
                       bool usePlasticity);

    void fixateBoundaryLtsTree();
    /**
     * Set up the internal structure.
     **/
    void initializeMemoryLayout();

    /**
     * Gets global data on the host.
     **/
    GlobalData* getGlobalDataOnHost() {
      return &m_globalDataOnHost;
    }

    /**
     * Gets the global data on device.
     **/
    GlobalData* getGlobalDataOnDevice() {
      assert(seissol::isDeviceOn() && "application is not compiled for acceleration device");
      return &m_globalDataOnDevice;
    }

    /**
     * Gets the global data on both host and device.
    **/
    CompoundGlobalData getGlobalData() {
      CompoundGlobalData global{};
      global.onHost = &m_globalDataOnHost;
      global.onDevice = nullptr;
      if constexpr (seissol::isDeviceOn()) {
        global.onDevice = &m_globalDataOnDevice;
      }
      return global;
    }

    /**
     * Gets the memory layout of a time cluster.
     *
     * @param i_cluster local id of the time cluster.
     * @param o_meshStructure mesh structure.
     * @param o_globalData global data.
     * @param o_globalDataCopies several copies of global data
     **/
    std::pair<MeshStructure*, CompoundGlobalData>
    getMemoryLayout(unsigned int i_cluster);
                          
    inline LTSTree* getLtsTree() {
      return &m_ltsTree;
    }
                          
    inline LTS* getLts() {
      return &m_lts;
    }

    inline Lut* getLtsLut() {
      return &m_ltsLut;
    }

    // TODO(David): remove again (this method is merely a temporary construction to transition from C++ to FORTRAN and should be removed in the next refactoring step)
    inline Lut& getLtsLutUnsafe() {
      return m_ltsLut;
    }

    inline std::array<LTSTree*, MULTIPLE_SIMULATIONS> getDynamicRuptureTree() {
      return m_dynRupTree;
    }
                          
    inline std::array<std::shared_ptr<DynamicRupture>, MULTIPLE_SIMULATIONS> getDynamicRupture() {
      return m_dynRup;
    }

    inline LTSTree* getBoundaryTree() {
      return &m_boundaryTree;
    }

    inline Boundary* getBoundary() {
      return &m_boundary;
    }

    inline void setInitialConditions(std::vector<std::unique_ptr<physics::InitialField>>&& iniConds) {
      m_iniConds = std::move(iniConds);
    }

    inline const std::vector<std::unique_ptr<physics::InitialField>>& getInitialConditions() {
      return m_iniConds;
    }

    inline void setLtsToFace(unsigned int* ptr) {
      ltsToFace = ptr;
    }

    inline unsigned int* ltsToFaceMap() const {
      return ltsToFace;
    }

    void initializeEasiBoundaryReader(const char* fileName);

    inline EasiBoundary* getEasiBoundaryReader() {
      return &m_easiBoundary;
    }

    inline std::array<std::shared_ptr<dr::friction_law::FrictionSolver>, MULTIPLE_SIMULATIONS> getFrictionLaw() {
        return m_FrictionLaw;
    }

    inline std::array<std::shared_ptr<dr::initializer::BaseDRInitializer>, MULTIPLE_SIMULATIONS> getDRInitializer() {
        return m_DRInitializer;
    }
    
    inline std::array<std::shared_ptr<dr::friction_law::FrictionSolver>, MULTIPLE_SIMULATIONS> getFrictionLawDevice() {
        return m_FrictionLawDevice;
    }
    
    inline std::array<std::shared_ptr<seissol::dr::output::OutputManager>, MULTIPLE_SIMULATIONS> getFaultOutputManager() {
        return m_faultOutputManager;
    }
    // inline seissol::initializer::parameters::DRParameters* getDRParameters() {
    //     return m_seissolParams->drParameters.data();
    // }

    inline std::array<std::shared_ptr<seissol::initializer::parameters::DRParameters>, MULTIPLE_SIMULATIONS> getDRParameters(){
      return m_seissolParams->drParameters;
    }

    inline seissol::initializer::parameters::LtsParameters* getLtsParameters() {
        return &(m_seissolParams->timeStepping.lts);
    };

    void setInputParams(std::shared_ptr<seissol::initializer::parameters::SeisSolParameters> params) {
      m_seissolParams = params;
    }

    std::string getOutputPrefix() const {
      return m_seissolParams->output.prefix;
    }

    bool isLoopStatisticsNetcdfOutputOn() const {
      return m_seissolParams->output.loopStatisticsNetcdfOutput;
    }

#ifdef ACL_DEVICE
  void recordExecutionPaths(bool usePlasticity);

  /**
   * Derives sizes of scratch memory required during computations of Wave Propagation solver
   **/
  static void deriveRequiredScratchpadMemoryForWp(LTSTree &ltsTree, LTS& lts);

  /**
   * Derives sizes of scratch memory required during computations of Dynamic Rupture solver
   **/
  static void deriveRequiredScratchpadMemoryForDr(LTSTree &ltsTree, DynamicRupture& dynRup);
#endif

  void initializeFrictionLaw();
  void initFaultOutputManager(const std::string& backupTimeStamp);
  void initFrictionData();
  void synchronizeTo(seissol::initializer::AllocationPlace place);
};


    bool isAcousticSideOfElasticAcousticInterface(CellMaterialData &material,
                                                  unsigned int face);
    bool isElasticSideOfElasticAcousticInterface(CellMaterialData &material,
                                                 unsigned int face);
    bool isAtElasticAcousticInterface(CellMaterialData &material, unsigned int face);

    bool requiresDisplacement(CellLocalInformation cellLocalInformation,
                              CellMaterialData &material,
                              unsigned int face);
    bool requiresNodalFlux(FaceType f);
    } // namespace initializer
} // namespace seissol

#endif
