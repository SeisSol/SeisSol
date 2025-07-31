// SPDX-FileCopyrightText: 2013 SeisSol Group
// SPDX-FileCopyrightText: 2015 Intel Corporation
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Alexander Breuer
// SPDX-FileContributor: Alexander Heinecke (Intel Corp.)

#ifndef SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_
#define SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_

#include "Memory/Tree/Layer.h"
#include "Initializer/Parameters/SeisSolParameters.h"
#include <Memory/Descriptor/Surface.h>
#include <Initializer/TimeStepping/ClusterLayout.h>
#include <Memory/Tree/Backmap.h>
#include <mpi.h>

#include <utils/logger.h>

#include "Initializer/Typedefs.h"
#include "Memory/MemoryAllocator.h"

#include "Memory/Descriptor/LTS.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Initializer/InputAux.h"
#include "Memory/Descriptor/Boundary.h"
#include "Initializer/ParameterDB.h"

#include "Physics/InitialField.h"

#include <utility>
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

    /*
     * Cross-cluster
     */
    //! global data
    GlobalData            m_globalDataOnHost;
    GlobalData            m_globalDataOnDevice;

    //! Memory organization storage
    LTS::Storage               ltsStorage;
    LTS::Backmap backmap;

    std::optional<ClusterLayout> clusterLayout;

    std::vector<std::unique_ptr<physics::InitialField>> m_iniConds;

    DynamicRupture::Storage drStorage;
    std::unique_ptr<DynamicRupture> m_dynRup = nullptr;
    std::unique_ptr<dr::initializer::BaseDRInitializer> m_DRInitializer = nullptr;
    std::unique_ptr<dr::friction_law::FrictionSolver> m_FrictionLaw = nullptr;
    std::unique_ptr<dr::friction_law::FrictionSolver> m_FrictionLawDevice = nullptr;
    std::unique_ptr<dr::output::OutputManager> m_faultOutputManager = nullptr;

    Boundary::Storage m_boundaryTree;

    SurfaceLTS::Storage surfaceStorage;

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
                                  LTS::Layer& layer);

    /**
     * Initializes the pointers of the internal state.
     **/
    void initializeBuffersDerivatives();

    /**
     * Derives the size of the displacement accumulation buffer.
     */
    void deriveDisplacementsBucket();

    /**
     * Initializes the displacement accumulation buffer.
     */
    void initializeDisplacements();

    /**
     * Initializes the communication structure.
     **/
    void initializeCommunicationStructure();

  public:
    /**
     * Constructor
     **/
    MemoryManager(seissol::SeisSol& instance) : seissolInstance(instance) {}

    /**
     * Destructor, memory is freed by managed allocator
     **/
    ~MemoryManager() = default;
    
    /**
     * Initialization function, which allocates memory for the global matrices and initializes them.
     **/
    void initialize();
    
    /**
     * Sets the number of cells in each leaf of the lts storage, fixates the variables, and allocates memory.
     * Afterwards the storage cannot be changed anymore.
     *
     * @param i_meshStructrue mesh structure.
     **/
    void fixateLtsStorage(struct ClusterLayout& clusterLayout,
                       struct MeshStructure* meshStructure,
                       const std::vector<std::size_t>& volumeSizes,
                       const std::vector<std::size_t>& drSizes,
                       bool usePlasticity);

    void fixateBoundaryStorage();
    /**
     * Set up the internal structure.
     **/
    void initializeMemoryLayout();

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
                          
    LTS::Storage& getLtsStorage() {
      return ltsStorage;
    }

    LTS::Backmap& getBackmap() {
      return backmap;
    }

    DynamicRupture::Storage& getDRStorage() {
      return drStorage;
    }

    DynamicRupture& getDynamicRupture() {
      return *m_dynRup;
    }

    SurfaceLTS::Storage& getSurfaceStorage() {
      return surfaceStorage;
    }

    void setInitialConditions(std::vector<std::unique_ptr<physics::InitialField>>&& iniConds) {
      m_iniConds = std::move(iniConds);
    }

    const std::vector<std::unique_ptr<physics::InitialField>>& getInitialConditions() {
      return m_iniConds;
    }

    void setLtsToFace(unsigned int* ptr) {
      ltsToFace = ptr;
    }

    unsigned int* ltsToFaceMap() const {
      return ltsToFace;
    }

    void initializeEasiBoundaryReader(const char* fileName);

    EasiBoundary* getEasiBoundaryReader() {
      return &m_easiBoundary;
    }

    dr::friction_law::FrictionSolver* getFrictionLaw() {
        return m_FrictionLaw.get();
    }
    dr::friction_law::FrictionSolver* getFrictionLawDevice() {
        return m_FrictionLawDevice.get();
    }
    seissol::dr::output::OutputManager* getFaultOutputManager() {
        return m_faultOutputManager.get();
    }

#ifdef ACL_DEVICE
  void recordExecutionPaths(bool usePlasticity);

  /**
   * Derives sizes of scratch memory required during computations of Wave Propagation solver
   **/
  static void deriveRequiredScratchpadMemoryForWp(bool plasticity, LTS::Storage& ltsStorage);

  /**
   * Derives sizes of scratch memory required during computations of Dynamic Rupture solver
   **/
  static void deriveRequiredScratchpadMemoryForDr(DynamicRupture::Storage& drStorage);
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
    }
}


#endif // SEISSOL_SRC_INITIALIZER_MEMORYMANAGER_H_
