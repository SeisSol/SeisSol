/*
Copyright (c) 2015, Intel Corporation

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

    * Redistributions of source code must retain the above copyright notice,
      this list of conditions and the following disclaimer.
    * Redistributions in binary form must reproduce the above copyright
      notice, this list of conditions and the following disclaimer in the
      documentation and/or other materials provided with the distribution.
    * Neither the name of Intel Corporation nor the names of its contributors
      may be used to endorse or promote products derived from this software
      without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/

/*
 * Copyright (c) 2013-2014, SeisSol Group
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
 **/
 
struct ProxyCells {
  unsigned int numberOfCells;
  real (*dofs)[NUMBER_OF_ALIGNED_DOFS];
  real **buffers;
  real **derivatives;
  real *(*faceNeighbors)[4];
};

GlobalData **m_globalDataArray;
GlobalData *m_globalData; 
CellLocalInformation *m_cellInformation;
CellData *m_cellData;
ProxyCells *m_cells;
LocalIntegrationData *m_localIntegration;
NeighboringIntegrationData * m_neighboringIntegration;

seissol::kernels::Time      m_timeKernel;
seissol::kernels::Local     m_localKernel;
seissol::kernels::Neighbor  m_neighborKernel;

seissol::memory::ManagedAllocator m_allocator;


/* This option is needed to avoid polution of low-level caches */
#define NUMBER_OF_THREADS_PER_GLOBALDATA_COPY 4
#ifndef NUMBER_OF_THREADS_PER_GLOBALDATA_COPY 
#define NUMBER_OF_THREADS_PER_GLOBALDATA_COPY 16383
#endif

#define HBM_DOFS seissol::memory::Standard
#define HBM_TDOFS seissol::memory::Standard
#define HBM_DERS seissol::memory::Standard
#define HBM_CELLLOCAL_LOCAL seissol::memory::Standard
#define HBM_CELLLOCAL_NEIGH seissol::memory::Standard
#define HBM_GLOBALDATA seissol::memory::Standard

real m_timeStepWidthSimulation = (real)1.0;
real* m_dofs;
real* m_tdofs;
#ifdef __USE_DERS
real* m_ders;
#endif
real** m_ptdofs;
real** m_pder;
real* m_faceNeighbors;
real** m_globalPointerArray;
real* m_globalPointer;

unsigned int readScenarioSize(std::string const& scenario)
{
  unsigned cells;
  unsigned reads;
  std::string file;
  std::ifstream data;

  file = scenario + ".size";
  data.open(file.c_str());
  if (!data) {
    printf("size of scenario couldn't be read!\n");
    exit(-1);
  }
  reads = 0;
  while (data >> cells) {
      printf("scenario name is: %s\n", scenario.c_str());
      printf("scenario has %i cells\n", cells);
      ++reads;
  }
  if (reads != 1) {
    printf("wrong number of sizes (%i) in scenario were read!\n", reads);
    exit(-1);
  }
  
  return cells;
}

void readQuadruples(std::string const& fileName, unsigned cells, unsigned (*target)[4])
{
  std::ifstream data;
  unsigned reads;
  unsigned value;
  
  target = new unsigned[cells][4];
  
  data.open(fileName.c_str());
  if (!data) {
    printf("scenario file %s couldn't be read!\n", fileName.c_str());
    exit(-1);
  }
    
  reads = 0;
  while (data >> value) {
    target[reads/4][reads%4] = value;
    reads++;
  }
  data.close();
  if (reads != cells*4) {
    printf("wrong (%i) in scenario file %s were read!\n", reads, fileName.c_str());
    exit(-1);
  }
}

void allocateEverything(unsigned i_cells, unsigned globalDataCopies)
{
  m_cellInformation = (CellLocalInformation*) m_allocator.allocateMemory(i_cells * sizeof(CellLocalInformation));
  m_dofs = (real*) m_allocator.allocateMemory(i_cells * NUMBER_OF_ALIGNED_DOFS * sizeof(real), PAGESIZE_HEAP, HBM_DOFS);
  m_tdofs = (real*) m_allocator.allocateMemory(i_cells * NUMBER_OF_ALIGNED_DOFS * sizeof(real), PAGESIZE_HEAP, HBM_TDOFS);
#ifdef __USE_DERS
  m_ders = (real*) m_allocator.allocateMemory(i_cells * NUMBER_OF_ALIGNED_DOFS * sizeof(real), PAGESIZE_HEAP, HBM_DERS);
#endif
  m_ptdofs = (real**) m_allocator.allocateMemory(i_cells * sizeof(real*));
  m_pder = (real**) m_allocator.allocateMemory(i_cells * sizeof(real*));
  m_faceNeighbors = (real*) m_allocator.allocateMemory(i_cells * sizeof(real*[4]));
  m_cells = (ProxyCells*) m_allocator.allocateMemory(sizeof(ProxyCells));
  m_cellData = (CellData*) m_allocator.allocateMemory(sizeof(CellData));
  m_localIntegration = (LocalIntegrationData*) m_allocator.allocateMemory(i_cells * sizeof(LocalIntegrationData), PAGESIZE_HEAP, HBM_CELLLOCAL_LOCAL);
  m_neighboringIntegration = (NeighboringIntegrationData*) m_allocator.allocateMemory(i_cells * sizeof(NeighboringIntegrationData), PAGESIZE_HEAP, HBM_CELLLOCAL_NEIGH);
  
  // global data
  m_globalPointerArray = (real**) m_allocator.allocateMemory(globalDataCopies * sizeof(real*));
  m_globalDataArray = (GlobalData**) m_allocator.allocateMemory(globalDataCopies * sizeof(GlobalData*));
  for (unsigned int global = 0; global < globalDataCopies; ++global) {
    m_globalPointerArray[global] = (real*) m_allocator.allocateMemory(  seissol::model::globalMatrixOffsets[seissol::model::numGlobalMatrices] * sizeof(real),
                                                                        PAGESIZE_HEAP,
                                                                        HBM_GLOBALDATA  );
    m_globalDataArray[global] = (GlobalData*) m_allocator.allocateMemory(sizeof(GlobalData));
  }
}

unsigned int init_data_structures(unsigned int i_cells)
{
  char* pScenario;
  std::string s_scenario;
  bool bUseScenario = false;
  unsigned int (*scenario_faceType)[4];
  unsigned int (*scenario_neighbor)[4];
  unsigned int (*scenario_side)[4];
  unsigned int (*scenario_orientation)[4];

  pScenario = getenv ("SEISSOL_PROXY_SCENARIO");
  if (pScenario !=NULL ) {
    bUseScenario = true;
    s_scenario.assign(pScenario);
    
    i_cells = readScenarioSize(s_scenario);

    readQuadruples(s_scenario + ".neigh", i_cells, scenario_neighbor);
    readQuadruples(s_scenario + ".bound", i_cells, scenario_faceType);
    readQuadruples(s_scenario + ".sides", i_cells, scenario_side);
    readQuadruples(s_scenario + ".orient", i_cells, scenario_orientation);
  }

  // init RNG
  srand48(i_cells);  
  
  // determine number of global data copies
  unsigned int l_numberOfThreads = 1;
#ifdef _OPENMP
  #pragma omp parallel
  {
    #pragma omp master
    {
      l_numberOfThreads = omp_get_num_threads();
    }
  }
#endif
  unsigned int l_numberOfCopiesCeil = (l_numberOfThreads%NUMBER_OF_THREADS_PER_GLOBALDATA_COPY == 0) ? 0 : 1;
  unsigned int l_numberOfCopies = (l_numberOfThreads/NUMBER_OF_THREADS_PER_GLOBALDATA_COPY) + l_numberOfCopiesCeil;

  allocateEverything(i_cells, l_numberOfCopies);

  /* cell information */
  for (unsigned int l_cell = 0; l_cell < i_cells; l_cell++) {
    for (unsigned int f = 0; f < 4; f++) {
      if (bUseScenario == true ) {
        switch (scenario_faceType[l_cell][f]) {
          case 0: 
            m_cellInformation[l_cell].faceTypes[f] = regular;
            break;
          case 1:
            m_cellInformation[l_cell].faceTypes[f] = freeSurface;
            break;
          case 3:
            m_cellInformation[l_cell].faceTypes[f] = dynamicRupture;
            break;
          case 5:
            m_cellInformation[l_cell].faceTypes[f] = outflow;
            break;
          case 6:
            m_cellInformation[l_cell].faceTypes[f] = periodic;
            break;
          default:
            printf("unsupported faceType (%i)!\n", scenario_faceType[l_cell][f]); 
            exit(-1);
            break;
        }
        m_cellInformation[l_cell].faceRelations[f][0] = scenario_side[l_cell][f];
        m_cellInformation[l_cell].faceRelations[f][1] = scenario_orientation[l_cell][f];
        m_cellInformation[l_cell].faceNeighborIds[f] =  scenario_neighbor[l_cell][f];
      } else {
        m_cellInformation[l_cell].faceTypes[f] = regular;
        m_cellInformation[l_cell].faceRelations[f][0] = ((unsigned int)lrand48() % 4);
        m_cellInformation[l_cell].faceRelations[f][1] = ((unsigned int)lrand48() % 3);
        m_cellInformation[l_cell].faceNeighborIds[f] = ((unsigned int)lrand48() % i_cells);
      }
    }
#ifdef __USE_DERS
    m_cellInformation[l_cell].ltsSetup = 4095;
#else
    m_cellInformation[l_cell].ltsSetup = 0;
#endif
    m_cellInformation[l_cell].timeStepWidth = 1.0;
  }

  /* init dofs */
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned int l_cell = 0; l_cell < i_cells; l_cell++) {
    for (unsigned int i = 0; i < NUMBER_OF_ALIGNED_DOFS; i++) {
      m_dofs[(l_cell*NUMBER_OF_ALIGNED_DOFS)+i] = (real)drand48();
    }
    for (unsigned int i = 0; i < NUMBER_OF_ALIGNED_DOFS; i++) {
      m_tdofs[(l_cell*NUMBER_OF_ALIGNED_DOFS)+i] = (real)drand48();
    }
  }

  /* init ders */
#ifdef __USE_DERS
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned int l_cell = 0; l_cell < i_cells; l_cell++) {
    for (unsigned int i = 0; i < NUMBER_OF_ALIGNED_DERS; i++) {
      m_ders[(l_cell*NUMBER_OF_ALIGNED_DERS)+i] = (real)drand48();
    }
  }
#endif

  /* init dof pointers */
  for (unsigned int l_cell = 0; l_cell < i_cells; l_cell++) {
    m_ptdofs[l_cell] = &(m_tdofs[(l_cell*NUMBER_OF_ALIGNED_DOFS)]);
#ifdef __USE_DERS
    m_pder[l_cell] = &(m_ders[(l_cell*NUMBER_OF_ALIGNED_DERS)]);
#else
    m_pder[l_cell] = NULL; 
#endif
  }

  /* init cell local information */
  m_cells->numberOfCells = i_cells;
  m_cells->dofs = (real(*)[NUMBER_OF_ALIGNED_DOFS])m_dofs;
  m_cells->buffers = m_ptdofs;
  m_cells->derivatives = m_pder;
  m_cells->faceNeighbors = (real*(*)[4])m_faceNeighbors;

  for (unsigned int l_cell = 0; l_cell < i_cells; l_cell++) {
    for (unsigned int f = 0; f < 4; f++) {
      if (m_cellInformation[l_cell].faceTypes[f] == outflow) {
        m_cells->faceNeighbors[l_cell][f] = NULL;
      } else if (m_cellInformation[l_cell].faceTypes[f] == freeSurface) {
#ifdef __USE_DERS
        m_cells->faceNeighbors[l_cell][f] = m_cells->derivatives[l_cell];
#else
        m_cells->faceNeighbors[l_cell][f] = m_cells->buffers[l_cell];
#endif
      } else if (m_cellInformation[l_cell].faceTypes[f] == periodic || m_cellInformation[l_cell].faceTypes[f] == regular) {
#ifdef __USE_DERS
        m_cells->faceNeighbors[l_cell][f] = m_cells->derivatives[m_cellInformation[l_cell].faceNeighborIds[f]];
#else
        m_cells->faceNeighbors[l_cell][f] = m_cells->buffers[m_cellInformation[l_cell].faceNeighborIds[f]];
#endif
      } else {
        printf("unsupported boundary type -> exit\n");
        exit(-1);
      }
    }
  }

  /* init local integration data */
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned int l_cell = 0; l_cell < i_cells; l_cell++) {
    // init star matrices
    for (size_t m = 0; m < 3; m++) { 
      for (size_t j = 0; j < seissol::model::AstarT::reals; j++) { 
        m_localIntegration[l_cell].starMatrices[m][j] = (real)drand48();
      }
    }
    // init flux solver
    for (size_t m = 0; m < 4; m++) { 
      for (size_t j = 0; j < seissol::model::AplusT::reals; j++) { 
        m_localIntegration[l_cell].nApNm1[m][j] = (real)drand48();
      }
    }
#ifdef EQUATIONS_VISCOELASTIC
    // init source matrix
    for (size_t j = 0; j < seissol::model::source::reals; j++) { 
      m_localIntegration[l_cell].specific.sourceMatrix[j] = (real)drand48();
    }
#endif
#ifdef EQUATIONS_VISCOELASTIC2
    // init source matrix
    for (size_t j = 0; j < NUMBER_OF_RELAXATION_MECHANISMS * seissol::model::ET::reals; j++) { 
      m_localIntegration[l_cell].specific.ET[j] = (real)drand48();
    }
    // init omega
    for (size_t j = 0; j < NUMBER_OF_RELAXATION_MECHANISMS; j++) { 
      m_localIntegration[l_cell].specific.omega[j] = (real)drand48();
    }
#endif
  }

  /* init neighbor integration data */
#ifdef _OPENMP
  #pragma omp parallel for schedule(static)
#endif
  for (unsigned int l_cell = 0; l_cell < i_cells; l_cell++) {
    // init flux solver
    for (size_t m = 0; m < 4; m++) { 
      for (size_t j = 0; j < seissol::model::AminusT::reals; j++) { 
        m_neighboringIntegration[l_cell].nAmNm1[m][j] = (real)drand48();
      }      
#ifdef EQUATIONS_VISCOELASTIC2
      // init omega
      for (size_t j = 0; j < NUMBER_OF_RELAXATION_MECHANISMS; j++) { 
        m_neighboringIntegration[l_cell].specific.omega[j] = (real)drand48();
      }
#endif
    }
  }

  // CellData
  m_cellData->localIntegration = m_localIntegration;
  m_cellData->neighboringIntegration = m_neighboringIntegration;

  // Global matrices
  unsigned int l_globalMatrices  = seissol::model::globalMatrixOffsets[seissol::model::numGlobalMatrices] * sizeof(real);
  // @TODO: for NUMA we need to bind this
  for (unsigned int l_globalDataCount = 0; l_globalDataCount < l_numberOfCopies; l_globalDataCount++) {
    m_globalPointer = m_globalPointerArray[l_globalDataCount];
    m_globalData =  m_globalDataArray[l_globalDataCount];
   
    for (unsigned int i = 0; i < (l_globalMatrices/sizeof(real)); i++) {
      m_globalPointer[i] = (real)drand48();
    }

    // stiffnes for time integration
    for( unsigned int l_transposedStiffnessMatrix = 0; l_transposedStiffnessMatrix < 3; l_transposedStiffnessMatrix++ ) {
      m_globalData->stiffnessMatricesTransposed[l_transposedStiffnessMatrix] = m_globalPointer + seissol::model::globalMatrixOffsets[l_transposedStiffnessMatrix];
    }

    // stiffnes for volume integration
    for( unsigned int l_stiffnessMatrix = 0; l_stiffnessMatrix < 3; l_stiffnessMatrix++ ) {
      m_globalData->stiffnessMatrices[l_stiffnessMatrix] = m_globalPointer + seissol::model::globalMatrixOffsets[3 + l_stiffnessMatrix];
    }

    // flux matrices for boundary integration
    for( unsigned int l_fluxMatrix = 0; l_fluxMatrix < 52; l_fluxMatrix++ ) {
      m_globalData->fluxMatrices[l_fluxMatrix] = m_globalPointer + seissol::model::globalMatrixOffsets[6 + l_fluxMatrix];
    }
  }

  // set default to first chunk
  m_globalPointer = m_globalPointerArray[0];
  m_globalData = m_globalDataArray[0];

  if (bUseScenario == true ) {
    delete[] scenario_faceType;
    delete[] scenario_neighbor;
    delete[] scenario_side;
    delete[] scenario_orientation;
  }
  
  return i_cells;
}
