/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 * 
 * @section LICENSE
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
 *
 * @section DESCRIPTION
 * Test suite, which tests the time manager.
 **/

// setup path of the XML-file, which contains information about the matrices.
// TODO: Add this to input-parameters of SeisSol.
#ifdef __MIC__
#define MATRIXXMLFILE "matrices_" STR(NUMBEROFBASISFUNCTIONS) ".mic.xml"
#else
#define MATRIXXMLFILE "matrices_" STR(NUMBEROFBASISFUNCTIONS) ".xml"
#endif

#include <limits>
#include <algorithm>

#include <cxxtest/TestSuite.h>
#include <utils/logger.h>

#include <Initializer/typedefs.hpp>
#include <Solver/time_stepping/TimeManager.h>
#include <seissol_kernels/unit_tests/SimpleTimeIntegrator.hpp>
#include <seissol_kernels/unit_tests/SimpleVolumeIntegrator.hpp>
#include <seissol_kernels/unit_tests/SimpleBoundaryIntegrator.hpp>

struct cellLocalDataDense {
  real aStar[NUMBEROFVARIABLESSQUARED];
  real bStar[NUMBEROFVARIABLESSQUARED];
  real cStar[NUMBEROFVARIABLESSQUARED];
};

namespace unit_test {
  class TimeManagerTestSuite;
}

class unit_test::TimeManagerTestSuite: public CxxTest::TestSuite {
  //private:
    //! number of cells
    int m_numberOfCells;

    //! data structure for computations
    real  *m_degreesOfFreedomTimeDerivatives;
    real  *m_degreesOfFreedomUnitTests;
    real  *m_degreesOfFreedomTimeDerivativesUnitTests;
    real  *m_degreesOfFreedomTimeIntegratedUnitTests;

    //! cell local data: star matrices, flux solvers
    struct CellLocalInformation *m_cellInformation;
    struct CellLocalData      *m_cellData;
    struct CellLocalDataDense *m_cellDataDense;
    struct degreesOfFreedom m_degreesOfFreedom;
 
    //! simple element-based time integreator
    unit_test::SimpleTimeIntegrator m_simpleTimeIntegrator;

    //! simple element-based volume integrator
    unit_test::SimpleVolumeIntegrator m_simpleVolumeIntegrator;

    //! simple element-based boundary integrator
    unit_test::SimpleBoundaryIntegrator m_simpleBoundaryIntegrator;

    //! dense matrix functionality
    unit_test::DenseMatrix m_denseMatrix;

  public:
    TimeManagerTestSuite() {
      m_numberOfCells = 253;

      // allocate memory
      m_degreesOfFreedom.regular                 = (real (*)[NUMBEROFUNKNOWNS]) _mm_malloc( NUMBEROFUNKNOWNS*m_numberOfCells*sizeof(real), 64                              );
      m_degreesOfFreedomTimeDerivatives          = (real*)                      _mm_malloc( ORDEROFTAYLORSERIESEXPANSION*NUMBEROFUNKNOWNS*m_numberOfCells*sizeof(real), 64 ); 
      m_degreesOfFreedom.timeIntegrated          = (real (*)[NUMBEROFUNKNOWNS]) _mm_malloc( NUMBEROFUNKNOWNS*m_numberOfCells*sizeof(real), 64                              ); 
      m_degreesOfFreedomUnitTests                = (real*)                      _mm_malloc( NUMBEROFUNKNOWNS*m_numberOfCells*sizeof(real), 64                              );
      m_degreesOfFreedomTimeDerivativesUnitTests = (real*)                      _mm_malloc( ORDEROFTAYLORSERIESEXPANSION*NUMBEROFUNKNOWNS*m_numberOfCells*sizeof(real), 64 ); 
      m_degreesOfFreedomTimeIntegratedUnitTests  = (real*)                      _mm_malloc( NUMBEROFUNKNOWNS*m_numberOfCells*sizeof(real), 64                              ); 
      
      // allocate memory for the star matrices
      m_cellInformation = new struct CellLocalInformation[ m_numberOfCells ];
      m_cellData        = new struct CellLocalData[        m_numberOfCells ];
      m_cellDataDense   = new struct cellLocalDataDense[   m_numberOfCells ];
      
      // path to matrices
      std::string l_matricesPath = MATRIXXMLFILE;

      // read matrix structures for our unit tests
      m_denseMatrix.readMatrices( l_matricesPath );
     
      // initialize the integrators
      m_simpleTimeIntegrator.initialize(     l_matricesPath );
      m_simpleVolumeIntegrator.initialize(   l_matricesPath );
      m_simpleBoundaryIntegrator.initialize( l_matricesPath ); 

    }

    ~TimeManagerTestSuite() {
      // free memory
      _mm_free( m_degreesOfFreedom.regular );
      _mm_free( m_degreesOfFreedomTimeDerivatives );
      _mm_free( m_degreesOfFreedom.timeIntegrated );
      _mm_free( m_degreesOfFreedomUnitTests );
      _mm_free( m_degreesOfFreedomTimeDerivativesUnitTests );
      _mm_free( m_degreesOfFreedomTimeIntegratedUnitTests );

      delete[] m_cellInformation;
      delete[] m_cellData;
      delete[] m_cellDataDense;
    }

    void setUp() {
      // initialize random seed
      srand (time(NULL));

      // initialize dofs and time integrated dofs
      m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS*m_numberOfCells, &m_degreesOfFreedom.regular[0][0]        ); 
      m_denseMatrix.setRandomValues( NUMBEROFUNKNOWNS*m_numberOfCells, &m_degreesOfFreedom.timeIntegrated[0][0] );
      std::copy( &m_degreesOfFreedom.regular[0][0],        &m_degreesOfFreedom.regular[       m_numberOfCells][0], m_degreesOfFreedomUnitTests               );
      std::copy( &m_degreesOfFreedom.timeIntegrated[0][0], &m_degreesOfFreedom.timeIntegrated[m_numberOfCells][0], m_degreesOfFreedomTimeIntegratedUnitTests );
      
      // initialize a regular computational domain
      for( int l_cell = 0; l_cell < m_numberOfCells; l_cell++ ) {
        for( int l_face = 0; l_face < 4; l_face++ ) {
          m_cellInformation[l_cell].faceTypes[l_face] = regular;
          m_cellInformation[l_cell].faceNeighborIds[l_face]  = rand()%m_numberOfCells; 
          m_cellInformation[l_cell].faceRelations[l_face][0] = rand()%4;
          m_cellInformation[l_cell].faceRelations[l_face][1] = rand()%3;
        }
      }
      
      // iterate over the three coordinates \f$ \xi, \eta, \zeta \f$
      for( unsigned int l_coordinate = 0; l_coordinate < 3; l_coordinate++ ) { 
        for( int l_cell = 0; l_cell < m_numberOfCells; l_cell++ ) {
          // initialize dense star matrices
          m_denseMatrix.initializeMatrix(  59,
                                           NUMBEROFVARIABLES,
                                           NUMBEROFVARIABLES,
                                           m_cellDataDense[l_cell].aStar + (l_coordinate*NUMBEROFVARIABLESSQUARED) );

          // copy non-zeros to sparse star matrix
          m_denseMatrix.copyDenseToSparse( NUMBEROFVARIABLES,
                                           NUMBEROFVARIABLES,
                                           m_cellDataDense[l_cell].aStar + (l_coordinate*NUMBEROFVARIABLESSQUARED),
                                           m_cellData[l_cell].aStar + (STARMATRIX_NUMBEROFNONZEROS*l_coordinate) );
        }
      }
    }

    /**
     * Tests the time step width management, which operates on time step widths and synchronization points.
     **/
    void testGlobalTimeStepWidthManagement() {
      // time manager
      seissol::time_stepping::TimeManager l_timeManager;
      
      // maximum time step width
      real l_maximumTimeStepWidth;
      
      // verify the default (infinity), where no cluster is present
      l_timeManager.getMaximumTimeStepWidth( l_maximumTimeStepWidth );
      TS_ASSERT_DELTA( l_maximumTimeStepWidth, std::numeric_limits<real>::max(), s_zeroTolerance );

      // add a single cluster 
      l_timeManager.addCluster( m_numberOfCells,
                                0.9,
                                m_cellInformation,
                                m_cellData,
                                m_degreesOfFreedom );
      
      // set next synchronization time
      l_timeManager.setSynchronizationPoint( 0.5 );

      // assert time step width matches the synchronization point
      l_timeManager.getMaximumTimeStepWidth( l_maximumTimeStepWidth );
      TS_ASSERT_DELTA( l_maximumTimeStepWidth, 0.5, s_zeroTolerance );

      // advance in time an check that the synchronization point was reached
      l_timeManager.advanceInTime();
      l_timeManager.getMaximumTimeStepWidth( l_maximumTimeStepWidth );
      TS_ASSERT_DELTA( l_maximumTimeStepWidth, 0.0, s_zeroTolerance );

      // set a new synchronization point far away in time
      l_timeManager.setSynchronizationPoint( 360.6 );
      for( int l_timeStep = 0; l_timeStep < 400; l_timeStep++ ) {
        l_timeManager.getMaximumTimeStepWidth( l_maximumTimeStepWidth );
        TS_ASSERT_DELTA( l_maximumTimeStepWidth, 0.9, s_zeroTolerance ); 
        
        l_timeManager.advanceInTime();
      }

      // assert the time step width of the last update before the synchroniaztion point matches
      l_timeManager.getMaximumTimeStepWidth( l_maximumTimeStepWidth );
      TS_ASSERT_DELTA( l_maximumTimeStepWidth, 0.1, s_zeroTolerance ); 

    }

    /**
     * Tests the time step width computation
     **/
    void testTimeStepWidthComputation() {
/*      logDebug() << "testing the time step width calculation";

      const int l_numberOfElements = 3;
      int l_numberOfWaves = 5;

      real l_waveSpeeds[l_numberOfElements] = {1,2,3};
      real l_inRadii[l_numberOfElements]    = {1,2,3};

      seissol::TimeManager l_timeManager;

      l_timeManager.computeMaximumTimeStepWidth( l_waveSpeeds,
                                                 l_inRadii,
                                                 l_numberOfWaves,
                                                 l_numberOfElements );
*/
    }

    /**
     * Test the face setup building an internal pointer structure base on the face type and neighbors.
     **/
    void testFaceSetup() {
      /*
       * Example Setup:
       *
       *  9 Tetrahedrons, where tetrahedron 0 to 7 are part of a cube.
       *  Tetrahedron 8 is attached to 6 -> 6 lies completely inside
       *  the computational domain
       *
       *
       *  Sketch of a cut through the cube and tetrahedron 8 (+++),
       *    all edges for tetrahedron 5 are shown:
       *  
       *                       *
       *                     * *    *
       *                   *   *         *
       *                 *     *             *
       *               ******* * ****************+++++++++++
       *              * *      *  5/1        *  *        +
       *             *    *    *         *     *       +
       *            *       *  *     *        *  8   +
       *           *   0/4    ** *           *     +
       *          *          *  *     6/3   *    +
       *         *       *        *        *   +
       *        *    *       2/7    *     *  +
       *       * *                     * * +
       *      ***************************
       *
       *  This leads to the following regular face neighbors with arbitrary chosen face numbers:
       *   -- format: (cell, face)  --
       *
       *   1)  (0,1) - (4,1)
       *   2)  (0,2) - (2,0)
       *   3)  (0,0) - (5,1)
       *   4)  (1,3) - (3,3)
       *   5)  (1,0) - (4,3)
       *   6)  (1,1) - (5,3)
       *   7)  (2,1) - (6,1)
       *   8)  (2,2) - (7,0)
       *   9)  (3,1) - (7,2)
       *   10) (3,2) - (6,0)
       *   11) (4,2) - (7,3)
       *   12) (5,0) - (6,2)
       *   13) (6,3) - (8,3)
       *
       *
       *  Except 6, all of the cube-tetrahedrons have a unique face (second index) not part of
       *  the computational domain. The type these faces is defined as:
       *   (0,3) periodic -> tetrahedron 7
       *   (1,2) MPI-1
       *   (2,3) free surface 
       *   (3,0) MPI (shared with 9)-0
       *   (4,0) outflow 
       *   (5,2) free surface
       *   (7,1) periodic -> tetrahedron 0
       *
       *  Tetrahedron 8 is connected to the computational domain only via
       *  face 3, which leads to three face-types:
       *   (8,0) outflow
       *   (8,1) MPI-2
       *   (8,2) MPI (shared with 3)-0
       *
       *  Tetrahedron 2 and 7 share a dynamic rupture face:
       *  Faces (2,2) and (7,0) are dynamic rupture faces.
       */

      // allocate memory 
      cellLocalInformation *l_cellInformation = new CellLocalInformation[9];
      CellLocalData        *l_cellData        = new CellLocalData[9];
      
      degreesOfFreedom l_degreesOfFreedom;
      l_degreesOfFreedom.regular           = new real[9][NUMBEROFUNKNOWNS];
      l_degreesOfFreedom.timeIntegrated    = new real[9][NUMBEROFUNKNOWNS];
      l_degreesOfFreedom.timeIntegratedMpi = new real[3][NUMBEROFUNKNOWNS];

      // set up face types
      for( int l_tetrahedron = 0; l_tetrahedron < 9; l_tetrahedron++ ) {
        for( int l_face = 0; l_face < 4; l_face++ ) {
          l_cellInformation[l_tetrahedron].faceTypes[l_face] = regular;
        }
      }
      l_cellInformation[0].faceTypes[3] = periodic;
      
      l_cellInformation[1].faceTypes[2] = mpi;
      
      l_cellInformation[2].faceTypes[2] = dynamicRupture;
      l_cellInformation[2].faceTypes[3] = freeSurface;
      
      l_cellInformation[3].faceTypes[0] = mpi;

      l_cellInformation[4].faceTypes[0] = outflow;
      
      l_cellInformation[4].faceTypes[0] = outflow;

      l_cellInformation[5].faceTypes[2] = freeSurface;

      l_cellInformation[7].faceTypes[1] = periodic;
      
      l_cellInformation[8].faceTypes[0] = outflow;
      l_cellInformation[8].faceTypes[1] = mpi;
      l_cellInformation[8].faceTypes[2] = mpi; 
      
      l_cellInformation[2].faceTypes[2] = dynamicRupture;
      l_cellInformation[7].faceTypes[0] = dynamicRupture; 

      // setup internal structure of faces
      l_cellInformation[0].faceNeighborIds[0] =  5;
      l_cellInformation[0].faceNeighborIds[1] =  4;
      l_cellInformation[0].faceNeighborIds[2] =  2;
      l_cellInformation[0].faceNeighborIds[3] =  7;
       
      l_cellInformation[1].faceNeighborIds[0] =  4;
      l_cellInformation[1].faceNeighborIds[1] =  5;
      l_cellInformation[1].faceNeighborIds[2] =  1;
      l_cellInformation[1].faceNeighborIds[3] =  3;
      
      l_cellInformation[2].faceNeighborIds[0] =  0;
      l_cellInformation[2].faceNeighborIds[1] =  6;
      l_cellInformation[2].faceNeighborIds[2] =  7;
      l_cellInformation[2].faceNeighborIds[3] = -1;
      
      l_cellInformation[3].faceNeighborIds[0] =  0;
      l_cellInformation[3].faceNeighborIds[1] =  7;
      l_cellInformation[3].faceNeighborIds[2] =  6;
      l_cellInformation[3].faceNeighborIds[3] =  1;
      
      l_cellInformation[4].faceNeighborIds[0] =  -1;
      l_cellInformation[4].faceNeighborIds[1] =  0;
      l_cellInformation[4].faceNeighborIds[2] =  7;
      l_cellInformation[4].faceNeighborIds[3] =  1;
      
      l_cellInformation[5].faceNeighborIds[0] =  6;
      l_cellInformation[5].faceNeighborIds[1] =  0;
      l_cellInformation[5].faceNeighborIds[2] =  -1;
      l_cellInformation[5].faceNeighborIds[3] =  1;
      
      l_cellInformation[6].faceNeighborIds[0] =  3;
      l_cellInformation[6].faceNeighborIds[1] =  2;
      l_cellInformation[6].faceNeighborIds[2] =  5;
      l_cellInformation[6].faceNeighborIds[3] =  8;
      
      l_cellInformation[7].faceNeighborIds[0] =  2;
      l_cellInformation[7].faceNeighborIds[1] =  0;
      l_cellInformation[7].faceNeighborIds[2] =  3;
      l_cellInformation[7].faceNeighborIds[3] =  4;
      
      l_cellInformation[8].faceNeighborIds[0] =  -1;
      l_cellInformation[8].faceNeighborIds[1] =  2;
      l_cellInformation[8].faceNeighborIds[2] =  0;
      l_cellInformation[8].faceNeighborIds[3] =  6;

      // add the setup to a time manager and compute pointers to face-neigbors
      seissol::time_stepping::TimeManager l_timeManager;
      l_timeManager.addCluster( 9,
                                0.9,
                                l_cellInformation,
                                l_cellData,
                                l_degreesOfFreedom );

      // check pointer structure to face neighbors
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[0][0] == &(l_degreesOfFreedom.timeIntegrated[   5][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[0][1] == &(l_degreesOfFreedom.timeIntegrated[   4][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[0][2] == &(l_degreesOfFreedom.timeIntegrated[   2][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[0][3] == &(l_degreesOfFreedom.timeIntegrated[   7][0]) );

      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[1][0] == &(l_degreesOfFreedom.timeIntegrated[   4][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[1][1] == &(l_degreesOfFreedom.timeIntegrated[   5][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[1][2] == &(l_degreesOfFreedom.timeIntegratedMpi[1][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[1][3] == &(l_degreesOfFreedom.timeIntegrated[   3][0]) );

      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[2][0] == &(l_degreesOfFreedom.timeIntegrated[   0][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[2][1] == &(l_degreesOfFreedom.timeIntegrated[   6][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[2][2] == NULL                                          ); // TODO: dynamic rupture
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[2][3] == &(l_degreesOfFreedom.timeIntegrated[   2][0]) );
      
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[3][0] == &(l_degreesOfFreedom.timeIntegratedMpi[0][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[3][1] == &(l_degreesOfFreedom.timeIntegrated[   7][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[3][2] == &(l_degreesOfFreedom.timeIntegrated[   6][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[3][3] == &(l_degreesOfFreedom.timeIntegrated[   1][0]) );
      
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[4][0] == NULL                                          );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[4][1] == &(l_degreesOfFreedom.timeIntegrated[   0][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[4][2] == &(l_degreesOfFreedom.timeIntegrated[   7][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[4][3] == &(l_degreesOfFreedom.timeIntegrated[   1][0]) );

      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[5][0] == &(l_degreesOfFreedom.timeIntegrated[   6][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[5][1] == &(l_degreesOfFreedom.timeIntegrated[   0][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[5][2] == &(l_degreesOfFreedom.timeIntegrated[   5][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[5][3] == &(l_degreesOfFreedom.timeIntegrated[   1][0]) );
      
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[6][0] == &(l_degreesOfFreedom.timeIntegrated[   3][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[6][1] == &(l_degreesOfFreedom.timeIntegrated[   2][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[6][2] == &(l_degreesOfFreedom.timeIntegrated[   5][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[6][3] == &(l_degreesOfFreedom.timeIntegrated[   8][0]) );
      
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[7][0] == NULL                                          ); // TODO: dynamic rupture
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[7][1] == &(l_degreesOfFreedom.timeIntegrated[   0][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[7][2] == &(l_degreesOfFreedom.timeIntegrated[   3][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[7][3] == &(l_degreesOfFreedom.timeIntegrated[   4][0]) );

      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[8][0] == NULL                                          );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[8][1] == &(l_degreesOfFreedom.timeIntegratedMpi[2][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[8][2] == &(l_degreesOfFreedom.timeIntegratedMpi[0][0]) );
      TS_ASSERT( l_degreesOfFreedom.faceNeighbors[8][3] == &(l_degreesOfFreedom.timeIntegrated[   6][0]) );

      // free memory
      delete[] l_cellInformation;
      delete[] l_cellData;
      delete[] l_degreesOfFreedom.regular;
      delete[] l_degreesOfFreedom.timeIntegrated;
      delete[] l_degreesOfFreedom.timeIntegratedMpi;
    }

    /**
     * Tests time, volume and boundary integration.
     **/
     void testIntegrations() {
       /*
        * Initialization: Dofs, time integrated dofs, stiffness, star matrices
        */
       real l_maximumTimeStepWidth = 0.1;

       // initialize flux solvers to random values
       for( int l_cell = 0; l_cell < m_numberOfCells; l_cell++ ) {
         m_denseMatrix.setRandomValues( 4*NUMBEROFVARIABLESSQUARED, &m_cellData[l_cell].nApNm1[0][0] );
         m_denseMatrix.setRandomValues( 4*NUMBEROFVARIABLESSQUARED, &m_cellData[l_cell].nAmNm1[0][0] );
       }

       // set LTS integration to GTS time step width
       real l_deltaT[2];
       l_deltaT[0] = 0;
       l_deltaT[1] = l_maximumTimeStepWidth;

       // setup time manager
       seissol::time_stepping::TimeManager l_timeManager;
       l_timeManager.addCluster( m_numberOfCells,
                                 l_maximumTimeStepWidth,
                                 m_cellInformation,
                                 m_cellData,
                                 m_degreesOfFreedom );
       l_timeManager.setSynchronizationPoint( l_maximumTimeStepWidth );

       /*
        * Compute integrations.
        */
       // compute time and volume integration for verification
       for( int l_cell = 0; l_cell < m_numberOfCells; l_cell++ ) {
         // time integration
         m_simpleTimeIntegrator.computeTimeDerivation(  m_degreesOfFreedomUnitTests                + (NUMBEROFUNKNOWNS                             *l_cell),
                                                        m_cellDataDense[l_cell].aStar,                                                        
                                                        m_cellDataDense[l_cell].bStar,                                                        
                                                        m_cellDataDense[l_cell].cStar,                                                        
                         (real (*)[NUMBEROFUNKNOWNS]) ( m_degreesOfFreedomTimeDerivativesUnitTests + (ORDEROFTAYLORSERIESEXPANSION*NUMBEROFUNKNOWNS*l_cell) ) );

         m_simpleTimeIntegrator.computeTimeIntegration(
                         (real (*)[NUMBEROFUNKNOWNS]) ( m_degreesOfFreedomTimeDerivativesUnitTests + (ORDEROFTAYLORSERIESEXPANSION*NUMBEROFUNKNOWNS*l_cell) ),
                                                        l_deltaT,
                                                        m_degreesOfFreedomTimeIntegratedUnitTests  + (NUMBEROFUNKNOWNS                             *l_cell) );
         // volume integration
         m_simpleVolumeIntegrator.computeVolumeIntegration( m_degreesOfFreedomTimeIntegratedUnitTests  + (NUMBEROFUNKNOWNS*l_cell),
                       (real (*)[NUMBEROFVARIABLESSQUARED]) m_cellDataDense[l_cell].aStar,
                                                            m_degreesOfFreedomUnitTests                + (NUMBEROFUNKNOWNS*l_cell) );
         
       }

       // compute boundary integration for verification
       for( int l_cell = 0; l_cell < m_numberOfCells; l_cell++ ) {
         // collect neighboring time integrated dofs
         real l_neighboringDegreesOfFreedomTimeIntegrated[4][NUMBEROFUNKNOWNS];
         for( int l_face = 0; l_face < 4; l_face ++ ) {
           int l_faceNeighborId = m_cellInformation[l_cell].faceNeighborIds[l_face];
           for( int l_value = 0; l_value < NUMBEROFUNKNOWNS; l_value++ ) { 
             l_neighboringDegreesOfFreedomTimeIntegrated[l_face][l_value] = m_degreesOfFreedomTimeIntegratedUnitTests[(NUMBEROFUNKNOWNS*l_faceNeighborId) + l_value];
           }
         }

         m_simpleBoundaryIntegrator.computeBoundaryIntegration( m_degreesOfFreedomTimeIntegratedUnitTests  + (NUMBEROFUNKNOWNS*l_cell),
                                                                l_neighboringDegreesOfFreedomTimeIntegrated,
                                                         (int*) m_cellInformation[l_cell].faceTypes,
                                                                m_cellInformation[l_cell].faceRelations,
                                                                m_cellData[l_cell].nApNm1,
                                                                m_cellData[l_cell].nAmNm1,
                                                                m_degreesOfFreedomUnitTests + (NUMBEROFUNKNOWNS*l_cell)
                                                              );
       }

       // compute integrations using the time manager
       l_timeManager.advanceInTime();

       /*
        * Check time integration and volume+boundary integration
        */
       m_denseMatrix.checkResult( NUMBEROFUNKNOWNS*m_numberOfCells, m_degreesOfFreedomTimeIntegratedUnitTests, ( (real*) m_degreesOfFreedom.timeIntegrated ) ); 
       m_denseMatrix.checkResult( NUMBEROFUNKNOWNS*m_numberOfCells, m_degreesOfFreedomUnitTests, &m_degreesOfFreedom.regular[0][0] ); 
        
    }

};
