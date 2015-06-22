/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2014-2015, SeisSol Group
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
 * IMPLIED WARRANTIES OF  MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE  USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 *
 * @section DESCRIPTION
 * Test suite, which tests the setup of the internal state.
 **/
#include <cxxtest/TestSuite.h>
#include <cstdlib>
#include <ctime>
#include <limits>
#include <utils/logger.h>

// overwrite const qualifiers for testing
#define const 
#include <Initializer/typedefs.hpp>
#undef const

#include <Initializer/InternalState.h>
#include <Initializer/time_stepping/common.hpp>

namespace unit_tests {
  class InternalStateTestSuite;
}

class unit_tests::InternalStateTestSuite: public CxxTest::TestSuite {
  public:
    void setUp() {
      srand( time(NULL) );
    }

    /**
     * Tests the derivation of communication layer layouts in global time stepping.
     **/
    void testGtsCommunicationLayerLayout() {
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        /*
         * set up communication structure to test again
         */
        struct CommunicationStructure l_communicationStructure;
        struct CellLocalInformation *l_cellLocalInformation;

        l_communicationStructure.numberOfRegions    = rand() % 1000;
        l_communicationStructure.numberOfGhostCells = 0;
        l_communicationStructure.numberOfGhostRegionCells = (unsigned int*) malloc(l_communicationStructure.numberOfRegions*sizeof(unsigned int) );
        l_communicationStructure.numberOfCopyCells  = 0;
        l_communicationStructure.numberOfCopyRegionCells = (unsigned int*) malloc(l_communicationStructure.numberOfRegions*sizeof(unsigned int) );

        // iterate over regions and set the number ghost and copy region cells
        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          l_communicationStructure.numberOfGhostRegionCells[l_region] = rand() % 10000;
          l_communicationStructure.numberOfGhostCells += l_communicationStructure.numberOfGhostRegionCells[l_region];         

          l_communicationStructure.numberOfCopyRegionCells[l_region] = rand() % 10000;
          l_communicationStructure.numberOfCopyCells += l_communicationStructure.numberOfCopyRegionCells[l_region];
        }

        /*
         * set up cell GTS-non DR cell information
         */
        l_cellLocalInformation = (struct CellLocalInformation*) malloc( ( l_communicationStructure.numberOfGhostCells +
                                                                          l_communicationStructure.numberOfCopyCells ) * sizeof( struct CellLocalInformation )
                                                                      );
        // iterate of ghost and copy layer and set face information
        for( unsigned int l_cell = 0; l_cell < l_communicationStructure.numberOfGhostCells + l_communicationStructure.numberOfCopyCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
            // set faces pointing to neighboring ranks
            if( l_cell < l_communicationStructure.numberOfGhostCells && rand()%3 == 0 ) {
              l_cellLocalInformation[l_cell].faceNeighborIds[l_face] = std::numeric_limits<unsigned int>::max();
            }
            else {
              l_cellLocalInformation[l_cell].faceNeighborIds[l_face] = rand()%10000;
            }

            // init face types
            unsigned int l_faceType = rand() % 10;

            if( l_faceType == 0 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = freeSurface;
            }
            else if( l_faceType == 1 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = outflow;
            }
            else if( l_faceType == 2 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = periodic;
            }
            else {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
            }
          }

          unsigned int l_mpiFace = rand()%3;
          // enforce a mpiface
          l_cellLocalInformation[l_cell].faceNeighborIds[l_mpiFace] = std::numeric_limits<unsigned int>::max();

          // set at least one true face neighbor
          l_cellLocalInformation[l_cell].faceNeighborIds[(rand()+l_mpiFace+1)%3] = rand()%10000;

          // set gts
          l_cellLocalInformation[l_cell].ltsSetup = (1 << 4);
        }

        /*
         * Check that the InternalState produces #(ghost cells) and #(copy cells) time buffers and 0 time derivatives
         */
        unsigned int *l_numberOfBuffers     = (unsigned int*) malloc( l_communicationStructure.numberOfRegions * sizeof(unsigned int) );
        unsigned int *l_numberOfDerivatives = (unsigned int*) malloc( l_communicationStructure.numberOfRegions * sizeof(unsigned int) );

        seissol::initializers::InternalState::deriveLayerLayout(  Layer::ghost,
                                                                  1,
                                                                 &l_communicationStructure.numberOfRegions,
                                                                 &l_communicationStructure.numberOfGhostRegionCells,
                                                                  l_cellLocalInformation,
                                                                 &l_numberOfBuffers,
                                                                 &l_numberOfDerivatives );

        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          TS_ASSERT_EQUALS( l_numberOfBuffers[l_region],     l_communicationStructure.numberOfGhostRegionCells[l_region] );
          TS_ASSERT_EQUALS( l_numberOfDerivatives[l_region], 0                                                           );
        }

        seissol::initializers::InternalState::deriveLayerLayout(  Layer::copy,
                                                                  1,
                                                                 &l_communicationStructure.numberOfRegions,
                                                                 &l_communicationStructure.numberOfCopyRegionCells,
                                                                 &l_cellLocalInformation[l_communicationStructure.numberOfGhostCells],
                                                                 &l_numberOfBuffers,
                                                                 &l_numberOfDerivatives );

        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          TS_ASSERT_EQUALS( l_numberOfBuffers[l_region],     l_communicationStructure.numberOfCopyRegionCells[l_region] );
          TS_ASSERT_EQUALS( l_numberOfDerivatives[l_region], 0                                                          );
        }

        /*
         * set up cell GTS with DR cell information
         */
        unsigned int *l_numberOfBuffersUT[2];     // 0: ghost, 1: copy
        unsigned int *l_numberOfDerivativesUT[2]; // 0: ghost, 1: copy
        for( unsigned int l_i = 0; l_i < 2; l_i++ ) {
          l_numberOfBuffersUT[l_i]     = (unsigned int*) malloc( l_communicationStructure.numberOfRegions * sizeof(unsigned int) );
          l_numberOfDerivativesUT[l_i] = (unsigned int*) malloc( l_communicationStructure.numberOfRegions * sizeof(unsigned int) );
        }

        // iterate over ghost layer and set face information
        unsigned int l_firstGhostRegionCell = 0;
        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          l_numberOfBuffersUT[0][l_region]     = 0; l_numberOfDerivativesUT[0][l_region] = 0;

          for( unsigned int l_ghost = 0; l_ghost < l_communicationStructure.numberOfGhostRegionCells[l_region]; l_ghost++ ) {
            unsigned int l_cell = l_firstGhostRegionCell+l_ghost;
            bool l_dynamicRupture = false;

            for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
              unsigned int l_faceType = rand() % 10;

              if( l_faceType == 0 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = freeSurface;
              }
              else if( l_faceType == 1 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = outflow;
              }
              else if( l_faceType == 2 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = periodic;
              }
              else if( l_faceType == 3 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = dynamicRupture;
                if( l_cellLocalInformation[l_cell].faceNeighborIds[l_face] != std::numeric_limits<unsigned int>::max() ) l_dynamicRupture = true;
              }
              else {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
              }
            }
            if( !l_dynamicRupture ) l_numberOfBuffersUT[0][l_region]++;
            else {
              l_cellLocalInformation[l_cell].ltsSetup ^= (1 << 4);
              l_cellLocalInformation[l_cell].ltsSetup |= (1 << 5);
              l_numberOfDerivativesUT[0][l_region]++;
            }
          }

          l_firstGhostRegionCell += l_communicationStructure.numberOfGhostRegionCells[l_region];
        }

        // iterate over copy layer and set face information
        unsigned int l_firstCopyRegionCell = l_communicationStructure.numberOfGhostCells;
        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          l_numberOfBuffersUT[1][l_region] = 0; l_numberOfDerivativesUT[1][l_region] = 0;

          for( unsigned int l_copy = 0; l_copy < l_communicationStructure.numberOfCopyRegionCells[l_region]; l_copy++ ) {
            unsigned int l_cell = l_firstCopyRegionCell+l_copy;
            unsigned int l_dynamicRupture = 0;

            for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
              unsigned int l_faceType = rand() % 10;

              if( l_faceType == 0 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = freeSurface;
              }
              else if( l_faceType == 1 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = outflow;
              }
              else if( l_faceType == 2 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = periodic;
              }
              else if( l_faceType == 3 ) {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = dynamicRupture;
                l_dynamicRupture++;
                l_cellLocalInformation[l_cell].ltsSetup |= (1 << 5);
              }
              else {
                l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
              }
            }
            if( l_dynamicRupture > 0  ) l_numberOfDerivativesUT[1][l_region]++;
            if( l_dynamicRupture != 4 ) l_numberOfBuffersUT[1][l_region]++;

            // no buffer for DR only
            if( l_dynamicRupture == 4 ) l_cellLocalInformation[l_cell].ltsSetup ^= (1 << 4);
          }

          l_firstCopyRegionCell += l_communicationStructure.numberOfCopyRegionCells[l_region];
        }

        /*
         * check that the InternalState produces the right values:
         * ghost layer:
         *   #derivatives = #(cells who share dynamic rupture faces with the copy layer)
         *   #(time buffers) = #(ghost layer cells) - #(derivatives)
         * copy layer:
         *   #derivatives = #(cells with dynamic rupture faces)
         *   #(time buffers) = #(copy layer cells)       
         */
        seissol::initializers::InternalState::deriveLayerLayout(  Layer::ghost,
                                                                  1,
                                                                 &l_communicationStructure.numberOfRegions,
                                                                 &l_communicationStructure.numberOfGhostRegionCells,
                                                                  l_cellLocalInformation,
                                                                 &l_numberOfBuffers,
                                                                 &l_numberOfDerivatives );

        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          TS_ASSERT_EQUALS( l_numberOfBuffers[l_region],     l_numberOfBuffersUT[0][l_region]     );
          TS_ASSERT_EQUALS( l_numberOfDerivatives[l_region], l_numberOfDerivativesUT[0][l_region] );
        }

        seissol::initializers::InternalState::deriveLayerLayout(  Layer::copy,
                                                                  1,
                                                                 &l_communicationStructure.numberOfRegions,
                                                                 &l_communicationStructure.numberOfCopyRegionCells,
                                                                 &l_cellLocalInformation[l_communicationStructure.numberOfGhostCells],
                                                                 &l_numberOfBuffers,
                                                                 &l_numberOfDerivatives );

        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          TS_ASSERT_EQUALS( l_numberOfBuffers[l_region],     l_numberOfBuffersUT[1][l_region]     );
          TS_ASSERT_EQUALS( l_numberOfDerivatives[l_region], l_numberOfDerivativesUT[1][l_region] );
        }

        /*
         * free memory
         */
        free( l_numberOfBuffersUT[0] );
        free( l_numberOfDerivativesUT[0] );
        free( l_numberOfBuffersUT[1] );
        free( l_numberOfDerivativesUT[1] );

        free( l_numberOfBuffers );
        free( l_numberOfDerivatives );

        free( (void*) l_communicationStructure.numberOfGhostRegionCells );
        free( (void*) l_communicationStructure.numberOfCopyRegionCells );
        free(         l_cellLocalInformation );
      }
    }

    /**
     * Test the derivation of the interior layout in global time stepping.
     **/
    void testGtsInteriorLayout() {
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        /*
         * set up random communication and interior structure.
         */
        struct CommunicationStructure l_communicationStructure;
        if( l_repeat % 4 == 0 ) {
          l_communicationStructure.numberOfGhostCells = 0;
          l_communicationStructure.numberOfCopyCells = 0;
        }
        else { 
          l_communicationStructure.numberOfGhostCells = rand() % 1000;
          l_communicationStructure.numberOfCopyCells  = rand() % 1000;
        }

        unsigned int l_numberOfInteriorCells;
        if( l_repeat < 10 ) l_numberOfInteriorCells = l_repeat;
        else                l_numberOfInteriorCells = rand() % 100000;


        unsigned int l_numberOfCells = l_communicationStructure.numberOfGhostCells
                                     + l_communicationStructure.numberOfCopyCells
                                     + l_numberOfInteriorCells;

        struct CellLocalInformation* l_cellLocalInformation = (CellLocalInformation*) malloc( l_numberOfCells * sizeof(CellLocalInformation) );

        /*
         * set up non-dr cell information
         */
        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
            // init face types
            unsigned int l_faceType = rand() % 10;

            if( l_faceType == 0 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = freeSurface;
            }
            else if( l_faceType == 1 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = outflow;
            }
            else if( l_faceType == 2 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = periodic;
            }
            else {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
            }

            // set face neighbor id
            l_cellLocalInformation[l_cell].faceNeighborIds[l_face] = l_communicationStructure.numberOfGhostCells;
          }

          // set gts
          l_cellLocalInformation[l_cell].ltsSetup = (1 << 4);
        }

        /*
         * check non-dr
         */
        unsigned int l_numberOfBuffers, l_numberOfDerivatives;

        unsigned int l_numberOfCommunicationCells = l_communicationStructure.numberOfGhostCells + l_communicationStructure.numberOfCopyCells;

        seissol::initializers::InternalState::deriveInteriorLayout(  1,
                                                                    &l_numberOfInteriorCells,
                                                                    &l_cellLocalInformation[l_numberOfCommunicationCells],
                                                                    &l_numberOfBuffers,
                                                                    &l_numberOfDerivatives );

        TS_ASSERT_EQUALS( l_numberOfBuffers,     l_numberOfInteriorCells );
        TS_ASSERT_EQUALS( l_numberOfDerivatives, 0                       );

        /*
         * set up dr cell information
         */
        unsigned int l_numberOfBuffersUT = 0;
        unsigned int l_numberOfDerivativesUT = 0;

        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          unsigned int l_dynamicRupture = 0;

          for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
            // init face types
            unsigned int l_faceType = rand() % 10;

            if( l_faceType == 0 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = freeSurface;
            }
            else if( l_faceType == 1 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = outflow;
            }
            else if( l_faceType == 2 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = periodic;
            }
            else if( l_faceType == 3) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = dynamicRupture;
              // adopt lts information
              l_cellLocalInformation[l_cell].ltsSetup |= (1 << l_face);
              l_cellLocalInformation[l_cell].ltsSetup |= (1 << 5);
              l_dynamicRupture++;
            }
            else {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
            }
          }

          if( l_cell >= l_numberOfCells - l_numberOfInteriorCells) {
            if( l_dynamicRupture > 0 )  l_numberOfDerivativesUT++;
            if( l_dynamicRupture != 4 ) l_numberOfBuffersUT++;
          }
          if( l_dynamicRupture == 4 ) l_cellLocalInformation[l_cell].ltsSetup ^= (1 << 4);
        }

        /*
         * check dr
         */
        seissol::initializers::InternalState::deriveInteriorLayout(  1,
                                                                    &l_numberOfInteriorCells,
                                                                    &l_cellLocalInformation[l_numberOfCommunicationCells],
                                                                    &l_numberOfBuffers,
                                                                    &l_numberOfDerivatives );

        TS_ASSERT_EQUALS( l_numberOfBuffers,     l_numberOfBuffersUT     );
        TS_ASSERT_EQUALS( l_numberOfDerivatives, l_numberOfDerivativesUT );


        free( l_cellLocalInformation );
      }
    }

    /**
     * Test that the pointers in the GTS case are set correctly.
     **/
    void testGtsGhostPointers() {
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        /*
         * set up communication structure to test again
         */
        struct CommunicationStructure l_communicationStructure;

        l_communicationStructure.numberOfRegions    = rand() % 1000;
        l_communicationStructure.numberOfGhostCells = 0;
        l_communicationStructure.numberOfGhostRegionCells = new unsigned int [l_communicationStructure.numberOfRegions];

        // iterate over regions and set the number ghost and copy region cells
        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          l_communicationStructure.numberOfGhostRegionCells[l_region] = rand() % 10000;
          l_communicationStructure.numberOfGhostCells += l_communicationStructure.numberOfGhostRegionCells[l_region];
        }

        /*
         * setup cell local non-dr information
         */
        struct CellLocalInformation* l_cellLocalInformation = new CellLocalInformation[ l_communicationStructure.numberOfGhostCells ];
        for( unsigned int l_cell = 0; l_cell < l_communicationStructure.numberOfGhostCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
          l_cellLocalInformation[l_cell].ltsSetup = (1 << 4);
        }

        real*  l_ghostLayer =  new real [ 1 ]; // actual memory is never used
        real** l_buffers     = new real*[ l_communicationStructure.numberOfGhostCells ];
        real** l_derivatives = new real*[ l_communicationStructure.numberOfGhostCells ];

        unsigned int *l_numberOfBuffers     = new unsigned int[l_communicationStructure.numberOfRegions];
        unsigned int *l_numberOfDerivatives = new unsigned int[l_communicationStructure.numberOfRegions];

        seissol::initializers::InternalState::deriveLayerLayout(  Layer::ghost,
                                                                  1,
                                                                 &l_communicationStructure.numberOfRegions,
                                                                 &l_communicationStructure.numberOfGhostRegionCells,
                                                                  l_cellLocalInformation,
                                                                 &l_numberOfBuffers,
                                                                 &l_numberOfDerivatives );

        /*
         * check the non-dr pointers
         */
        seissol::initializers::InternalState::setUpLayerPointers( Layer::ghost,
                                                                  l_communicationStructure.numberOfRegions,
                                                                  l_communicationStructure.numberOfGhostRegionCells,
                                                                  l_cellLocalInformation,
                                                                  l_numberOfBuffers,
                                                                  l_numberOfDerivatives,
                                                                  l_ghostLayer,
                                                                  l_buffers,
                                                                  l_derivatives );

        for( unsigned int l_ghost = 0; l_ghost < l_communicationStructure.numberOfGhostCells; l_ghost++ ) {
          TS_ASSERT_EQUALS( l_buffers[l_ghost],     l_ghostLayer + (l_ghost*NUMBER_OF_ALIGNED_DOFS) );
          TS_ASSERT_EQUALS( l_derivatives[l_ghost], (real*)0                                        );
        }

        /*
         * setup cell local dr information
         */
        for( unsigned int l_cell = 0; l_cell < l_communicationStructure.numberOfGhostCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
            if( rand()%4 == 0 ) l_cellLocalInformation[l_cell].faceTypes[l_face] = dynamicRupture;
            else l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;

            if( rand()%3 == 0 ) l_cellLocalInformation[l_cell].faceNeighborIds[l_face] = std::numeric_limits<unsigned int>::max();
            else l_cellLocalInformation[l_cell].faceNeighborIds[l_face] = rand()%10000;
          }
        }

          seissol::initializers::InternalState::deriveLayerLayout(  Layer::ghost,
                                                                    1,
                                                                   &l_communicationStructure.numberOfRegions,
                                                                   &l_communicationStructure.numberOfGhostRegionCells,
                                                                    l_cellLocalInformation,
                                                                   &l_numberOfBuffers,
                                                                   &l_numberOfDerivatives );
        /*
         * check the dr pointers
         */
        seissol::initializers::InternalState::setUpLayerPointers(  Layer::ghost,
                                                                   l_communicationStructure.numberOfRegions,
                                                                   l_communicationStructure.numberOfGhostRegionCells,
                                                                   l_cellLocalInformation,
                                                                   l_numberOfBuffers,
                                                                   l_numberOfDerivatives,
                                                                   l_ghostLayer,
                                                                   l_buffers,
                                                                   l_derivatives );
        // limits for memory
        real* l_limits[3]; // 0: first time integrated position, 1: first derivatives position, 2: first non-derivatives position
        l_limits[0] = l_limits[1] = l_limits[2] = l_ghostLayer;

        unsigned int l_firstGhostRegionCell = 0;
        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          unsigned int l_firstNonGhostRegionCell = l_firstGhostRegionCell + l_communicationStructure.numberOfGhostRegionCells[l_region];

          unsigned int l_numberOfBuffersUT = 0;
          unsigned int l_numberOfDerivativesUT = 0;

          l_limits[0] = l_limits[2];
          l_limits[1] = l_limits[0] + l_numberOfBuffers[l_region]     * NUMBER_OF_ALIGNED_DOFS;
          l_limits[2] = l_limits[1] + l_numberOfDerivatives[l_region] * NUMBER_OF_ALIGNED_DERS;

          for( unsigned int l_ghost = l_firstGhostRegionCell; l_ghost < l_firstNonGhostRegionCell;  l_ghost++ ) {
            if( l_buffers[l_ghost] != NULL ) {
              l_numberOfBuffersUT++;

              // assert that the buffer is in bounds
              TS_ASSERT_LESS_THAN_EQUALS( l_limits[0],        l_buffers[l_ghost] );
              TS_ASSERT_LESS_THAN(        l_buffers[l_ghost], l_limits[1]        );
            }
            if( l_derivatives[l_ghost] != NULL ) {
              l_numberOfDerivativesUT++;

              // assert that the derivative is in bounds
              TS_ASSERT_LESS_THAN(        l_derivatives[l_ghost] , l_limits[2]            );
              TS_ASSERT_LESS_THAN_EQUALS( l_limits[1],             l_derivatives[l_ghost] );
            }
          }

          l_firstGhostRegionCell = l_firstNonGhostRegionCell;

          // check the number of non-NULL pointers
          TS_ASSERT_EQUALS( l_numberOfBuffers[l_region],     l_numberOfBuffersUT     );
          TS_ASSERT_EQUALS( l_numberOfDerivatives[l_region], l_numberOfDerivativesUT );
        }

        /*
         * free memory
         */
        delete[] l_numberOfBuffers;
        delete[] l_numberOfDerivatives;

        delete[] l_ghostLayer;

        delete[] l_derivatives;
        delete[] l_buffers;

        delete[] l_cellLocalInformation;
        delete[] l_communicationStructure.numberOfGhostRegionCells;
      }
    }

    void testGtsCopyPointers() {
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        /*
         * set up communication structure to test again
         */
        struct CommunicationStructure l_communicationStructure;

        l_communicationStructure.numberOfRegions    = rand() % 1000;
        l_communicationStructure.numberOfCopyCells = 0;
        l_communicationStructure.numberOfCopyRegionCells = new unsigned int [l_communicationStructure.numberOfRegions];

        // iterate over regions and set the number copy region cells
        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          l_communicationStructure.numberOfCopyRegionCells[l_region] = rand() % 10000;
          l_communicationStructure.numberOfCopyCells += l_communicationStructure.numberOfCopyRegionCells[l_region];
        }

        /*
         * setup cell local non-dr information
         */
        struct CellLocalInformation* l_cellLocalInformation = new CellLocalInformation[ l_communicationStructure.numberOfCopyCells ];
        for( unsigned int l_cell = 0; l_cell < l_communicationStructure.numberOfCopyCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
          l_cellLocalInformation[l_cell].ltsSetup = (1 << 4);
        }

        real*  l_copyLayer =  new real [ 1 ]; // actual memory is never used
        real** l_buffers     = new real*[ l_communicationStructure.numberOfCopyCells ];
        real** l_derivatives = new real*[ l_communicationStructure.numberOfCopyCells ];

        unsigned int *l_numberOfBuffers     = new unsigned int[l_communicationStructure.numberOfRegions];
        unsigned int *l_numberOfDerivatives = new unsigned int[l_communicationStructure.numberOfRegions];

        seissol::initializers::InternalState::deriveLayerLayout(  Layer::copy,
                                                                  1,
                                                                 &l_communicationStructure.numberOfRegions,
                                                                 &l_communicationStructure.numberOfCopyRegionCells,
                                                                  l_cellLocalInformation,
                                                                 &l_numberOfBuffers,
                                                                 &l_numberOfDerivatives );

        /*
         * check the non-dr pointers
         */
        seissol::initializers::InternalState::setUpLayerPointers( Layer::copy,
                                                                  l_communicationStructure.numberOfRegions,
                                                                  l_communicationStructure.numberOfCopyRegionCells,
                                                                  l_cellLocalInformation,
                                                                  l_numberOfBuffers,
                                                                  l_numberOfDerivatives,
                                                                  l_copyLayer,
                                                                  l_buffers,
                                                                  l_derivatives );

        for( unsigned int l_copy = 0; l_copy < l_communicationStructure.numberOfCopyCells; l_copy++ ) {
          TS_ASSERT_EQUALS( l_buffers[l_copy],     l_copyLayer + (l_copy*NUMBER_OF_ALIGNED_DOFS) );
          TS_ASSERT_EQUALS( l_derivatives[l_copy], (real*)0                                        );
        }

        /*
         * setup cell local dr information
         */
        for( unsigned int l_cell = 0; l_cell < l_communicationStructure.numberOfCopyCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
            if( rand()%4 == 0 ) l_cellLocalInformation[l_cell].faceTypes[l_face] = dynamicRupture;
            else l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;

            l_cellLocalInformation[l_cell].faceNeighborIds[l_face] = rand()%10000;
          }
        }

        seissol::initializers::InternalState::deriveLayerLayout(  Layer::copy,
                                                                  1,
                                                                 &l_communicationStructure.numberOfRegions,
                                                                 &l_communicationStructure.numberOfCopyRegionCells,
                                                                  l_cellLocalInformation,
                                                                 &l_numberOfBuffers,
                                                                 &l_numberOfDerivatives );
        /*
         * check the dr pointers
         */
        seissol::initializers::InternalState::setUpLayerPointers(  Layer::copy,
                                                                   l_communicationStructure.numberOfRegions,
                                                                   l_communicationStructure.numberOfCopyRegionCells,
                                                                   l_cellLocalInformation,
                                                                   l_numberOfBuffers,
                                                                   l_numberOfDerivatives,
                                                                   l_copyLayer,
                                                                   l_buffers,
                                                                   l_derivatives );
        // limits for memory
        real* l_limits[3]; // 0: first time integrated position, 1: first derivatives position, 2: first non-derivatives position
        l_limits[0] = l_limits[1] = l_limits[2] = l_copyLayer;

        unsigned int l_firstCopyRegionCell = 0;
        for( unsigned int l_region = 0; l_region < l_communicationStructure.numberOfRegions; l_region++ ) {
          unsigned int l_firstNonCopyRegionCell = l_firstCopyRegionCell + l_communicationStructure.numberOfCopyRegionCells[l_region];

          unsigned int l_numberOfBuffersUT = 0;
          unsigned int l_numberOfDerivativesUT = 0;

          l_limits[0] = l_limits[2];
          l_limits[1] = l_limits[0] + l_numberOfBuffers[l_region]     * NUMBER_OF_ALIGNED_DOFS;
          l_limits[2] = l_limits[1] + l_numberOfDerivatives[l_region] * NUMBER_OF_ALIGNED_DERS;

          unsigned int l_copyRegionCell = 0;

          for( unsigned int l_copy = l_firstCopyRegionCell; l_copy < l_firstNonCopyRegionCell;  l_copy++ ) {
            // the GTS-copy layer forces time integrated DOFs to be present
            TS_ASSERT_DIFFERS( l_buffers[l_copy], (real*) 0 );
            l_numberOfBuffersUT++;

            // assert linear pointers in the time integrated DOFs
            TS_ASSERT_EQUALS( l_buffers[l_copy], l_limits[0] + (l_copyRegionCell*NUMBER_OF_ALIGNED_DOFS) );
            l_copyRegionCell++;

            // assert that the buffer is in bounds
            TS_ASSERT_LESS_THAN_EQUALS( l_limits[0],       l_buffers[l_copy] );
            TS_ASSERT_LESS_THAN(        l_buffers[l_copy], l_limits[1]       );

            if( l_derivatives[l_copy] != NULL ) {
              l_numberOfDerivativesUT++;

              // assert that the derivative is in bounds
              TS_ASSERT_LESS_THAN(        l_derivatives[l_copy] , l_limits[2]            );
              TS_ASSERT_LESS_THAN_EQUALS( l_limits[1],            l_derivatives[l_copy] );
            }
          }

          l_firstCopyRegionCell = l_firstNonCopyRegionCell;

          // check the number of non-NULL pointers
          TS_ASSERT_EQUALS( l_numberOfBuffers[l_region],     l_numberOfBuffersUT     );
          TS_ASSERT_EQUALS( l_numberOfDerivatives[l_region], l_numberOfDerivativesUT );
        }

        /*
         * free memory
         */
        delete[] l_numberOfBuffers;
        delete[] l_numberOfDerivatives;

        delete[] l_copyLayer;

        delete[] l_derivatives;
        delete[] l_buffers;

        delete[] l_cellLocalInformation;
        delete[] l_communicationStructure.numberOfCopyRegionCells;
      }
    }

    /*
     * Tests the setup of the pointers in the interior.
     */
    void testGtsInteriorPointers() {
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        unsigned int l_numberOfCells = rand()%10000;

        /*
         * set up non-dr cell information
         */
        struct CellLocalInformation* l_cellLocalInformation = new CellLocalInformation[ l_numberOfCells ];
        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;
          l_cellLocalInformation[l_cell].ltsSetup = ( 1 << 4 );
        }

        real*  l_interiorMemory = new real [ 1 ]; // actual memory is never used
        real** l_buffers        = new real*[ l_numberOfCells ];
        real** l_derivatives    = new real*[ l_numberOfCells ];

        /*
         * test GTS without DR
         */
        seissol::initializers::InternalState::setUpInteriorPointers( l_numberOfCells,
                                                                     l_cellLocalInformation,
                                                                     l_numberOfCells,
                                                                     0,
                                                                     l_interiorMemory,
                                                                     l_buffers,
                                                                     l_derivatives );
        // check the pointers
        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          TS_ASSERT_DIFFERS( l_buffers[l_cell],     (real*) NULL );
          TS_ASSERT_EQUALS(  l_derivatives[l_cell], (real*) NULL );

          TS_ASSERT_EQUALS(  l_buffers[l_cell], l_interiorMemory + l_cell * NUMBER_OF_ALIGNED_DOFS );
        }

        /*
         * set up some DR information
         */
        unsigned int l_numberOfDerivativesUT = 0;
        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          bool l_derivative = false;
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) {
            if( rand()%10 == 0 ) {
              l_cellLocalInformation[l_cell].faceTypes[l_face] = dynamicRupture;
              l_derivative = true;
            }
          }
          if( l_derivative ) l_numberOfDerivativesUT++;
        }

        /*
         * test GTS with DR
         */
        unsigned int l_numberOfBuffers, l_numberOfDerivatives;
        seissol::initializers::InternalState::deriveInteriorLayout(  1,
                                                                    &l_numberOfCells,
                                                                     l_cellLocalInformation,
                                                                    &l_numberOfBuffers,
                                                                    &l_numberOfDerivatives );

        seissol::initializers::InternalState::setUpInteriorPointers( l_numberOfCells,
                                                                     l_cellLocalInformation,
                                                                     l_numberOfBuffers,
                                                                     l_numberOfDerivatives,
                                                                     l_interiorMemory,
                                                                     l_buffers,
                                                                     l_derivatives );

        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          TS_ASSERT_DIFFERS( l_buffers[l_cell],     (real*) NULL );
          TS_ASSERT_EQUALS(  l_buffers[l_cell], l_interiorMemory + l_cell * NUMBER_OF_ALIGNED_DOFS );

          if( l_derivatives[l_cell] != NULL ) {
            l_numberOfDerivatives--;
            TS_ASSERT_LESS_THAN( l_interiorMemory + ( l_numberOfBuffers * NUMBER_OF_ALIGNED_DOFS ) - 1, l_derivatives[l_cell] );
            TS_ASSERT_LESS_THAN( l_derivatives[l_cell], l_interiorMemory + ( l_numberOfBuffers       * NUMBER_OF_ALIGNED_DOFS )
                                                                         + ( l_numberOfDerivativesUT * NUMBER_OF_ALIGNED_DERS ) );
          }
        }
        TS_ASSERT_EQUALS( l_numberOfDerivatives, 0 );
        TS_ASSERT( l_derivatives[l_numberOfCells-1] == NULL ||
                   l_derivatives[l_numberOfCells-1] == l_interiorMemory + (  l_numberOfBuffers          * NUMBER_OF_ALIGNED_DOFS )
                                                                        + ( (l_numberOfDerivativesUT-1) * NUMBER_OF_ALIGNED_DERS ) );

        delete[] l_cellLocalInformation;
        delete[] l_buffers;
        delete[] l_derivatives;
      }
    }

    /**
     * Test the derivation of the interior layout in local time stepping.
     **/
    void testLtsInteriorLayout() {
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        /*
         * setup lts information
         */
        unsigned int  l_numberOfTimeClusters = rand()%5+1;
        unsigned int *l_numberOfClusterCells = new unsigned int[l_numberOfTimeClusters];
        unsigned int  l_numberOfCells = 0;
        for( unsigned int l_cluster = 0; l_cluster < l_numberOfTimeClusters; l_cluster++ ) {
          l_numberOfClusterCells[l_cluster] = rand()%500+1;
          l_numberOfCells += l_numberOfClusterCells[l_cluster];
        }

       CellLocalInformation *l_cellLocalInformation = new CellLocalInformation[l_numberOfCells];

       unsigned int *l_numberOfBuffers       = new unsigned int[l_numberOfTimeClusters];
       unsigned int *l_numberOfBuffersUT     = new unsigned int[l_numberOfTimeClusters];
       unsigned int *l_numberOfDerivatives   = new unsigned int[l_numberOfTimeClusters];
       unsigned int *l_numberOfDerivativesUT = new unsigned int[l_numberOfTimeClusters];

       unsigned int l_firstClusterCell = 0;

       // iterate over clusters
       for( unsigned int l_cluster = 0; l_cluster < l_numberOfTimeClusters; l_cluster++ ) {
         l_numberOfBuffersUT[l_cluster] = 0;
         l_numberOfDerivativesUT[l_cluster] = 0;

         // iterate over cells
         for( unsigned int l_cell = 0; l_cell < l_numberOfClusterCells[l_cluster]; l_cell++ ) {
           unsigned int l_ltsRandom = rand()%3;

           if( l_ltsRandom == 0 ) {
             l_cellLocalInformation[l_cell+l_firstClusterCell].ltsSetup = (1 << 4);
             l_numberOfBuffersUT[l_cluster]++;
           }
           else if( l_ltsRandom == 1 ) {
             l_cellLocalInformation[l_cell+l_firstClusterCell].ltsSetup  = (1 << 4);
             l_cellLocalInformation[l_cell+l_firstClusterCell].ltsSetup |= (1 << 5);
             l_numberOfBuffersUT[l_cluster]++;
             l_numberOfDerivativesUT[l_cluster]++;
           }
           else if( l_ltsRandom == 2 ) {
             l_cellLocalInformation[l_cell+l_firstClusterCell].ltsSetup  = (1 << 5);
             l_numberOfDerivativesUT[l_cluster]++;
           }
         }
         l_firstClusterCell += l_numberOfClusterCells[l_cluster];
       }

       /*
        * derive the interior layout
        */
       seissol::initializers::InternalState::deriveInteriorLayout(  l_numberOfTimeClusters,
                                                                    l_numberOfClusterCells,
                                                                    l_cellLocalInformation,
                                                                    l_numberOfBuffers,
                                                                    l_numberOfDerivatives );

       /*
        * check solution
        */
       for( unsigned int l_cluster = 0; l_cluster < l_numberOfTimeClusters; l_cluster++ ) {
         TS_ASSERT_EQUALS( l_numberOfBuffers[l_cluster],     l_numberOfBuffersUT[l_cluster]     );
         TS_ASSERT_EQUALS( l_numberOfDerivatives[l_cluster], l_numberOfDerivativesUT[l_cluster] );
       }

       /*
        * free memory
        */
       delete[] l_numberOfClusterCells;
       delete[] l_cellLocalInformation;
       delete[] l_numberOfBuffers;     delete[] l_numberOfBuffersUT;
       delete[] l_numberOfDerivatives; delete[] l_numberOfDerivativesUT;
      }
    }

    /*
     * Tests the derivation of the interior pointers in local times stepping.
     */
    void testLtsInteriorPointers() {
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        unsigned int l_numberOfCells = rand()%10000;
        unsigned int l_numberOfBuffers = 0;
        unsigned int l_numberOfDerivatives = 0;

        /*
         * set up cell information
         */
        struct CellLocalInformation* l_cellLocalInformation = new CellLocalInformation[ l_numberOfCells ];
        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          for( unsigned int l_face = 0; l_face < 4; l_face++ ) l_cellLocalInformation[l_cell].faceTypes[l_face] = regular;

          // random lts
          unsigned int l_ltsRandom = rand()%3;

          if( l_ltsRandom == 0 ) {
            l_cellLocalInformation[l_cell].ltsSetup  = (1 << 4);
            l_numberOfBuffers++;
          }
          else if( l_ltsRandom == 1 ) {
            l_cellLocalInformation[l_cell].ltsSetup  = (1 << 4);
            l_numberOfBuffers++;
            l_cellLocalInformation[l_cell].ltsSetup |= (1 << 5);
            l_numberOfDerivatives++;
          }
          else if( l_ltsRandom == 2 ) {
            l_cellLocalInformation[l_cell].ltsSetup  = (1 << 5);
            l_numberOfDerivatives++;
          }
        }

        real*  l_interiorMemory = new real [ 1 ]; // actual memory is never used
        real** l_buffers        = new real*[ l_numberOfCells ];
        real** l_derivatives    = new real*[ l_numberOfCells ];

        /*
         * run setup
         */
        seissol::initializers::InternalState::setUpInteriorPointers( l_numberOfCells,
                                                                     l_cellLocalInformation,
                                                                     l_numberOfBuffers,
                                                                     l_numberOfDerivatives,
                                                                     l_interiorMemory,
                                                                     l_buffers,
                                                                     l_derivatives );
        /*
         * check the pointers
         */
        real *l_buffer     = l_interiorMemory;
        real *l_derivative = l_interiorMemory + l_numberOfBuffers*NUMBER_OF_ALIGNED_DOFS;

        for( unsigned int l_cell = 0; l_cell < l_numberOfCells; l_cell++ ) {
          if( (l_cellLocalInformation[l_cell].ltsSetup >> 4) % 2 ) {
            TS_ASSERT_DIFFERS( l_buffers[l_cell],     (real*) NULL );
            TS_ASSERT_EQUALS(  l_buffers[l_cell],     l_buffer     );
            l_buffer += NUMBER_OF_ALIGNED_DOFS;
          }
          else TS_ASSERT_EQUALS( l_buffers[l_cell], (real*) NULL );

          if( (l_cellLocalInformation[l_cell].ltsSetup >> 5) % 2 ) {
            TS_ASSERT_DIFFERS( l_derivatives[l_cell], (real*) NULL );
            TS_ASSERT_EQUALS(  l_derivatives[l_cell], l_derivative );
            l_derivative += NUMBER_OF_ALIGNED_DERS;
          }
          else TS_ASSERT_EQUALS( l_derivatives[l_cell], (real*) NULL );
        }

        /*
         * free memory
         */
        delete[] l_cellLocalInformation;
        delete[] l_interiorMemory;
        delete[] l_buffers;
        delete[] l_derivatives;
      }
    }

    /*
     * Tests the setup of the interior layer.
     */
    void testSetUpInteriorLayer() {
      // repeat the test
      for( unsigned int l_repeat = 0; l_repeat < 100; l_repeat++ ) {
        /*
         * Setup the inverse problem.
         */
        unsigned int  l_numberOfClusters = rand()%25 + 1;
        unsigned int  l_numberOfUniqueCells = 0;
        unsigned int *l_clusterIds             = new unsigned int[l_numberOfClusters];
        l_clusterIds[0] = rand()%2;
        unsigned int *l_numberOfClusterCellsUT = new unsigned int[l_numberOfClusters];

        for( unsigned int l_cluster = 0; l_cluster < l_numberOfClusters; l_cluster++ ) {
          // generate cells for this cluster and add to total number of cells
          l_numberOfClusterCellsUT[l_cluster] = rand()%10000 + 1;
          l_numberOfUniqueCells += l_numberOfClusterCellsUT[l_cluster];
          if( l_cluster > 0 )
            l_clusterIds[l_cluster] = l_clusterIds[l_cluster-1] + 1 + rand()%2;
        }

        CellLocalInformation *l_cellInformationUT = new CellLocalInformation[l_numberOfUniqueCells];
        unsigned int l_firstClusterCell = 0;
        for( unsigned int l_cluster = 0; l_cluster < l_numberOfClusters; l_cluster++ ) {
          for( unsigned int l_cell = 0; l_cell < l_numberOfClusterCellsUT[l_cluster]; l_cell++ ) {
            l_cellInformationUT[l_firstClusterCell+l_cell].clusterId = l_clusterIds[l_cluster];
          }
          l_firstClusterCell += l_numberOfClusterCellsUT[l_cluster];
        }

        /*
         * Run the setup of the layer based on the mesh.
         */
        unsigned int *l_numberOfClusterCells;
        CellLocalInformation *l_cellLocalInformation;

        seissol::initializers::InternalState::setUpLayer(  interior,
                                                           l_numberOfClusters,
                                                           l_numberOfUniqueCells,
                                                           l_cellInformationUT,
                                                          &l_numberOfClusterCells,
                                                          &l_cellLocalInformation );

        /*
         * Check the results.
         */
        TS_ASSERT_EQUALS( l_cellInformationUT, l_cellLocalInformation );

        for( unsigned int l_cluster = 0; l_cluster < l_numberOfClusters; l_cluster++ ) {
          TS_ASSERT_EQUALS( l_numberOfClusterCellsUT[l_cluster], l_numberOfClusterCells[l_cluster] );
        }

        /*
         * Free memory.
         */
        delete[] l_clusterIds;
        delete[] l_numberOfClusterCellsUT;
        delete[] l_cellInformationUT;
        delete[] l_numberOfClusterCells;
      }
    }
};
