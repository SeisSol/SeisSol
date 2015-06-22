/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2015, SeisSol Group
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
 * Test suite, which tests the common functions of the LTS setup.
 **/

#include <Initializer/time_stepping/common.hpp>

namespace unit_tests {
  class commonTestSuite;
}

class unit_tests::commonTestSuite: public CxxTest::TestSuite {
  public:
    void testNormalize() {
      /*
       * Set up a simple clustering requiring normalization.
       * cell:      old / new cell id
       * TC:        associated time cluster
       * FT0 - FT3: face types
       * FN0 - FN3: face neighbors; *x* marks a critical face neighbor
       *  ___________________________________________________________
       * | cell | TC | FT0 | FT1 | FT2 | FT3 | FN0 | FN1 | FN2 | FN3 |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       * |  0   | 2  | reg | reg | reg | reg | *1* |  3  |  4  |  2  |
       * |  1   | 5  | DR  | out | reg | per |  3  |  -  |  0  |  2  |
       * |  2   | 2  | reg | reg | per | reg |  0  |  4  | *1* |  3  |
       * |  3   | 1  | DR  | DR  | reg | reg |  2  | *4* | *1* |  0  |
       * |  4   | 3  | reg | reg | reg | reg |  2  |  1  |  3  |  0  |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       * |               -----> normalization ------>                |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       * |  0   | 2  | reg | reg | reg | reg |  1  |  3  |  4  |  2  |
       * |  1   | 2  | DR  | out | reg | per |  3  |  -  |  0  |  2  |
       * |  2   | 2  | reg | reg | per | reg |  0  |  4  |  1  |  3  |
       * |  3   | 1  | DR  | DR  | reg | reg |  2  |  4  |  1  |  0  |
       * |  4   | 2  | reg | reg | reg | reg |  2  |  1  |  3  |  0  |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       */
       CellLocalInformation l_cellLocalInformation[5];

       // set cluster ids
       l_cellLocalInformation[0].clusterId = 2;
       l_cellLocalInformation[1].clusterId = 5;
       l_cellLocalInformation[2].clusterId = 2;
       l_cellLocalInformation[3].clusterId = 1;
       l_cellLocalInformation[4].clusterId = 3;

       // first face face types
       l_cellLocalInformation[0].faceTypes[0] = regular;
       l_cellLocalInformation[1].faceTypes[0] = dynamicRupture;
       l_cellLocalInformation[2].faceTypes[0] = regular;
       l_cellLocalInformation[3].faceTypes[0] = dynamicRupture;
       l_cellLocalInformation[4].faceTypes[0] = regular;

       // second face face types
       l_cellLocalInformation[0].faceTypes[1] = regular;
       l_cellLocalInformation[1].faceTypes[1] = outflow;
       l_cellLocalInformation[2].faceTypes[1] = regular;
       l_cellLocalInformation[3].faceTypes[1] = dynamicRupture;
       l_cellLocalInformation[4].faceTypes[1] = regular;

       // third face face types
       l_cellLocalInformation[0].faceTypes[2] = regular;
       l_cellLocalInformation[1].faceTypes[2] = regular;
       l_cellLocalInformation[2].faceTypes[2] = periodic;
       l_cellLocalInformation[3].faceTypes[2] = regular;
       l_cellLocalInformation[4].faceTypes[2] = regular;

       // fourth face face types
       l_cellLocalInformation[0].faceTypes[3] = regular;
       l_cellLocalInformation[1].faceTypes[3] = periodic;
       l_cellLocalInformation[2].faceTypes[3] = regular;
       l_cellLocalInformation[3].faceTypes[3] = regular;
       l_cellLocalInformation[4].faceTypes[3] = regular;


       // first face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[0] =  1;
       l_cellLocalInformation[1].faceNeighborIds[0] =  3;
       l_cellLocalInformation[2].faceNeighborIds[0] =  0;
       l_cellLocalInformation[3].faceNeighborIds[0] =  2;
       l_cellLocalInformation[4].faceNeighborIds[0] =  2;

       // second face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[1] =  3;
       l_cellLocalInformation[1].faceNeighborIds[1] = -1;
       l_cellLocalInformation[2].faceNeighborIds[1] =  4;
       l_cellLocalInformation[3].faceNeighborIds[1] =  4;
       l_cellLocalInformation[4].faceNeighborIds[1] =  1;

       // third face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[2] =  4;
       l_cellLocalInformation[1].faceNeighborIds[2] =  0;
       l_cellLocalInformation[2].faceNeighborIds[2] =  1;
       l_cellLocalInformation[3].faceNeighborIds[2] =  1;
       l_cellLocalInformation[4].faceNeighborIds[2] =  3;

       // fourth face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[3] =  2;
       l_cellLocalInformation[1].faceNeighborIds[3] =  2;
       l_cellLocalInformation[2].faceNeighborIds[3] =  3;
       l_cellLocalInformation[3].faceNeighborIds[3] =  0;
       l_cellLocalInformation[4].faceNeighborIds[3] =  0;


       // do the normalization
       unsigned int l_numberOfNormalizations = seissol::initializers::time_stepping::normalize( 4,
                                                                                                l_cellLocalInformation );

       // check the normalized cluster ids
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].clusterId, 2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].clusterId, 2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].clusterId, 2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].clusterId, 1 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].clusterId, 2 );
    }

    void testReorder() {
      /*
       * Set up a simple mesh requiring reordering:
       * cell:      old / new cell id
       * TC:        associated time cluster
       * FT0 - FT3: face types
       * FN0 - FN3: face neighbors
       *  ___________________________________________________________
       * | cell | TC | FT0 | FT1 | FT2 | FT3 | FN0 | FN1 | FN2 | FN3 |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       * |  0   | 2  | reg | reg | reg | reg |  1  |  3  |  4  |  2  |
       * |  1   | 5  | DR  | out | reg | per |  3  |  -  |  0  |  2  |
       * |  2   | 2  | reg | reg | per | reg |  0  |  4  |  1  |  3  |
       * |  3   | 1  | ff  | DR  | out | reg |  -  |  4  |  -  |  0  |
       * |  4   | 3  | reg | reg | reg | reg |  2  |  1  |  3  |  0  |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       * |               -----> reordering ------>                   |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       * |  3   | 1  | ff  | DR  | out | reg |  -  |  3  |  -  |  2  |
       * |  2   | 2  | reg | reg | per | reg |  2  |  3  |  4  |  0  |
       * |  0   | 2  | reg | reg | reg | reg |  4  |  0  |  3  |  1  |
       * |  4   | 3  | reg | reg | reg | reg |  1  |  4  |  0  |  2  |
       * |  1   | 5  | DR  | out | reg | per |  0  |  -  |  2  |  1  |
       * |------|----|-----|-----|-----|-----|-----|-----|-----|-----|
       */
       CellLocalInformation l_cellLocalInformation[5];

       // set cluster ids
       l_cellLocalInformation[0].clusterId = 2;
       l_cellLocalInformation[1].clusterId = 5;
       l_cellLocalInformation[2].clusterId = 2;
       l_cellLocalInformation[3].clusterId = 1;
       l_cellLocalInformation[4].clusterId = 3;

       // first face face types
       l_cellLocalInformation[0].faceTypes[0] = regular;
       l_cellLocalInformation[1].faceTypes[0] = dynamicRupture;
       l_cellLocalInformation[2].faceTypes[0] = regular;
       l_cellLocalInformation[3].faceTypes[0] = freeSurface;
       l_cellLocalInformation[4].faceTypes[0] = regular;

       // second face face types
       l_cellLocalInformation[0].faceTypes[1] = regular;
       l_cellLocalInformation[1].faceTypes[1] = outflow;
       l_cellLocalInformation[2].faceTypes[1] = regular;
       l_cellLocalInformation[3].faceTypes[1] = dynamicRupture;
       l_cellLocalInformation[4].faceTypes[1] = regular;

       // third face face types
       l_cellLocalInformation[0].faceTypes[2] = regular;
       l_cellLocalInformation[1].faceTypes[2] = regular;
       l_cellLocalInformation[2].faceTypes[2] = periodic;
       l_cellLocalInformation[3].faceTypes[2] = outflow;
       l_cellLocalInformation[4].faceTypes[2] = regular;

       // fourth face face types
       l_cellLocalInformation[0].faceTypes[3] = regular;
       l_cellLocalInformation[1].faceTypes[3] = periodic;
       l_cellLocalInformation[2].faceTypes[3] = regular;
       l_cellLocalInformation[3].faceTypes[3] = regular;
       l_cellLocalInformation[4].faceTypes[3] = regular;


       // first face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[0] =  1;
       l_cellLocalInformation[1].faceNeighborIds[0] =  3;
       l_cellLocalInformation[2].faceNeighborIds[0] =  0;
       l_cellLocalInformation[3].faceNeighborIds[0] = -1;
       l_cellLocalInformation[4].faceNeighborIds[0] =  2;

       // second face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[1] =  3;
       l_cellLocalInformation[1].faceNeighborIds[1] = -1;
       l_cellLocalInformation[2].faceNeighborIds[1] =  4;
       l_cellLocalInformation[3].faceNeighborIds[1] =  4;
       l_cellLocalInformation[4].faceNeighborIds[1] =  1;

       // third face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[2] =  4;
       l_cellLocalInformation[1].faceNeighborIds[2] =  0;
       l_cellLocalInformation[2].faceNeighborIds[2] =  1;
       l_cellLocalInformation[3].faceNeighborIds[2] = -1;
       l_cellLocalInformation[4].faceNeighborIds[2] =  3;

       // fourth face face neighbors
       l_cellLocalInformation[0].faceNeighborIds[3] =  2;
       l_cellLocalInformation[1].faceNeighborIds[3] =  2;
       l_cellLocalInformation[2].faceNeighborIds[3] =  3;
       l_cellLocalInformation[3].faceNeighborIds[3] =  0;
       l_cellLocalInformation[4].faceNeighborIds[3] =  0;

       // local cluster
       unsigned int l_localClusterIds[4];
       l_localClusterIds[0] = 1; l_localClusterIds[1] = 2;
       l_localClusterIds[2] = 3; l_localClusterIds[3] = 5;

       // forward and backward mapping of the mesh ids to lts ids
       unsigned int l_meshToLts[5];
       unsigned int l_ltsToMesh[5];

       // do the reordering
       seissol::initializers::time_stepping::reorder( 4,
                                                      l_localClusterIds,
                                                      5,
                                                      l_cellLocalInformation,
                                                      l_meshToLts,
                                                      l_ltsToMesh );

       // check the mappings
       TS_ASSERT_EQUALS( l_meshToLts[0], 2 );
       TS_ASSERT_EQUALS( l_meshToLts[1], 4 );
       TS_ASSERT_EQUALS( l_meshToLts[2], 1 );
       TS_ASSERT_EQUALS( l_meshToLts[3], 0 );
       TS_ASSERT_EQUALS( l_meshToLts[4], 3 );

       TS_ASSERT_EQUALS( l_ltsToMesh[0], 3 );
       TS_ASSERT_EQUALS( l_ltsToMesh[1], 2 );
       TS_ASSERT_EQUALS( l_ltsToMesh[2], 0 );
       TS_ASSERT_EQUALS( l_ltsToMesh[3], 4 );
       TS_ASSERT_EQUALS( l_ltsToMesh[4], 1 );

       // check the time clusters
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].clusterId, 1 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].clusterId, 2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].clusterId, 2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].clusterId, 3 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].clusterId, 5 );

       // check the face types of face 1
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceTypes[0], freeSurface    );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceTypes[0], regular        );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceTypes[0], regular        );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceTypes[0], regular        );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceTypes[0], dynamicRupture );

       // check the face types of face 2
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceTypes[1], dynamicRupture );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceTypes[1], regular        );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceTypes[1], regular        );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceTypes[1], regular        );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceTypes[1], outflow        );

       // check the face types of face 3
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceTypes[2], outflow  );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceTypes[2], periodic );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceTypes[2], regular  );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceTypes[2], regular  );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceTypes[2], regular  );

       // check the face types of face 4
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceTypes[3], regular  );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceTypes[3], regular  );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceTypes[3], regular  );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceTypes[3], regular  );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceTypes[3], periodic );

       // check the face neighbor ids of face 1
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceNeighborIds[0], -1 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceNeighborIds[0],  2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceNeighborIds[0],  4 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceNeighborIds[0],  1 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceNeighborIds[0],  0 );

       // check the face neighbor ids of face 2
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceNeighborIds[1],  3 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceNeighborIds[1],  3 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceNeighborIds[1],  0 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceNeighborIds[1],  4 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceNeighborIds[1], -1 );

       // check the face neighbor ids of face 3
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceNeighborIds[2], -1 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceNeighborIds[2],  4 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceNeighborIds[2],  3 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceNeighborIds[2],  0 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceNeighborIds[2],  2 );

       // check the face neighbor ids of face 4
       TS_ASSERT_EQUALS( l_cellLocalInformation[0].faceNeighborIds[3],  2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[1].faceNeighborIds[3],  0 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[2].faceNeighborIds[3],  1 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[3].faceNeighborIds[3],  2 );
       TS_ASSERT_EQUALS( l_cellLocalInformation[4].faceNeighborIds[3],  1 );
    }
};
