/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
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
 */

#include "Parallel/MPI.h"

#include <cxxtest/TestSuite.h>

#include "Geometry/GambitReader.h"
#include "Geometry/NetcdfReader.h"
#include "Geometry/MeshReaderFBinding.h"

#include <algorithm>
#include <iostream>

#define CUBE_SIZE 4

std::ostream& operator<<(std::ostream& o, const Fault& f)
{
	o << f.element << ' ' << f.side << ' ' << f.neighborElement << ' '  << f.neighborSide;
	return o;
}

class TestGambitReader : public CxxTest::TestSuite
{

private:
	/**
	 * Checks whether two vertices (one new and one old) are the same
	 */
	bool isSameVertex(const VrtxCoords &coords1, const double* coords2)
	{
		for (int i = 0; i < 3; i++) {
			if (coords1[i] != coords2[i])
				return false;
		}

		return true;
	}

	/**
	 * Tests a mesh reader for rank 0
	 */
	void testMeshReader0(MeshReader &meshReader, bool checkFault = false)
	{

		int nVertices;
		double* verticesXY;
		getverticesxy(&nVertices, &verticesXY);

		const std::vector<Element>& elementsNew = meshReader.getElements();
		const std::vector<Vertex>& verticesNew = meshReader.getVertices();

		int nElements;
		int* elementVertices;
		getelementvertices(&nElements, &elementVertices);

		const std::map<int, MPINeighbor>& mpiNeighborsNew = meshReader.getMPINeighbors();

		TS_ASSERT_EQUALS(nVertices, verticesNew.size());
		TS_ASSERT_EQUALS(nElements, elementsNew.size());

		// Check coordinates
		// Does not work, because they are ordered different
		//for (int i = 0; i < nVertices; i++) {
		//	for (int j = 0; j < 3; j++)
		//		TS_ASSERT_DELTA(verticesXY[i*3+j], verticesNew[i].coords[j], 1e-20);
		//}

		// Check coordinates
		for (int i = 0; i < nElements; i++) {
			for (int j = 0; j < 4; j++) {
				for (int k = 0; k < 3; k++)
					TS_ASSERT_DELTA(verticesXY[(elementVertices[i*4+j]-1)*3+k],
							verticesNew[elementsNew[i].vertices[j]].coords[k], 1e-20);
			}
		}

		// Check mesh reconstruction
		int size;

		int* vrtxnelements;
		getverticesnelements(&size, &vrtxnelements);

		int* vrtxelements;
		getverticeselements(&size, &vrtxelements);

		int* reference;
		getreference(&size, &reference);

		int* mpireference;
		getmpireference(&size, &mpireference);

		int* mpinumber;
		getmpinumber(&size, &mpinumber);

		int* boundarytoobject;
		getboundarytoobject(&size, &boundarytoobject);

		int* sideneighbor;
		getsideneighbor(&size, &sideneighbor);

		int* localneighborside;
		getlocalneighborside(&size, &localneighborside);

		int* localneighborvrtx;
		getlocalneighborvrtx(&size, &localneighborvrtx);

		for (int i = 0; i < nVertices; i++) {
			for (int j = 0; j < nVertices; j++) {
				if (isSameVertex(verticesNew[i].coords, &verticesXY[j*3])) {
					TS_ASSERT_EQUALS(verticesNew[i].elements.size(), vrtxnelements[j]);

					for (unsigned int k = 0; k < verticesNew[i].elements.size(); k++) {
						TS_ASSERT_EQUALS(verticesNew[i].elements[k], vrtxelements[j+k*nVertices]-1);
					}
				}
			}
		}

		// References -> boundary condition, MPI References -> neighbor, sideneighbor -> element
		//for (int i = 0; i < size; i++) {
		//	for (int j = 0; j < 4; j++) {
		//		if (mpireferences[i*5+j+1] == 1 || references[i*5+j+1] == 5) {
		//			TS_ASSERT_EQUALS(sideneighbor[i*4+j], CUBE_SIZE*CUBE_SIZE*CUBE_SIZE*5/2+1);
		//		} else {
		//			TS_ASSERT_LESS_THAN(sideneighbor[i*4+j], CUBE_SIZE*CUBE_SIZE*CUBE_SIZE*5/2+1);
		//		}
		//	}
		//}

		for (int i = 0; i < nElements; i++) {
			for (int j = 0; j < 4; j++) {
				if (mpireference[i*5+j+1] == 1) {
					TS_ASSERT_EQUALS(elementsNew[i].neighborRanks[j], 1);
					TS_ASSERT_EQUALS(mpiNeighborsNew.at(elementsNew[i].neighborRanks[j]).localID, boundarytoobject[i*4+j]-1)

					// Does not work because the ordering inside one element may be different
					//TS_ASSERT_EQUALS(elementsNew[i].mpiIndices[j], mpinumber[i*4+j]-1);
					int k;
					for (k = 0; k < 4; k++) {
						if (elementsNew[i].mpiIndices[j] == mpinumber[i*4+k]-1)
							break;
					}
					TS_ASSERT_LESS_THAN(k, 4);
				} else {
					TS_ASSERT_EQUALS(elementsNew[i].neighborRanks[j], 0);

					if (elementsNew[i].boundaries[j] == 3) {
						TS_ASSERT_EQUALS(elementsNew[elementsNew[i].neighbors[j]].boundaries[elementsNew[i].neighborSides[j]], 3);
					}
				}

				TS_ASSERT_EQUALS(elementsNew[i].neighbors[j], sideneighbor[i*4+j]-1);
				TS_ASSERT_EQUALS(elementsNew[i].boundaries[j], reference[i*5+j+1]);

				if (reference[i*5+j+1] == 0) {
					TS_ASSERT_EQUALS(elementsNew[i].neighborSides[j], localneighborside[i*4+j]-1);
					TS_ASSERT_EQUALS(elementsNew[i].sideOrientations[j], localneighborvrtx[i*4+j]-1);
				}
			}
		}

		// Check MPI boundary structures
		int mpisize;

		getbndsize(&mpisize);
		TS_ASSERT_EQUALS(mpiNeighborsNew.size(), mpisize);

		int rank;
		int bndnelem;
		int* domainelements;
		for (int i = 1; i <= mpisize; i++) {
			getbndrank(i, &rank);

			getbndnelem(i, &bndnelem);
			TS_ASSERT_EQUALS(mpiNeighborsNew.at(rank).elements.size(), bndnelem);

			getbnddomainelements(i, &size, &domainelements);
			for (int j = 0; j < size; j++) {
				TS_ASSERT_EQUALS(mpiNeighborsNew.at(rank).elements[j].localElement, domainelements[j]-1);
			}
		}

		if (checkFault) {
			VrtxCoords center;
                        int refPointMethod;
			getfaultreferencepoint(&center[0], &center[1], &center[2], &refPointMethod);
			meshReader.findFault(center, refPointMethod);

			int* mpinumberdr;
			getmpinumberdr(&size, &mpinumberdr);

			for (int i = 0; i < nElements; i++) {
				for (int j = 0; j < 4; j++) {
					if (elementsNew[i].mpiFaultIndices[j] >= 0) {
						// Does not work because the ordering inside one element may be different
						//TS_ASSERT_EQUALS(elementsNew[i].mpiFaultIndices[j], mpifaultnumber[i*4+j]-1);
						int k;
						for (k = 0; k < 4; k++) {
							if (elementsNew[i].mpiFaultIndices[j] == mpinumberdr[i*4+k]-1)
								break;
						}
						TS_ASSERT_LESS_THAN(k, 4);
					}
				}
			}

			const std::vector<Fault>& faultNew = meshReader.getFault();

			int* faultface;
			getfaultface(&size, &faultface);

			double* faultnormals;
			getfaultnormals(&size, &faultnormals);

			double* faulttangent1;
			getfaulttangent1(&size, &faulttangent1);

			double* faulttangent2;
			getfaulttangent2(&size, &faulttangent2);

			TS_ASSERT_EQUALS(faultNew.size(), size);

			for (unsigned int i = 0; i < faultNew.size(); i++) {
				unsigned int j;
				for (j = 0; j < faultNew.size(); j++) {
					if (faultNew[i].element == faultface[j]-1 && faultNew[i].neighborElement == faultface[j + size*2]-1) {
						//TS_ASSERT_EQUALS(faultNew[i].element, faultface[i]-1);
						TS_ASSERT_EQUALS(faultNew[i].side, faultface[j + size]-1);
						//TS_ASSERT_EQUALS(faultNew[i].neighborElement, faultface[i + size*2]-1);
						TS_ASSERT_EQUALS(faultNew[i].neighborSide, faultface[j + size*3]-1);

						TS_ASSERT_DELTA(faultNew[i].normal[0], faultnormals[j*3], 10e-15);
						TS_ASSERT_DELTA(faultNew[i].normal[1], faultnormals[j*3+1], 10e-15);
						TS_ASSERT_DELTA(faultNew[i].normal[2], faultnormals[j*3+2], 10e-15);

						TS_ASSERT_DELTA(faultNew[i].tangent1[0], faulttangent1[j*3], 10e-15);
						TS_ASSERT_DELTA(faultNew[i].tangent1[1], faulttangent1[j*3+1], 10e-15);
						TS_ASSERT_DELTA(faultNew[i].tangent1[2], faulttangent1[j*3+2], 10e-15);

						TS_ASSERT_DELTA(faultNew[i].tangent2[0], faulttangent2[j*3], 10e-15);
						TS_ASSERT_DELTA(faultNew[i].tangent2[1], faulttangent2[j*3+1], 10e-15);
						TS_ASSERT_DELTA(faultNew[i].tangent2[2], faulttangent2[j*3+2], 10e-15);
						break;
					}
				}
				if (j == faultNew.size())
					TS_FAIL("No matching fault face found");
			}
		}
	}

	/**
	 * Tests MPI boundaries for a mesh on rank 1
	 */
	void testMeshReader1(MeshReader &meshReader)
	{
		int size;

		const std::map<int, MPINeighbor>& mpiNeighborsNew = meshReader.getMPINeighbors();

		int mpisize;
		getbndsize(&mpisize);
		TS_ASSERT_EQUALS(mpiNeighborsNew.size(), mpisize);

		int rank;
		int bndnelem;
		int* domainelements;
		for (int i = 1; i <= mpisize; i++) {
			getbndrank(i, &rank);

			getbndnelem(i, &bndnelem);
			TS_ASSERT_EQUALS(mpiNeighborsNew.at(rank).elements.size(), bndnelem);

			getbnddomainelements(i, &size, &domainelements);
			for (int j = 0; j < size; j++) {
				TS_ASSERT_EQUALS(mpiNeighborsNew.at(rank).elements[j].localElement, domainelements[j]-1);
			}
		}
	}

public:
	void testGambitReader()
	{
		GambitReader meshReader(seissol::MPI::mpi.rank(),
				SEISSOL_TESTS "Geometry/cube4.neu", SEISSOL_TESTS "Geometry/cube4.met.epart.2");

		if (seissol::MPI::mpi.rank() == 0)
			testMeshReader0(meshReader);
		else
			testMeshReader1(meshReader);
	}

	void testNetcdfReader()
	{
#ifdef USE_NETCDF
		NetcdfReader meshReader(seissol::MPI::mpi.rank(), 2, SEISSOL_TESTS "Geometry/cube4.nc");

		if (seissol::MPI::mpi.rank() == 0)
			testMeshReader0(meshReader);
		else
			testMeshReader1(meshReader);
#endif // USE_NETCDF
	}

	void testDynamicRupture()
	{
		GambitReader meshReader0(0, SEISSOL_TESTS "Geometry/scec.neu", SEISSOL_TESTS "Geometry/scec.met.epart.2");
		//GambitReader meshReader1(1, "src/tests/Geometry/scec.neu", "src/tests/Geometry/scec.met.epart.2");

		testMeshReader0(meshReader0, true);
		//readdrboxold(1);
		// Does not work for this scenario because one element can have multiple neighbors in the same domain
		//testMeshReader1(meshReader1);
	}

};
