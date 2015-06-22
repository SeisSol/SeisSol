/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2012-2013, SeisSol Group
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
 * @DESCRIPTION
 * Reorder a mesh using Zoltan SFC
 */

#include <mpi.h>
#include <zoltan_cpp.h>

#include <cstring>
#include <iostream>
#include <sstream>

typedef double coord_t[3];

/**
 * Zoltan callback function
 *
 * @return The number of objects (elements)
 */
int numberOfObjects(void *data, int *ierr)
{
	*ierr = ZOLTAN_OK;

	return *static_cast<int*>(data);
}

/**
 * Zoltan callback function
 * Sets the object IDs; object IDs are numbered continuously from
 * 0 to (n-1)
 */
void getObjectList(void *data, int sizeGID, int sizeLID,
	ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
	int wgt_dim, float *obj_wgts, int *ierr)
{
	if (wgt_dim != 0) {
		std::cerr << "Zoltan needs weights, but we do not have any" << std::endl;
		*ierr = ZOLTAN_FATAL;
		return;
	}

	if (sizeGID != 1 || sizeLID != 1) {
		std::cerr << "Local or global ID is not 1; makes me helpless" << std::endl;
		*ierr = ZOLTAN_FATAL;
		return;
	}

	*ierr = ZOLTAN_OK;

	int size = *static_cast<int*>(data);

	for (int i = 0; i < size; i++) {
		globalID[i] = i;
		localID[i] = i;
	}
}

/**
 * Zoltan callback function
 *
 * @return The number of dimension. This is currently always 3.
 */
int numGeom(void *data, int *ierr)
{
	*ierr = ZOLTAN_OK;

	return 3;
}

/**
 * Zoltan callback function
 * Provides the geometry infomation required by Zoltan
 */
void geom(void *data, int sizeGID, int sizeLID,
	int numObj, ZOLTAN_ID_PTR globalID, ZOLTAN_ID_PTR localID,
	int numDim, double *geom, int *ierr)
{
	if (sizeGID != 1 || sizeLID != 1) {
		std::cerr << "Local or global ID is not 1; makes me helpless" << std::endl;
		*ierr = ZOLTAN_FATAL;
		return;
	}

	if (numDim != 3) {
		std::cerr << "This is only designed for 3 dimensions" << std::endl;
		*ierr = ZOLTAN_FATAL;
		return;
	}

	*ierr = ZOLTAN_OK;

	memcpy(geom, data, numObj*numDim*sizeof(double));
}

int main(int argc, char *argv[])
{
	MPI_Init(&argc, &argv);

	// Check number of MPI processes
	int mpiSize;
	MPI_Comm_size(MPI_COMM_WORLD, &mpiSize);
	if (mpiSize > 1) {
		std::cerr << "This program currently supports only one MPI process" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Check input parameter
	if (argc != 2) {
		std::cerr << "Number of elements not provides" << std::endl;
		std::cerr << "Usage: zoltan <number of elements>" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 1);
	}

	// Initialize Zoltan
	Zoltan *zz = new Zoltan(MPI_COMM_WORLD);

	zz->Set_Param("DEBUG_LEVEL", "0");      /* No debug output */
	zz->Set_Param("ORDER_METHOD", "LOCAL_HSFC");    /* Zoltan method */
	zz->Set_Param("NUM_GID_ENTRIES", "1");  /* global ID is 1 integer */
	zz->Set_Param("NUM_LID_ENTRIES", "1");  /* local ID is 1 integer */

	// Get number of elements
	int nElements;
	std::istringstream ss(argv[1]);
	ss >> nElements;

	zz->Set_Num_Obj_Fn(numberOfObjects, &nElements);
	zz->Set_Obj_List_Fn(getObjectList, &nElements);

	zz->Set_Num_Geom_Fn(numGeom, 0L);

	// Read the elements coordinates from stdin
	coord_t *coords = new coord_t[nElements];

	for (int i = 0; i < nElements; i++) {
		std::cin >> coords[i][0] >> coords[i][1] >> coords[i][2];
		std::cout << coords[i][0] << ' ' << coords[i][1] << ' ' << coords[i][2] << std::endl;
	}

	zz->Set_Geom_Multi_Fn(geom, coords);

	// The actual Zoltan call
	ZOLTAN_ID_PTR globalID = new ZOLTAN_ID_TYPE[nElements];
	ZOLTAN_ID_PTR permutedID = new ZOLTAN_ID_TYPE[nElements];

	for (int i = 0; i < nElements; i++)
		globalID[i] = i;

	int rc = zz->Order(1, nElements, globalID, permutedID);

	if (rc != ZOLTAN_OK) {
		std::cerr << "Could not order elements" << std::endl;
		MPI_Abort(MPI_COMM_WORLD, 2);
	}

	// Print the permutation to stdout
	for (int i = 0; i < nElements; i++)
		std::cout << globalID[i] << ' ' << permutedID[i] << std::endl;

	// Cleanup Zoltan
	delete zz;

	// Free memory
	delete [] coords;
	delete [] globalID;
	delete [] permutedID;

	MPI_Finalize();

	return 0;
}
