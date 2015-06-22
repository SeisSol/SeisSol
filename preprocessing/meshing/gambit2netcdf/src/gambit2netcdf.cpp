/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 *
 * @section LICENSE
 * Copyright (c) 2013, SeisSol Group
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
 */

#include <mpi.h>

#include <algorithm>
#include <cassert>
#include <cstring>
#include <fstream>
#include <vector>

#include <netcdf.h>
#include <netcdf_par.h>

#include "mesh/GambitReader.h"
#include "utils/logger.h"

void checkNcError(int error)
{
	if (error != NC_NOERR)
		logError() << "Error while writing netCDF file:" << nc_strerror(error);
}

int main(int argc, char *argv[]) {
	// Initialize MPI
	MPI_Init(&argc, &argv);

	int rank, procs;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &procs);

	// Parse cmd line arguments
	// TODO Use a more advanced tool here
	if (argc != 4)
		logError() << "Usage:" << argv[0] << "<mesh file> <partition file> <output netcdf>";

	// Parse the partition once on rank 0 to get the number of partitions
	int partInfo[2]; // Partition information (#partitions, max partition size)
	if (rank == 0) {
		std::ifstream partitionFile(argv[2]);
		if (!partitionFile)
			logError() << "Could not open partition file" << argv[2];

		int p;
		std::vector<int> partitionSizes;

		while (!partitionFile.eof()) {
			partitionFile >> p;

			int oldSize = partitionSizes.size();
			if (p+1 >= partitionSizes.size())
				partitionSizes.resize(p+1);
			for (int i = oldSize; i < partitionSizes.size(); i++)
				partitionSizes[i] = 0;

			partitionSizes[p]++;

			// Skip white spaces
			partitionFile >> std::ws;
		}

		logInfo() << "Found" << partitionSizes.size() << "partitions in" << argv[2];

		partInfo[0] = partitionSizes.size();
		partInfo[1] = *std::max_element(partitionSizes.begin(), partitionSizes.end());
	}

	// Broadcast partition information
	MPI_Bcast(&partInfo, 2, MPI_INT, 0, MPI_COMM_WORLD);

	// Create netcdf file
	int ncFile;
	checkNcError(nc_create_par(argv[3], NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncFile));

	// Create netcdf dimensions
	int ncDimDimension;
	nc_def_dim(ncFile, "dimension", 3, &ncDimDimension);

	int ncDimPart;
	nc_def_dim(ncFile, "partitions", partInfo[0], &ncDimPart);

	int ncDimElem, ncDimElemSides, ncDimElemVertices;
	nc_def_dim(ncFile, "elements", partInfo[1], &ncDimElem);
	nc_def_dim(ncFile, "element_sides", 4, &ncDimElemSides);
	nc_def_dim(ncFile, "element_vertices", 4, &ncDimElemVertices);

	int ncDimVrtx;
	nc_def_dim(ncFile, "vertices", NC_UNLIMITED, &ncDimVrtx);

	int ncDimBnd, ncDimBndElem;
	nc_def_dim(ncFile, "boundaries", NC_UNLIMITED, &ncDimBnd);
	nc_def_dim(ncFile, "boundary_elements", NC_UNLIMITED, &ncDimBndElem);

	// Create netcdf variables
	int ncVarElemSize;
	checkNcError(nc_def_var(ncFile, "element_size", NC_INT, 1, &ncDimPart, &ncVarElemSize));
	checkNcError(nc_var_par_access(ncFile, ncVarElemSize, NC_COLLECTIVE));

	int ncVarElemVertices;
	int dimsElemVertices[] = {ncDimPart, ncDimElem, ncDimElemVertices};
	checkNcError(nc_def_var(ncFile, "element_vertices", NC_INT, 3, dimsElemVertices, &ncVarElemVertices));
	checkNcError(nc_var_par_access(ncFile, ncVarElemVertices, NC_COLLECTIVE));

	int ncVarElemNeighbors;
	int dimsElemSides[] = {ncDimPart, ncDimElem, ncDimElemSides};
	checkNcError(nc_def_var(ncFile, "element_neighbors", NC_INT, 3, dimsElemSides, &ncVarElemNeighbors));
	checkNcError(nc_var_par_access(ncFile, ncVarElemNeighbors, NC_COLLECTIVE));

	int ncVarElemBoundaries;
	checkNcError(nc_def_var(ncFile, "element_boundaries", NC_INT, 3, dimsElemSides, &ncVarElemBoundaries));
	checkNcError(nc_var_par_access(ncFile, ncVarElemBoundaries, NC_COLLECTIVE));

	int ncVarElemNeighborSides;
	checkNcError(nc_def_var(ncFile, "element_neighbor_sides", NC_INT, 3, dimsElemSides, &ncVarElemNeighborSides));
	checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborSides, NC_COLLECTIVE));

	int ncVarElemSideOrientations;
	checkNcError(nc_def_var(ncFile, "element_side_orientations", NC_INT, 3, dimsElemSides, &ncVarElemSideOrientations));
	checkNcError(nc_var_par_access(ncFile, ncVarElemSideOrientations, NC_COLLECTIVE));

	int ncVarElemNeighborRanks;
	checkNcError(nc_def_var(ncFile, "element_neighbor_ranks", NC_INT, 3, dimsElemSides, &ncVarElemNeighborRanks));
	checkNcError(nc_var_par_access(ncFile, ncVarElemNeighborRanks, NC_COLLECTIVE));

	int ncVarElemMPIIndices;
	checkNcError(nc_def_var(ncFile, "element_mpi_indices", NC_INT, 3, dimsElemSides, &ncVarElemMPIIndices));
	checkNcError(nc_var_par_access(ncFile, ncVarElemMPIIndices, NC_COLLECTIVE));

	int ncVarVrtxSize;
	checkNcError(nc_def_var(ncFile, "vertex_size", NC_INT, 1, &ncDimPart, &ncVarVrtxSize));
	checkNcError(nc_var_par_access(ncFile, ncVarVrtxSize, NC_COLLECTIVE));

	int ncVarVrtxCoords;
	int dimsVrtxCoords[] = {ncDimPart, ncDimVrtx, ncDimDimension};
	checkNcError(nc_def_var(ncFile, "vertex_coordinates", NC_DOUBLE, 3, dimsVrtxCoords, &ncVarVrtxCoords));
	checkNcError(nc_var_par_access(ncFile, ncVarVrtxCoords, NC_COLLECTIVE));

	int ncVarBndSize;
	checkNcError(nc_def_var(ncFile, "boundary_size", NC_INT, 1, &ncDimPart, &ncVarBndSize));
	checkNcError(nc_var_par_access(ncFile, ncVarBndSize, NC_COLLECTIVE));

	int ncVarBndElemSize;
	int dimsBndElemSize[] = {ncDimPart, ncDimBnd};
	checkNcError(nc_def_var(ncFile, "boundary_element_size", NC_INT, 2, dimsBndElemSize, &ncVarBndElemSize));
	checkNcError(nc_var_par_access(ncFile, ncVarBndElemSize, NC_COLLECTIVE));

	int ncVarBndElemRank;
	checkNcError(nc_def_var(ncFile, "boundary_element_rank", NC_INT, 2, dimsBndElemSize, &ncVarBndElemRank));
	checkNcError(nc_var_par_access(ncFile, ncVarBndElemRank, NC_COLLECTIVE));

	int ncVarBndElemLocalIds;
	int dimsBndElemLocalIds[] = {ncDimPart, ncDimBnd, ncDimBndElem};
	checkNcError(nc_def_var(ncFile, "boundary_element_localids", NC_INT, 3, dimsBndElemLocalIds, &ncVarBndElemLocalIds));
	checkNcError(nc_var_par_access(ncFile, ncVarBndElemLocalIds, NC_COLLECTIVE));

	// Start writing data
	checkNcError(nc_enddef(ncFile));

	// Variable buffers
	ElemVertices* elemVertices = new ElemVertices[partInfo[1]];
	ElemNeighbors* elemNeighbors = new ElemNeighbors[partInfo[1]];
	ElemNeighborSides* elemNeighborSides = new ElemNeighborSides[partInfo[1]];
	ElemSideOrientations* elemSideOrientations = new ElemSideOrientations[partInfo[1]];
	ElemBoundaries* elemBoundaries = new ElemBoundaries[partInfo[1]];
	ElemNeighborRanks* elemNeighborRanks = new ElemNeighborRanks[partInfo[1]];
	ElemMPIIndices* elemMPIIndices = new ElemMPIIndices[partInfo[1]];

	int i;
	for (i = rank; i < partInfo[0]; i += procs) {
		logInfo() << "Start mesh reconstruction for partition" << i;
		GambitReader gambitReader(i, argv[1], argv[2]);

		logInfo() << "Start writing mesh information for partition" << i;

		// Elements
		const std::vector<Element> elements = gambitReader.getElements();
		assert(elements.size() <= partInfo[1]);

		// Copy element values into buffers
		for (int j = 0; j < elements.size(); j++) {
			memcpy(&elemVertices[j], elements[j].vertices, sizeof(ElemVertices));
			memcpy(&elemNeighbors[j], elements[j].neighbors, sizeof(ElemNeighbors));
			memcpy(&elemNeighborSides[j], elements[j].neighborSides, sizeof(ElemNeighborSides));
			memcpy(&elemSideOrientations[j], elements[j].sideOrientations, sizeof(ElemSideOrientations));
			memcpy(&elemBoundaries[j], elements[j].boundaries, sizeof(ElemBoundaries));
			memcpy(&elemNeighborRanks[j], elements[j].neighborRanks, sizeof(ElemNeighborRanks));
			memcpy(&elemMPIIndices[j], elements[j].mpiIndices, sizeof(ElemMPIIndices));
		}

		// Write element buffers to netcdf
		size_t start[3] = {i, 0, 0};
		int size = elements.size();
		checkNcError(nc_put_var1_int(ncFile, ncVarElemSize, start, &size));
		size_t count[3] = {1, size, 4};
		checkNcError(nc_put_vara_int(ncFile, ncVarElemVertices, start, count, reinterpret_cast<int*>(elemVertices)));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighbors, start, count, reinterpret_cast<int*>(elemNeighbors)));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborSides, start, count, reinterpret_cast<int*>(elemNeighborSides)));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemSideOrientations, start, count, reinterpret_cast<int*>(elemSideOrientations)));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemBoundaries, start, count, reinterpret_cast<int*>(elemBoundaries)));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborRanks, start, count, reinterpret_cast<int*>(elemNeighborRanks)));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemMPIIndices, start, count, reinterpret_cast<int*>(elemMPIIndices)));

		// Vertices
		const std::vector<Vertex> vertices = gambitReader.getVertices();
		VrtxCoords* vrtxCoords = new VrtxCoords[vertices.size()];

		// Copy vertex values into buffers
		for (int j = 0; j < vertices.size(); j++) {
			memcpy(&vrtxCoords[j], vertices[j].coords, sizeof(VrtxCoords));
		}

		// Write vertex buffer to netcdf
		size = vertices.size();
		checkNcError(nc_put_var1_int(ncFile, ncVarVrtxSize, start, &size));
		count[0] = 1; count[1] = size; count[2] = 3;
		checkNcError(nc_put_vara_double(ncFile, ncVarVrtxCoords, start, count, reinterpret_cast<double*>(vrtxCoords)));

		delete [] vrtxCoords;

		// Boundaries (MPI neighbors)
		const std::map<int, MPINeighbor> mpiNeighbors = gambitReader.getMPINeighbors();

		// Get maximum number of neighbors (required to get collective MPI-IO right)
		int maxNeighbors = mpiNeighbors.size();
		MPI_Allreduce(MPI_IN_PLACE, &maxNeighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		// Write number of boundaries to netcdf
		size = mpiNeighbors.size();
		checkNcError(nc_put_var1_int(ncFile, ncVarBndSize, start, &size));

		int j = 0;
		for (std::map<int, MPINeighbor>::const_iterator iter = mpiNeighbors.begin();
			iter != mpiNeighbors.end(); iter++, j++) {
			size_t bndStart[3] = {i, j, 0};

			// Write size of this boundary to netcdf
			int elemSize = iter->second.elements.size();
			checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, bndStart, &elemSize));

			// Write neighbor rank to netcdf
			int bndRank = iter->first;
			checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, bndStart, &bndRank));

			// Max sure rank resizes the netCDF dimension
			// -> everybody has to write maxBndElements
			int maxBndElements = iter->second.elements.size();
			MPI_Allreduce(MPI_IN_PLACE, &maxBndElements, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			// Copy local element ids to buffer
			int* bndElemLocalIds = new int[maxBndElements];
			for (int k = 0; k < iter->second.elements.size(); k++) {
				bndElemLocalIds[k] = iter->second.elements[k].localElement;
			}

			// Write local element ids to netcdf
			size_t bndCount[3] = {1, 1, maxBndElements};
			checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, bndStart, bndCount, bndElemLocalIds));

			delete [] bndElemLocalIds;
		}

		for (; j < maxNeighbors; j++) {
			size_t bndStart[3] = {i, j, 0};

			// Write some dummy values so netCDF run into a deadlock
			int bndSize = 0;
			checkNcError(nc_put_var1_int(ncFile, ncVarBndElemSize, bndStart, &bndSize));
			int bndRank = 0;
			checkNcError(nc_put_var1_int(ncFile, ncVarBndElemRank, bndStart, &bndRank));

			// -> everybody has to write maxBndElements
			int maxBndElements = 0;
			MPI_Allreduce(MPI_IN_PLACE, &maxBndElements, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			int* bndElemLocalIds = new int[maxBndElements];

			size_t bndCount[3] = {1, 1, maxBndElements};
			checkNcError(nc_put_vara(ncFile, ncVarBndElemLocalIds, bndStart, bndCount, bndElemLocalIds));

			delete [] bndElemLocalIds;
		}
	}

	// Some processors may idle during the last iteration
	if (i - rank < partInfo[0] && i >= partInfo[0]) {
		size_t start[3] = {0, 0, 0};
		size_t count[3] = {0, 0, 0};
		checkNcError(nc_put_vara_int(ncFile, ncVarElemSize, start, count, 0L));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemVertices, start, count, 0L));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighbors, start, count, 0L));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborSides, start, count, 0L));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemSideOrientations, start, count, 0L));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemBoundaries, start, count, 0L));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemNeighborRanks, start, count, 0L));
		checkNcError(nc_put_vara_int(ncFile, ncVarElemMPIIndices, start, count, 0L));

		checkNcError(nc_put_vara_int(ncFile, ncVarVrtxSize, start, count, 0L));
		checkNcError(nc_put_vara_double(ncFile, ncVarVrtxCoords, start, count, 0L));

		int maxNeighbors = 0;
		MPI_Allreduce(MPI_IN_PLACE, &maxNeighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		checkNcError(nc_put_vara_int(ncFile, ncVarBndSize, start, count, 0L));
		for (int j = 0; j < maxNeighbors; j++) {
			checkNcError(nc_put_vara_int(ncFile, ncVarBndElemSize, start, count, 0L));
			checkNcError(nc_put_vara_int(ncFile, ncVarBndElemRank, start, count, 0L));

			int maxBndElements = 0;
			MPI_Allreduce(MPI_IN_PLACE, &maxBndElements, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

			checkNcError(nc_put_vara_int(ncFile, ncVarBndElemLocalIds, start, count, 0L));
		}
	}

	// Free buffers
	delete [] elemVertices;
	delete [] elemNeighbors;
	delete [] elemNeighborSides;
	delete [] elemSideOrientations;
	delete [] elemBoundaries;
	delete [] elemNeighborRanks;
	delete [] elemMPIIndices;

	// Close netcdf file
	checkNcError(nc_close(ncFile));

	// Cleanup MPI
	MPI_Finalize();
}
