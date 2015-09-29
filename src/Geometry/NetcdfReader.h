/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
 * Read Gambit Mesh and Metis Partition in memory efficient way
 **/

#ifndef NETCDF_READER_H
#define NETCDF_READER_H

#ifdef USE_NETCDF

#include "MeshReader.h"

#include <cassert>
#include <cstring>

#ifdef USE_MPI
#include <netcdf_par.h>
#endif // USE_MPI
#include <netcdf.h>

#include "Initializer/preProcessorMacros.fpp"

#include "utils/logger.h"

class NetcdfReader : public MeshReader
{
public:
	NetcdfReader(int rank, int nProcs, const char* meshFile)
		: MeshReader(rank)
	{
//		EPIK_USER_REG(r_init_netcdf, "init_netcdf");
//		EPIK_USER_START(r_init_netcdf);
		// Open nc file
		int ncFile;
#ifdef USE_MPI
		checkNcError(nc_open_par(meshFile, NC_NETCDF4 | NC_MPIIO, MPI_COMM_WORLD, MPI_INFO_NULL, &ncFile));
#else // USE_MPI
		checkNcError(nc_open(meshFile, NC_NETCDF4, &ncFile));
#endif // USE_MPI

		// Get important dimensions
		int ncDimPart;
		checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
		size_t partitions;
		checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &partitions));

		if (partitions != static_cast<unsigned int>(nProcs))
			logError() << "Number of partitions in netCDF file does not match number of MPI ranks.";

		int ncDimBndSize;
		checkNcError(nc_inq_dimid(ncFile, "boundaries", &ncDimBndSize));
		size_t bndSize;
		checkNcError(nc_inq_dimlen(ncFile, ncDimBndSize, &bndSize));

		int ncDimBndElem;
		checkNcError(nc_inq_dimid(ncFile, "boundary_elements", &ncDimBndElem));
		size_t bndElemSize;
		checkNcError(nc_inq_dimlen(ncFile, ncDimBndElem, &bndElemSize));

		// Create netcdf variables
		int ncVarElemSize;
		checkNcError(nc_inq_varid(ncFile, "element_size", &ncVarElemSize));
		collectiveAccess(ncFile, ncVarElemSize);

		int ncVarElemVertices;
		checkNcError(nc_inq_varid(ncFile, "element_vertices", &ncVarElemVertices));
		collectiveAccess(ncFile, ncVarElemVertices);

		int ncVarElemNeighbors;
		checkNcError(nc_inq_varid(ncFile, "element_neighbors", &ncVarElemNeighbors));
		collectiveAccess(ncFile, ncVarElemNeighbors);

		int ncVarElemBoundaries;
		checkNcError(nc_inq_varid(ncFile, "element_boundaries", &ncVarElemBoundaries));
		collectiveAccess(ncFile, ncVarElemBoundaries);

		int ncVarElemNeighborSides;
		checkNcError(nc_inq_varid(ncFile, "element_neighbor_sides", &ncVarElemNeighborSides));
		collectiveAccess(ncFile, ncVarElemNeighborSides);

		int ncVarElemSideOrientations;
		checkNcError(nc_inq_varid(ncFile, "element_side_orientations", &ncVarElemSideOrientations));
		collectiveAccess(ncFile, ncVarElemSideOrientations);

		int ncVarElemNeighborRanks;
		checkNcError(nc_inq_varid(ncFile, "element_neighbor_ranks", &ncVarElemNeighborRanks));
		collectiveAccess(ncFile, ncVarElemNeighborRanks);

		int ncVarElemMPIIndices;
		checkNcError(nc_inq_varid(ncFile, "element_mpi_indices", &ncVarElemMPIIndices));
		collectiveAccess(ncFile, ncVarElemMPIIndices);

		int ncVarElemGroup;
		bool hasGroup = false;
		int ncResult = nc_inq_varid(ncFile, "element_group", &ncVarElemGroup);
		if (ncResult != NC_ENOTVAR) {
			checkNcError(ncResult);
			hasGroup = true;
			collectiveAccess(ncFile, ncVarElemGroup);
		}

		int ncVarVrtxSize;
		checkNcError(nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize));
		collectiveAccess(ncFile, ncVarVrtxSize);

		int ncVarVrtxCoords;
		checkNcError(nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords));
		collectiveAccess(ncFile, ncVarVrtxCoords);

		int ncVarBndSize;
		checkNcError(nc_inq_varid(ncFile, "boundary_size", &ncVarBndSize));
		collectiveAccess(ncFile, ncVarBndSize);

		int ncVarBndElemSize;
		checkNcError(nc_inq_varid(ncFile, "boundary_element_size", &ncVarBndElemSize));
		collectiveAccess(ncFile, ncVarBndElemSize);

		int ncVarBndElemRank;
		checkNcError(nc_inq_varid(ncFile, "boundary_element_rank", &ncVarBndElemRank));
		collectiveAccess(ncFile, ncVarBndElemRank);

		int ncVarBndElemLocalIds;
		checkNcError(nc_inq_varid(ncFile, "boundary_element_localids", &ncVarBndElemLocalIds));
		collectiveAccess(ncFile, ncVarBndElemLocalIds);
//		EPIK_USER_END(r_init_netcdf);

		logInfo(rank) << "Start reading mesh from netCDF file";

		// Elements
		size_t start[3] = {static_cast<size_t>(rank), 0, 0};
		int size;
		checkNcError(nc_get_var1_int(ncFile, ncVarElemSize, start, &size));

		ElemVertices* elemVertices = new ElemVertices[size];
		ElemNeighbors* elemNeighbors = new ElemNeighbors[size];
		ElemNeighborSides* elemNeighborSides = new ElemNeighborSides[size];
		ElemSideOrientations* elemSideOrientations = new ElemSideOrientations[size];
		ElemBoundaries* elemBoundaries = new ElemBoundaries[size];
		ElemNeighborRanks* elemNeighborRanks = new ElemNeighborRanks[size];
		ElemMPIIndices* elemMPIIndices = new ElemMPIIndices[size];
		ElemMaterial* elemMaterial = new ElemMaterial[size];

		m_elements.resize(size);

		EPIK_USER_REG(r_read_elements, "read_elements");
		SCOREP_USER_REGION_DEFINE( r_read_elements )
		EPIK_USER_START(r_read_elements);
		SCOREP_USER_REGION_BEGIN( r_read_elements, "read_elements", SCOREP_USER_REGION_TYPE_COMMON )
		// Read element buffers from netcdf
		size_t count[3] = {1, static_cast<size_t>(size), 4};
		checkNcError(nc_get_vara_int(ncFile, ncVarElemVertices, start, count, reinterpret_cast<int*>(elemVertices)));
		checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighbors, start, count, reinterpret_cast<int*>(elemNeighbors)));
		checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborSides, start, count, reinterpret_cast<int*>(elemNeighborSides)));
		checkNcError(nc_get_vara_int(ncFile, ncVarElemSideOrientations, start, count, reinterpret_cast<int*>(elemSideOrientations)));
		checkNcError(nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count, reinterpret_cast<int*>(elemBoundaries)));
		checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborRanks, start, count, reinterpret_cast<int*>(elemNeighborRanks)));
		checkNcError(nc_get_vara_int(ncFile, ncVarElemMPIIndices, start, count, reinterpret_cast<int*>(elemMPIIndices)));
		if (hasGroup)
			checkNcError(nc_get_vara_int(ncFile, ncVarElemGroup, start, count, reinterpret_cast<int*>(elemMaterial)));
		EPIK_USER_END(r_read_elements);
		SCOREP_USER_REGION_END( r_read_elements )

		// Copy buffers to elements
		for (int i = 0; i < size; i++) {
			m_elements[i].localId = i;

			memcpy(m_elements[i].vertices, &elemVertices[i], sizeof(ElemVertices));
			memcpy(m_elements[i].neighbors, &elemNeighbors[i], sizeof(ElemNeighbors));
			memcpy(m_elements[i].neighborSides, &elemNeighborSides[i], sizeof(ElemNeighborSides));
			memcpy(m_elements[i].sideOrientations, &elemSideOrientations[i], sizeof(ElemSideOrientations));
			memcpy(m_elements[i].boundaries, &elemBoundaries[i], sizeof(ElemBoundaries));
			memcpy(m_elements[i].neighborRanks, &elemNeighborRanks[i], sizeof(ElemNeighborRanks));
			memcpy(m_elements[i].mpiIndices, &elemMPIIndices[i], sizeof(ElemMPIIndices));
			m_elements[i].material = elemMaterial[i];
		}

		delete [] elemVertices;
		delete [] elemNeighbors;
		delete [] elemNeighborSides;
		delete [] elemSideOrientations;
		delete [] elemBoundaries;
		delete [] elemNeighborRanks;
		delete [] elemMPIIndices;

		// Vertices
		checkNcError(nc_get_var1_int(ncFile, ncVarVrtxSize, start, &size));

		VrtxCoords* vrtxCoords = new VrtxCoords[size];

		m_vertices.resize(size);

		EPIK_USER_REG(r_read_vertices, "read_vertices");
		SCOREP_USER_REGION_DEFINE( r_read_vertices )
		EPIK_USER_START(r_read_vertices);
		SCOREP_USER_REGION_BEGIN( r_read_vertices, "read_vertices", SCOREP_USER_REGION_TYPE_COMMON )
		// Read vertex buffer from netcdf
		count[0] = 1; count[1] = size; count[2] = 3;
		checkNcError(nc_get_vara_double(ncFile, ncVarVrtxCoords, start, count, reinterpret_cast<double*>(vrtxCoords)));
		EPIK_USER_END(r_read_vertices);
		SCOREP_USER_REGION_END( r_read_vertices )

		// Copy buffers to vertices
		for (int i = 0; i < size; i++) {
			memcpy(m_vertices[i].coords, &vrtxCoords[i], sizeof(VrtxCoords));
		}

		delete [] vrtxCoords;

		// Boundaries (MPI neighbors)
		checkNcError(nc_get_var1_int(ncFile, ncVarBndSize, start, &size));

		// Get maximum number of neighbors (required to get collective MPI-IO right)
		int maxNeighbors = bndSize;
		//MPI_Allreduce(MPI_IN_PLACE, &maxNeighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);

		int* bndElemLocalIds = new int[bndElemSize];

		EPIK_USER_REG(r_read_boundaries, "read_boundaries");
  		SCOREP_USER_REGION_DEFINE( r_read_boundaries )
		EPIK_USER_START(r_read_boundaries);
		SCOREP_USER_REGION_BEGIN( r_read_boundaries, "read_boundaries", SCOREP_USER_REGION_TYPE_COMMON )
		int i = 0;
		for (; i < maxNeighbors; i++) {
			MPINeighbor neighbor;
			neighbor.localID = i;

			size_t bndStart[3] = {static_cast<size_t>(rank), static_cast<size_t>(i), 0};

			// Get neighbor rank from netcdf
			int bndRank;
			checkNcError(nc_get_var1_int(ncFile, ncVarBndElemRank, bndStart, &bndRank));

			// Read size of this boundary from netcdf
			int elemSize;
			checkNcError(nc_get_var1_int(ncFile, ncVarBndElemSize, bndStart, &elemSize));

			//int* bndElemLocalIds = new int[bndElemSize];

			if (i < size) {
				neighbor.elements.resize(elemSize);
			}

			// Read local element ids from netcdf
			size_t bndCount[3] = {1, 1, bndElemSize};
			checkNcError(nc_get_vara_int(ncFile, ncVarBndElemLocalIds, bndStart, bndCount, bndElemLocalIds));

			if (i < size) {
				// Copy buffer to boundary
				for (int j = 0; j < elemSize; j++) {
					neighbor.elements[j].localElement = bndElemLocalIds[j];
				}
			}

			//delete [] bndElemLocalIds;

			if (i < size) {
				m_MPINeighbors[bndRank] = neighbor;
			}
		}

		delete [] bndElemLocalIds;

		EPIK_USER_END(r_read_boundaries);
		SCOREP_USER_REGION_END( r_read_boundaries )

		logInfo(rank) << "Finished reading mesh";

		// Close netcdf file
		checkNcError(nc_close(ncFile));

		// Recompute additional information
		findElementsPerVertex();
	}

private:
	/**
	 * Finds all locals elements for each vertex
	 */
	void findElementsPerVertex()
	{
		for (std::vector<Element>::const_iterator i = m_elements.begin();
				i != m_elements.end(); i++) {
			for (int j = 0; j < 4; j++) {
				assert(i->vertices[j] < m_vertices.size());
				m_vertices[i->vertices[j]].elements.push_back(i->localId);
			}
		}
	}

private:
	/**
	 * Switch to collective access for a netCDf variable
	 */
	static void collectiveAccess(int ncFile, int ncVar)
	{
#ifdef USE_MPI
		checkNcError(nc_var_par_access(ncFile, ncVar, NC_COLLECTIVE));
#endif // USE_MPI
	}

	static void checkNcError(int error)
	{
		if (error != NC_NOERR)
			logError() << "Error while reading netCDF file:" << nc_strerror(error);
	}
};

#endif // USE_NETCDF

#endif // NETCDF_READER_H
