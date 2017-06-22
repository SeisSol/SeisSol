/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2013-2016, SeisSol Group
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

#include "Parallel/MPI.h"

#include "MeshReader.h"

#include <cassert>
#include <cstring>


#ifndef NETCDF_PASSIVE

#ifdef USE_MPI
#include <netcdf_par.h>
#endif // USE_MPI
#include <netcdf.h>

#endif // NETCDF_PASSIVE

#include "MeshReader.h"
#include "Initializer/preProcessorMacros.fpp"

#include "utils/logger.h"

class NetcdfReader : public MeshReader
{
public:
	NetcdfReader(int rank, int nProcs, const char* meshFile)
		: MeshReader(rank)
	{
		// Open nc file
		int ncFile;
		int masterRank;
		unsigned int groupSize = 1;
#ifdef USE_MPI
		groupSize = utils::Env::get<unsigned int>("SEISSOL_NETCDF_GROUP_SIZE", 1);
		if (nProcs % groupSize != 0)
			logError() << "#Processes must be a multiple of the group size" << groupSize;

		int master = (rank / groupSize) * groupSize;
		MPI_Comm commMaster;
		MPI_Comm_split(seissol::MPI::mpi.comm(), rank % groupSize == 0 ? 1 : MPI_UNDEFINED,
			rank, &commMaster);

		masterRank = -1;

		if (commMaster != MPI_COMM_NULL) {
#ifdef NETCDF_PASSIVE
			logError() << "netCDF master found with netCDF passive support only";
#else // NETCDF_PASSIVE
			MPI_Comm_rank(commMaster, &masterRank);
			checkNcError(nc_open_par(meshFile, NC_NETCDF4 | NC_MPIIO, commMaster, MPI_INFO_NULL, &ncFile));
#endif // NETCDF_PASSIVE
		}
#else // USE_MPI
		masterRank = rank; // = 0;
		checkNcError(nc_open(meshFile, NC_NETCDF4, &ncFile));
#endif // USE_MPI

		size_t bndSize = -1;
		size_t bndElemSize = -1;

		int ncVarElemSize = -1;
		int ncVarElemVertices = -1;
		int ncVarElemNeighbors = -1;
		int ncVarElemBoundaries = -1;
		int ncVarElemNeighborSides = -1;
		int ncVarElemSideOrientations = -1;
		int ncVarElemNeighborRanks = -1;
		int ncVarElemMPIIndices = -1;
		int ncVarElemGroup = -1;
		bool hasGroup = false;
		int ncVarVrtxSize = -1;
		int ncVarVrtxCoords = -1;
		int ncVarBndSize = -1;
		int ncVarBndElemSize = -1;
		int ncVarBndElemRank = -1;
		int ncVarBndElemLocalIds = -1;

		int* sizes = 0L;
		int maxSize = 0;
		if (masterRank == 0) {
#ifdef NETCDF_PASSIVE
			assert(false);
#else // NETCDF_PASSIVE
			// Get important dimensions
			int ncDimPart;
			checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
			size_t partitions;
			checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &partitions));

			if (partitions != static_cast<unsigned int>(nProcs))
				logError() << "Number of partitions in netCDF file does not match number of MPI ranks.";

			int ncDimBndSize;
			checkNcError(nc_inq_dimid(ncFile, "boundaries", &ncDimBndSize));
			checkNcError(nc_inq_dimlen(ncFile, ncDimBndSize, &bndSize));

			int ncDimBndElem;
			checkNcError(nc_inq_dimid(ncFile, "boundary_elements", &ncDimBndElem));
			checkNcError(nc_inq_dimlen(ncFile, ncDimBndElem, &bndElemSize));
#endif // NETCDF_PASSIVE
		}

#ifdef USE_MPI
		// Broadcast boundary sizes
		unsigned long buf[2] = {bndSize, bndElemSize};
		MPI_Bcast(buf, 2, MPI_UNSIGNED_LONG, 0, seissol::MPI::mpi.comm());
		bndSize = buf[0];
		bndElemSize = buf[1];
#endif // USE_MPI

		if (masterRank >= 0) {
#ifdef NETCDF_PASSIVE
			assert(false);
#else // NETCDF_PASSIVE
			// Create netcdf variables
			checkNcError(nc_inq_varid(ncFile, "element_size", &ncVarElemSize));
			collectiveAccess(ncFile, ncVarElemSize);

			checkNcError(nc_inq_varid(ncFile, "element_vertices", &ncVarElemVertices));
			collectiveAccess(ncFile, ncVarElemVertices);

			checkNcError(nc_inq_varid(ncFile, "element_neighbors", &ncVarElemNeighbors));
			collectiveAccess(ncFile, ncVarElemNeighbors);

			checkNcError(nc_inq_varid(ncFile, "element_boundaries", &ncVarElemBoundaries));
			collectiveAccess(ncFile, ncVarElemBoundaries);

			checkNcError(nc_inq_varid(ncFile, "element_neighbor_sides", &ncVarElemNeighborSides));
			collectiveAccess(ncFile, ncVarElemNeighborSides);

			checkNcError(nc_inq_varid(ncFile, "element_side_orientations", &ncVarElemSideOrientations));
			collectiveAccess(ncFile, ncVarElemSideOrientations);

			checkNcError(nc_inq_varid(ncFile, "element_neighbor_ranks", &ncVarElemNeighborRanks));
			collectiveAccess(ncFile, ncVarElemNeighborRanks);

			checkNcError(nc_inq_varid(ncFile, "element_mpi_indices", &ncVarElemMPIIndices));
			collectiveAccess(ncFile, ncVarElemMPIIndices);

			int ncResult = nc_inq_varid(ncFile, "element_group", &ncVarElemGroup);

			if (ncResult != NC_ENOTVAR) {
				checkNcError(ncResult);
				hasGroup = true;
				collectiveAccess(ncFile, ncVarElemGroup);
			}

			checkNcError(nc_inq_varid(ncFile, "vertex_size", &ncVarVrtxSize));
			collectiveAccess(ncFile, ncVarVrtxSize);

			checkNcError(nc_inq_varid(ncFile, "vertex_coordinates", &ncVarVrtxCoords));
			collectiveAccess(ncFile, ncVarVrtxCoords);

			checkNcError(nc_inq_varid(ncFile, "boundary_size", &ncVarBndSize));
			collectiveAccess(ncFile, ncVarBndSize);

			checkNcError(nc_inq_varid(ncFile, "boundary_element_size", &ncVarBndElemSize));
			collectiveAccess(ncFile, ncVarBndElemSize);

			checkNcError(nc_inq_varid(ncFile, "boundary_element_rank", &ncVarBndElemRank));
			collectiveAccess(ncFile, ncVarBndElemRank);

			checkNcError(nc_inq_varid(ncFile, "boundary_element_localids", &ncVarBndElemLocalIds));
			collectiveAccess(ncFile, ncVarBndElemLocalIds);

			logInfo(rank) << "Start reading mesh from netCDF file";

			// Elements
			sizes = new int[groupSize];

			for (int i = groupSize-1; i >= 0; i--){
				size_t start = static_cast<size_t>(i+rank);

				int size;
				checkNcError(nc_get_var1_int(ncFile, ncVarElemSize, &start, &size));
				sizes[i] = size;
				if (i != 0) {
#ifdef USE_MPI
					// Forward size (except the own)
					MPI_Send(&size, 1, MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
#else // USE_MPI
					assert(false);
#endif // USE_MPI
				}

				maxSize = std::max(maxSize, size);
			}
#endif
		} else {
#ifdef USE_MPI
			sizes = new int[1];
			MPI_Recv(sizes, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);

			maxSize = sizes[0];
#else // USE_MPI
			assert(false);
#endif // USE_MPI
		}

#ifdef USE_MPI
		// Broadcast group information
		int iHasGroup = hasGroup;
		MPI_Bcast(&iHasGroup, 1, MPI_INT, 0, seissol::MPI::mpi.comm());
		hasGroup = iHasGroup != 0;
#endif // USE_MPI

		m_elements.resize(sizes[0]);

		ElemVertices* elemVertices = new ElemVertices[maxSize];
		ElemNeighbors* elemNeighbors = new ElemNeighbors[maxSize];
		ElemNeighborSides* elemNeighborSides = new ElemNeighborSides[maxSize];
		ElemSideOrientations* elemSideOrientations = new ElemSideOrientations[maxSize];
		ElemBoundaries* elemBoundaries = new ElemBoundaries[maxSize];
		ElemNeighborRanks* elemNeighborRanks = new ElemNeighborRanks[maxSize];
		ElemMPIIndices* elemMPIIndices = new ElemMPIIndices[maxSize];
		ElemMaterial* elemMaterial = new ElemMaterial[maxSize];

//		SCOREP_USER_REGION_DEFINE( r_read_elements )
//		SCOREP_USER_REGION_BEGIN( r_read_elements, "read_elements", SCOREP_USER_REGION_TYPE_COMMON )

		// Read element buffers from netcdf
		if (masterRank >= 0) {
#ifdef NETCDF_PASSIVE
			assert(false);
#else // NETCDF_PASSIVE
			size_t start[3] = {0, 0, 0};
			size_t count[3] = {1, 0, 4};

			for(int i = groupSize - 1; i >= 0; i--){
				start[0] = static_cast<size_t>(i+rank);
				count[1] = static_cast<size_t>(sizes[i]);

				checkNcError(nc_get_vara_int(ncFile, ncVarElemVertices, start, count, reinterpret_cast<int*>(elemVertices)));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighbors, start, count, reinterpret_cast<int*>(elemNeighbors)));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborSides, start, count, reinterpret_cast<int*>(elemNeighborSides)));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemSideOrientations, start, count, reinterpret_cast<int*>(elemSideOrientations)));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count, reinterpret_cast<int*>(elemBoundaries)));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborRanks, start, count, reinterpret_cast<int*>(elemNeighborRanks)));
				checkNcError(nc_get_vara_int(ncFile, ncVarElemMPIIndices, start, count, reinterpret_cast<int*>(elemMPIIndices)));
				if (hasGroup)
					checkNcError(nc_get_vara_int(ncFile, ncVarElemGroup, start, count, reinterpret_cast<int*>(elemMaterial)));

				if (i != 0) {
#ifdef USE_MPI
					MPI_Send(elemVertices, 4*sizes[i], MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
					MPI_Send(elemNeighbors, 4*sizes[i], MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
					MPI_Send(elemNeighborSides, 4*sizes[i], MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
					MPI_Send(elemSideOrientations, 4*sizes[i], MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
					MPI_Send(elemBoundaries, 4*sizes[i], MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
					MPI_Send(elemNeighborRanks, 4*sizes[i], MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
					MPI_Send(elemMPIIndices, 4*sizes[i], MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
					MPI_Send(elemMaterial, sizes[i], MPI_INT, i+rank ,0 , seissol::MPI::mpi.comm());
#else // USE_MPI
					assert(false);
#endif // USE_MPI
				}
			}
#endif // NETCDF_PASSIVE
		} else {
#ifdef USE_MPI
			MPI_Recv(elemVertices, 4*sizes[0], MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			MPI_Recv(elemNeighbors, 4*sizes[0] ,MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			MPI_Recv(elemNeighborSides, 4*sizes[0], MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			MPI_Recv(elemSideOrientations, 4*sizes[0], MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			MPI_Recv(elemBoundaries, 4*sizes[0], MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			MPI_Recv(elemNeighborRanks, 4*sizes[0], MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			MPI_Recv(elemMPIIndices, 4*sizes[0], MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			MPI_Recv(elemMaterial, sizes[0], MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
#else // USE_MPI
			assert(false);
#endif // USE_MPI
		}

		// Copy buffers to elements
		for (int i = 0; i < sizes[0]; i++) {
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

//		SCOREP_USER_REGION_END( r_read_elements )

		delete [] elemVertices;
		delete [] elemNeighbors;
		delete [] elemNeighborSides;
		delete [] elemSideOrientations;
		delete [] elemBoundaries;
		delete [] elemNeighborRanks;
		delete [] elemMPIIndices;

		// Vertices
		maxSize = 0;
		if (masterRank >= 0) {
#ifdef NETCDF_PASSIVE
			assert(false);
#else // NETCDF_PASSIVE
			for (int i = groupSize - 1; i >= 0; i--) {
				size_t start = static_cast<size_t>(i+rank);

				int size;
				checkNcError(nc_get_var1_int(ncFile, ncVarVrtxSize, &start, &size));
				sizes[i] = size;
				if (i != 0) {
#ifdef USE_MPI
					MPI_Send(&size, 1, MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
#else // USE_MPI
					assert(false);
#endif // USE_MPI
				}

				maxSize = std::max(maxSize,size);
			}
#endif // NETCDF_PASSIVE
		} else {
#ifdef USE_MPI
			MPI_Recv(sizes, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
			maxSize = sizes[0];
#else // USE_MPI
			assert(false);
#endif // USE_MPi
		}

		m_vertices.resize(sizes[0]);
		VrtxCoords* vrtxCoords = new VrtxCoords[maxSize];

//		SCOREP_USER_REGION_DEFINE( r_read_vertices )
//		SCOREP_USER_REGION_BEGIN( r_read_vertices, "read_vertices", SCOREP_USER_REGION_TYPE_COMMON )

		// Read vertex buffer from netcdf
		if (masterRank >= 0) {
#ifdef NETCDF_PASSIVE
			assert(false);
#else // NETCDF_PASSIVE
			size_t start[3] = {0, 0, 0};
			size_t count[3] = {1, 0, 3};

			for (int i = groupSize - 1; i >= 0; i--) {
				start[0] = static_cast<size_t>(i+rank);
				count[1] = static_cast<size_t>(sizes[i]);

				checkNcError(nc_get_vara_double(ncFile, ncVarVrtxCoords, start, count, reinterpret_cast<double*>(vrtxCoords)));
				if (i != 0) {
#ifdef USE_MPI
					MPI_Send(vrtxCoords, 3*sizes[i], MPI_DOUBLE, i+rank, 0, seissol::MPI::mpi.comm());
#else // USE_MPI
					assert(false);
#endif // USE_MPI
				}
			}
#endif // NETCDF_PASSIVE
		} else {
#ifdef USE_MPI
			MPI_Recv(vrtxCoords, 3*sizes[0], MPI_DOUBLE, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
#else // USE_MPI
			assert(false);
#endif // USE_MPI
		}

//		SCOREP_USER_REGION_END( r_read_vertices )

		// Copy buffers to vertices
		for (int i = 0; i < sizes[0]; i++) {
			memcpy(m_vertices[i].coords, &vrtxCoords[i], sizeof(VrtxCoords));
		}

		delete [] vrtxCoords;

		// Boundaries (MPI neighbors)
		if (masterRank >= 0) {
#ifdef NETCDF_PASSIVE
			assert(false);
#else // NETCDF_PASSIVE
			for (int i = groupSize - 1; i >= 0; i--) {
				size_t start = static_cast<size_t>(i+rank);

				int size;
				checkNcError(nc_get_var1_int(ncFile, ncVarBndSize, &start, &size));
				sizes[i] = size;

				if (i != 0) {
#ifdef USE_MPI
					MPI_Send(&size, 1, MPI_INT, i+rank, 0, seissol::MPI::mpi.comm());
#else // USE_MPI
					assert(false);
#endif // USE_MPI
				}
			}
#endif // NETCDF_PASSIVE
		} else {
#ifdef USE_MPI
			MPI_Recv(sizes, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
#else // USE_MPI
			assert(false);
#endif // USE_MPI
		}


		// Get maximum number of neighbors (required to get collective MPI-IO right)
		int maxNeighbors = bndSize;
		//MPI_Allreduce(MPI_IN_PLACE, &maxNeighbors, 1, MPI_INT, MPI_MAX, seissol::MPI::mpi.comm());
		int* bndElemLocalIds = new int[bndElemSize];

//		SCOREP_USER_REGION_DEFINE( r_read_boundaries );
//		SCOREP_USER_REGION_BEGIN( r_read_boundaries, "read_boundaries", SCOREP_USER_REGION_TYPE_COMMON );

		size_t bndStart[3] = {0, 0, 0};
		for (int i = 0; i < maxNeighbors; i++) {
			bndStart[1] = static_cast<size_t>(i);

			if (masterRank >= 0) {
#ifdef NETCDF_PASSIVE
				assert(false);
#else // NETCDF_PASSIVE
				for (int j = groupSize - 1; j >= 0; j--) {
					bndStart[0] = static_cast<size_t>(j+rank);

					// Get neighbor rank from netcdf
					int bndRank;
					checkNcError(nc_get_var1_int(ncFile, ncVarBndElemRank, bndStart, &bndRank));

					// Read size of this boundary from netcdf
					int elemSize;
					checkNcError(nc_get_var1_int(ncFile, ncVarBndElemSize, bndStart, &elemSize));

					size_t bndCount[3] = {1, 1, bndElemSize};
					checkNcError(nc_get_vara_int(ncFile, ncVarBndElemLocalIds, bndStart, bndCount, bndElemLocalIds));

					if (i < sizes[j]) {
						if (j == 0) {
							addMPINeighbor(i, bndRank, elemSize, bndElemLocalIds);
						} else {
#ifdef USE_MPI
							MPI_Send(&bndRank, 1, MPI_INT, j+rank, 0, seissol::MPI::mpi.comm());
							MPI_Send(&elemSize, 1, MPI_INT, j+rank, 0, seissol::MPI::mpi.comm());
							MPI_Send(bndElemLocalIds, elemSize, MPI_INT, j+rank, 0, seissol::MPI::mpi.comm());
#else // USE_MPI
							assert(false);
#endif // USE_MPI
						}
					}
				}
#endif // NETCDF_PASSIVE
			} else {
				if (i < sizes[0]) {
#ifdef USE_MPI
					int bndRank;
					MPI_Recv(&bndRank, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);
					int elemSize;
					MPI_Recv(&elemSize, 1, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);

					MPI_Recv(bndElemLocalIds, elemSize, MPI_INT, master, 0, seissol::MPI::mpi.comm(), MPI_STATUS_IGNORE);

					addMPINeighbor(i, bndRank, elemSize, bndElemLocalIds);
#else // USE_MPI
					assert(false);
#endif // USE_MPI
				}
			}
		}

		delete [] bndElemLocalIds;

//		SCOREP_USER_REGION_END( r_read_boundaries )

		delete [] sizes;

		logInfo(rank) << "Finished reading mesh";

		// Close netcdf file
		if (masterRank >= 0) {
#ifndef NETCDF_PASSIVE
			checkNcError(nc_close(ncFile));
#ifdef USE_MPI
			MPI_Comm_free(&commMaster);
#endif // USE_MPI
#endif // NETCDF_PASSIVE
		}

		// Recompute additional information
		findElementsPerVertex();
	}

private:
	/**
	 * Adds the information about an MPI neighbor
	 *
	 * @localID The local id of the neighborhood
	 * @bndRank The neighbor rank
	 * @elemSize Number of boundary elements
	 * @bndElemLocalIds List of boundary ids
	 */
	void addMPINeighbor(int localID, int bndRank, int elemSize, const int* bndElemLocalIds)
	{
		MPINeighbor neighbor;
		neighbor.localID = localID;

		neighbor.elements.resize(elemSize);
		for (int i = 0; i < elemSize; i++)
			neighbor.elements[i].localElement = bndElemLocalIds[i];

		m_MPINeighbors[bndRank] = neighbor;
	}

	/**
	 * Finds all locals elements for each vertex
	 */
	void findElementsPerVertex()
	{
		for (std::vector<Element>::const_iterator i = m_elements.begin();
				i != m_elements.end(); i++) {
			for (int j = 0; j < 4; j++) {
				assert(i->vertices[j] < static_cast<int>(m_vertices.size()));
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
#ifndef NETCDF_PASSIVE
		checkNcError(nc_var_par_access(ncFile, ncVar, NC_COLLECTIVE));
#else // NETCDF_PASSIVE
		assert(false);
#endif // NETCDF_PASSIVE
#endif // USE_MPI
	}

	static void checkNcError(int error)
	{
#ifndef NETCDF_PASSIVE
		if (error != NC_NOERR)
			logError() << "Error while reading netCDF file:" << nc_strerror(error);
#endif // NETCDF_PASSIVE
	}
};

#endif // USE_NETCDF

#endif // NETCDF_READER_H
