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

#ifndef NETCDF_COPROCESSOR

#ifdef USE_MPI

#include <netcdf_par.h>
#endif // USE_MPI
#include <netcdf.h>

#endif //NETCDF_COPROCESSOR

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
		int modulo=utils::Env::get<int>("SEISSOL_NETCDF_HOST_RANKS",1);

		int host= rank/modulo*modulo;

		int color= rank % modulo == 0 ? 1 : MPI_UNDEFINED;

		MPI_Comm comm_host;
		MPI_Comm_split(MPI_COMM_WORLD,color,rank,&comm_host);

		int host_rank = -1;


		if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		  MPI_Comm_rank(comm_host,&host_rank);
		  checkNcError(nc_open_par(meshFile, NC_NETCDF4 | NC_MPIIO, comm_host, MPI_INFO_NULL, &ncFile));
#endif //NETCDF_COPROCESSOR
		}

#else // USE_MPI
		  checkNcError(nc_open(meshFile, NC_NETCDF4, &ncFile));

#endif // USE_MPI




		int size;
		size_t start[3] = {static_cast<size_t>(rank), 0, 0};
		MPI_Status status;

		unsigned long bndElemSize_buffer =0;
		unsigned long bndSize_buffer =0;
		size_t bndSize=-1;
		size_t bndElemSize=-1;
		int ncVarElemVertices=-1;
		int ncVarElemNeighbors=-1;
		int ncVarElemNeighborSides=-1;
		int ncVarElemSideOrientations=-1;
		int ncVarElemBoundaries=-1;
		int ncVarElemNeighborRanks=-1;
		int ncVarElemMPIIndices=-1;
		int ncVarElemGroup=-1;
		int ncVarElemSize=-1;
		bool hasGroup = false;
		int ncVarVrtxSize=-1;
		int ncVarVrtxCoords=-1;
		int ncVarBndSize=-1;
		int ncVarBndElemSize=-1;
		int ncVarBndElemRank=-1;

		int ncVarBndElemLocalIds=-1;
		int ncDimBndSize;
		int ncDimBndElem;
				
		printf("World_Rank: %d Host_Rank: %d\n",rank,host_rank) ;
		int* size_ar = (int*)malloc(sizeof(int)*modulo);		
		int max_size=0;

		if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		  // Get important dimensions
		  int ncDimPart;
		  checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));
		  size_t partitions;
		  checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &partitions));
		  
		  if (partitions != static_cast<unsigned int>(nProcs))
		    logError() << "Number of partitions in netCDF file does not match number of MPI ranks.";
		  checkNcError(nc_inq_dimid(ncFile, "boundaries", &ncDimBndSize));
		  checkNcError(nc_inq_dimlen(ncFile, ncDimBndSize, &bndSize));
		  checkNcError(nc_inq_dimid(ncFile, "boundary_elements", &ncDimBndElem));
		  checkNcError(nc_inq_dimlen(ncFile, ncDimBndElem, &bndElemSize));
		  
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
		  //		pEPIK_USER_END(r_init_netcdf);
		
		  logInfo(rank) << "Start reading mesh from netCDF file";

		  // Elements

		  bndElemSize_buffer=static_cast<unsigned long>(bndElemSize);
		  bndSize_buffer=static_cast<unsigned long>(bndSize);

		  for(int i = rank+(modulo -1) ; i >= rank ; i--){
		    start[0] = static_cast<size_t>(i);
		    checkNcError(nc_get_var1_int(ncFile, ncVarElemSize, start, &size));
		    size_ar[i-rank]=size;
		    max_size = std::max(max_size,size);
		    if(rank != i) {
		      //		      printf("Send size %d bndElemSize %d and bndSize %d to %d\n",size,static_cast<int>(bndElemSize),static_cast<int>(bndSize),i); 
		      MPI_Send(&size,1,MPI_INT,i,0,MPI_COMM_WORLD);
		      MPI_Send(&bndElemSize_buffer,1,MPI_UNSIGNED_LONG,i,0,MPI_COMM_WORLD);
		      MPI_Send(&bndSize_buffer,1,MPI_UNSIGNED_LONG,i,0,MPI_COMM_WORLD);
		    }
		  }
#endif		 
		}else{
		  MPI_Recv(&size,1,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(&bndElemSize_buffer,1,MPI_UNSIGNED_LONG,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(&bndSize_buffer,1,MPI_UNSIGNED_LONG,host,0,MPI_COMM_WORLD,&status);
		  bndSize=static_cast<size_t>(bndSize_buffer);
		  bndElemSize=static_cast<size_t>(bndElemSize_buffer); 
		  //	  printf("Recv size %d bndElemSize %d and bndSize %d from %d\n",size,static_cast<int>(bndElemSize),static_cast<int>(bndSize),host);  
		  max_size=size;
		}

		m_elements.resize(size);

		ElemVertices* elemVertices = new ElemVertices[max_size];
		ElemNeighbors* elemNeighbors = new ElemNeighbors[max_size];
		ElemNeighborSides* elemNeighborSides = new ElemNeighborSides[max_size];
		ElemSideOrientations* elemSideOrientations = new ElemSideOrientations[max_size];
		ElemBoundaries* elemBoundaries = new ElemBoundaries[max_size];
		ElemNeighborRanks* elemNeighborRanks = new ElemNeighborRanks[max_size];
		ElemMPIIndices* elemMPIIndices = new ElemMPIIndices[max_size];
		ElemMaterial* elemMaterial = new ElemMaterial[max_size];
		

		    
		EPIK_USER_REG(r_read_elements, "read_elements");
		SCOREP_USER_REGION_DEFINE( r_read_elements )
		EPIK_USER_START(r_read_elements);
		SCOREP_USER_REGION_BEGIN( r_read_elements, "read_elements", SCOREP_USER_REGION_TYPE_COMMON )
		// Read element buffers from netcdf
		size_t count[3] = {1, static_cast<size_t>(size), 4};
		if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		  for(int i = rank+(modulo -1) ; i >= rank ; i--){
		    start[0] = static_cast<size_t>(i);
		    count[1]=static_cast<size_t>(size_ar[i-rank]);
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
		  if(rank != i){
		    MPI_Send(elemVertices,4*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    MPI_Send(elemNeighbors,4*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    MPI_Send(elemNeighborSides,4*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    MPI_Send(elemSideOrientations,4*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    MPI_Send(elemBoundaries,4*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    MPI_Send(elemNeighborRanks,4*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    MPI_Send(elemMPIIndices,4*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    MPI_Send(elemMaterial,1*size_ar[i-rank],MPI_INT,i,0,MPI_COMM_WORLD);
		    //		    printf("Send elements to %d\n",i);
		  }
    		  
		  }
#endif
		}else{
		  MPI_Recv(elemVertices,4*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(elemNeighbors,4*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(elemNeighborSides,4*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(elemSideOrientations,4*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(elemBoundaries,4*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(elemNeighborRanks,4*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(elemMPIIndices,4*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  MPI_Recv(elemMaterial,1*size,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  //		  printf("Recv elements from %d\n",host); 
		}
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
		max_size=0;
		if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		  for(int i = rank+(modulo -1) ; i >= rank ; i--){
		    start[0] = static_cast<size_t>(i);
		    checkNcError(nc_get_var1_int(ncFile, ncVarVrtxSize, start, &size));
		    size_ar[i-rank]=size;
		    max_size=std::max(max_size,size);
		    if(rank != i) {
		      //		      printf("2. Send size %d to %d\n",size,i); 
		      MPI_Send(&size,1,MPI_INT,i,0,MPI_COMM_WORLD);
		    }
		  }
#endif
		}else{
		  MPI_Recv(&size,1,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  max_size=size;
		  //		  printf("2. Recv size %d from %d\n",size, host); 
		}
	
		VrtxCoords* vrtxCoords = new VrtxCoords[max_size];
		m_vertices.resize(size);
	
		EPIK_USER_REG(r_read_vertices, "read_vertices");
		SCOREP_USER_REGION_DEFINE( r_read_vertices )
		EPIK_USER_START(r_read_vertices);
		SCOREP_USER_REGION_BEGIN( r_read_vertices, "read_vertices", SCOREP_USER_REGION_TYPE_COMMON )
		  // Read vertex buffer from netcdf
		  if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		    count[0] = 1;
		    count[2] = 3;
		    for(int i = rank+(modulo -1) ; i >= rank ; i--){
		      start[0] = static_cast<size_t>(i);
		      count[1] = static_cast<size_t>(size_ar[i-rank]);
		      checkNcError(nc_get_vara_double(ncFile, ncVarVrtxCoords, start, count, reinterpret_cast<double*>(vrtxCoords)));
		      if(rank != i){
			//			printf("3.Send vertexCoords to %d\n",i); 
			MPI_Send(vrtxCoords,3*size_ar[i-rank],MPI_DOUBLE,i,0,MPI_COMM_WORLD);
		      }
		    }
#endif
		  }else{
		    MPI_Recv(vrtxCoords,3*size,MPI_DOUBLE,host,0,MPI_COMM_WORLD,&status); 
		    //		    printf("3.Recv vertexCoords from %d\n",host); 
		  }
		EPIK_USER_END(r_read_vertices);
		SCOREP_USER_REGION_END( r_read_vertices )
		  
		// Copy buffers to vertices
		for (int i = 0; i < size; i++) {
		  //  printf("%d copy vertex Coords %d \t",rank,i);
		  for(int j = 0 ; j< 3 ; j++){
		    //		    printf("%f \t",vrtxCoords[i][j]);
		    //  printf("%f \t",m_vertices[i].coords[j]);
		  }
		  //		  printf("\n");
		  
		  memcpy(m_vertices[i].coords, &vrtxCoords[i], sizeof(VrtxCoords));
		  // 		  printf("%d copied vertex Coords %d \n", rank,i);
		}
		//		printf("%d wrote all vertex Coords\n", rank);
		delete [] vrtxCoords;

		// Boundaries (MPI neighbors)
		if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		  for(int i = rank+(modulo -1) ; i >= rank ; i--){
		    start[0] = static_cast<size_t>(i);
		    checkNcError(nc_get_var1_int(ncFile, ncVarBndSize, start, &size));
		    if(rank != i) {
		      //		      printf("4.Send size to %d\n",i); 
		      MPI_Send(&size,1,MPI_INT,i,0,MPI_COMM_WORLD);
		    }
		  }
#endif
		}else{
		  MPI_Recv(&size,1,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		  //		  printf("4.Recv size from %d\n",host); 
		}
		
		
		// Get maximum number of neighbors (required to get collective MPI-IO right)
		int maxNeighbors = bndSize;
		//MPI_Allreduce(MPI_IN_PLACE, &maxNeighbors, 1, MPI_INT, MPI_MAX, MPI_COMM_WORLD);
		//		printf("5.1 %d bndElemSize %lu \n",rank, bndElemSize); 
		int* bndElemLocalIds = new int[bndElemSize];
		
		EPIK_USER_REG(r_read_boundaries, "read_boundaries");
  		SCOREP_USER_REGION_DEFINE( r_read_boundaries )
		EPIK_USER_START(r_read_boundaries);
		SCOREP_USER_REGION_BEGIN( r_read_boundaries, "read_boundaries", SCOREP_USER_REGION_TYPE_COMMON )
		int i = 0;
		int elemSize;
		size_t bndStart[3] = {0,0, 0};	
		int bndRank = -1;
		//		printf("5.2 %d\n",rank); 		
		for (; i < maxNeighbors; i++) {
		  bndStart[1]=static_cast<size_t>(i);
		  MPINeighbor neighbor;
		  neighbor.localID = i;
		  if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		    for(int j = rank+(modulo -1) ; j >= rank ; j--){
		      bndStart[0]=static_cast<size_t>(j); 
		      
		      // Get neighbor rank from netcdf
		      checkNcError(nc_get_var1_int(ncFile, ncVarBndElemRank, bndStart, &bndRank));
		      
		      // Read size of this boundary from netcdf
		      checkNcError(nc_get_var1_int(ncFile, ncVarBndElemSize, bndStart, &elemSize));
		      
		      if(rank != j){
			MPI_Send(&bndRank,1,MPI_INT,j,0,MPI_COMM_WORLD);
			MPI_Send(&elemSize,1,MPI_INT,j,0,MPI_COMM_WORLD);
			//			printf("5.Send Elem size to %d\n",j); 
		      }
		      //int* bndElemLocalIds = new int[bndElemSize];
		    }
#endif
		  }else{
		    //		    printf("5.3 Recv Elem size from %d\n",host); 
		    MPI_Recv(&bndRank,1,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		    MPI_Recv(&elemSize,1,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		    //		    printf("5.Recv Elem size from %d\n",host); 
		  }
		  
		  if (i < size) {
		    neighbor.elements.resize(elemSize);
		  }
		  
		  if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		    for(int j = rank+(modulo -1) ; j >= rank ; j--){
		      bndStart[0]=static_cast<size_t>(j); 
		      // Read local element ids from netcdf
		      size_t bndCount[3] = {1, 1, bndElemSize};
		      checkNcError(nc_get_vara_int(ncFile, ncVarBndElemLocalIds, bndStart, bndCount, bndElemLocalIds));
		      if(rank != j){
			MPI_Send(bndElemLocalIds,bndElemSize,MPI_INT,j,0,MPI_COMM_WORLD);
			//			printf("6.Send LocalIds to %d\n",j); 
		      }
		    }
#endif
		  }else{
		    MPI_Recv(bndElemLocalIds,bndElemSize,MPI_INT,host,0,MPI_COMM_WORLD,&status);
		    //		    printf("6.Recv LocalIds from %d\n",host); 
		  }

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
		free(size_ar);
		delete [] bndElemLocalIds;

		EPIK_USER_END(r_read_boundaries);
		SCOREP_USER_REGION_END( r_read_boundaries )

		logInfo(rank) << "Finished reading mesh";

		// Close netcdf file
		if(MPI_COMM_NULL != comm_host){
#ifndef NETCDF_COPROCESSOR
		  checkNcError(nc_close(ncFile));
#endif
		  MPI_Comm_free(&comm_host);
		}

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
#ifndef NETCDF_COPROCESSOR
		checkNcError(nc_var_par_access(ncFile, ncVar, NC_COLLECTIVE));
#endif
#endif // USE_MPI
	}

	static void checkNcError(int error)
	{
#ifndef NETCDF_COPROCESSOR
		if (error != NC_NOERR)
			logError() << "Error while reading netCDF file:" << nc_strerror(error);
#endif
	}
};

#endif // USE_NETCDF

#endif // NETCDF_READER_H
