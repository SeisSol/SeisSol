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

/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (rettenbs AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
 * @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
 *
 * @section LICENSE
 * Copyright (c) 2014, SeisSol Group
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

#include "utils/args.h"
#include "utils/logger.h"

#include <netcdf.h>
#include <list>

#include <omp.h>

#include <cstring>
#include <iostream>
#include <cassert>

class RankDescriptor {
public:
  unsigned long m_numberOfElements;
  unsigned long m_numberOfNeighboringMpiRanks;
  unsigned long m_maximumNumberOfMpiFacesPerRank;
  unsigned long m_numberOfInteriorElements;
  unsigned long m_numberOfMpiElements;
  unsigned long m_numberOfMpiCopyLayerElements;
  unsigned long m_numberOfMpiCopyLayerDrElements;
  unsigned long m_numberOfMpiFaces;
  unsigned long m_numberOfMpiDynamicRuptureFaces;
  unsigned long m_numberOfDynamicRuptureFaces;
  unsigned long m_numberOfFreeSurfaceFaces;
  unsigned long m_numberOfTransparentFaces;
  std::map<int, int> m_numberOfElementsInCopyLayerPerRank;
  std::map<int, int> m_numberOfDrElementsInCopyLayerPerRank;

  RankDescriptor() : m_numberOfElements(0),
                     m_numberOfNeighboringMpiRanks(0),
                     m_maximumNumberOfMpiFacesPerRank(0),
                     m_numberOfInteriorElements(0), 
                     m_numberOfMpiElements(0),
                     m_numberOfMpiCopyLayerElements(0),
                     m_numberOfMpiCopyLayerDrElements(0),
                     m_numberOfMpiFaces(0), 
                     m_numberOfMpiDynamicRuptureFaces(0),
                     m_numberOfDynamicRuptureFaces(0),
                     m_numberOfFreeSurfaceFaces(0),
                     m_numberOfTransparentFaces(0)
  {
    m_numberOfElementsInCopyLayerPerRank.clear();
    m_numberOfDrElementsInCopyLayerPerRank.clear();
  }
};

static void getMinMaxAvg(int* values, unsigned int count, int &min, int &max, double &avg)
{
  min = values[0];
  max = values[0];
  unsigned long sum = values[0];
  for (unsigned int i = 1; i < count; i++) {
    if (values[i] < min)
      min = values[i];
    if (values[i] > max)
      max = values[i];
    sum += values[i];
  }
  avg = (static_cast<double>(sum)/static_cast<double>(count));
}

void checkNcError(int error)
{
  if (error != NC_NOERR)
    logError() << "Error while writing netCDF file:" << nc_strerror(error);
}

int main(int argc, char* argv[])
{
  // Parse command line arguments
  utils::Args args;
  args.addOption("fault", 'f', "Include fault statistics", utils::Args::No, false);
  args.addOption("matrices", 'm', "Include matrix statistics", utils::Args::No, false);
  args.addAdditionalOption("input", "netCDF mesh file");
        args.addAdditionalOption("partitions", "prints the information in csv-format per partition.", false);

  // Parse/check command line arguments
  if (args.parse(argc, argv) != utils::Args::Success)
    return 1;
  
  bool includeFault = args.getArgument("fault", false);
  bool includeMatrices = args.getArgument("matrices", false);

  logInfo() << "Using" << omp_get_max_threads() << "threads";
  if (includeFault)
    logInfo() << "Including fault statistics";
  if (includeMatrices)
    logInfo() << "Including matrix statistics";

  // Open the netcdf file
  int ncFile;
  checkNcError(nc_open(args.getAdditionalArgument<std::string>("input").c_str(), 0, &ncFile));

  // Get dimension ids
  int ncDimPart;
  checkNcError(nc_inq_dimid(ncFile, "partitions", &ncDimPart));

  int ncDimElem;
  checkNcError(nc_inq_dimid(ncFile, "elements", &ncDimElem));

  int ncDimBnd;
  checkNcError(nc_inq_dimid(ncFile, "boundaries", &ncDimBnd));

  // Get Variable ids
  int ncVarElemSize;
  checkNcError(nc_inq_varid(ncFile, "element_size", &ncVarElemSize));

  int ncVarBndSize;
  checkNcError(nc_inq_varid(ncFile, "boundary_size", &ncVarBndSize));
  int ncVarBndElemSize;
  checkNcError(nc_inq_varid(ncFile, "boundary_element_size", &ncVarBndElemSize));

  int ncVarElemBoundaries;
  checkNcError(nc_inq_varid(ncFile, "element_boundaries", &ncVarElemBoundaries));
  int ncVarElemNeighorSides;
  checkNcError(nc_inq_varid(ncFile, "element_neighbor_sides", &ncVarElemNeighorSides));
  int ncVarElemSideOrientations;
  checkNcError(nc_inq_varid(ncFile, "element_side_orientations", &ncVarElemSideOrientations));
  int ncVarElemNeighborRanks;
  checkNcError(nc_inq_varid(ncFile, "element_neighbor_ranks", &ncVarElemNeighborRanks));

  // Get partition information
  size_t partitions;
  checkNcError(nc_inq_dimlen(ncFile, ncDimPart, &partitions));

  logInfo() << "Mesh contains" << partitions << "partitions";

  // Get max number of elements
  size_t maxElements;
  checkNcError(nc_inq_dimlen(ncFile, ncDimElem, &maxElements));

  // Get element statistics
  int* element_size = new int[partitions];
  checkNcError(nc_get_var_int(ncFile, ncVarElemSize, element_size));

  int min, max;
  double avg;
  getMinMaxAvg(element_size, partitions, min, max, avg);

  logInfo() << "Min number of elements per partition:" << min;
  logInfo() << "Avg number of elements per partition" << avg;
  logInfo() << "Max number of elements per partition" << max;

  // Get boundary statistics
  int* boundary_size = new int[partitions];
  checkNcError(nc_get_var_int(ncFile, ncVarBndSize, boundary_size));

  getMinMaxAvg(boundary_size, partitions, min, max, avg);

  logInfo() << "Min number of MPI boundaries per partition:" << min;
  logInfo() << "Avg number of MPI boundaries per partition" << avg;
  logInfo() << "Max number of MPI boundaries per partition" << max;

  size_t max_boundaries;
  checkNcError(nc_inq_dimlen(ncFile, ncDimBnd, &max_boundaries));

  int* boundary_element_size = new int[partitions*max_boundaries];
  checkNcError(nc_get_var_int(ncFile, ncVarBndElemSize, boundary_element_size));

  int* boundary_element_sum = new int[partitions];
  for (unsigned int i = 0; i < partitions; i++) {
    boundary_element_sum[i] = 0;
    for (int j = 0; j < boundary_size[i]; j++)
      boundary_element_sum[i] += boundary_element_size[i*max_boundaries + j];
  }

  getMinMaxAvg(boundary_element_sum, partitions, min, max, avg);

  logInfo() << "Min number of MPI boundary faces per partition:" << min;
  logInfo() << "Avg number of MPI boundary faces per partition" << avg;
  logInfo() << "Max number of MPI boundary faces per partition" << max;


  delete [] boundary_size;
  delete [] boundary_element_size;
  delete [] boundary_element_sum;

  if (includeFault) {
    int* elemBoundaries = new int[maxElements*4];
    int* elemNeighborRanks = new int[maxElements*4];
    
    int* dr_sum = new int[partitions];
    int* dr_mpi_sum = new int[partitions];
    
    for (size_t i = 0; i < partitions; i++) {
      size_t start[3] = {i, 0, 0};
      size_t count[3] = {1, static_cast<size_t>(element_size[i]), 4};
      checkNcError(nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count, elemBoundaries));
      checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborRanks, start, count, elemNeighborRanks));
      
      int dr = 0;
      int dr_mpi = 0;
      
      #pragma omp parallel for reduction(+:dr,dr_mpi)
      for (int j = 0; j < element_size[i]*4; j++) {
        if (elemBoundaries[j] == 3) {
          dr++;
          if (static_cast<size_t>(elemNeighborRanks[j]) != i)
            dr_mpi++;
        }
      }
      
      dr_sum[i] = dr;
      dr_mpi_sum[i] = dr_mpi;
    }
    
    delete [] elemBoundaries;
    delete [] elemNeighborRanks;
    
    getMinMaxAvg(dr_sum, partitions, min, max, avg);
    
    logInfo() << "Min number of DR boundary faces per partition:" << min;
    logInfo() << "Avg number of DR boundary faces per partition" << avg;
    logInfo() << "Max number of DR boundary faces per partition" << max;
    
    getMinMaxAvg(dr_mpi_sum, partitions, min, max, avg);
    
    logInfo() << "Min number of DR MPI boundary faces per partition:" << min;
    logInfo() << "Avg number of DR MPI boundary faces per partition" << avg;
    logInfo() << "Max number of DR MPI boundary faces per partition" << max;
    
    delete [] dr_sum;
    delete [] dr_mpi_sum;
  }


  if( args.isSetAdditional("partitions") &&
    args.getAdditionalArgument<std::string>("partitions") == "yes" ) {
    logInfo() << "printing statistics per partition.";
    logInfo() << "dynamic rupture information is counted twice: as part of mpi_dr_faces and dynamic_rupture_faces.";
    std::cout << "partition,"
              << "elements,"
              << "neighboring_ranks,"
              << "max_mpi_faces_per_neighboring_rank,"
              << "interior_elements,"
              << "mpi_elements,"
              << "mpi_faces,"
              << "mpi_dr_faces,"
              << "dynamic_rupture_faces,"
              << "free_surface_faces,"
              << "transparent_faces,"
              << "mpi_copy_elements,"
              << "mpi_ghost_elements,"
              << "mpi_copy_dr_elements,"
              << "mpi_ghost_dr_elements"
              << std::endl;
    int* elemBoundaries = new int[maxElements*4];
    int* elemNeighborRanks = new int[maxElements*4];
    std::vector<RankDescriptor> rankInfo;

    for (size_t i = 0; i < partitions; i++) {
      RankDescriptor l_myRank;
      size_t start[3] = {i, 0, 0};
      size_t count[3] = {1, element_size[i], 4};
      checkNcError(nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count, elemBoundaries));
      checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborRanks, start, count, elemNeighborRanks));

      l_myRank.m_numberOfElements = element_size[i];

      // get information about neighboring partitions
      std::list<size_t> l_neighboringRanks( elemNeighborRanks,
                                            elemNeighborRanks+(element_size[i]*4) );
      // remove interior faces
      l_neighboringRanks.remove( i );

      // sort the list
      l_neighboringRanks.sort();

      unsigned long l_numberOfNeighboringMpiRanks = 0;
      unsigned long l_maximumNumberOfMpiFacesPerRank = 0;
      unsigned long l_currentNumberOfMpiFaces = 0;
      unsigned long l_currentRank = i;

      // get number of neighboring ranks and maximum number of faces to a single rank
      for ( std::list<size_t>::const_iterator l_iterator = l_neighboringRanks.begin(), l_end = l_neighboringRanks.end(); l_iterator != l_end; ++l_iterator) {
        if( l_currentRank != *l_iterator ) {
          // new rank
          l_currentRank = *l_iterator;
          l_currentNumberOfMpiFaces = 1;
          l_numberOfNeighboringMpiRanks++;
        }
        else {
          l_currentNumberOfMpiFaces++;
        }
      
        l_maximumNumberOfMpiFacesPerRank = std::max( l_currentNumberOfMpiFaces,
                                                     l_maximumNumberOfMpiFacesPerRank);
      }

      l_myRank.m_numberOfNeighboringMpiRanks = l_numberOfNeighboringMpiRanks;
      l_myRank.m_maximumNumberOfMpiFacesPerRank = l_maximumNumberOfMpiFacesPerRank;

      unsigned long l_numberOfInteriorElements = 0;
      unsigned long l_numberOfMpiElements = 0;
      unsigned long l_numberOfMpiCopyLayerElements = 0;
      unsigned long l_numberOfMpiCopyLayerDrElements = 0;
      unsigned long l_numberOfMpiFaces = 0;
      unsigned long l_numberOfMpiDynamicRuptureFaces = 0;
      unsigned long l_numberOfDynamicRuptureFaces = 0;
      unsigned long l_numberOfFreeSurfaceFaces = 0;
      unsigned long l_numberOfTransparentFaces = 0;
      long l_lastMpiNeighbor = -1;
      bool l_mpiElement = false;
      std::map<int, int> l_numberOfElementsInCopyLayerPerRank;
      std::map<int, int> l_currentNeighbors;
      l_numberOfElementsInCopyLayerPerRank.clear();
      l_currentNeighbors.clear();
      std::map<int, int> l_numberOfDrElementsInCopyLayerPerRank;
      std::map<int, int> l_currentDrNeighbors;
      l_numberOfDrElementsInCopyLayerPerRank.clear();
      l_currentDrNeighbors.clear();

      // iterate over faces
      for( int l_face = 0; l_face < element_size[i]*4; l_face++ ) {
        // count mpi-parts
        size_t l_currentNeighbor = static_cast<size_t>(elemNeighborRanks[l_face]);

        if (l_currentNeighbor != i) {
          l_mpiElement = true;
          l_numberOfMpiFaces++;

          if (elemBoundaries[l_face] == 3) {
            l_currentDrNeighbors[l_currentNeighbor] = 1;
            l_numberOfMpiDynamicRuptureFaces++;
          } else {
            l_currentNeighbors[l_currentNeighbor] = 1;
          }
        }

        // count boundary conditions
        if( elemBoundaries[l_face] == 3 ) {
          l_numberOfDynamicRuptureFaces++;
        }
        else if( elemBoundaries[l_face] == 1 ) {
          l_numberOfFreeSurfaceFaces++;
        }
        else if( elemBoundaries[l_face] == 5 ) {
          l_numberOfTransparentFaces++;
        }
        else {
          // unsupported face-type
          assert( false );
        }
        
        // update mpi-element count & reset mpi-element switch
        if( (l_face+1)%4 == 0 ) {
          if( l_mpiElement == false ) {
            l_numberOfInteriorElements++;
          }
          else {
            l_numberOfMpiElements++;
            // lets determine how often this element is in the copy layer
            l_numberOfMpiCopyLayerElements += l_currentNeighbors.size();
            for (std::map<int,int>::iterator it = l_currentNeighbors.begin(); it != l_currentNeighbors.end(); it++) {
              if (l_numberOfElementsInCopyLayerPerRank.find(it->first) != l_numberOfElementsInCopyLayerPerRank.end()) {
                l_numberOfElementsInCopyLayerPerRank[it->first]++;
              } else {
                l_numberOfElementsInCopyLayerPerRank[it->first] = 1;
              }
            }
            // lets determine how often this element is in the copy layer (DR)
            l_numberOfMpiCopyLayerDrElements += l_currentDrNeighbors.size();
            for (std::map<int,int>::iterator it = l_currentDrNeighbors.begin(); it != l_currentDrNeighbors.end(); it++) {
              if (l_numberOfDrElementsInCopyLayerPerRank.find(it->first) != l_numberOfDrElementsInCopyLayerPerRank.end()) {
                l_numberOfDrElementsInCopyLayerPerRank[it->first]++;
              } else {
                l_numberOfDrElementsInCopyLayerPerRank[it->first] = 1;
              }
            }
          }
          l_mpiElement = false;
          l_currentNeighbors.clear();
          l_currentDrNeighbors.clear();
        }
      }

      l_myRank.m_numberOfInteriorElements = l_numberOfInteriorElements;
      l_myRank.m_numberOfMpiElements = l_numberOfMpiElements;
      l_myRank.m_numberOfMpiFaces = l_numberOfMpiFaces;
      l_myRank.m_numberOfMpiCopyLayerElements = l_numberOfMpiCopyLayerElements;
      l_myRank.m_numberOfMpiCopyLayerDrElements = l_numberOfMpiCopyLayerDrElements;
      l_myRank.m_numberOfMpiDynamicRuptureFaces = l_numberOfMpiDynamicRuptureFaces;
      l_myRank.m_numberOfDynamicRuptureFaces = l_numberOfDynamicRuptureFaces;
      l_myRank.m_numberOfFreeSurfaceFaces = l_numberOfFreeSurfaceFaces;
      l_myRank.m_numberOfTransparentFaces = l_numberOfTransparentFaces;
      l_myRank.m_numberOfElementsInCopyLayerPerRank = l_numberOfElementsInCopyLayerPerRank;
      l_myRank.m_numberOfDrElementsInCopyLayerPerRank = l_numberOfDrElementsInCopyLayerPerRank;

      rankInfo.push_back(l_myRank);
    }

    for (size_t i = 0; i < partitions; i++) {
      RankDescriptor l_myRank = rankInfo[i];

      unsigned long l_numberOfMpiGhostLayerElements = 0;
      for (std::map<int,int>::iterator it = l_myRank.m_numberOfElementsInCopyLayerPerRank.begin(); it != l_myRank.m_numberOfElementsInCopyLayerPerRank.end(); it++) {
         l_numberOfMpiGhostLayerElements += rankInfo[it->first].m_numberOfElementsInCopyLayerPerRank[i];
      }

      unsigned long l_numberOfMpiGhostLayerDrElements = 0;
      for (std::map<int,int>::iterator it = l_myRank.m_numberOfDrElementsInCopyLayerPerRank.begin(); it != l_myRank.m_numberOfDrElementsInCopyLayerPerRank.end(); it++) {
         l_numberOfMpiGhostLayerDrElements += rankInfo[it->first].m_numberOfDrElementsInCopyLayerPerRank[i];
      }

      // print statistics for this partition
      std::cout << i << ","
                << l_myRank.m_numberOfElements << ","
                << l_myRank.m_numberOfNeighboringMpiRanks << ","
                << l_myRank.m_maximumNumberOfMpiFacesPerRank << "," 
                << l_myRank.m_numberOfInteriorElements << ","
                << l_myRank.m_numberOfMpiElements << ","
                << l_myRank.m_numberOfMpiFaces << ","
                << l_myRank.m_numberOfMpiDynamicRuptureFaces << ","
                << l_myRank.m_numberOfDynamicRuptureFaces << ","
                << l_myRank.m_numberOfFreeSurfaceFaces << ","
                << l_myRank.m_numberOfTransparentFaces << ","
                << l_myRank.m_numberOfMpiCopyLayerElements << ","
                << l_numberOfMpiGhostLayerElements << ","
                << l_myRank.m_numberOfMpiCopyLayerDrElements << ","
                << l_numberOfMpiGhostLayerDrElements
                << std::endl;
    }
  }
  
  if (includeMatrices) {
    int* elemBoundaries = new int[maxElements*4];
    int* elemNeighborSides = new int[maxElements*4];
    int* elemSideOrientations = new int[maxElements*4];
    int* elemNeighborRanks = new int[maxElements*4];
    
    unsigned long int* localMatrices = new unsigned long int[omp_get_max_threads()*4];
    memset(localMatrices, 0, sizeof(unsigned long int)*4*omp_get_max_threads());
    unsigned long int* fluxMatrices = new unsigned long int[omp_get_max_threads()*4*4*3];
    memset(fluxMatrices, 0, sizeof(unsigned long int)*4*4*3*omp_get_max_threads());
    unsigned long int* boundaries = new unsigned long int[omp_get_max_threads()*3];
    memset(boundaries, 0, sizeof(unsigned long int)*3*omp_get_max_threads());
    unsigned long int mpi_dr_boundary = 0;
    
    for (size_t i = 0; i < partitions; i++) {
      size_t start[3] = {i, 0, 0};
      size_t count[3] = {1, static_cast<size_t>(element_size[i]), 4};
      
      checkNcError(nc_get_vara_int(ncFile, ncVarElemBoundaries, start, count, elemBoundaries));
      checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighorSides, start, count, elemNeighborSides));
      checkNcError(nc_get_vara_int(ncFile, ncVarElemSideOrientations, start, count, elemSideOrientations));
      checkNcError(nc_get_vara_int(ncFile, ncVarElemNeighborRanks, start, count, elemNeighborRanks));
      
      #pragma omp parallel for reduction(+:mpi_dr_boundary)
      for (int j = 0; j < element_size[i]; j++) {
        for (int k = 0; k < 4; k++) {
          switch (elemBoundaries[j*4 + k]) {
          case 0:
            fluxMatrices[((omp_get_thread_num()*4 + k)*4 + elemNeighborSides[j*4 + k])*3 + elemSideOrientations[j*4 + k]]++;
            localMatrices[omp_get_thread_num()*4 + k]++;
            break;
          case 1:
            boundaries[omp_get_thread_num()*3]++;
            localMatrices[omp_get_thread_num()*4 + k] += 2;
            break;
          case 3:
            boundaries[omp_get_thread_num()*3 + 1]++;
            if (elemNeighborRanks[j*4 + k] != static_cast<int>(i))
              mpi_dr_boundary++;
            break;
          case 5:
            boundaries[omp_get_thread_num()*3 + 2]++;
            localMatrices[omp_get_thread_num()*4 + k]++;
            break;
          case 6:
            fluxMatrices[((omp_get_thread_num()*4 + k)*4 + elemNeighborSides[j*4 + k])*3 + elemSideOrientations[j*4 + k]]++;
            localMatrices[omp_get_thread_num()*4 + k]++;
            break;
          default:
            logError() << "Unknown boundaries condition" << elemBoundaries[j*4 + k] << "found";
          }
        }
      }
    }
    
    for (int i = 1; i < omp_get_max_threads(); i++) {
      for (int j = 0; j < 3; j++) {
        boundaries[j] += boundaries[i*3 + j];
      }
      
      for (int j = 0; j < 4; j++) {
        localMatrices[j] += localMatrices[i*4 + j];
      }
      
      for (int j = 0; j < 4*4*3; j++) {
        fluxMatrices[j] += fluxMatrices[i*4*4*3 + j];
      }
    }
    
    logInfo() << "# Boundary faces:";
    std::cout << "1:" << boundaries[0] << std::endl;
    std::cout << "3:" << boundaries[1] << std::endl;
    std::cout << "5:" << boundaries[2] << std::endl;
    
    logInfo() << "# DR MPI Boundary faces:";
    std::cout << mpi_dr_boundary << std::endl;
    
    logInfo() << "Local Matrices:";
    for (int i = 0; i < 4; i++)
      std::cout << i << ',' << localMatrices[i] << std::endl;
    
    logInfo() << "Local Face,Neighbor Face,Orientation";
    for (int i = 0; i < 4; i++)
      for (int j = 0; j < 4; j++)
        for (int k = 0; k < 3; k++)
          std::cout << i << ',' << j << ',' << k << ',' << fluxMatrices[(i*4 + j)*3 + k] << std::endl;
    
    
    delete [] localMatrices;
  }
  
  delete [] element_size;

  checkNcError(nc_close(ncFile));
}
