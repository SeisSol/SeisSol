/**
 * @file
 * This file is part of SeisSol.
 *
 * @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
 *
 * @section LICENSE
 * Copyright (c) 2017, SeisSol Group
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

#include <algorithm>
#include <cassert>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <array> // TODO(Lukas) Remove

#include "PUML/PUML.h"
#include "PUML/PartitionMetis.h"
#include "PUML/Downward.h"
#include "PUML/Neighbor.h"

#include "PUMLReader.h"
#include "Monitoring/instrumentation.fpp"

#include "Initializer/time_stepping/LtsWeights.h"

#include <hdf5.h>
#include <sstream>
#include <fstream>


class GlobalFaceSorter
{
private:
	const PUML::TETPUML &m_puml;

public:
	GlobalFaceSorter(const PUML::TETPUML &puml)
		: m_puml(puml)
	{ }

	bool operator()(unsigned int a, unsigned int b) const
	{
		return m_puml.faces()[a].gid() < m_puml.faces()[b].gid();
	}
};

/**
 * @todo Cleanup this code
 */
seissol::PUMLReader::PUMLReader(const char *meshFile, initializers::time_stepping::LtsWeights* ltsWeights, double tpwgt, bool readPartitionFromFile)
	: MeshReader(MPI::mpi.rank())
{
	PUML::TETPUML puml;
	puml.setComm(MPI::mpi.comm());

	read(puml, meshFile);
	logInfo(MPI::mpi.rank()) << "Read mesh.";
  
	if (ltsWeights != nullptr) {
		generatePUML(puml);
		ltsWeights->computeWeights(puml);
	}
	partition(puml, ltsWeights, tpwgt, meshFile, readPartitionFromFile);

	generatePUML(puml);

	getMesh(puml);
}

void seissol::PUMLReader::read(PUML::TETPUML &puml, const char* meshFile)
{
	SCOREP_USER_REGION("PUMLReader_read", SCOREP_USER_REGION_TYPE_FUNCTION);

	std::string file(meshFile);

	puml.open((file + ":/connect").c_str(), (file + ":/geometry").c_str());
	puml.addData((file + ":/group").c_str(), PUML::CELL);
	puml.addData((file + ":/boundary").c_str(), PUML::CELL);
}

int seissol::PUMLReader::readPartition(PUML::TETPUML &puml, int* partition, const char* meshFile)
{
	/*
	write the partionning array to an hdf5 file using parallel access
	see https://support.hdfgroup.org/ftp/HDF5/examples/parallel/coll_test.c for more info about the hdf5 functions
	*/
	SCOREP_USER_REGION("PUMLReader_readPartition", SCOREP_USER_REGION_TYPE_FUNCTION);
	const int rank = seissol::MPI::mpi.rank();
	const int nrank = seissol::MPI::mpi.size();
	int nPartitionCells = puml.numOriginalCells();

	/* 
	 Gather number of cells in each nodes. This is necessary to be able to write the data in the correct location
	*/
	int *num_cells = new int[nrank];
	int *offsets = new int[nrank];
	MPI_Allgather(&nPartitionCells, 1, MPI_INT, num_cells, 1, MPI_INT, MPI_COMM_WORLD);

	offsets[0] = 0;
	for (int rk = 1; rk < nrank; ++rk) {
		offsets[rk] = offsets[rk-1] + num_cells[rk-1];
	}
	const hsize_t dimMem[] = {static_cast<hsize_t>(nPartitionCells)};

	/* 
	 Open file and dataset 
	*/
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, seissol::MPI::mpi.comm(), info);

	std::ostringstream os;
	os << meshFile<<"_partitions_o" <<CONVERGENCE_ORDER<<"_n"<< nrank << ".h5";
	std::string fname = os.str();

	std::ifstream ifile(fname.c_str());
	if (!ifile) { 
		logInfo(rank) <<fname.c_str()<<"does not exist";
		return -1;
	}

	hid_t file = H5Fopen(fname.c_str(), H5F_ACC_RDONLY, plist_id);
	H5Pclose(plist_id);

	hid_t dataset = H5Dopen2(file, "/partition", H5P_DEFAULT);
	/* 
	 Create memspace (portion of filespace) and read collectively the data
	*/
	hid_t memspace = H5Screate_simple(1, dimMem, NULL);
	hid_t filespace = H5Dget_space(dataset);

	hsize_t start[] = {static_cast<hsize_t>(offsets[rank])};
	hsize_t count[] = {static_cast<hsize_t>(nPartitionCells)};
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0L, count, 0L);

	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	int status = H5Dread(dataset, H5T_NATIVE_INT, memspace, filespace, plist_id, partition);

	if (status<0)
		logError() << "An error occured when reading the partitionning with HDF5";
	H5Dclose(dataset);
	H5Fclose(file);

	logInfo(rank)<<"partitionning was read successfully from "<<fname.c_str();
	return 0;
}


void seissol::PUMLReader::writePartition(PUML::TETPUML &puml, int* partition, const char *meshFile)
{
	/*
	write the partionning array to an hdf5 file using parallel access
	see https://support.hdfgroup.org/ftp/HDF5/examples/parallel/coll_test.c for more info about the hdf5 functions
	*/
	SCOREP_USER_REGION("PUMLReader_writePartition", SCOREP_USER_REGION_TYPE_FUNCTION);
	const int rank = seissol::MPI::mpi.rank();
	const int nrank = seissol::MPI::mpi.size();
	int nPartitionCells = puml.numOriginalCells();

	/* 
	 Gather number of cells in each nodes. This is necessary to be able to write the data in the correct location
	*/
	int *num_cells = new int[nrank];
	int *offsets = new int[nrank];
	MPI_Allgather(&nPartitionCells, 1, MPI_INT, num_cells, 1, MPI_INT, MPI_COMM_WORLD);

	offsets[0] = 0;
	for (int rk = 1; rk < nrank; ++rk) {
		offsets[rk] = offsets[rk-1] + num_cells[rk-1];
	}
	int nCells = offsets[nrank-1]+num_cells[nrank-1];

	const hsize_t dim[] = {static_cast<hsize_t>(nCells)};
	const hsize_t dimMem[] = {static_cast<hsize_t>(nPartitionCells)};

	/* 
	 Create file and file space 
	*/
	MPI_Info info  = MPI_INFO_NULL;
	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, seissol::MPI::mpi.comm(), info);

	std::ostringstream os;
	os << meshFile<<"_partitions_o" <<CONVERGENCE_ORDER<<"_n"<< nrank << ".h5";
	std::string fname = os.str();

	hid_t file = H5Fcreate(fname.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);
	H5Pclose(plist_id);

	hid_t filespace = H5Screate_simple(1, dim, NULL);
	hid_t dataset = H5Dcreate(file, "/partition", H5T_NATIVE_INT, filespace, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	H5Sclose(filespace);

	/* 
	 Create memspace (portion of filespace) and write collectively the data
	*/
	hid_t memspace = H5Screate_simple(1, dimMem, NULL);
	filespace = H5Dget_space(dataset);

	hsize_t start[] = {static_cast<hsize_t>(offsets[rank])};
	hsize_t count[] = {static_cast<hsize_t>(nPartitionCells)};
	H5Sselect_hyperslab(filespace, H5S_SELECT_SET, start, 0L, count, 0L);

	plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);

	int status = H5Dwrite(dataset, H5T_NATIVE_INT, memspace, filespace, plist_id, partition);

	if (status<0)
		logError() << "An error occured when writing the partitionning with HDF5";
	H5Dclose(dataset);
	H5Fclose(file);
}

void seissol::PUMLReader::partition(  PUML::TETPUML &puml,
                                      initializers::time_stepping::LtsWeights* ltsWeights,
                                      double tpwgt,
                                      const char *meshFile,
                                      bool readPartitionFromFile  )
{
	SCOREP_USER_REGION("PUMLReader_partition", SCOREP_USER_REGION_TYPE_FUNCTION);

	int* partition = new int[puml.numOriginalCells()];

  auto partitionMetis = [&] {
    PUML::TETPartitionMetis metis(puml.originalCells(), puml.numOriginalCells());
#ifdef USE_MPI
    double* nodeWeights = new double[seissol::MPI::mpi.size()];
    MPI_Allgather(&tpwgt, 1, MPI_DOUBLE, nodeWeights, 1, MPI_DOUBLE, seissol::MPI::mpi.comm());
    double sum = 0.0;
    for (int rk = 0; rk < seissol::MPI::mpi.size(); ++rk) {
     sum += nodeWeights[rk];
    }
    for (int rk = 0; rk < seissol::MPI::mpi.size(); ++rk) {
     nodeWeights[rk] /= sum;
    }
#else
    tpwgt = 1.0;
    double* nodeWeights = &tpwgt;
#endif

    metis.partition(partition, ltsWeights->vertexWeights(), ltsWeights->nWeightsPerVertex(), nodeWeights, 1.01);

#ifdef USE_MPI
    delete[] nodeWeights;
#endif
  };

  if (readPartitionFromFile) {
    int status = readPartition(puml, &partition[0], meshFile);
    if (status < 0) {
      partitionMetis();
      writePartition(puml, partition, meshFile);
    }
  } else {
    partitionMetis();
  }

	puml.partition(partition);
	delete [] partition;
}

void seissol::PUMLReader::generatePUML(PUML::TETPUML &puml)
{
	SCOREP_USER_REGION("PUMLReader_generate", SCOREP_USER_REGION_TYPE_FUNCTION);

	puml.generateMesh();
}

void seissol::PUMLReader::getMesh(const PUML::TETPUML &puml)
{
	SCOREP_USER_REGION("PUMLReader_getmesh", SCOREP_USER_REGION_TYPE_FUNCTION);

	const int rank = MPI::mpi.rank();

	const std::vector<PUML::TETPUML::cell_t> &cells = puml.cells();
	const std::vector<PUML::TETPUML::face_t> &faces = puml.faces();
	const std::vector<PUML::TETPUML::vertex_t> &vertices = puml.vertices();

	const int* material = puml.cellData(0);
	const int* boundaryCond = puml.cellData(1);

	std::unordered_map<int, std::vector<unsigned int> > neighborInfo; // List of shared local face ids

	// Compute everything local
	m_elements.resize(cells.size());

	// First check, whether we have a periodic boundary condition.
	// In this loop, we also do some preliminary work
	// such as setting the local cell ids and storing the vertices information.
	bool foundPeriodicBc = false;
	for (unsigned int i = 0; i < cells.size(); ++i) {
	  m_elements[i].localId = i;

	  // Vertices
	  PUML::Downward::vertices(puml, cells[i], reinterpret_cast<unsigned int*>(m_elements[i].vertices));

	  for (unsigned int j = 0; j < 4; ++j) {
	    int bcCurrentFace = (boundaryCond[i] >> (j*8)) & 0xFF;
	    foundPeriodicBc |= bcCurrentFace == 6;
	  }
	}
	if (MPI::mpi.size() > 1 && foundPeriodicBc) {
	  logError() << "Periodic boundary conditions are currently not supported for MPI+HDF5.";
	}
	
	std::vector<Plane> planes;
	std::vector<std::unordered_set<PeriodicVertex>> planes2Vertices;
	// TODO(Lukas) Those should be arrays.
	std::unordered_map<size_t, size_t> faces2Planes;
	std::unordered_map<size_t, std::unordered_set<int>> vertices2Cells;
	// TODO(Lukas) Name for this... Vertex -> Mirrored vertex
	std::unordered_map<size_t, std::unordered_set<size_t>> vertices2Vertices;


	if (foundPeriodicBc) {
	  // We need to match vertices of opposite boundaries.
	  // For simplification, we sort the boundaries by planes.
	  for (unsigned int i = 0; i < cells.size(); i++) {
	    //m_elements[i].localId = i;

	    // Vertices
	    //PUML::Downward::vertices(puml, cells[i], reinterpret_cast<unsigned int*>(m_elements[i].vertices));

	    // Neighbor information
	    unsigned int faceids[4];
	    PUML::Downward::faces(puml, cells[i], faceids);
	    int neighbors[4];
	    PUML::Neighbor::face(puml, i, neighbors);
	    for (unsigned int j = 0; j < 4; j++) {
	      int bcCurrentFace = (boundaryCond[i] >> (j*8)) & 0xFF;
	      if (bcCurrentFace == 6) {
		auto &faceNodesIds = MeshTools::FACE2NODES[FACE_PUML2SEISSOL[j]];
		auto faceNodesCoords = std::array<const double*, 3>{};

		for (int k = 0; k < faceNodesCoords.size(); ++k) {
		  const auto curVertexId = m_elements[i].vertices[faceNodesIds[k]];
		  vertices2Cells[curVertexId].insert(i);
		  faceNodesCoords[k] = vertices[curVertexId].coordinate();
		}
		bool found = false;
		size_t idx = 0;
		for (const auto& plane : planes) {
		  ++idx;

		  // We need to check if all nodes of a face are contained in the plane.
		  // This is needed because some vertices can be at two planes at the same
		  // time but no face can be in two planes. 
		  bool containsAllPoints = true;
		  for (int z = 0; z < faceNodesCoords.size(); ++z) {
		    containsAllPoints &=
		      plane.containsPoint(faceNodesCoords[z]);
		  }
		  if (containsAllPoints) {
		    found = true;
		    break;
		  }
		}
		if (!found) {
		  // If the node isn't on an already existing plane, we need to create
		  // a new one.
		  Plane plane;
		  // Compute an orthonormal coordinate system for each plane.
		  // Basis is [normal, tangent1, tangent2].
		  VrtxCoords ac;
		  MeshTools::sub(faceNodesCoords[1], faceNodesCoords[0], plane.tangent1);
		  MeshTools::sub(faceNodesCoords[2], faceNodesCoords[0], ac);
		  MeshTools::cross(plane.tangent1, ac, plane.normal);
		  // normalize normal
		  MeshTools::mul(plane.normal,
				 1/(MeshTools::norm(plane.normal)),
				 plane.normal);
		  // take absolut value of normal to standarize direction
		  // This way, opposite planes have the same coordinate system.
		  for (int k = 0; k < 3; ++k) {
		    plane.normal[k] = std::abs(plane.normal[k]);
		  }
		  if (std::abs(plane.normal[0]) > std::abs(plane.normal[1])) {
		    plane.tangent1[0] = plane.normal[2];
		    plane.tangent1[1] = 0.0;
		    plane.tangent1[2] = -plane.normal[0];
		  } else {
		    plane.tangent1[0] = 0.0;
		    plane.tangent1[1] = -plane.normal[2];
		    plane.tangent1[2] = plane.normal[1];
		  }
		  MeshTools::cross(plane.normal, plane.tangent1, plane.tangent2);

		  for (int k = 0; k < 3; ++k) {
		    // A point on the plane is needed for the plane equations.
		    // This also differentiates opposite planes.
		    plane.point[k] = faceNodesCoords[0][k];
		  }
		  planes.push_back(plane);
		  planes2Vertices.push_back(std::unordered_set<PeriodicVertex>());
		  ++idx;
		}
		//logInfo(0) << idx << planes.size() << planes2Vertices.size() << faceNodesCoords.size();
		for (int k = 0; k < faceNodesCoords.size(); ++k) {
		  const auto curVertexId = m_elements[i].vertices[faceNodesIds[k]];
		  auto pVert = PeriodicVertex{faceNodesCoords[k],
					      planes[idx-1]};
		  pVert.vertexId = curVertexId;
		  planes2Vertices[idx-1].insert(pVert);
		  faces2Planes.insert({faceids[j], idx-1});
		}
	      }
	    }
	  }

	}

	logInfo(rank) << "Number of planes =" << planes.size();
	for (const auto& plane : planes) {
	  const auto &n = plane.normal;
	  const auto &p = plane.point;
	  const auto &t1 = plane.tangent1;
	  const auto &t2 = plane.tangent2;
	  std::cout << "Normal:\t" << n[0]
		    << ",\t" << n[1]
		    << ",\t" << n[2] << std::endl;
	  std::cout << "Tangent1:\t" << t1[0]
		    << ",\t" << t1[1]
		    << ",\t" << t1[2] << std::endl;
	  std::cout << "Tangent2: " << t2[0]
		    << ",\t" << t2[1]
		    << ",\t" << t2[2] << std::endl;
	  std::cout << "Point:\t" << p[0]
		    << ",\t" << p[1]
		    << ",\t" << p[2] << std::endl;
	  std::cout << std::endl;
	}

	for (unsigned int i = 0; i < cells.size(); i++) {
	  m_elements[i].localId = i;

	  // Neighbor information
	  unsigned int faceids[4];
	  PUML::Downward::faces(puml, cells[i], faceids);
	  int neighbors[4];
	  PUML::Neighbor::face(puml, i, neighbors);
	  for (unsigned int j = 0; j < 4; j++) {
	    auto &faceNodesIds = MeshTools::FACE2NODES[FACE_PUML2SEISSOL[j]];
	    int bcCurrentFace = (boundaryCond[i] >> (j*8)) & 0xFF;
	    if (bcCurrentFace == 6) {
	      int result = -1;
	      std::array<const double*, 3> faceNodesCoords;
	      getFaceNodes(m_elements[i],
			   vertices,
			   j,
			   faceNodesCoords);
	      const auto facePlaneIdx = faces2Planes[faceids[j]];
	      const auto &facePlane = planes[facePlaneIdx];

	      for (int k = 0; k < planes.size(); ++k) {
		if (k == facePlaneIdx) continue;
			    
		// Move current plane s.t. it cuts the other plane at its point.
		// Remember: opposite planes have different offsets but same coordinate system!
		const auto &curPlane = planes[k];
		auto movedPlane = planes[k];
		std::copy_n(facePlane.point, 3, movedPlane.point);
		bool containsAllPoints = true;
		for (int l = 0; l < 3; ++l) {
		  containsAllPoints &=
		    movedPlane.containsPoint(faceNodesCoords[l]);
		}
		if (containsAllPoints) {
		  // We count how to which face most points belong.
		  auto faceCounter = std::unordered_map<size_t, int>{};
		  // Iterate over all points of the edge.
		  for (int l = 0; l < 3; ++l) {
		    const auto* curCoords = faceNodesCoords[l];
		    const auto curVertex = PeriodicVertex{faceNodesCoords[l],
							  planes[k]};
		    const auto curVertexId = m_elements[i].vertices[faceNodesIds[l]];
		    // We compare with all other points in the plane.
		    // Note that this is O(N^2); this doesn't really matter,
		    // because normally N << 100 for a single node.
		    for (const auto &other : planes2Vertices[k]) {
		      // for (int m = 0; m < planes2Vertices[k].size(); ++m) {
		      //const auto &other = planes2Vertices[k][m];
		      if (curVertex.isApproxEq(other)) {
			vertices2Vertices[curVertexId].insert(other.vertexId);
			assert(vertices2Cells.find(curVertexId) !=
			       vertices2Cells.end());
			for (const auto otherCellId : vertices2Cells[other.vertexId]) {
			  faceCounter[otherCellId] += 1;
			}
		      }
		    }
		  }
		  // Check which cell was mentioned most often
		  // This is the matching cell, as most vertices belong to it.
		  auto bestCellId = -1;
		  auto maxCounter = -1;
		  for (const auto& it : faceCounter) {
		    if (it.first != i && it.second > maxCounter) {
		      maxCounter = it.second;
		      bestCellId = it.first;
		    }
		  }
		  //logInfo(0) << faceCounter.size()
		  //   << vertices2Cells.size()
		  //   << bestCellId
		  //   << maxCounter;
		  result = bestCellId;
		}

	      }

	      if (result >= 0) {
		neighbors[j] = result;
		//m_elements[i].neighbors[FACE_PUML2SEISSOL[j]] = neighbors[j];
	      } else {
		logError() << "Did not find periodic counterpart!";
	      }
	    }
	    // Set the neighbor directly in seissol data structure.
	    m_elements[i].neighbors[FACE_PUML2SEISSOL[j]] = neighbors[j];
	  }
	  m_elements[i].material = material[i];
	}

	for (unsigned int i = 0; i < cells.size(); i++) {
		// Neighbor information
		unsigned int faceids[4];
		PUML::Downward::faces(puml, cells[i], faceids);
		int neighbors[4];
		PUML::Neighbor::face(puml, i, neighbors);
		for (unsigned int j = 0; j < 4; j++) {
		  int bcCurrentFace = (boundaryCond[i] >> (j*8)) & 0xFF;
		  if (bcCurrentFace == 6) {
		    // For periodic bcs, the array is already correct!
		    neighbors[j] = m_elements[i].neighbors[FACE_PUML2SEISSOL[j]];
		  }

		  if (neighbors[j] < 0) {
		    // No local neighbour found.
		    m_elements[i].neighbors[FACE_PUML2SEISSOL[j]] = cells.size();

		    // isShared is true if face belongs to more than one rank.
		    if (!faces[faceids[j]].isShared()) {
		      // Boundary sides
		      m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = rank;
		    } else {
		      // MPI Boundary
		      neighborInfo[faces[faceids[j]].shared()[0]].push_back(faceids[j]);

		      m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = faces[faceids[j]].shared()[0];
		    }
		  } else {
		    assert(neighbors[j] >= 0 && static_cast<unsigned>(neighbors[j]) < cells.size());

		    m_elements[i].neighbors[FACE_PUML2SEISSOL[j]] = neighbors[j];

		    if (bcCurrentFace == 6) {
		      // TODO(Lukas) Always 
		      const auto& neighbor = m_elements[neighbors[j]];
		      const auto* back = std::find(neighbor.neighbors, neighbor.neighbors + 4, i);
		      assert(back < neighbor.neighbors + 4);
		      /*
		      for (int z = 0; z < 4; ++z) {
			if (neighbor.neighbors[z] == i) {
			  //m_elements[i].neighborSides[FACE_PUML2SEISSOL[j]] = z; //FACE_PUML2SEISSOL[z];
			  //logInfo(0) << "Matched face " << j
			  //	     << "with face " << z;
			  break;
			}
		      }
		      */
		      // The rhs of this assignment is already a seissol face id
		      m_elements[i].neighborSides[FACE_PUML2SEISSOL[j]] = back - neighbor.neighbors;
		    } else {
		      int nfaces[4];
		      PUML::Neighbor::face(puml, neighbors[j], nfaces);
		      int* back = std::find(nfaces, nfaces+4, i);
		      assert(back < nfaces+4);
		      
		      m_elements[i].neighborSides[FACE_PUML2SEISSOL[j]] = FACE_PUML2SEISSOL[back-nfaces];
		    }

		    const unsigned int firstVertex = m_elements[i].vertices[FIRST_FACE_VERTEX[FACE_PUML2SEISSOL[j]]];

		    unsigned int nvertices[4];
		    PUML::Downward::vertices(puml, cells[neighbors[j]], nvertices);
		    
		    if (bcCurrentFace == 6) {
		      assert(vertices2Vertices.find(firstVertex) !=
			     vertices2Vertices.end());
		      // Find all vertices that are matched with ours.
		      auto& matchedVertices = vertices2Vertices[firstVertex];
		      assert(matchedVertices.size() > 0);
		      int matchedVertex = -1;
		      for (auto m : matchedVertices) {
			const auto &curCells = vertices2Cells[m];
			// Check if vertex is in considered cell.
			auto it = curCells.find(neighbors[j]);
			if (it != curCells.end()) {
			  //logInfo(0) << "Matched with vId ="  << matchedVertex;
			  matchedVertex = m;
			  break;
			}
		      }
		      assert(matchedVertex >= 0);
		      for (int z = 0; z < 4; ++z) {
			//logInfo(0) << "Back vertex search (matched):"
			//   << matchedVertex << j << z << nvertices[z];
			if (nvertices[z] == matchedVertex) {
			  // TODO(Lukas) Do i need puml2seissol for z as well?
			  m_elements[i].sideOrientations[FACE_PUML2SEISSOL[j]]
			    = FACEVERTEX2ORIENTATION[m_elements[i].neighborSides[FACE_PUML2SEISSOL[j]]][z];
			  break;
			}
		      }
		      //logInfo(0) << "Search resulted in orientation of"
		      //	 << m_elements[i].sideOrientations[FACE_PUML2SEISSOL[j]];
		    } else {
		      for (int z = 0; z < 4; ++z) {
			//logInfo(0) << "Back vertex search:"
			//   << firstVertex << j << z << nvertices[z];
		      }
		      unsigned int* neighborFirstVertex = std::find(nvertices, nvertices+4, firstVertex);
		      assert(neighborFirstVertex < nvertices+4);
		      m_elements[i].sideOrientations[FACE_PUML2SEISSOL[j]]
			= FACEVERTEX2ORIENTATION[m_elements[i].neighborSides[FACE_PUML2SEISSOL[j]]][neighborFirstVertex-nvertices];
		    }

		    assert(m_elements[i].sideOrientations[FACE_PUML2SEISSOL[j]] >= 0);

		    m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = rank;
		  }

		  int faultTag = bcCurrentFace;
		  // Boundary values higher than 6 indicate dynamic rupture.
		  if (bcCurrentFace > 6) {
		    // Dynamic rupture is tagged as 3.
		    bcCurrentFace = 3;
		  }
			
		  m_elements[i].boundaries[FACE_PUML2SEISSOL[j]] = bcCurrentFace;
		  m_elements[i].faultTags[FACE_PUML2SEISSOL[j]] = faultTag;
		  m_elements[i].mpiIndices[FACE_PUML2SEISSOL[j]] = 0;
		}
	}

	// Exchange ghost layer information and generate neighbor list
	char** copySide = new char*[neighborInfo.size()];
	char** ghostSide = new char*[neighborInfo.size()];
	unsigned long** copyFirstVertex = new unsigned long*[neighborInfo.size()];
	unsigned long** ghostFirstVertex = new unsigned long*[neighborInfo.size()];

	MPI_Request* requests = new MPI_Request[neighborInfo.size() * 4];

	std::unordered_set<unsigned int> t;
	unsigned int sum = 0;
	unsigned int k = 0;
	for (std::unordered_map<int, std::vector<unsigned int> >::iterator it = neighborInfo.begin();
			it != neighborInfo.end(); ++it, k++) {
		// Need to sort the neighborInfo vectors once
		GlobalFaceSorter faceSorter(puml);
		std::sort(it->second.begin(), it->second.end(), faceSorter);

		t.insert(it->second.begin(), it->second.end());
		sum += it->second.size();

		// Create MPI neighbor list
		addMPINeighor(puml, it->first, it->second);

		copySide[k] = new char[it->second.size()];
		ghostSide[k] = new char[it->second.size()];
		copyFirstVertex[k] = new unsigned long[it->second.size()];
		ghostFirstVertex[k] = new unsigned long[it->second.size()];

		MPI_Irecv(ghostSide[k], it->second.size(), MPI_CHAR, it->first, 0, MPI::mpi.comm(), &requests[k]);
		MPI_Irecv(ghostFirstVertex[k], it->second.size(), MPI_UNSIGNED_LONG, it->first, 0, MPI::mpi.comm(),
			&requests[neighborInfo.size() + k]);

		// Neighbor side
		for (unsigned int i = 0; i < it->second.size(); i++) {
			// The side of boundary
			int cellIds[2];
			PUML::Upward::cells(puml, faces[it->second[i]], cellIds);
			int side = PUML::Downward::faceSide(puml, cells[cellIds[0]], it->second[i]);
			assert(side >= 0 && side < 4);
			copySide[k][i] = side;

			// First vertex of the face on the boundary
			const unsigned int firstVertex
				= m_elements[cellIds[0]].vertices[FIRST_FACE_VERTEX[FACE_PUML2SEISSOL[side]]];
			copyFirstVertex[k][i] = vertices[firstVertex].gid();

			// Set the MPI index
			assert(m_elements[cellIds[0]].mpiIndices[FACE_PUML2SEISSOL[side]] == 0);
			m_elements[cellIds[0]].mpiIndices[FACE_PUML2SEISSOL[side]] = i;
		}

		MPI_Isend(copySide[k], it->second.size(), MPI_CHAR, it->first, 0, MPI::mpi.comm(),
			&requests[neighborInfo.size() * 2 + k]);
		MPI_Isend(copyFirstVertex[k], it->second.size(), MPI_UNSIGNED_LONG, it->first, 0, MPI::mpi.comm(),
			&requests[neighborInfo.size() * 3 + k]);
	}
	assert(t.size() == sum);

	MPI_Waitall(neighborInfo.size() * 4, requests, MPI_STATUSES_IGNORE);

	k = 0;
	for (std::unordered_map<int, std::vector<unsigned int> >::const_iterator it = neighborInfo.begin();
			it != neighborInfo.end(); ++it, k++) {
		for (unsigned int i = 0; i < it->second.size(); i++) {
			// Set neighbor side
			int cellIds[2];
			PUML::Upward::cells(puml, faces[it->second[i]], cellIds);
			assert(cellIds[1] < 0);

			int side = copySide[k][i];
			int gSide = ghostSide[k][i];
			m_elements[cellIds[0]].neighborSides[FACE_PUML2SEISSOL[side]] = FACE_PUML2SEISSOL[gSide];

			// Set side Orientation
			unsigned long nvertices[4];
			PUML::Downward::gvertices(puml, cells[cellIds[0]], nvertices);

			unsigned long* localFirstVertex = std::find(nvertices, nvertices+4, ghostFirstVertex[k][i]);
			assert(localFirstVertex != nvertices+4);

			m_elements[cellIds[0]].sideOrientations[FACE_PUML2SEISSOL[side]]
				= FACEVERTEX2ORIENTATION[FACE_PUML2SEISSOL[side]][localFirstVertex-nvertices];
			assert(m_elements[cellIds[0]].sideOrientations[FACE_PUML2SEISSOL[side]] >= 0);
		}

		delete [] copySide[k];
		delete [] ghostSide[k];
		delete [] copyFirstVertex[k];
		delete [] ghostFirstVertex[k];
	}

	delete [] copySide;
	delete [] ghostSide;
	delete [] copyFirstVertex;
	delete [] ghostFirstVertex;
	delete [] requests;

	// Set vertices
	m_vertices.resize(vertices.size());
	for (unsigned int i = 0; i < vertices.size(); i++) {
		memcpy(m_vertices[i].coords, vertices[i].coordinate(), 3*sizeof(double));

		PUML::Upward::cells(puml, vertices[i], m_vertices[i].elements);
	}

	//std::abort();

}

void seissol::PUMLReader::addMPINeighor(const PUML::TETPUML &puml, int rank, const std::vector<unsigned int> &faces)
{
	unsigned int id = m_MPINeighbors.size();
	MPINeighbor& neighbor = m_MPINeighbors[rank];

	neighbor.localID = id;
	neighbor.elements.resize(faces.size());

	for (unsigned int i = 0; i < faces.size(); i++) {
		int cellIds[2];
		PUML::Upward::cells(puml, puml.faces()[faces[i]], cellIds);

		neighbor.elements[i].localElement = cellIds[0];
	}
}


int seissol::PUMLReader::FACE_PUML2SEISSOL[4] = {0, 1, 3, 2};

int seissol::PUMLReader::FACEVERTEX2ORIENTATION[4][4] = {
	{0, 2, 1, -1}, {0, 1, -1, 2}, {0, -1, 2, 1}, {-1, 0, 1, 2}
};

int seissol::PUMLReader::FIRST_FACE_VERTEX[4] = {0, 0, 0, 1};
