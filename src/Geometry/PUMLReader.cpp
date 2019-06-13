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
	  //logError() << "Periodic boundary conditions are currently not supported for MPI+HDF5.";
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
		  planes.push_back(computePlane(faceNodesCoords));
		  planes2Vertices.push_back(std::unordered_set<PeriodicVertex>());
		  ++idx;
		}
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

	auto planeFragments = std::vector<PlaneFragment>();
	for (int z = 0; z < planes.size(); ++z) {
	  const auto& plane = planes[z];
	  planeFragments.push_back(PlaneFragment(planes2Vertices[z],
						 plane));
	}
	
	
	logInfo(rank) << "Number of planes =" << planeFragments.size();
	for (const auto& planeFragment : planeFragments) {
	  const auto &plane = planeFragment.plane;
	  const auto &n = plane.normal;
	  const auto &p = plane.point;
	  const auto &t1 = plane.tangent1;
	  const auto &t2 = plane.tangent2;
	  const auto constDim = planeFragment.constDim;
	  const auto &convexHull = planeFragment.convexHull;
	  std::cout << "Constant dimension"
		    << constDim << std::endl;
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
	  std::cout << "Convex hull:\n";
	  for (const auto &p : convexHull) {
	    std::cout << p.coords[0] << ",\t"
		      << p.coords[1] << ",\t"
		      << std::endl;
	  }
	    
	  std::cout << std::endl;
	}
	assert(planes.size() == planeFragments.size()); // TODO(Lukas): Communicate plane frag ids!

	std::vector<std::array<int,3>> periodicSendList; // [element id, face idx (0..3), planeIdx]
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
	      std::array<const double*, 3> faceNodesCoords;
	      getFaceNodes(m_elements[i],
			   vertices,
			   j,
			   faceNodesCoords);
	      assert(faces2Planes.find(faceids[j]) != faces2Planes.end());
	      const auto facePlaneIdx = faces2Planes[faceids[j]];
	      const auto facePlane = planes[facePlaneIdx];
	      
	      const auto result = findMatchingCell(i,
						   j,
						   facePlane,
						   faceNodesCoords,
						   planes,
						   planes2Vertices,
						   vertices2Cells,
						   vertices2Vertices);
	      if (result >= 0) {
		neighbors[j] = result;
	      } else {
		periodicSendList.push_back({i, j, facePlaneIdx});
	      }
	    }
	    // Set the neighbor directly in seissol data structure.
	    m_elements[i].neighbors[FACE_PUML2SEISSOL[j]] = neighbors[j];
	  }
	  m_elements[i].material = material[i];
	}

	MPI_Barrier(MPI::mpi.comm());
	logInfo() << "Rank" << rank << "needs to send" << periodicSendList.size() << "faces.";

	// Plane fragments
	auto planeFragmentsDataBuffer = std::vector<double>();
	auto planeFragmentsLengthsBuffer = std::vector<int>();
	for (const auto &frag : planeFragments) {
	  encodePlaneFragment(frag,
			      planeFragmentsDataBuffer,
			      planeFragmentsLengthsBuffer);
	}
	auto serializedPlaneFragments = PlaneFragmentsSerialized(
								 planeFragmentsDataBuffer,
								 planeFragmentsLengthsBuffer);

	// Communicate size of plane fragment buffers.
	const auto planesDataSizeSend = serializedPlaneFragments.dataSize;
	const auto planesLengthsSizeSend = serializedPlaneFragments.lengthsSize;
	auto planesDataSizePerNode = std::vector<int>(MPI::mpi.size());
	auto planesLengthsSizePerNode = std::vector<int>(MPI::mpi.size());
	MPI_Allgather(&planesDataSizeSend, 1, MPI_INT, planesDataSizePerNode.data(),
		      1, MPI_INT, MPI::mpi.comm());
	MPI_Allgather(&planesLengthsSizeSend, 1, MPI_INT, planesLengthsSizePerNode.data(),
		      1, MPI_INT, MPI::mpi.comm());

	for (int i = 0; i < MPI::mpi.size(); ++i) {
	  std::cout << "Plane fragments from rank " << i
		     << " with data size" << planesDataSizePerNode[i]
		     << " and lenghts size" << planesLengthsSizePerNode[i]
		     << std::endl;
	}

	auto serializedPlaneFragmentsDispl = std::vector<int>(MPI::mpi.size());
	for (int z = 0; z < serializedPlaneFragmentsDispl.size(); ++z) {
	  serializedPlaneFragmentsDispl[z] = (z > 0) ?
	    serializedPlaneFragmentsDispl[z-1] + planesLengthsSizePerNode[z]
	    : 0;
	  std::cout << "planeFragmentsDispl z =" << z << ": "
		    << serializedPlaneFragmentsDispl[z] << std::endl;
	}
	MPI_Barrier(MPI::mpi.comm());

	// Initialize receive buffers and their displacements.
	// Plane fragments
	auto planesDataRecvDispl = std::vector<int>(MPI::mpi.size());
	auto planesLengthsRecvDispl = std::vector<int>(MPI::mpi.size());
	for (int i = 0; i < MPI::mpi.size(); ++i) {
	  planesDataRecvDispl[i] = (i > 0) ?
	    (planesDataRecvDispl[i-1] + planesDataSizePerNode[i-1]) : 0; 
	  planesLengthsRecvDispl[i] = (i > 0) ?
	    (planesLengthsRecvDispl[i-1] + planesLengthsSizePerNode[i-1]) : 0; 

	}
	const auto planesDataRecvSize = planesDataRecvDispl[planesDataRecvDispl.size() - 1]
	  + planesDataSizePerNode[planesDataSizePerNode.size() - 1];
	auto planesDataRecvBuffer = std::vector<double>(planesDataRecvSize);
	const auto planesLengthsRecvSize = planesLengthsRecvDispl[planesLengthsRecvDispl.size() - 1]
	  + planesLengthsSizePerNode[planesLengthsSizePerNode.size() - 1];
	auto planesLengthsRecvBuffer = std::vector<int>(planesLengthsRecvSize);

	// Communicate data + lengths
	MPI_Allgatherv(serializedPlaneFragments.data,
		       planesDataSizePerNode[rank],
		       MPI_DOUBLE,
		       planesDataRecvBuffer.data(),
		       planesDataSizePerNode.data(),
		       planesDataRecvDispl.data(),
		       MPI_DOUBLE,
		       MPI::mpi.comm()
		       );

	MPI_Allgatherv(serializedPlaneFragments.lengths,
		       planesLengthsSizePerNode[rank],
		       MPI_INT,
		       planesLengthsRecvBuffer.data(),
		       planesLengthsSizePerNode.data(),
		       planesLengthsRecvDispl.data(),
		       MPI_INT,
		       MPI::mpi.comm()
		       );

	const auto serializedPlaneFragmentsRecv =
	  PlaneFragmentsSerialized(planesDataRecvBuffer,
				   planesLengthsRecvBuffer);

	const auto planeFragmentsRecv = decodePlaneFragments(serializedPlaneFragmentsRecv,
							     planesLengthsSizePerNode);

	auto facesSendBuffers = std::unordered_map<int, std::vector<FaceMPI>>{};
	// This map maps ranks to vector of array indices of periodicSendList
	// can be used to identify vertices.
	auto facesSendMap = std::unordered_map<int, std::vector<int>>{};
	// Iterate over all vertices and check which nodes could hold the mirrored vertices.
	for (int z = 0; z < periodicSendList.size(); ++z) {
	  const auto curPlaneIdx = periodicSendList[z][2];

	  const auto& facePlane = planes[curPlaneIdx];
	  // Check all plane fragments
	  for (const auto& frag: planeFragmentsRecv) {
	    const auto& plane = frag.plane;
	    // TODO: In this way, we always match on own rank
	    // with same plane as the vertex that should be mirrored
	    // stems from.
	    
	    // if (plane == facePlane) continue;
	    auto movedPlane = frag.plane;
	    std::copy_n(facePlane.point, 3, movedPlane.point);

	    std::array<const double*, 3> faceNodesCoords;
	    getFaceNodes(m_elements[periodicSendList[z][0]],
			 vertices,
			 periodicSendList[z][1],
			 faceNodesCoords);

	    // Check if vertex is in moved plane.
	    bool containsAllPoints = true;
	    for (int y = 0; y < 3; ++y) {
	      containsAllPoints &=
		movedPlane.containsPoint(faceNodesCoords[y]);
	    }
	    if (!containsAllPoints) continue;
	    
	    //std::cout << "Found plane." << std::endl;
	      
	    // Check if vertex is inside plane fragment.
	    int pointsInside = 0;
	    for (int y = 0; y < 3; ++y) {
	      const auto* curCoords = faceNodesCoords[y];
	      PointPlane pointPlane;
	      const auto planeConstDim = plane.getConstantDim();
	      int planeOtherDims[2];
	      plane.getOtherDims(planeOtherDims);
	      pointPlane.coords[planeOtherDims[0]] = curCoords[planeOtherDims[0]];
	      pointPlane.coords[planeOtherDims[1]] = curCoords[planeOtherDims[1]];
	    
	      bool isInside = isInsideConvexPolygon(pointPlane,
						    frag.convexHull);
	      if (isInside) ++pointsInside;

	    }
	    if (true || pointsInside >= 0) {
	      FaceMPI faceMpi;
	      // Copy coordinates to buffer
	      for (int y = 0; y < 3; ++y) {
		std::copy_n(faceNodesCoords[y],
			    3,
			    faceMpi.coords[y]);
	      }
	      faceMpi.planeId = periodicSendList[z][2];
	      facesSendBuffers[frag.rank].push_back(faceMpi);
	      facesSendMap[frag.rank].push_back(z);
	    }
	  }
	}

	auto sizeRequests = std::vector<MPI_Request>(2 * MPI::mpi.size());
	auto facesSendSizePerNode = std::vector<int>(MPI::mpi.size(), 0);
	for (const auto& s : facesSendBuffers) {
	  facesSendSizePerNode[s.first] = s.second.size();
	  std::cout << "Send " << s.second.size() << " to rank " << s.first
		    << " from rank " << rank << std::endl;
	}


	// TODO(Lukas) Switch to allgather
	auto facesRecvSize = std::vector<int>(MPI::mpi.size());
	for (int otherRank = 0; otherRank < MPI::mpi.size(); ++otherRank) {
	  // Get size from other ranks
	  MPI_Irecv(&facesRecvSize[otherRank],
		    1,
		    MPI_INT,
		    otherRank,
		    0,
		    MPI::mpi.comm(),
		    &sizeRequests[otherRank]);

	  MPI_Isend(&facesSendSizePerNode[otherRank],
		    1,
		    MPI_INT,
		    otherRank,
		    0,
		    MPI::mpi.comm(),
		    &sizeRequests[MPI::mpi.size() + otherRank]);
	}
	MPI_Waitall(sizeRequests.size(), sizeRequests.data(), MPI_STATUSES_IGNORE);

	// Initialize faces recv buffers
	// TODO(Lukas) Change type to correct one
	auto faceMpi_t = registerFaceMPI();
	auto facesRecvBuffers = std::unordered_map<int, std::vector<FaceMPI>>{};
	for (int otherRank = 0; otherRank < facesRecvSize.size(); ++otherRank) {
	  const auto otherRankSendSize = facesRecvSize[otherRank];
	  // TODO(Lukas) Maybe exclude own rank?
	  if (otherRankSendSize == 0) continue; // no empty messages
	  facesRecvBuffers[otherRank] = std::vector<FaceMPI>(otherRankSendSize);
	}
	for (int y = 0; y < MPI::mpi.size(); ++y) {
	    std::cout << "Sending " << facesRecvSize[y]
		      << " to rank " << rank
		      << " faces from rank " << y << std::endl;
	}
	MPI_Barrier(MPI::mpi.comm());
	//throw -1;
	
	std::cout << "Receiving faces from " << facesRecvBuffers.size()
		  << " ranks" << std::endl;
	// Recv/Send faces
	auto faceRequests = std::vector<MPI_Request>(2 * facesRecvBuffers.size());
	auto faceRequestsIdx = 0;
	//for (int otherRank = 0; otherRank < facesRecvBuffers.size(); ++otherRank) {
	// TODO(Lukas) Switch to explicit all-to-all comm?
	for (const auto& buf : facesRecvBuffers) {
	  MPI_Irecv(facesRecvBuffers[buf.first].data(),
		    facesRecvBuffers[buf.first].size(),
		    faceMpi_t,
		    buf.first,
		    0,
		    MPI::mpi.comm(),
		    &faceRequests[faceRequestsIdx++]);
	}
	//for (int otherRank = 0; otherRank < facesSendBuffers.size(); ++otherRank) {
	for (const auto& buf : facesSendBuffers) {
	  MPI_Isend(facesSendBuffers[buf.first].data(),
		    facesSendBuffers[buf.first].size(),
		    faceMpi_t,
		    buf.first,
		    0,
		    MPI::mpi.comm(),
		    &faceRequests[faceRequestsIdx++]);

	}
	faceRequests.resize(faceRequestsIdx); // TODO(Lukas) Is this resize correct?
	std::cout << "faceRequestsIdx, faceRequests.size() = " << faceRequestsIdx
		  << ", " << faceRequests.size() << std::endl;
	MPI_Waitall(faceRequests.size(), faceRequests.data(), MPI_STATUSES_IGNORE);
	std::cout << "Send/Recv faces of size " << faceRequests.size() << std::endl;
	
	// Initialise buffers for the matching matching.
	// TODO(Lukas) Are the sizes of the buffers correct or exchanged?
	auto matchingRequests = std::vector<MPI_Request>(2 * MPI::mpi.size());
	auto matchingRequestsIdx = 0;
	auto matchingSendBuffers = std::unordered_map<int, std::vector<int>>{};
	auto matchingRecvBuffers = std::unordered_map<int, std::vector<int>>{};

	// TODO(Lukas) Switch to all-to-all communication?
	for (const auto &buf : facesSendBuffers) {
	  matchingRecvBuffers[buf.first] = std::vector<int>(facesSendSizePerNode[buf.first], -1);
	  MPI_Irecv(matchingRecvBuffers[buf.first].data(),
		    matchingRecvBuffers[buf.first].size(),
		    MPI_INT,
		    buf.first,
		    0,
		    MPI::mpi.comm(),
		    &matchingRequests[matchingRequestsIdx++]);
	}
	for (const auto &buf : facesRecvBuffers) {
	  matchingSendBuffers[buf.first] = std::vector<int>(facesRecvSize[buf.first], -1);
	}

	for (const auto &buf : facesRecvBuffers) {
	  for (int z = 0; z < buf.second.size(); ++z) {
	    // TODO(Lukas) Matching is completely wrong probably.
	    const auto &faceMpi = buf.second[z];
	    auto faceNodesCoords = std::array<const double*, 3>{
	      faceMpi.coords[0],
	      faceMpi.coords[1],
	      faceMpi.coords[2]};
	    // TODO: Is offset correct?
	    const auto facePlaneIdx = serializedPlaneFragmentsDispl[buf.first] + faceMpi.planeId;
	    assert(facePlaneIdx < planeFragmentsRecv.size());
	    // std::cout << "Displ = " << serializedPlaneFragmentsDispl[buf.first]
	    // 	      << "facePlIdx = " << faceMpi.planeId
	    // 	      << "Size = " << planeFragmentsRecv.size()
	    // 	      << std::endl;
	    const auto &facePlane = planeFragmentsRecv[facePlaneIdx].plane;
	    //logInfo(rank) << "Try matching";
	    //MPI_Barrier(MPI::mpi.comm());
	    const auto result = findMatchingCell(-1, // TODO
						 0, // TODO
						 facePlane,
						 faceNodesCoords,
						 planes,
						 // TODO: Shouldnt modify these maps
						 planes2Vertices,
						 vertices2Cells,
						 vertices2Vertices);
	    matchingSendBuffers[buf.first][z] = result;
	    if (buf.first == rank) continue;
	    if (result > 0) {
	      std::cout << "Matching after MPI" << std::endl;
	    } else {
	      std::cout << "No matching after MPI" << std::endl;
	    }

	  }
	}
	
	for (const auto &buf : matchingSendBuffers) {
	  MPI_Isend(matchingSendBuffers[buf.first].data(),
		    matchingSendBuffers[buf.first].size(),
		    MPI_INT,
		    buf.first,
		    0,
		    MPI::mpi.comm(),
		    &matchingRequests[matchingRequestsIdx++]);

	}

	matchingRequests.resize(matchingRequestsIdx);
	MPI_Waitall(matchingRequests.size(), matchingRequests.data(), MPI_STATUSES_IGNORE);
	std::cout << "Send/Recv matchings of size " << matchingRequests.size() << std::endl;

	auto combinedMatching = std::unordered_map<std::pair<int,int>,
						   std::pair<int,int>,
						   pair_hash>{};
	// Vertices id in periodicSendList
	// Matching in matchingRecvBuffers
	for (const auto& recvMatching : matchingRecvBuffers) {
	  const auto fromRank = recvMatching.first;
	  if (fromRank == rank) continue; // Ignore matchings from own rank.
	  const auto& curMatchings = recvMatching.second;
	  std::cout << "curMatchings.size(), periodicSendList.size() "
		    << curMatchings.size() << ", "
		    << periodicSendList.size()
		    << std::endl;
	  assert(curMatchings.size() == facesSendMap[fromRank].size());
	  for (int matchingId = 0; matchingId < curMatchings.size(); ++matchingId) {
	    const auto sendListIdx = facesSendMap[fromRank][matchingId];
	    const auto curCellId = periodicSendList[sendListIdx][0];
	    const auto curFaceId = periodicSendList[sendListIdx][1];
	    const auto key = std::make_pair(curCellId, curFaceId);
	    const auto matching = curMatchings[matchingId];
	    if (matching < 0) continue; // no matching
	    if (combinedMatching.find(key) != combinedMatching.end()
		&& combinedMatching[key].first != matching) {
	      std::cout << "\nDuplicate, different matching" << std::endl;
	      std::cout << "Old: " << combinedMatching[key].first << std::endl;
	    }
	    // Make sure we didn't match the cell with itself.
	    assert(matching != curCellId);
	    // std::cout << "cell idx, face idx, matching "
	    // 	      << curCellId << ", "
	    // 	      << curFaceId << ", "
	    // 	      << matching << std::endl;

	    // Just overwrite matching...
	    combinedMatching[key] =  {matching, fromRank};
	    // TODO(Lukas) Save rank.
	  }
	}
	
	std::cout << combinedMatching.size() << " matchings found. "
		  << "We sent " << periodicSendList.size() << " faces"
		  << std::endl;
	
	MPI_Barrier(MPI::mpi.comm());
	//throw -1;

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

		    if (bcCurrentFace == 6) {
		      // Periodic boundary with MPI
		      const auto key = std::make_pair(i, j);
		      const auto ownerRankIt = combinedMatching.find(key);
		      if (ownerRankIt != combinedMatching.end()) {
			//const auto cellId = ownerRankIt->second; // Matched
			const auto ownerRank = ownerRankIt->second.second; // TODO(Lukas) Set owner rank
			neighborInfo[ownerRank].push_back(faceids[j]);
			m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = ownerRank;
			//std::cout << "Matching for i, j: ";
		      } else {
			//std::cout << "No matching for i, j: ";
		      }
		      std::cout << i << ", "
				<< j << ", "
				<< std::endl;
		      
		    } else if (!faces[faceids[j]].isShared()) {
		      // isShared is true if face belongs to more than one rank.
		      // Boundary sides
		      m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = rank;
		    } else {
		      // MPI Boundary
		      neighborInfo[faces[faceids[j]].shared()[0]].push_back(faceids[j]);

		      m_elements[i].neighborRanks[FACE_PUML2SEISSOL[j]] = faces[faceids[j]].shared()[0];
		    }
		  } else {
		    continue; // TODO(Lukas)
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

	MPI_Barrier(MPI::mpi.comm());
	//throw -1;

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
