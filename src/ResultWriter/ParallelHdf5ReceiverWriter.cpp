// ParallelHdf5ReceiverWriter.cpp

#include "ParallelHdf5ReceiverWriter.h"
#include <iostream>

ParallelHdf5ReceiverWriter::ParallelHdf5ReceiverWriter(MPI_Comm comm,
                                                       const std::string& filename,
                                                       hsize_t totalReceivers,
                                                       hsize_t numVariables,
                                                       hsize_t totalTimeSteps)
    : comm_(comm)
{


  dims_[0] = totalTimeSteps;
  dims_[1] = totalReceivers;
  dims_[2] = numVariables;

  // Create property list for parallel I/O
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plistId, comm_, MPI_INFO_NULL);

  // Create file
  fileId_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
  H5Pclose(plistId);

  // Create dataspace
  //hsize_t initialDims[RANK] = {0, totalReceivers, numVariables};
  // hsize_t maxDims[RANK] = {H5S_UNLIMITED, totalReceivers, numVariables};
  //filespaceId_ = H5Screate_simple(RANK, initialDims, maxDims);

  hsize_t datasetDims[RANK] = {totalTimeSteps, totalReceivers, numVariables};
  filespaceId_ = H5Screate_simple(RANK, datasetDims, nullptr); // No maxDims needed

  // Create a dataset creation property list for "ReceiverData".
  // hid_t dcplId = H5Pcreate(H5P_DATASET_CREATE);
  // hsize_t chunkDims[RANK] = {128, 4096, 10};

  // 4. Enable chunking
  // H5Pset_chunk(dcplId, 3, chunkDims);

  // Create dataset "ReceiverData"
  dsetId_ = H5Dcreate(fileId_,
                      "ReceiverData",
                      H5T_NATIVE_DOUBLE,
                      filespaceId_,
                      H5P_DEFAULT,
                      H5P_DEFAULT,
                      H5P_DEFAULT);

  // H5Pclose(dcplId);

  // Create dataspace for pointIds
  hsize_t pointIdDims_[1] = {totalReceivers}; // 1D array for pointIds
  pointIdSpaceId_ = H5Screate_simple(1, pointIdDims_, nullptr);
  pointIdDsetId_ = H5Dcreate(fileId_,
                             "PointIds",
                             H5T_NATIVE_ULLONG, // Use unsigned long long for pointIds
                             pointIdSpaceId_,
                             H5P_DEFAULT,
                             H5P_DEFAULT,
                             H5P_DEFAULT);
}

void ParallelHdf5ReceiverWriter::writeChunk(hsize_t timeOffset,
                                            hsize_t timeCount,
                                            hsize_t receiverOffset,
                                            hsize_t localReceiverCount,
                                            const std::vector<double>& data) {
  // Get my MPI rank for debugging
  int myRank = -1;
  MPI_Comm_rank(comm_, &myRank); // comm_ is presumably stored from the constructor

  std::cout << "[Rank " << myRank << "] Entering writeChunk()"
            << " with timeOffset=" << timeOffset << " timeCount=" << timeCount
            << " receiverOffset=" << receiverOffset << " localReceiverCount=" << localReceiverCount
            << " dataSize=" << data.size() << std::endl;
  std::fflush(stdout); // Force immediate flush to see this output in real time

  // Get the filespace directly
  hid_t filespaceId = H5Dget_space(dsetId_);
  checkStatus(filespaceId, "H5Dget_space in writeChunk");

  // ----------------------------------------------------
  // 3) Define the hyperslab in the file where we'll write
  // ----------------------------------------------------
  hsize_t start[RANK] = {timeOffset, receiverOffset, 0};
  hsize_t count[RANK] = {timeCount, localReceiverCount, dims_[2]};

  std::cout << "[Rank " << myRank << "] Selecting hyperslab with start=(" << start[0] << ", "
            << start[1] << ", " << start[2] << ") "
            << " and count=(" << count[0] << ", " << count[1] << ", " << count[2] << ")"
            << std::endl;
  std::fflush(stdout);

  herr_t status = H5Sselect_hyperslab(filespaceId, H5S_SELECT_SET, start, nullptr, count, nullptr);
  checkStatus(status, "H5Sselect_hyperslab in writeChunk");

  std::cout << "[Rank " << myRank << "] Hyperslab selection successful." << std::endl;
  std::fflush(stdout);

  // ----------------------------------------------------
  // 4) Create a matching memory dataspace for our data
  // ----------------------------------------------------
  hid_t memspaceId = H5Screate_simple(RANK, count, nullptr);
  checkStatus(memspaceId, "H5Screate_simple in writeChunk");

  std::cout << "[Rank " << myRank << "] Created memory dataspace." << std::endl;
  std::fflush(stdout);

  // ----------------------------------------------------
  // 5) Perform the write
  // ----------------------------------------------------
  hid_t xferPlist = H5Pcreate(H5P_DATASET_XFER);
  checkStatus(xferPlist, "H5Pcreate(H5P_DATASET_XFER) in writeChunk");

  status = H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE);
  checkStatus(status, "H5Pset_dxpl_mpio in writeChunk");

  std::cout << "[Rank " << myRank << "] About to call H5Dwrite (independent mode)..." << std::endl;
  std::fflush(stdout);

  status = H5Dwrite(dsetId_, H5T_NATIVE_DOUBLE, memspaceId, filespaceId, xferPlist, data.data());
  checkStatus(status, "H5Dwrite in writeChunk");

  std::cout << "[Rank " << myRank << "] H5Dwrite call completed successfully." << std::endl;
  std::fflush(stdout);

  // ----------------------------------------------------
  // 6) (Optional) Flush to ensure data is fully written
  // ----------------------------------------------------
  std::cout << "[Rank " << myRank << "] Flushing file..." << std::endl;
  std::fflush(stdout);

  status = H5Fflush(fileId_, H5F_SCOPE_GLOBAL);
  checkStatus(status, "H5Fflush in writeChunk");

  std::cout << "[Rank " << myRank << "] Flush completed." << std::endl;
  std::fflush(stdout);

  // ----------------------------------------------------
  // 7) Cleanup
  // ----------------------------------------------------
  status = H5Pclose(xferPlist);
  checkStatus(status, "H5Pclose(xferPlist) in writeChunk");

  status = H5Sclose(memspaceId);
  checkStatus(status, "H5Sclose(memspaceId) in writeChunk");

  status = H5Sclose(filespaceId);
  checkStatus(status, "H5Sclose(filespaceId) in writeChunk");

  std::cout << "[Rank " << myRank << "] Exiting writeChunk()." << std::endl;
  std::fflush(stdout);
}

void ParallelHdf5ReceiverWriter::writePointIds(hsize_t receiverOffset,
                                               hsize_t localReceiverCount,
                                               const std::vector<hsize_t>& pointIds) {

  std::cout << "localReceiverCount: " << localReceiverCount
            << ", pointIds.size(): " << pointIds.size() << std::endl;

  if (pointIds.size() != localReceiverCount) {
    throw std::runtime_error("writePointIds() size mismatch!");
  }

  // Define the hyperslab
  hsize_t start[1] = {receiverOffset};
  hsize_t count[1] = {localReceiverCount};

  H5Sselect_hyperslab(pointIdSpaceId_, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Create memory dataspace
  hid_t memspaceId = H5Screate_simple(1, count, nullptr);

  // Collective write for pointIds
  hid_t xferPlist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE);

  herr_t status = H5Dwrite(
      pointIdDsetId_, H5T_NATIVE_ULLONG, memspaceId, pointIdSpaceId_, xferPlist, pointIds.data());
  checkStatus(status, "H5Dwrite");

  // Flush
  status = H5Fflush(fileId_, H5F_SCOPE_GLOBAL);
  checkStatus(status, "H5Fflush");

  // Cleanup
  H5Pclose(xferPlist);
  H5Sclose(memspaceId);
}

void ParallelHdf5ReceiverWriter::writeCoordinates(std::vector<Eigen::Vector3d> points) {
  // Convert Eigen vectors into a contiguous buffer of doubles
  std::vector<double> coords(points.size() * 3);
  for (size_t i = 0; i < points.size(); ++i) {
    coords[i * 3 + 0] = points[i](0);
    coords[i * 3 + 1] = points[i](1);
    coords[i * 3 + 2] = points[i](2);
  }

  // Create a 2D dataspace: (nPoints x 3)
  hsize_t dims[2];
  dims[0] = static_cast<hsize_t>(points.size()); // number of points
  dims[1] = 3;                                   // x,y,z
  hid_t dataspaceId = H5Screate_simple(2, dims, nullptr);

  // Create the dataset.
  hid_t datasetId = H5Dcreate(fileId_,           // existing HDF5 file handle (from the class)
                              "Coordinates",     // name of the new dataset
                              H5T_NATIVE_DOUBLE, // store as double precision
                              dataspaceId,
                              H5P_DEFAULT, // link creation property list
                              H5P_DEFAULT, // dataset creation property list
                              H5P_DEFAULT  // dataset access property list
  );

  // Write the entire coords array in one shot
  herr_t status = H5Dwrite(datasetId,
                           H5T_NATIVE_DOUBLE,
                           H5S_ALL, // memory dataspace = entire buffer
                           H5S_ALL, // file dataspace = entire dataset
                           H5P_DEFAULT,
                           coords.data());
  checkStatus(status, "H5Dwrite in writeCoordinates");

  // Close HDF5 handles
  H5Dclose(datasetId);
  H5Sclose(dataspaceId);
}

ParallelHdf5ReceiverWriter::~ParallelHdf5ReceiverWriter() {
  H5Dclose(dsetId_);         // Close the main dataset
  H5Dclose(pointIdDsetId_);  // Close the point ID dataset
  H5Sclose(filespaceId_);    // Close the main dataspace
  H5Sclose(pointIdSpaceId_); // Close the point ID dataspace
  H5Fclose(fileId_);         // Finally, close the file
  }
