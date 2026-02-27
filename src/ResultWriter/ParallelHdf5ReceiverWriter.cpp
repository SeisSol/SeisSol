// SPDX-FileCopyrightText: 2024-2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause

#include "ParallelHdf5ReceiverWriter.h"

#include <array>
#include <utils/logger.h>

ParallelHdf5ReceiverWriter::ParallelHdf5ReceiverWriter(MPI_Comm comm,
                                                       const std::string& filename,
                                                       hsize_t totalReceivers,
                                                       hsize_t numVariables,
                                                       hsize_t totalTimeSteps)
    : comm_(comm), dims_{totalTimeSteps, totalReceivers, numVariables} {
  // Create property list for parallel I/O
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plistId, comm_, MPI_INFO_NULL);

  // Create file
  fileId_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
  H5Pclose(plistId);

  // Create dataspace for ReceiverData
  std::array<hsize_t, Rank> datasetDims = {totalTimeSteps, totalReceivers, numVariables};
  filespaceId_ = H5Screate_simple(Rank, datasetDims.data(), nullptr);

  // Create dataset "ReceiverData"
  dsetId_ = H5Dcreate(fileId_,
                      "ReceiverData",
                      H5T_NATIVE_DOUBLE,
                      filespaceId_,
                      H5P_DEFAULT,
                      H5P_DEFAULT,
                      H5P_DEFAULT);

  // Create dataspace for pointIds
  std::array<hsize_t, 1> pointIdDims = {totalReceivers};
  pointIdSpaceId_ = H5Screate_simple(1, pointIdDims.data(), nullptr);
  pointIdDsetId_ = H5Dcreate(fileId_,
                             "PointIds",
                             H5T_NATIVE_ULLONG,
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
  logDebug() << "writeChunk: timeOffset=" << timeOffset << " timeCount=" << timeCount
             << " receiverOffset=" << receiverOffset << " localReceiverCount=" << localReceiverCount
             << " dataSize=" << data.size();

  // Get the filespace directly
  hid_t filespaceId = H5Dget_space(dsetId_);
  checkStatus(filespaceId, "H5Dget_space in writeChunk");

  // Define the hyperslab in the file where we'll write
  std::array<hsize_t, Rank> start = {timeOffset, receiverOffset, 0};
  std::array<hsize_t, Rank> count = {timeCount, localReceiverCount, dims_[2]};

  logDebug() << "writeChunk: hyperslab start=(" << start[0] << ", " << start[1] << ", " << start[2]
             << ") count=(" << count[0] << ", " << count[1] << ", " << count[2] << ")";

  herr_t status = H5Sselect_hyperslab(
      filespaceId, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);
  checkStatus(status, "H5Sselect_hyperslab in writeChunk");

  // Create a matching memory dataspace for our data
  hid_t memspaceId = H5Screate_simple(Rank, count.data(), nullptr);
  checkStatus(memspaceId, "H5Screate_simple in writeChunk");

  // Perform the collective write
  hid_t xferPlist = H5Pcreate(H5P_DATASET_XFER);
  checkStatus(xferPlist, "H5Pcreate(H5P_DATASET_XFER) in writeChunk");

  status = H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE);
  checkStatus(status, "H5Pset_dxpl_mpio in writeChunk");

  status = H5Dwrite(dsetId_, H5T_NATIVE_DOUBLE, memspaceId, filespaceId, xferPlist, data.data());
  checkStatus(status, "H5Dwrite in writeChunk");

  // Flush to ensure data is fully written
  status = H5Fflush(fileId_, H5F_SCOPE_GLOBAL);
  checkStatus(status, "H5Fflush in writeChunk");

  // Cleanup
  H5Pclose(xferPlist);
  H5Sclose(memspaceId);
  H5Sclose(filespaceId);
}

void ParallelHdf5ReceiverWriter::writePointIds(hsize_t receiverOffset,
                                               hsize_t localReceiverCount,
                                               const std::vector<hsize_t>& pointIds) {
  logDebug() << "writePointIds: localReceiverCount=" << localReceiverCount
             << " pointIds.size()=" << pointIds.size();

  if (pointIds.size() != localReceiverCount) {
    throw std::runtime_error("writePointIds() size mismatch!");
  }

  // Define the hyperslab
  std::array<hsize_t, 1> start = {receiverOffset};
  std::array<hsize_t, 1> count = {localReceiverCount};

  H5Sselect_hyperslab(
      pointIdSpaceId_, H5S_SELECT_SET, start.data(), nullptr, count.data(), nullptr);

  // Create memory dataspace
  hid_t memspaceId = H5Screate_simple(1, count.data(), nullptr);

  // Collective write for pointIds
  hid_t xferPlist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE);

  herr_t status = H5Dwrite(
      pointIdDsetId_, H5T_NATIVE_ULLONG, memspaceId, pointIdSpaceId_, xferPlist, pointIds.data());
  checkStatus(status, "H5Dwrite in writePointIds");

  // Flush
  status = H5Fflush(fileId_, H5F_SCOPE_GLOBAL);
  checkStatus(status, "H5Fflush in writePointIds");

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
  std::array<hsize_t, 2> dims = {static_cast<hsize_t>(points.size()), 3};
  hid_t dataspaceId = H5Screate_simple(2, dims.data(), nullptr);

  // Create the dataset
  hid_t datasetId = H5Dcreate(fileId_,
                              "Coordinates",
                              H5T_NATIVE_DOUBLE,
                              dataspaceId,
                              H5P_DEFAULT,
                              H5P_DEFAULT,
                              H5P_DEFAULT);

  // Write the entire coords array in one shot
  herr_t status =
      H5Dwrite(datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, coords.data());
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
