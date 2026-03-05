// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "ParallelHdf5ReceiverWriter.h"

#include "Parallel/MPI.h"

#include <array>
#include <utils/logger.h>

namespace {
#define _eh(x) _ehh((x), __FILE__, __LINE__)

hid_t _ehh(hid_t data, const char* file, int line) {
  if (data < 0) {
    H5Eprint(H5Eget_current_stack(), stdout);
    throw std::runtime_error("HDF5 error at " + std::string(file) + ":" + std::to_string(line));
  }
  return data;
}
} // namespace

namespace seissol::writer {

ParallelHdf5ReceiverWriter::ParallelHdf5ReceiverWriter(MPI_Comm comm,
                                                       const std::string& filename,
                                                       hsize_t totalReceivers,
                                                       hsize_t numVariables)
    : comm_(comm), dims_{0, totalReceivers, numVariables} {
  // Create property list for parallel I/O
  hid_t plistId = _eh(H5Pcreate(H5P_FILE_ACCESS));
  _eh(H5Pset_fapl_mpio(plistId, comm_, MPI_INFO_NULL));

  // Create file
  fileId_ = _eh(H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId));
  _eh(H5Pclose(plistId));

  // Create dataspace for ReceiverData with an unlimited time dimension
  std::array<hsize_t, Rank> datasetDims = {0, totalReceivers, numVariables};
  std::array<hsize_t, Rank> maxDims = {H5S_UNLIMITED, totalReceivers, numVariables};
  filespaceId_ = _eh(H5Screate_simple(Rank, datasetDims.data(), maxDims.data()));

  // Create dataset creation property list and set chunking
  hid_t dcplId = _eh(H5Pcreate(H5P_DATASET_CREATE));
  // Chunk over chunk-sized time steps, but span all receivers and variables
  // (adjust chunk sizing according to typical sync frequency and memory limits, using 100
  // arbitrarily here)
  std::array<hsize_t, Rank> chunkDims = {100, totalReceivers, numVariables};
  _eh(H5Pset_chunk(dcplId, Rank, chunkDims.data()));

  // Set default fill value to 0 for unlimited expansion
  double fillValue = 0.0;
  _eh(H5Pset_fill_value(dcplId, H5T_NATIVE_DOUBLE, &fillValue));

  // Create dataset "ReceiverData"
  dsetId_ = _eh(H5Dcreate(
      fileId_, "ReceiverData", H5T_NATIVE_DOUBLE, filespaceId_, H5P_DEFAULT, dcplId, H5P_DEFAULT));

  _eh(H5Pclose(dcplId));

  // Create dataspace for pointIds (remains fixed size)
  pointIdSpaceId_ = _eh(H5Screate_simple(1, &totalReceivers, nullptr));
  pointIdDsetId_ = _eh(H5Dcreate(fileId_,
                                 "PointIds",
                                 H5T_NATIVE_UINT64,
                                 pointIdSpaceId_,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT,
                                 H5P_DEFAULT));
}

void ParallelHdf5ReceiverWriter::writeChunk(hsize_t timeOffset,
                                            hsize_t timeCount,
                                            const std::vector<std::uint64_t>& pointIds,
                                            const std::vector<double>& data) {
  hsize_t localReceiverCount = pointIds.size();

  logDebug() << "writeChunk: timeOffset=" << timeOffset << " timeCount=" << timeCount
             << " localReceiverCount=" << localReceiverCount
             << " dataSize=" << data.size();

  // Find the maximum time extent needed across all MPI ranks
  hsize_t localMaxTime = timeOffset + timeCount;
  hsize_t globalMaxTime = 0;
  MPI_Allreduce(
      &localMaxTime, &globalMaxTime, 1, seissol::Mpi::castToMpiType<hsize_t>(), MPI_MAX, comm_);

  // Get the filespace directly
  hid_t filespaceId = _eh(H5Dget_space(dsetId_));

  // Extend the dataset collectively if necessary.
  // All ranks must do this simultaneously with the same global dimensions.
  if (globalMaxTime > dims_[0]) {
    dims_[0] = globalMaxTime;
    _eh(H5Dset_extent(dsetId_, dims_.data()));
    // Refresh the filespaceId since the extent just changed
    _eh(H5Sclose(filespaceId));
    filespaceId = _eh(H5Dget_space(dsetId_));
  }

  // Handle case where this rank has no data to write or there's no data to be written.
  // With collective MPI I/O, all ranks must participate in the call regardless.
  if (timeCount == 0 || localReceiverCount == 0) {
    _eh(H5Sselect_none(filespaceId));

    // Create an empty memory dataspace
    hid_t memspaceId = _eh(H5Screate(H5S_NULL));

    hid_t xferPlist = _eh(H5Pcreate(H5P_DATASET_XFER));
    _eh(H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE));

    _eh(H5Dwrite(dsetId_, H5T_NATIVE_DOUBLE, memspaceId, filespaceId, xferPlist, data.data()));

    _eh(H5Pclose(xferPlist));
    _eh(H5Sclose(memspaceId));
    _eh(H5Sclose(filespaceId));
    return;
  }

  // Define the hyperslabs in the file where we'll write scattered across receivers
  _eh(H5Sselect_none(filespaceId));
  for (size_t r = 0; r < localReceiverCount; ++r) {
    if (pointIds[r] >= dims_[1]) {
      throw std::runtime_error("writeChunk(): receiver range exceeds allocated HDF5 dataset extent");
    }
    std::array<hsize_t, Rank> start = {timeOffset, static_cast<hsize_t>(pointIds[r]), 0};
    std::array<hsize_t, Rank> count = {timeCount, 1, dims_[2]};
    _eh(H5Sselect_hyperslab(filespaceId, H5S_SELECT_OR, start.data(), nullptr, count.data(), nullptr));
  }

  // Create a matching memory dataspace for our data
  // Even though it's a contiguous block of memory, mapping blocks will follow monotonically
  // increasing constraints since pointIds is sorted.
  std::array<hsize_t, Rank> memCount = {timeCount, localReceiverCount, dims_[2]};
  hid_t memspaceId = _eh(H5Screate_simple(Rank, memCount.data(), nullptr));

  // Perform the collective write
  hid_t xferPlist = _eh(H5Pcreate(H5P_DATASET_XFER));
  _eh(H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE));

  _eh(H5Dwrite(dsetId_, H5T_NATIVE_DOUBLE, memspaceId, filespaceId, xferPlist, data.data()));

  // Cleanup
  _eh(H5Pclose(xferPlist));
  _eh(H5Sclose(memspaceId));
  _eh(H5Sclose(filespaceId));
}

void ParallelHdf5ReceiverWriter::writePointIds(const std::vector<std::uint64_t>& pointIds) {
  hsize_t localReceiverCount = pointIds.size();
  logDebug() << "writePointIds: localReceiverCount=" << localReceiverCount
             << " pointIds.size()=" << pointIds.size();

  if (localReceiverCount == 0) {
    _eh(H5Sselect_none(pointIdSpaceId_));

    hid_t memspaceId = _eh(H5Screate(H5S_NULL));

    hid_t xferPlist = _eh(H5Pcreate(H5P_DATASET_XFER));
    _eh(H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE));

    _eh(H5Dwrite(pointIdDsetId_,
                 H5T_NATIVE_UINT64,
                 memspaceId,
                 pointIdSpaceId_,
                 xferPlist,
                 pointIds.data()));

    _eh(H5Pclose(xferPlist));
    _eh(H5Sclose(memspaceId));
    return;
  }

  // Define the scattered hyperslabs for pointIds
  _eh(H5Sselect_none(pointIdSpaceId_));
  for (size_t r = 0; r < localReceiverCount; ++r) {
    if (pointIds[r] >= dims_[1]) {
      throw std::runtime_error("writePointIds(): receiver range exceeds allocated HDF5 dataset extent");
    }
    std::array<hsize_t, 1> start = {static_cast<hsize_t>(pointIds[r])};
    std::array<hsize_t, 1> count = {1};
    _eh(H5Sselect_hyperslab(
        pointIdSpaceId_, H5S_SELECT_OR, start.data(), nullptr, count.data(), nullptr));
  }

  // Create memory dataspace
  std::array<hsize_t, 1> memCount = {localReceiverCount};
  hid_t memspaceId = _eh(H5Screate_simple(1, memCount.data(), nullptr));

  // Collective write for pointIds
  hid_t xferPlist = _eh(H5Pcreate(H5P_DATASET_XFER));
  _eh(H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE));

  _eh(H5Dwrite(
      pointIdDsetId_, H5T_NATIVE_UINT64, memspaceId, pointIdSpaceId_, xferPlist, pointIds.data()));

  // Cleanup
  _eh(H5Pclose(xferPlist));
  _eh(H5Sclose(memspaceId));
}

void ParallelHdf5ReceiverWriter::writeCoordinates(const std::vector<Eigen::Vector3d>& points) {
  int rank;
  MPI_Comm_rank(comm_, &rank);

  // Convert Eigen vectors into a contiguous buffer of doubles
  std::vector<double> coords(points.size() * 3);
  for (size_t i = 0; i < points.size(); ++i) {
    coords[i * 3 + 0] = points[i](0);
    coords[i * 3 + 1] = points[i](1);
    coords[i * 3 + 2] = points[i](2);
  }

  // Create a 2D dataspace: (nPoints x 3)
  std::array<hsize_t, 2> dims = {static_cast<hsize_t>(points.size()), 3};
  hid_t dataspaceId = _eh(H5Screate_simple(2, dims.data(), nullptr));

  // Create the dataset
  hid_t datasetId = _eh(H5Dcreate(fileId_,
                                  "Coordinates",
                                  H5T_NATIVE_DOUBLE,
                                  dataspaceId,
                                  H5P_DEFAULT,
                                  H5P_DEFAULT,
                                  H5P_DEFAULT));

  hid_t xferPlist = _eh(H5Pcreate(H5P_DATASET_XFER));
  _eh(H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE));

  // Only rank 0 actually writes the data, everyone else writes nothing.
  if (rank == 0) {
    _eh(H5Dwrite(datasetId, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xferPlist, coords.data()));
  } else {
    _eh(H5Sselect_none(dataspaceId));
    hid_t memspaceId = _eh(H5Screate(H5S_NULL));
    _eh(H5Dwrite(datasetId, H5T_NATIVE_DOUBLE, memspaceId, dataspaceId, xferPlist, nullptr));
    _eh(H5Sclose(memspaceId));
  }

  // Close HDF5 handles
  _eh(H5Pclose(xferPlist));
  _eh(H5Dclose(datasetId));
  _eh(H5Sclose(dataspaceId));
}

void ParallelHdf5ReceiverWriter::flush() { _eh(H5Fflush(fileId_, H5F_SCOPE_GLOBAL)); }

ParallelHdf5ReceiverWriter::~ParallelHdf5ReceiverWriter() {
  H5Dclose(dsetId_);         // Close the main dataset
  H5Dclose(pointIdDsetId_);  // Close the point ID dataset
  H5Sclose(filespaceId_);    // Close the main dataspace
  H5Sclose(pointIdSpaceId_); // Close the point ID dataspace
  H5Fclose(fileId_);         // Finally, close the file
}

} // namespace seissol::writer
