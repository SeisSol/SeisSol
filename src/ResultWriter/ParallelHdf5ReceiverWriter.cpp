// ParallelHdf5ReceiverWriter.cpp

#include "ParallelHdf5ReceiverWriter.h"
#include <iostream>

ParallelHdf5ReceiverWriter::ParallelHdf5ReceiverWriter(MPI_Comm comm,
                                                       const std::string& filename,
                                                       hsize_t totalReceivers,
                                                       hsize_t numVariables)
    : comm_(comm) // <-- fixes the bug
{

  dims_[0] = 1000000; // max time steps. TODO: calculate this
  dims_[1] = totalReceivers;
  dims_[2] = numVariables;

  // Create property list for parallel I/O
  hid_t plistId = H5Pcreate(H5P_FILE_ACCESS);
  H5Pset_fapl_mpio(plistId, comm_, MPI_INFO_NULL);

  // Create file
  fileId_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistId);
  H5Pclose(plistId);

  // Create dataspace
  filespaceId_ = H5Screate_simple(RANK, dims_, nullptr);

  // Create dataset "ReceiverData"
  dsetId_ = H5Dcreate(fileId_,
                      "ReceiverData",
                      H5T_NATIVE_DOUBLE,
                      filespaceId_,
                      H5P_DEFAULT,
                      H5P_DEFAULT,
                      H5P_DEFAULT);
}

void ParallelHdf5ReceiverWriter::writeChunk(hsize_t timeOffset,
                                            hsize_t receiverOffset,
                                            hsize_t timeCount,
                                            hsize_t localReceiverCount,
                                            const std::vector<double>& data) {

  // Check data size
  if (data.size() != static_cast<size_t>(timeCount * localReceiverCount * dims_[2])) {
    throw std::runtime_error("writeChunk() data size mismatch!");
  }

  // Define the hyperslab
  hsize_t start[RANK] = {timeOffset, receiverOffset, 0};
  hsize_t count[RANK] = {timeCount, localReceiverCount, dims_[2]};

  H5Sselect_hyperslab(filespaceId_, H5S_SELECT_SET, start, nullptr, count, nullptr);

  // Create a matching memory dataspace
  hid_t memspaceId = H5Screate_simple(RANK, count, nullptr);

  // Collective write
  hid_t xferPlist = H5Pcreate(H5P_DATASET_XFER);
  H5Pset_dxpl_mpio(xferPlist, H5FD_MPIO_COLLECTIVE);

  herr_t status =
      H5Dwrite(dsetId_, H5T_NATIVE_DOUBLE, memspaceId, filespaceId_, xferPlist, data.data());
  checkStatus(status, "H5Dwrite");

  // Flush
  status = H5Fflush(fileId_, H5F_SCOPE_GLOBAL);
  checkStatus(status, "H5Fflush");

  // Cleanup
  H5Pclose(xferPlist);
  H5Sclose(memspaceId);
}

ParallelHdf5ReceiverWriter::~ParallelHdf5ReceiverWriter() {
  // For example, close open HDF5 objects:
  H5Sclose(filespaceId_);
  H5Dclose(dsetId_);
  H5Fclose(fileId_);
}
