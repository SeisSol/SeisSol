#pragma once

#include <hdf5.h>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>

class ParallelHdf5ReceiverWriter {
  public:
  // Same constructor, but now we also keep a copy of comm
  ParallelHdf5ReceiverWriter(MPI_Comm comm,
                             const std::string& filename,
                             hsize_t totalReceivers,
                             hsize_t numVariables);

  void writeChunk(hsize_t timeOffset,
                  hsize_t receiverOffset,
                  hsize_t timeCount,
                  hsize_t localReceiverCount,
                  const std::vector<double>& data);

  ~ParallelHdf5ReceiverWriter();

  private:
  static constexpr int RANK = 3;

  // We store the communicator to print ranks
  MPI_Comm comm_;
  hid_t fileId_;
  hid_t dsetId_;
  hid_t filespaceId_;

  // dataset dimensions
  hsize_t dims_[RANK];

  void checkStatus(herr_t status, const std::string& msg) const {
    if (status < 0) {
      throw std::runtime_error("HDF5 error in: " + msg);
    }
  }
};
