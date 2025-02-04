#pragma once

#include <hdf5.h>
#include <mpi.h>
#include <stdexcept>
#include <string>
#include <vector>
#include <Eigen/Dense>

class ParallelHdf5ReceiverWriter {
  public:
  ParallelHdf5ReceiverWriter(MPI_Comm comm,
                             const std::string& filename,
                             hsize_t totalReceivers,
                             hsize_t numVariables);

  void writeChunk(hsize_t timeOffset,
                  hsize_t receiverOffset,
                  hsize_t timeCount,
                  hsize_t localReceiverCount,
                  const std::vector<double>& data);

  void writePointIds(hsize_t receiverOffset,
                     hsize_t localReceiverCount,
                     const std::vector<unsigned long long>& pointIds);

  void writeCoordinates(std::vector<Eigen::Vector3d> points);

  ~ParallelHdf5ReceiverWriter();

  private:
  static constexpr int RANK = 3;

  MPI_Comm comm_;
  hid_t fileId_;
  hid_t dsetId_;
  hid_t filespaceId_;
  hid_t pointIdSpaceId_;
  hid_t pointIdDsetId_;

  hsize_t dims_[RANK];

      void
      checkStatus(herr_t status, const std::string& msg) const {
    if (status < 0) {
      throw std::runtime_error("HDF5 error in: " + msg);
    }
  }
};
