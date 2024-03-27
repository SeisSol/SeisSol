#include "Hdf5Reader.hpp"

#include <H5Dpublic.h>
#include <H5Fpublic.h>
#include <H5Spublic.h>
#include <H5version.h>
#include <IO/Datatype/Datatype.hpp>
#include <IO/Datatype/HDF5Type.hpp>
#include <IO/Datatype/Inference.hpp>
#include <hdf5.h>
#include <memory>
#include <stack>
#include <string>
#include <vector>

#include "utils/logger.h"

namespace {
static hid_t _eh(hid_t data) {
  if (data < 0) {
    logError() << "HDF5 error:" << data;
  }
  return data;
}
} // namespace

namespace seissol::io::reader::file {
Hdf5Reader::Hdf5Reader(MPI_Comm comm) : comm(comm) {}
void Hdf5Reader::openFile(const std::string& name) {
  hid_t h5falist = _eh(H5Pcreate(H5P_FILE_ACCESS));
#ifdef H5F_LIBVER_V18
  _eh(H5Pset_libver_bounds(h5falist, H5F_LIBVER_V18, H5F_LIBVER_V18));
#else
  _eh(H5Pset_libver_bounds(h5falist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#endif
  _eh(H5Pset_fapl_mpio(h5falist, comm, MPI_INFO_NULL));
  hid_t file = _eh(H5Fopen(name.c_str(), H5F_ACC_RDONLY, h5falist));
  _eh(H5Pclose(h5falist));

  handles.push(file);
}
void Hdf5Reader::openGroup(const std::string& name) {
  hid_t handle = _eh(H5Gopen(handles.top(), name.c_str(), H5P_DEFAULT));
  handles.push(handle);
}
std::size_t Hdf5Reader::attributeCount(const std::string& name) {
  hid_t attr = _eh(H5Aopen(handles.top(), name.c_str(), H5P_DEFAULT));
  hid_t attrspace = _eh(H5Aget_space(attr));
  hid_t rank = _eh(H5Sget_simple_extent_ndims(attrspace));
  std::vector<hsize_t> dims(rank);
  _eh(H5Sget_simple_extent_dims(attrspace, dims.data(), nullptr));
  _eh(H5Sclose(attrspace));
  _eh(H5Aclose(attr));
  return dims[0];
}
void Hdf5Reader::readAttributeRaw(void* data, const std::string& name) {
  /*hid_t attr = _eh(H5Aopen(handles.top(), name.c_str(), H5P_DEFAULT));
  _eh(H5Aread(attr, TODO, data));
  _eh(H5Aclose(attr));*/
}
std::size_t Hdf5Reader::dataCount(const std::string& name) {
  hid_t dataset = _eh(H5Dopen(handles.top(), name.c_str(), H5P_DEFAULT));
  hid_t dataspace = _eh(H5Dget_space(dataset));
  hid_t rank = _eh(H5Sget_simple_extent_ndims(dataspace));
  std::vector<hsize_t> dims(rank);
  _eh(H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr));
  _eh(H5Sclose(dataspace));
  _eh(H5Dclose(dataset));
}
void Hdf5Reader::readDataRaw(void* data,
                             const std::string& name,
                             std::size_t count,
                             std::shared_ptr<datatype::Datatype> targetType) {
  /*hid_t h5alist = H5Pcreate(H5P_DATASET_XFER);
  _eh(h5alist);
#ifdef USE_MPI
  _eh(H5Pset_dxpl_mpio(h5alist, H5FD_MPIO_COLLECTIVE));
#endif // USE_MPI

  hid_t dataset = _eh(H5Dopen(handles.top(), name.c_str(), H5P_DEFAULT));
  hid_t dataspace = _eh(H5Dget_space(dataset));

  hid_t datatype = datatype::convertToHdf5(targetType);

  hid_t rank = _eh(H5Sget_simple_extent_ndims(dataspace));
  std::vector<hsize_t> dims(rank);
  _eh(H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr));

  hid_t memspace = _eh(H5Screate_simple(rank, dims.data(), nullptr));

  std::vector<hsize_t> nullstart(rank);
  std::vector<hsize_t> readcount(rank);
  std::vector<hsize_t> filepos(rank);

  for (std::size_t i = 0; i < rounds; ++i) {
    _eh(H5Sselect_hyperslab(
        dataspace, H5S_SELECT_SET, nullstart.data(), nullptr, readcount.data(), nullptr));
    _eh(H5Sselect_hyperslab(
        dataspace, H5S_SELECT_SET, const hsize_t* start, nullptr, readcount.data(), nullptr));
    _eh(H5Dread(dataset, datatype, memspace, dataspace, h5alist, data));
  }*/
}
void Hdf5Reader::closeGroup() {
  _eh(H5Gclose(handles.top()));
  handles.pop();
}
void Hdf5Reader::closeFile() {
  _eh(H5Fclose(handles.top()));
  handles.pop();
}
} // namespace seissol::io::reader::file
