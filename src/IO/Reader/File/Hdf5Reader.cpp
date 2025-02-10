// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Hdf5Reader.h"

#include <IO/Datatype/Datatype.h>
#include <IO/Datatype/HDF5Type.h>
#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <algorithm>
#include <cstddef>
#include <cstdio>
#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <vector>

#include "utils/logger.h"

namespace {
#define _eh(x) _ehh(x, __FILE__, __LINE__)

hid_t _ehh(hid_t data, const char* file, int line) {
  if (data < 0) {
    H5Eprint(H5Eget_current_stack(), stdout);
    logError() << "HDF5 error:" << data << "at" << file << ":" << line;
  }
  return data;
}
} // namespace

namespace seissol::io::reader::file {

Hdf5Reader::Hdf5Reader(MPI_Comm comm) : comm(comm) {}
void Hdf5Reader::checkExistence(const std::string& name, const std::string& type) {
  if (handles.empty()) {
    // we're opening a file, i.e. we'd need to check for the existence here
  } else if (type == "attribute") {
    // TODO: change type to enum
    if (_eh(H5Aexists(handles.top(), name.c_str())) == 0) {
      logError() << "The " << type.c_str() << name << "does not exist in the given Hdf5 file.";
    }
  } else {
    // groups + datasets
    if (_eh(H5Lexists(handles.top(), name.c_str(), H5P_DEFAULT)) == 0) {
      logError() << "The " << type.c_str() << name << "does not exist in the given Hdf5 file.";
    }
  }
}

void Hdf5Reader::openFile(const std::string& name) {
  checkExistence(name, "file");
  const hid_t h5falist = _eh(H5Pcreate(H5P_FILE_ACCESS));
#ifdef H5F_LIBVER_V18
  _eh(H5Pset_libver_bounds(h5falist, H5F_LIBVER_V18, H5F_LIBVER_V18));
#else
  _eh(H5Pset_libver_bounds(h5falist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#endif
  _eh(H5Pset_fapl_mpio(h5falist, comm, MPI_INFO_NULL));
  const hid_t file = _eh(H5Fopen(name.c_str(), H5F_ACC_RDONLY, h5falist));
  _eh(H5Pclose(h5falist));

  handles.push(file);
}
void Hdf5Reader::openGroup(const std::string& name) {
  checkExistence(name, "group");
  const hid_t handle = _eh(H5Gopen(handles.top(), name.c_str(), H5P_DEFAULT));
  handles.push(handle);
}
std::size_t Hdf5Reader::attributeCount(const std::string& name) {
  checkExistence(name, "attribute");
  const hid_t attr = _eh(H5Aopen(handles.top(), name.c_str(), H5P_DEFAULT));
  const hid_t attrspace = _eh(H5Aget_space(attr));
  const hid_t rank = _eh(H5Sget_simple_extent_ndims(attrspace));
  std::vector<hsize_t> dims(rank);
  _eh(H5Sget_simple_extent_dims(attrspace, dims.data(), nullptr));
  _eh(H5Sclose(attrspace));
  _eh(H5Aclose(attr));
  return dims[0];
}
void Hdf5Reader::readAttributeRaw(void* data,
                                  const std::string& name,
                                  const std::shared_ptr<datatype::Datatype>& type) {
  checkExistence(name, "attribute");
  const hid_t attr = _eh(H5Aopen(handles.top(), name.c_str(), H5P_DEFAULT));
  _eh(H5Aread(attr, datatype::convertToHdf5(type), data));
  _eh(H5Aclose(attr));
}
std::size_t Hdf5Reader::dataCount(const std::string& name) {
  checkExistence(name, "dataset");
  const hid_t dataset = _eh(H5Dopen(handles.top(), name.c_str(), H5P_DEFAULT));
  const hid_t dataspace = _eh(H5Dget_space(dataset));
  const hid_t rank = _eh(H5Sget_simple_extent_ndims(dataspace));
  std::vector<hsize_t> dims(rank);
  _eh(H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr));
  _eh(H5Sclose(dataspace));
  _eh(H5Dclose(dataset));

  int mpirank = 0;
  int mpisize = 0;
  MPI_Comm_size(comm, &mpisize);
  MPI_Comm_rank(comm, &mpirank);

  return dims[0] / mpisize + ((dims[0] % mpisize) > mpirank ? 1 : 0);
}
void Hdf5Reader::readDataRaw(void* data,
                             const std::string& name,
                             std::size_t count,
                             const std::shared_ptr<datatype::Datatype>& targetType) {
  checkExistence(name, "dataset");
  const hid_t h5alist = H5Pcreate(H5P_DATASET_XFER);
  _eh(h5alist);
#ifdef USE_MPI
  _eh(H5Pset_dxpl_mpio(h5alist, H5FD_MPIO_COLLECTIVE));
#endif // USE_MPI

  const hid_t dataset = _eh(H5Dopen(handles.top(), name.c_str(), H5P_DEFAULT));
  const hid_t dataspace = _eh(H5Dget_space(dataset));

  const hid_t datatype = datatype::convertToHdf5(targetType);

  const hid_t rank = _eh(H5Sget_simple_extent_ndims(dataspace));
  std::vector<hsize_t> dims(rank);
  _eh(H5Sget_simple_extent_dims(dataspace, dims.data(), nullptr));

  std::size_t dimprod = 1;
  for (std::size_t i = 1; i < dims.size(); ++i) {
    dimprod *= dims[i];
  }

  const std::size_t chunksize =
      std::max(std::size_t(1), std::size_t(2'000'000'000) / (targetType->size() * dimprod));
  std::size_t rounds = (count + chunksize - 1) / chunksize;
  std::size_t start = 0;
  MPI_Allreduce(MPI_IN_PLACE,
                &rounds,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_MAX,
                comm);
  MPI_Exscan(&count,
             &start,
             1,
             datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
             MPI_SUM,
             comm);

  std::vector<hsize_t> nullstart(rank);
  std::vector<hsize_t> readcount(rank);
  std::vector<hsize_t> filepos(rank);

  std::size_t subsize = 1;
  for (std::size_t i = 1; i < readcount.size(); ++i) {
    readcount[i] = dims[i];
    subsize *= dims[i];
  }
  readcount[0] = std::min(chunksize, count);

  filepos[0] = start;

  const hid_t memspace = _eh(H5Screate_simple(rank, readcount.data(), nullptr));

  std::size_t read = 0;

  auto* ptr = reinterpret_cast<unsigned char*>(data);

  for (std::size_t i = 0; i < rounds; ++i) {
    filepos[0] = start + read;
    readcount[0] = std::min(chunksize, count - read);

    _eh(H5Sselect_hyperslab(
        memspace, H5S_SELECT_SET, nullstart.data(), nullptr, readcount.data(), nullptr));
    _eh(H5Sselect_hyperslab(
        dataspace, H5S_SELECT_SET, filepos.data(), nullptr, readcount.data(), nullptr));
    _eh(H5Dread(dataset, datatype, memspace, dataspace, h5alist, ptr));

    read += readcount[0];
    ptr += readcount[0] * subsize * targetType->size();
  }

  _eh(H5Pclose(h5alist));
  _eh(H5Sclose(memspace));
  _eh(H5Sclose(dataspace));
  _eh(H5Dclose(dataset));
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
