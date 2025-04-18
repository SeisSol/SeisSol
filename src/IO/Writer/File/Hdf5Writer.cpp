// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Hdf5Writer.h"

#include <IO/Datatype/Datatype.h>
#include <IO/Datatype/HDF5Type.h>
#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <algorithm>
#include <async/ExecInfo.h>
#include <cassert>
#include <cstddef>
#include <hdf5.h>
#include <memory>
#include <mpi.h>
#include <stack>
#include <string>
#include <utility>
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

#define h5async(function, ...)                                                                     \
  [&]() {                                                                                          \
    if (asyncEnabled) {                                                                            \
      return _eh(function##_async(__VA_ARGS__, asyncHandle));                                      \
    } else {                                                                                       \
      return _eh(function(__VA_ARGS__));                                                           \
    }                                                                                              \
  }()
} // namespace

namespace seissol::io::writer::file {

Hdf5File::Hdf5File(MPI_Comm comm) : comm(comm) {}

void Hdf5File::openFile(const std::string& name) {
  const hid_t h5falist = _eh(H5Pcreate(H5P_FILE_ACCESS));
#ifdef H5F_LIBVER_V18
  _eh(H5Pset_libver_bounds(h5falist, H5F_LIBVER_V18, H5F_LIBVER_V18));
#else
  _eh(H5Pset_libver_bounds(h5falist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST));
#endif
  _eh(H5Pset_fapl_mpio(h5falist, comm, MPI_INFO_NULL));
  file = _eh(H5Fcreate(name.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, h5falist));
  _eh(H5Pclose(h5falist));

  handles.push(file);
}
void Hdf5File::openGroup(const std::string& name) {
  // cf. https://stackoverflow.com/a/18468735
  auto existenceTest = _eh(H5Lexists(handles.top(), name.c_str(), H5P_DEFAULT));
  hid_t handle = 0;
  if (existenceTest > 0) {
    handle = _eh(H5Gopen(handles.top(), name.c_str(), H5P_DEFAULT));
  } else {
    handle = _eh(H5Gcreate(handles.top(), name.c_str(), H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT));
  }
  handles.push(handle);
}
void Hdf5File::openDataset(const std::string& name) {
  // cf. https://stackoverflow.com/a/18468735
  auto existenceTest = _eh(H5Lexists(handles.top(), name.c_str(), H5P_DEFAULT));
  const hid_t handle = _eh(H5Dopen(handles.top(), name.c_str(), H5P_DEFAULT));
  handles.push(handle);
}
void Hdf5File::writeAttribute(const async::ExecInfo& info,
                              const std::string& name,
                              const std::shared_ptr<DataSource>& source) {
  // TODO: store datatype separately
  if (source->distributed()) {
    logError() << "Attempted to write an HDF5 attributed as a distributed source.";
  }
  hid_t h5space = 0;
  if (source->shape().empty()) {
    h5space = _eh(H5Screate(H5S_SCALAR));
  } else {
    const auto& shape = source->shape();
    std::vector<hsize_t> hshape(shape.begin(), shape.end());
    h5space = _eh(H5Screate_simple(shape.size(), hshape.data(), nullptr));
  }
  const hid_t h5type = datatype::convertToHdf5(source->datatype());
  const hid_t handle =
      _eh(H5Acreate(handles.top(), name.c_str(), h5type, h5space, H5P_DEFAULT, H5P_DEFAULT));
  _eh(H5Awrite(handle, h5type, source->getPointer(info)));
  _eh(H5Aclose(handle));
  _eh(H5Sclose(h5space));
}
void Hdf5File::writeData(const async::ExecInfo& info,
                         const std::string& name,
                         const std::shared_ptr<DataSource>& source,
                         const std::shared_ptr<datatype::Datatype>& targetType,
                         int compress) {
  MPI_Datatype sizetype = datatype::convertToMPI(datatype::inferDatatype<std::size_t>());

  std::size_t trueCount = source->count(info);
  std::size_t dimprod = 1;
  for (auto dimension : source->shape()) {
    assert(trueCount % dimension == 0);
    trueCount /= dimension;
    dimprod *= dimension;
  }

  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  // if we don't write distributed data, only one rank needs to do the work
  const std::size_t count = (source->distributed() || rank == 0) ? trueCount : 0;
  const auto& dimensions = source->shape();
  // TODO: adjust chunksize according to dimensions and datatype size
  const std::size_t chunksize =
      std::max(std::size_t(1), std::size_t(2'000'000'000) / (source->datatype()->size() * dimprod));

  const std::size_t actualDimensions =
      source->distributed() ? dimensions.size() + 1 : dimensions.size();

  std::size_t allcount = count;
  std::size_t offset = 0;

  std::size_t localRounds = (count + chunksize - 1) / chunksize;
  std::size_t rounds = localRounds;

  MPI_Allreduce(&localRounds, &rounds, 1, sizetype, MPI_MAX, comm);
  MPI_Allreduce(&count, &allcount, 1, sizetype, MPI_SUM, comm);
  MPI_Exscan(&count, &offset, 1, sizetype, MPI_SUM, comm);

  std::vector<hsize_t> globalSizes;
  std::vector<hsize_t> localSizes;

  const std::size_t chunkcount = std::min(chunksize, count);

  if (source->distributed()) {
    globalSizes.push_back(allcount);
    localSizes.push_back(chunkcount);
  }
  for (const auto& dim : dimensions) {
    globalSizes.push_back(dim);
    localSizes.push_back(dim);
  }
  const hid_t h5space = _eh(H5Screate_simple(globalSizes.size(), globalSizes.data(), nullptr));
  const hid_t h5memspace = _eh(H5Screate_simple(localSizes.size(), localSizes.data(), nullptr));

  // create empty spaces just in case (MPIO likes to get stuck otherwise)
  const hid_t h5spaceEmpty = _eh(H5Screate(H5S_NULL));
  const hid_t h5memspaceEmpty = _eh(H5Screate(H5S_NULL));

  std::vector<hsize_t> writeStart;
  std::vector<hsize_t> writeLength;

  if (source->distributed()) {
    writeStart.push_back(offset);
    writeLength.push_back(chunkcount);
  }
  for (const auto& dim : dimensions) {
    writeStart.push_back(0);
    writeLength.push_back(dim);
  }

  const hid_t h5dxlist = H5Pcreate(H5P_DATASET_XFER);
  _eh(h5dxlist);
#ifdef USE_MPI
  _eh(H5Pset_dxpl_mpio(h5dxlist, H5FD_MPIO_COLLECTIVE));
#endif

  const hid_t h5memtype = datatype::convertToHdf5(source->datatype());

  const hid_t preh5type = datatype::convertToHdf5(targetType);
  const hid_t h5type = _eh(H5Tcopy(preh5type));
  if (_eh(H5Tget_class(h5type)) == H5T_COMPOUND) {
    _eh(H5Tpack(h5type));
  }
  _eh(H5Tcommit(handles.top(),
                (name + std::string("_Type")).c_str(),
                h5type,
                H5P_DEFAULT,
                H5P_DEFAULT,
                H5P_DEFAULT));

  hid_t h5filter = H5P_DEFAULT;
  if (compress > 0) {
    h5filter = _eh(H5Pcreate(H5P_DATASET_CREATE));
    _eh(H5Pset_chunk(h5filter, actualDimensions, writeLength.data()));
    const int deflateStrength = compress;
    _eh(H5Pset_deflate(h5filter, deflateStrength));
  }

  const hid_t h5data =
      H5Dcreate(handles.top(), name.c_str(), h5type, h5space, H5P_DEFAULT, h5filter, H5P_DEFAULT);

  std::size_t written = 0;

  std::vector<hsize_t> nullstart(actualDimensions);

  const char* data = reinterpret_cast<const char*>(source->getPointer(info));

  // TODO(David): maybe always compute in the loop instead?
  std::size_t baseWriteSize = targetType->size();
  for (const auto& dim : dimensions) {
    baseWriteSize *= dim;
  }

  const char* dataloc = data;

  for (std::size_t i = 0; i < rounds; ++i) {
    if (source->distributed()) {
      writeStart[0] = offset + written;
      writeLength[0] = std::min(count - written, chunkcount);
    }

    const auto [memspace, space] = [&]() -> std::pair<hid_t, hid_t> {
      if (nullstart.empty()) {
        return {h5memspace, h5space};
      } else if (source->distributed() && writeLength[0] > 0) {
        _eh(H5Sselect_hyperslab(
            h5memspace, H5S_SELECT_SET, nullstart.data(), nullptr, writeLength.data(), nullptr));

        _eh(H5Sselect_hyperslab(
            h5space, H5S_SELECT_SET, writeStart.data(), nullptr, writeLength.data(), nullptr));

        return {h5memspace, h5space};
      } else {
        return {h5memspaceEmpty, h5spaceEmpty};
      }
    }();

    _eh(H5Dwrite(h5data, h5memtype, memspace, space, h5dxlist, dataloc));

    std::size_t writeSize = baseWriteSize;
    if (source->distributed()) {
      written += writeLength[0];
      writeSize *= writeLength[0];
    }
    dataloc += writeSize;
  }

  if (compress > 0) {
    _eh(H5Pclose(h5filter));
  }
  _eh(H5Tclose(h5type));
  _eh(H5Sclose(h5space));
  _eh(H5Sclose(h5memspace));
  _eh(H5Sclose(h5spaceEmpty));
  _eh(H5Sclose(h5memspaceEmpty));
  _eh(H5Dclose(h5data));
  _eh(H5Pclose(h5dxlist));
}
void Hdf5File::closeDataset() {
  _eh(H5Dclose(handles.top()));
  handles.pop();
}
void Hdf5File::closeGroup() {
  _eh(H5Gclose(handles.top()));
  handles.pop();
}
void Hdf5File::closeFile() {
  _eh(H5Fclose(file));
  handles.pop();
}

Hdf5Writer::Hdf5Writer(MPI_Comm comm) : comm(comm) {}

void Hdf5Writer::writeAttribute(const async::ExecInfo& info,
                                const instructions::Hdf5AttributeWrite& write) {
  Hdf5File file(comm);
  if (openFiles.find(write.location.file()) == openFiles.end()) {
    file.openFile(write.location.file());
    openFiles.insert({write.location.file(), file});
  }
  file = openFiles.at(write.location.file());
  for (const auto& groupname : write.location.groups()) {
    file.openGroup(groupname);
  }
  if (write.location.dataset().has_value()) {
    file.openDataset(write.location.dataset().value());
  }
  file.writeAttribute(info, write.name, write.dataSource);
  if (write.location.dataset().has_value()) {
    file.closeDataset();
  }
  for (auto _ : write.location.groups()) {
    file.closeGroup();
  }
}

void Hdf5Writer::writeData(const async::ExecInfo& info, const instructions::Hdf5DataWrite& write) {
  Hdf5File file(comm);
  if (openFiles.find(write.location.file()) == openFiles.end()) {
    file.openFile(write.location.file());
    openFiles.insert({write.location.file(), file});
  }
  file = openFiles.at(write.location.file());
  for (const auto& groupname : write.location.groups()) {
    file.openGroup(groupname);
  }
  if (write.location.dataset().has_value()) {
    file.openDataset(write.location.dataset().value());
  }
  file.writeData(info, write.name, write.dataSource, write.targetType, write.compress);
  if (write.location.dataset().has_value()) {
    file.closeDataset();
  }
  for (auto _ : write.location.groups()) {
    file.closeGroup();
  }
}

void Hdf5Writer::finalize() {
  for (auto [_, file] : openFiles) {
    file.closeFile();
  }
}

} // namespace seissol::io::writer::file
