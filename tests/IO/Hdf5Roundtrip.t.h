// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Datatype/Inference.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "Parallel/MPI.h"

#include <cstdio>
#include <cstring>
#include <hdf5.h>
#include <mpi.h>
#include <string>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::io;

// Helper RAII for temp HDF5 files
struct TempHdf5File {
  std::string path;
  TempHdf5File() {
    // Create a unique temp file name per rank
    int rank = 0;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    path = "/tmp/seissol_test_hdf5_r" + std::to_string(rank) + ".h5";
  }
  ~TempHdf5File() { std::remove(path.c_str()); }
};

// Helper: create a parallel HDF5 file with given MPI communicator
static hid_t createParallelFile(const std::string& path, MPI_Comm comm) {
  hid_t plist = H5Pcreate(H5P_FILE_ACCESS);
#ifdef H5F_LIBVER_V18
  H5Pset_libver_bounds(plist, H5F_LIBVER_V18, H5F_LIBVER_V18);
#else
  H5Pset_libver_bounds(plist, H5F_LIBVER_LATEST, H5F_LIBVER_LATEST);
#endif
  H5Pset_fapl_mpio(plist, comm, MPI_INFO_NULL);
  hid_t file = H5Fcreate(path.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist);
  H5Pclose(plist);
  return file;
}

// ---------------------------------------------------------------------------
// Roundtrip: write scalar attribute, read back with Hdf5Reader
// ---------------------------------------------------------------------------

TEST_CASE("HDF5 roundtrip scalar attribute" * doctest::test_suite("io")) {
  TempHdf5File tmp;
  const double writeValue = 3.14159;

  // Write with raw HDF5 API
  {
    hid_t file = createParallelFile(tmp.path, MPI_COMM_WORLD);
    REQUIRE(file >= 0);

    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(file, "pi", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_DOUBLE, &writeValue);
    H5Aclose(attr);
    H5Sclose(space);
    H5Fclose(file);
  }

  // Read back with Hdf5Reader
  {
    reader::file::Hdf5Reader reader(MPI_COMM_WORLD);
    reader.openFile(tmp.path);

    double readValue = reader.readAttributeScalar<double>("pi");
    CHECK(readValue == doctest::Approx(writeValue));

    reader.closeFile();
  }
}

// ---------------------------------------------------------------------------
// Roundtrip: write vector attribute, read back
// ---------------------------------------------------------------------------

TEST_CASE("HDF5 roundtrip vector attribute" * doctest::test_suite("io")) {
  TempHdf5File tmp;
  const std::vector<int> writeData = {10, 20, 30, 40, 50};

  {
    hid_t file = createParallelFile(tmp.path, MPI_COMM_WORLD);
    REQUIRE(file >= 0);

    hsize_t dims[1] = {writeData.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t attr = H5Acreate(file, "numbers", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    H5Awrite(attr, H5T_NATIVE_INT, writeData.data());
    H5Aclose(attr);
    H5Sclose(space);
    H5Fclose(file);
  }

  {
    reader::file::Hdf5Reader reader(MPI_COMM_WORLD);
    reader.openFile(tmp.path);

    auto readData = reader.readAttribute<int>("numbers");
    REQUIRE(readData.size() == writeData.size());
    for (std::size_t i = 0; i < writeData.size(); ++i) {
      CHECK(readData[i] == writeData[i]);
    }

    reader.closeFile();
  }
}

// ---------------------------------------------------------------------------
// Roundtrip: write 1D dataset, read back
// ---------------------------------------------------------------------------

TEST_CASE("HDF5 roundtrip 1D dataset" * doctest::test_suite("io")) {
  TempHdf5File tmp;
  const std::vector<double> writeData = {1.1, 2.2, 3.3, 4.4, 5.5, 6.6, 7.7, 8.8};

  {
    hid_t file = createParallelFile(tmp.path, MPI_COMM_WORLD);
    REQUIRE(file >= 0);

    hsize_t dims[1] = {writeData.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dset =
        H5Dcreate(file, "data", H5T_NATIVE_DOUBLE, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    // Write collectively
    hid_t xfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, xfer, writeData.data());
    H5Pclose(xfer);
    H5Dclose(dset);
    H5Sclose(space);
    H5Fclose(file);
  }

  {
    reader::file::Hdf5Reader reader(MPI_COMM_WORLD);
    reader.openFile(tmp.path);

    auto readData = reader.readData<double>("data");
    REQUIRE(readData.size() == writeData.size());
    for (std::size_t i = 0; i < writeData.size(); ++i) {
      CHECK(readData[i] == doctest::Approx(writeData[i]));
    }

    reader.closeFile();
  }
}

// ---------------------------------------------------------------------------
// Roundtrip: group navigation
// ---------------------------------------------------------------------------

TEST_CASE("HDF5 roundtrip with groups" * doctest::test_suite("io")) {
  TempHdf5File tmp;

  {
    hid_t file = createParallelFile(tmp.path, MPI_COMM_WORLD);
    REQUIRE(file >= 0);

    hid_t grp = H5Gcreate(file, "mygroup", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    hid_t space = H5Screate(H5S_SCALAR);
    hid_t attr = H5Acreate(grp, "answer", H5T_NATIVE_INT, space, H5P_DEFAULT, H5P_DEFAULT);
    int val = 42;
    H5Awrite(attr, H5T_NATIVE_INT, &val);
    H5Aclose(attr);
    H5Sclose(space);
    H5Gclose(grp);
    H5Fclose(file);
  }

  {
    reader::file::Hdf5Reader reader(MPI_COMM_WORLD);
    reader.openFile(tmp.path);
    reader.openGroup("mygroup");

    int readVal = reader.readAttributeScalar<int>("answer");
    CHECK(readVal == 42);

    reader.closeGroup();
    reader.closeFile();
  }
}

// ---------------------------------------------------------------------------
// Roundtrip: integer dataset
// ---------------------------------------------------------------------------

TEST_CASE("HDF5 roundtrip integer dataset" * doctest::test_suite("io")) {
  TempHdf5File tmp;
  const std::vector<long long> writeData = {-100, 0, 100, 999999999LL};

  {
    hid_t file = createParallelFile(tmp.path, MPI_COMM_WORLD);
    REQUIRE(file >= 0);

    hsize_t dims[1] = {writeData.size()};
    hid_t space = H5Screate_simple(1, dims, nullptr);
    hid_t dset =
        H5Dcreate(file, "integers", H5T_NATIVE_LLONG, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    hid_t xfer = H5Pcreate(H5P_DATASET_XFER);
    H5Pset_dxpl_mpio(xfer, H5FD_MPIO_COLLECTIVE);
    H5Dwrite(dset, H5T_NATIVE_LLONG, H5S_ALL, H5S_ALL, xfer, writeData.data());
    H5Pclose(xfer);
    H5Dclose(dset);
    H5Sclose(space);
    H5Fclose(file);
  }

  {
    reader::file::Hdf5Reader reader(MPI_COMM_WORLD);
    reader.openFile(tmp.path);

    auto readData = reader.readData<long long>("integers");
    REQUIRE(readData.size() == writeData.size());
    for (std::size_t i = 0; i < writeData.size(); ++i) {
      CHECK(readData[i] == writeData[i]);
    }

    reader.closeFile();
  }
}

} // namespace seissol::unit_test
