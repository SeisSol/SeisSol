// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Reader/Datafield/CellData.h"

#include "IO/Datatype/Inference.h"
#include "IO/Datatype/MPIType.h"
#include "IO/Reader/Distribution.h"
#include "IO/Reader/File/Hdf5Reader.h"
#include "Reader/Datafield/Grid.h"

#include <algorithm>
#include <cstddef>
#include <filesystem>
#include <functional>
#include <mpi.h>
#include <numeric>
#include <optional>
#include <stdexcept>
#include <string>
#include <unordered_set>
#include <utility>
#include <vector>

namespace seissol::reader::datafield {

namespace {

MPI_Datatype sizeType() {
  return io::datatype::convertToMPI(io::datatype::inferDatatype<std::size_t>());
}

std::string described(const GridDesc& desc) {
  return "the data on the mesh `" + desc.name + "` (the variable " + desc.variable + " of " +
         desc.path + ")";
}

} // namespace

std::string CellDataStore::keyOf(const GridDesc& desc) {
  return desc.path + '\x1d' + desc.variable;
}

void CellDataStore::load(const GridDesc& desc, const std::vector<std::size_t>& ids, MPI_Comm comm) {
  // the same on every rank, which all fail alike
  if (!std::filesystem::exists(desc.path)) {
    throw std::runtime_error(described(desc) + ": there is no such file");
  }
  io::reader::file::Hdf5Reader reader(comm);
  reader.openFile(desc.path);
  const auto fail = [&](const std::string& what) {
    reader.closeFile();
    throw std::runtime_error(described(desc) + ": " + what);
  };
  // the variable, and the groups on the way to it
  for (std::size_t slash = desc.variable.find('/', 1);;
       slash = desc.variable.find('/', slash + 1)) {
    const std::string prefix = desc.variable.substr(0, slash);
    if (!prefix.empty() && prefix != "/" && !reader.hasDataset(prefix)) {
      fail("the file has no " + prefix);
    }
    if (slash == std::string::npos) {
      break;
    }
  }
  const std::size_t width = reader.dataRowSize(desc.variable);
  if (width != desc.componentNames.size()) {
    fail("a cell has " + std::to_string(width) + " values in the file, but the declaration names " +
         std::to_string(desc.componentNames.size()) + " components");
  }

  // a row per cell of the mesh, of which every rank reads a part
  const std::size_t count = reader.dataCount(desc.variable);
  std::size_t rows = count;
  MPI_Allreduce(MPI_IN_PLACE, &rows, 1, sizeType(), MPI_SUM, comm);
  std::size_t first = 0;
  MPI_Exscan(&count, &first, 1, sizeType(), MPI_SUM, comm);
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    first = 0;
  }
  // the cells of the mesh: the greatest id of a cell, plus one, as every cell is on some rank
  std::size_t cells = 0;
  for (const auto id : ids) {
    cells = std::max(cells, id + 1);
  }
  MPI_Allreduce(MPI_IN_PLACE, &cells, 1, sizeType(), MPI_MAX, comm);
  if (cells != rows) {
    fail("it has " + std::to_string(rows) + " rows, but the mesh has " + std::to_string(cells) +
         " cells, a row for each in the order of the mesh file");
  }
  std::vector<double> block(count * width);
  reader.readDataRaw(block.data(), desc.variable, count, io::datatype::inferDatatype<double>());
  reader.closeFile();

  // the rows go to the ranks of their cells, once per rank
  std::vector<std::size_t> wanted;
  wanted.reserve(ids.size());
  {
    std::unordered_set<std::size_t> seen;
    seen.reserve(ids.size());
    for (const auto id : ids) {
      if (seen.insert(id).second) {
        wanted.push_back(id);
      }
    }
  }
  std::vector<std::size_t> source(count);
  std::iota(source.begin(), source.end(), first);
  io::reader::Distributor distributor(comm);
  distributor.setup(source, wanted);

  Rows target;
  target.width = width;
  target.values.resize(wanted.size() * width);
  MPI_Datatype row = MPI_DATATYPE_NULL;
  MPI_Type_contiguous(static_cast<int>(width), MPI_DOUBLE, &row);
  MPI_Type_commit(&row);
  distributor
      .distributeRaw(target.values.data(),
                     block.data(),
                     row,
                     row,
                     std::optional<std::function<void(void*, const void*)>>())
      .complete();
  MPI_Type_free(&row);
  target.position.reserve(wanted.size());
  for (std::size_t i = 0; i < wanted.size(); ++i) {
    target.position.emplace(wanted[i], i);
  }
  datasets_[keyOf(desc)] = std::move(target);
}

bool CellDataStore::loaded(const GridDesc& desc) const { return datasets_.count(keyOf(desc)) != 0; }

const double* CellDataStore::row(const GridDesc& desc, std::size_t id) const {
  const auto found = datasets_.find(keyOf(desc));
  if (found == datasets_.end()) {
    throw std::runtime_error(described(desc) +
                             " was not read: data on the mesh is read where SeisSol initializes "
                             "the materials, and bound to their points only");
  }
  const auto position = found->second.position.find(id);
  if (position == found->second.position.end()) {
    throw std::runtime_error(described(desc) + " was not read for the cell " + std::to_string(id));
  }
  return found->second.values.data() + position->second * found->second.width;
}

CellDataStore& sharedCellDataStore() {
  static CellDataStore store;
  return store;
}

} // namespace seissol::reader::datafield
