// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Grouping.h"

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Datatype/MPIType.h"
#include "IO/Instance/Point/TableWriter.h"

#include <algorithm>
#include <cstddef>
#include <map>
#include <mpi.h>
#include <set>
#include <sstream>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::instance::point {

std::string quantitySetKey(const std::vector<TableQuantity>& quantities) {
  YAML::Node node;
  for (const auto& quantity : quantities) {
    YAML::Node entry;
    entry["name"] = quantity.name;
    entry["datatype"] = quantity.datatype->serialize();
    node.push_back(entry);
  }
  std::ostringstream stream;
  YAML::Emitter emitter(stream);
  emitter << YAML::Flow << node;
  return stream.str();
}

std::vector<TableQuantity> quantitySetFromKey(const std::string& key) {
  const auto node = YAML::Load(key);
  std::vector<TableQuantity> quantities;
  quantities.reserve(node.size());
  for (const auto& entry : node) {
    quantities.push_back(TableQuantity{entry["name"].as<std::string>(),
                                       datatype::Datatype::deserialize(entry["datatype"])});
  }
  return quantities;
}

namespace {

//! Every key that occurs on any rank, in an order all ranks agree on.
std::vector<std::string> gatherKeys(const std::set<std::string>& local, MPI_Comm comm) {
  std::string packed;
  for (const auto& key : local) {
    packed += key;
    packed.push_back('\0');
  }

  int size = 0;
  MPI_Comm_size(comm, &size);
  std::vector<int> lengths(size);
  auto length = static_cast<int>(packed.size());
  MPI_Allgather(&length, 1, MPI_INT, lengths.data(), 1, MPI_INT, comm);

  std::vector<int> offsets(size, 0);
  for (int rank = 1; rank < size; ++rank) {
    offsets[rank] = offsets[rank - 1] + lengths[rank - 1];
  }
  std::vector<char> all(offsets.back() + lengths.back());
  MPI_Allgatherv(
      packed.data(), length, MPI_CHAR, all.data(), lengths.data(), offsets.data(), MPI_CHAR, comm);

  // sorted, so that the order does not depend on which rank saw a set first
  std::set<std::string> unique;
  std::size_t start = 0;
  for (std::size_t position = 0; position < all.size(); ++position) {
    if (all[position] == '\0') {
      unique.emplace(all.data() + start, position - start);
      start = position + 1;
    }
  }
  return {unique.begin(), unique.end()};
}

} // namespace

Grouping groupPoints(const std::vector<std::vector<TableQuantity>>& pointQuantities,
                     MPI_Comm comm) {
  std::vector<std::string> localKeys(pointQuantities.size());
  std::set<std::string> distinct;
  for (std::size_t point = 0; point < pointQuantities.size(); ++point) {
    localKeys[point] = quantitySetKey(pointQuantities[point]);
    distinct.emplace(localKeys[point]);
  }

  const auto keys = gatherKeys(distinct, comm);

  Grouping grouping;
  grouping.quantities.reserve(keys.size());
  std::map<std::string, std::size_t> groupOf;
  for (std::size_t index = 0; index < keys.size(); ++index) {
    groupOf.emplace(keys[index], index);
    grouping.quantities.push_back(quantitySetFromKey(keys[index]));
  }

  grouping.group.resize(pointQuantities.size());
  std::vector<std::size_t> localCount(keys.size(), 0);
  for (std::size_t point = 0; point < pointQuantities.size(); ++point) {
    const auto group = groupOf.at(localKeys[point]);
    grouping.group[point] = group;
    ++localCount[group];
  }

  MPI_Datatype sizetype = datatype::convertToMPI(datatype::inferDatatype<std::size_t>());
  std::vector<std::size_t> offset(keys.size(), 0);
  grouping.globalCount.assign(keys.size(), 0);
  if (!keys.empty()) {
    MPI_Exscan(
        localCount.data(), offset.data(), static_cast<int>(keys.size()), sizetype, MPI_SUM, comm);
    MPI_Allreduce(localCount.data(),
                  grouping.globalCount.data(),
                  static_cast<int>(keys.size()),
                  sizetype,
                  MPI_SUM,
                  comm);
  }
  int rank = 0;
  MPI_Comm_rank(comm, &rank);
  if (rank == 0) {
    // MPI_Exscan leaves the result untouched on the first rank
    std::fill(offset.begin(), offset.end(), 0);
  }

  grouping.index.resize(pointQuantities.size());
  std::vector<std::size_t> running(keys.size(), 0);
  for (std::size_t point = 0; point < pointQuantities.size(); ++point) {
    const auto group = grouping.group[point];
    grouping.index[point] = offset[group] + running[group];
    ++running[group];
  }

  return grouping;
}

} // namespace seissol::io::instance::point
