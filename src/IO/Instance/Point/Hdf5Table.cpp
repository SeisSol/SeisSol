// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Hdf5Table.h"

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Instance/Point/Grouping.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Dimension.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <string>
#include <utility>
#include <vector>

namespace seissol::io::instance::point {

namespace {
//! @brief The compound one sample of a point of this quantity set is written as.
std::shared_ptr<datatype::Datatype> sampleDatatype(const std::vector<TableQuantity>& quantities) {
  std::vector<datatype::StructDatatype::MemberInfo> members(quantities.size());
  std::size_t offset = 0;
  for (std::size_t i = 0; i < quantities.size(); ++i) {
    members[i] =
        datatype::StructDatatype::MemberInfo{quantities[i].name, offset, quantities[i].datatype};
    offset += quantities[i].datatype->size();
  }
  return std::make_shared<datatype::StructDatatype>(members);
}

std::string groupName(std::size_t group) { return "group" + std::to_string(group); }
} // namespace

Hdf5Table::Hdf5Table(std::string name,
                     const std::vector<std::vector<TableQuantity>>& pointQuantities,
                     MPI_Comm comm,
                     std::size_t sampleChunk)
    : name_(std::move(name)), grouping_(groupPoints(pointQuantities, comm)), comm_(comm),
      sampleChunk_(sampleChunk) {
  storage_.resize(grouping_.groupCount());
  samples_.resize(grouping_.groupCount(), 0);
  localPoints_.resize(grouping_.groupCount(), 0);
  for (const auto group : grouping_.group) {
    ++localPoints_[group];
  }
}

const Grouping& Hdf5Table::grouping() const { return grouping_; }

std::size_t Hdf5Table::sampleSize(std::size_t group) const {
  return sampleDatatype(grouping_.quantities.at(group))->size();
}

std::size_t Hdf5Table::localPointCount(std::size_t group) const { return localPoints_.at(group); }

char* Hdf5Table::prepare(std::size_t group, std::size_t samples) {
  samples_.at(group) = samples;
  storage_.at(group).assign(samples * localPoints_.at(group) * sampleSize(group), 0);
  return storage_.at(group).data();
}

void Hdf5Table::clear() {
  for (std::size_t group = 0; group < grouping_.groupCount(); ++group) {
    samples_[group] = 0;
    storage_[group].clear();
  }
}

std::function<writer::Writer(const std::string&, std::size_t, double)> Hdf5Table::makeWriter() {
  return [this](const std::string& prefix, std::size_t counter, double /*time*/) -> writer::Writer {
    const auto filename = prefix + "-" + name_ + ".h5";
    auto writer = writer::Writer();

    for (std::size_t group = 0; group < grouping_.groupCount(); ++group) {
      const auto& quantities = grouping_.quantities[group];
      const auto datatype = sampleDatatype(quantities);

      const std::vector<writer::Dimension> shape{
          writer::Dimension::appended(samples_[group], sampleChunk_),
          writer::Dimension::distributed()};

      // a row is one sample of one point, which is what the compound holds
      const auto rows = samples_[group] * localPoints_[group];
      const auto data =
          writer::WriteBuffer::createShaped<char>(storage_[group].data(), rows, shape, datatype);

      writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {name_}), groupName(group), data, datatype));

      if (counter == 0) {
        // what the compound of this group is, spelled out the way the grouping reads it back, and
        // how many points there are over all ranks -- neither follows from the dataset alone
        // before it has been written to for the first time
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {name_}, groupName(group)),
            "Quantities",
            writer::WriteInline::createString(quantitySetKey(quantities))));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5AttributeWrite>(
            writer::instructions::Hdf5Location(filename, {name_}, groupName(group)),
            "NumberOfPoints",
            writer::WriteInline::create<std::int64_t>(
                static_cast<std::int64_t>(grouping_.globalCount[group]))));
      }
    }

    if (counter == 0) {
      // Where a point of the caller's numbering ended up. The grouping renumbers the points so
      // that every rank owns one run of each group, so without this there is no way back from a
      // row of a dataset to the point it belongs to.
      index_.resize(grouping_.group.size() * 2);
      for (std::size_t point = 0; point < grouping_.group.size(); ++point) {
        index_[point * 2] = static_cast<std::uint64_t>(grouping_.group[point]);
        index_[point * 2 + 1] = static_cast<std::uint64_t>(grouping_.index[point]);
      }
      writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {name_}),
          "Index",
          writer::WriteBuffer::create(index_.data(), grouping_.group.size(), {2}),
          datatype::inferDatatype<std::uint64_t>()));

      for (const auto& entry : pointData_) {
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, {name_}),
            entry.name,
            writer::WriteBuffer::create(
                entry.bytes.data(), grouping_.group.size(), entry.shape, entry.datatype),
            entry.datatype));
      }
    }

    return writer;
  };
}

} // namespace seissol::io::instance::point
