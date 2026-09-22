// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "VtkHdf.h"

#include "Common/Filesystem.h"
#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Datatype/MPIType.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Metadata/Pvd.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"
#include "Parallel/MPI.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <optional>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace seissol::io::instance::mesh {
VtkHdfWriter::VtkHdfWriter(const std::string& name,
                           std::size_t localElementCount,
                           geometry::Shape shape,
                           std::size_t targetDegree,
                           bool temporal,
                           std::int32_t compress,
                           bool constFile,
                           std::optional<VertexMap> vertexMap)
    : name_(name), localElementCount_(localElementCount), globalElementCount_(localElementCount),
      pointsPerElement_(
          geometry::numPoints(std::max(targetDegree, static_cast<std::size_t>(1)), shape)),
      type_(geometry::vtkType(shape)), targetDegree_(targetDegree), constFile_(constFile),
      temporal_(temporal), compress_(compress) {
  MPI_Exscan(&localElementCount,
             &elementOffset_,
             1,
             datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
             MPI_SUM,
             seissol::Mpi::mpi.comm());
  MPI_Allreduce(&localElementCount,
                &globalElementCount_,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::Mpi::mpi.comm());
  connectivityOffset_ = elementOffset_ * pointsPerElement_;
  localPointCount_ =
      vertexMap.has_value() ? vertexMap->localPointCount : localElementCount * pointsPerElement_;
  if (vertexMap.has_value()) {
    // shared points, so how many this rank has is no longer a multiple of the cell count.
    // MPI_Exscan leaves the result untouched on rank 0, so it has to start at zero.
    pointOffset_ = 0;
    globalPointCount_ = localPointCount_;
    MPI_Exscan(&localPointCount_,
               &pointOffset_,
               1,
               datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
               MPI_SUM,
               seissol::Mpi::mpi.comm());
    MPI_Allreduce(&localPointCount_,
                  &globalPointCount_,
                  1,
                  datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                  MPI_SUM,
                  seissol::Mpi::mpi.comm());
    pointSourceCount_ = localPointCount_;
    pointsPerSource_ = 1;
  } else {
    pointOffset_ = connectivityOffset_;
    globalPointCount_ = globalElementCount_ * pointsPerElement_;
    pointSourceCount_ = localElementCount_;
    pointsPerSource_ = pointsPerElement_;
  }

  const auto version = temporal ? std::vector<int64_t>{2, 0} : std::vector<int64_t>{1, 0};

  addData("Type",
          {},
          temporal,
          writer::WriteInline::create("UnstructuredGrid",
                                      std::make_shared<datatype::StringDatatype>(16)),
          true);
  addData("Version",
          {},
          temporal,
          writer::WriteInline::createArray<int64_t>({version.size()}, version),
          true);

  // to capture by value
  const auto selfGlobalElementCount = globalElementCount_;
  const auto selfLocalElementCount = localElementCount_;
  const auto selfGlobalPointCount = globalPointCount_;
  const auto selfLocalPointCount = localPointCount_;
  const auto selfPointOffset = pointOffset_;
  const auto selfConnectivityOffset = connectivityOffset_;
  const auto selfPointsPerElement = pointsPerElement_;
  const auto selfType = type_;

  // TODO: auto-generate using a managed buffer maybe?

  addData("NumberOfCells",
          {},
          temporal,
          writer::WriteInline::createArray<int64_t>(
              {1}, {static_cast<int64_t>(selfGlobalElementCount)}));
  // one entry per corner of every cell, which stays the same when the points are shared
  addData("NumberOfConnectivityIds",
          {},
          temporal,
          writer::WriteInline::createArray<int64_t>(
              {1}, {static_cast<int64_t>(selfGlobalElementCount * selfPointsPerElement)}));
  addData(
      "NumberOfPoints",
      {},
      temporal,
      writer::WriteInline::createArray<int64_t>({1}, {static_cast<int64_t>(selfGlobalPointCount)}));

  const bool isLastRank = Mpi::mpi.size() == Mpi::mpi.rank() + 1;
  addData("Offsets",
          {},
          true,
          writer::GeneratedBuffer::createElementwise<int64_t>(
              selfLocalElementCount + (isLastRank ? 1 : 0),
              1,
              std::vector<std::size_t>(),
              [=](int64_t* target, std::size_t index) {
                target[0] = index * selfPointsPerElement + selfConnectivityOffset;
              }));
  addData("Types",
          {},
          true,
          writer::GeneratedBuffer::createElementwise<uint8_t>(
              selfLocalElementCount,
              1,
              std::vector<std::size_t>(),
              [=](uint8_t* target, std::size_t /*index*/) { target[0] = selfType; }));
  addData(
      "Connectivity",
      {},
      true,
      vertexMap.has_value()
          ? writer::GeneratedBuffer::createElementwise<int64_t>(
                selfLocalElementCount,
                selfPointsPerElement,
                std::vector<std::size_t>(),
                [=, map = std::move(vertexMap->connectivity)](int64_t* target, std::size_t index) {
                  for (std::size_t corner = 0; corner < selfPointsPerElement; ++corner) {
                    target[corner] = static_cast<int64_t>(
                        map[index * selfPointsPerElement + corner] + selfPointOffset);
                  }
                })
          : writer::GeneratedBuffer::createElementwise<int64_t>(
                selfLocalPointCount,
                1,
                std::vector<std::size_t>(),
                [=](int64_t* target, std::size_t index) {
                  target[0] = static_cast<int64_t>(index + selfPointOffset);
                }));

  if (temporal) {
    // https://docs.vtk.org/en/latest/vtk_file_formats/vtkhdf_file_format/vtkhdf_specifications.html
    // The mesh itself is written once and every step reads it again, so all of the geometry
    // offsets stay at zero ("Offset value can be repeated for static data"); only the attribute
    // data grows.
    instructions_.emplace_back(
        [](const std::string& filename, std::size_t counter, double /*time*/) {
          return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
              writer::instructions::Hdf5Location(filename, {GroupName, StepsName}),
              "NSteps",
              writer::WriteInline::createArray<int64_t>({}, {static_cast<int64_t>(counter) + 1}));
        });

    instructions_.emplace_back(
        [](const std::string& filename, std::size_t /*counter*/, double time) {
          const auto data =
              writer::WriteInline::createShaped<double>({writer::Dimension::appended(1)}, {time});
          return std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, {GroupName, StepsName}),
              "Values",
              data,
              data->datatype());
        });

    // NumberOfPoints and friends hold one entry for the single part the mesh was written as, not
    // one per step, so the part count has to be given explicitly rather than inferred.
    instructions_.emplace_back(
        [](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
          const auto data =
              writer::WriteInline::createShaped<int64_t>({writer::Dimension::appended(1)}, {1});
          return std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, {GroupName, StepsName}),
              "NumberOfParts",
              data,
              data->datatype());
        });

    addStepOffset("PartOffsets", {}, 0);
    addStepOffset("PointOffsets", {}, 0);
    // CellOffsets and ConnectivityIdOffsets are (NSteps, NTopologies), and an unstructured grid
    // has one topology
    addStepOffset("CellOffsets", {}, 0);
    addStepOffset("ConnectivityIdOffsets", {}, 0);

    // the two field data entries makeWriter always adds; both are single scalars
    addStepOffset("Time", FieldDataName + "Offsets", 1);
    addStepOffset("Index", FieldDataName + "Offsets", 1);
    addStepFieldDataSize("Time", 1, 1);
    addStepFieldDataSize("Index", 1, 1);
  }
}

void VtkHdfWriter::addStepOffset(const std::string& name,
                                 const std::optional<std::string>& group,
                                 std::size_t perStep) {
  std::vector<std::string> groups{GroupName, StepsName};
  if (group.has_value()) {
    groups.emplace_back(group.value());
  }
  // one entry per step, and CellOffsets and ConnectivityIdOffsets carry one per topology on top
  // of that, of which an unstructured grid has one
  std::vector<writer::Dimension> dimensions{writer::Dimension::appended(1)};
  if (name == "CellOffsets" || name == "ConnectivityIdOffsets") {
    dimensions.push_back(writer::Dimension::replicated(1));
  }

  instructions_.emplace_back(
      [=](const std::string& filename, std::size_t counter, double /*time*/) {
        const auto data = writer::WriteInline::createShaped<uint64_t>(
            dimensions, {static_cast<uint64_t>(counter * perStep)});
        return std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, groups), name, data, data->datatype());
      });
}

void VtkHdfWriter::addData(const std::string& name,
                           const std::optional<std::string>& group,
                           bool isConst,
                           const std::shared_ptr<writer::DataSource>& data,
                           bool attribute) {
  auto& instrarray = isConst && (constFile_ || temporal_) ? instructionsConst_ : instructions_;

  std::vector<std::string> groups{GroupName};
  if (group.has_value()) {
    groups.emplace_back(group.value());
  }

  const auto compress = this->compress_;

  if (attribute) {
    instrarray.emplace_back(
        [=](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
          return std::make_shared<writer::instructions::Hdf5AttributeWrite>(
              writer::instructions::Hdf5Location(filename, groups), name, data);
        });
  } else {
    instrarray.emplace_back(
        [=](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
          return std::make_shared<writer::instructions::Hdf5DataWrite>(
              writer::instructions::Hdf5Location(filename, groups),
              name,
              data,
              data->datatype(),
              compress);
        });

    if (isConst && constFile_ && !temporal_) {
      instructionsConstLink_.emplace_back(
          [=](const std::string& filename, const std::string& filenameConst) {
            return std::make_shared<writer::instructions::Hdf5LinkExternalWrite>(
                writer::instructions::Hdf5Location(filename, groups),
                name,
                writer::instructions::Hdf5Location(filenameConst, groups, name));
          });
    }
  }
}

void VtkHdfWriter::addStepFieldDataSize(const std::string& name,
                                        std::size_t components,
                                        std::size_t tuples) {
  const std::vector<std::string> groups{GroupName, StepsName, FieldDataName + "Sizes"};
  instructions_.emplace_back(
      [=](const std::string& filename, std::size_t /*counter*/, double /*time*/) {
        // (NSteps, 2), the appended dimension supplying the leading one
        const auto data = writer::WriteInline::createShaped<int64_t>(
            {writer::Dimension::appended(1), writer::Dimension::replicated(2)},
            {static_cast<int64_t>(components), static_cast<int64_t>(tuples)});
        return std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, groups), name, data, data->datatype());
      });
}

void VtkHdfWriter::addHook(const std::function<void(std::size_t, double)>& hook) {
  hooks_.push_back(hook);
}

std::function<writer::Writer(const std::string&, std::size_t, double)> VtkHdfWriter::makeWriter() {
  logInfo() << "Adding VTK writer" << name_ << "of order" << targetDegree_;
  const auto self = *this;
  return
      [self, pvu = std::vector<metadata::PvuEntry>(), constCounter = std::optional<std::size_t>()](
          const std::string& prefix, std::size_t counter, double time) mutable -> writer::Writer {
        // The unchanging data is written once per run, not once per counter: a run resuming from a
        // checkpoint starts at whatever counter the schedule gives it and would otherwise link to a
        // file nobody wrote. Its name carries that counter, so the const file of one run cannot
        // collide with the one of a previous run -- rewriting it in place would fail, since none of
        // its datasets are appendable.
        const auto fullWrite = !constCounter.has_value();
        if (fullWrite) {
          constCounter = counter;
        }
        for (const auto& hook : self.hooks_) {
          hook(counter, time);
        }

        // the .pvd and the links between the files name them relative to the directory they are
        // in, which is the one of the file that refers to them
        const auto inDirectory = [](const std::string& path) {
          return seissol::filesystem::path(path).filename().string();
        };
        // a time series is one file holding every step; a snapshot is one file per step
        const auto suffix = self.temporal_ ? std::string() : "-" + std::to_string(counter);
        const auto filename = prefix + "-" + self.name_ + suffix + ".vtkhdf";
        const auto filenameFile = inDirectory(filename);
        const auto constSuffix = "-const-" + std::to_string(constCounter.value()) + ".vtkhdf";
        const auto filenameConst = prefix + "-" + self.name_ + constSuffix;
        const auto filenameConstFile = inDirectory(filenameConst);
        const auto filenamePvu = prefix + "-" + self.name_ + ".pvd";
        if (!self.temporal_) {
          pvu.emplace_back(metadata::PvuEntry{filenameFile, time});
        }
        auto writer = writer::Writer();

        if (fullWrite) {
          // with a const file the unchanging data lives next to the snapshots; in a time series it
          // goes into the one file, once
          const auto& constTarget = self.constFile_ ? filenameConst : filename;
          for (const auto& instruction : self.instructionsConst_) {
            writer.addInstruction(instruction(constTarget, counter, time));
          }
        }
        for (const auto& instruction : self.instructionsConstLink_) {
          writer.addInstruction(instruction(filename, filenameConstFile));
        }
        for (const auto& instruction : self.instructions_) {
          writer.addInstruction(instruction(filename, counter, time));
        }
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
            "Time",
            self.temporal_ ? writer::WriteInline::createShaped<double>(
                                 {writer::Dimension::appended(1)}, {time})
                           : writer::WriteInline::createArray<double>({1}, {time}),
            datatype::inferDatatype<decltype(time)>()));
        writer.addInstruction(std::make_shared<writer::instructions::Hdf5DataWrite>(
            writer::instructions::Hdf5Location(filename, {GroupName, FieldDataName}),
            "Index",
            self.temporal_ ? writer::WriteInline::createShaped<std::size_t>(
                                 {writer::Dimension::appended(1)}, {counter})
                           : writer::WriteInline::createArray<std::size_t>({1}, {counter}),
            datatype::inferDatatype<decltype(counter)>()));
        // A time series carries its own step list, and every step would enter the collection under
        // the same file name, so there is nothing for a pvd to say.
        if (!self.temporal_) {
          writer.addInstructions(metadata::makePvu(pvu).instructions(filenamePvu));
        }
        return writer;
      };
}
} // namespace seissol::io::instance::mesh
