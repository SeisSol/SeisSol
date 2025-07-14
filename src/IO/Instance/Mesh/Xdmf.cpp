// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Xdmf.h"
#include <IO/Datatype/Datatype.h>
#include <IO/Datatype/Inference.h>
#include <IO/Datatype/MPIType.h>
#include <IO/Instance/Geometry/Typedefs.h>
#include <IO/Instance/Metadata/Xml.h>
#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <IO/Writer/Writer.h>
#include <Parallel/MPI.h>
#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <utils/logger.h>
#include <vector>

namespace {
using namespace seissol::io;
using namespace seissol::io::instance::metadata;

struct XdmfDataset {
  std::string name;
  std::string type;

  std::shared_ptr<datatype::Datatype> datatype;
  std::string format;
  std::vector<std::size_t> offset;
  std::vector<std::size_t> dimensions;
  std::string location;

  XdmfDataset(const instance::mesh::XdmfWriter::WriteResult& result)
      : name(result.name), type(result.type), datatype(result.datatype), format(result.format),
        offset(result.offset), dimensions(result.dimensions), location(result.location) {}

  [[nodiscard]] std::shared_ptr<XmlEntry> dataItem(const std::string& numberType,
                                                   int precision,
                                                   const std::string& format,
                                                   const std::vector<std::size_t>& dimensions,
                                                   const std::string& payload) const {
    std::ostringstream dimensionStream;
    dimensionStream << dimensions[0];
    for (std::size_t i = 1; i < dimensions.size(); ++i) {
      dimensionStream << " " << dimensions[i];
    }

    auto data = std::make_shared<XmlData>("DataItem");
    data->addAttribute(XmlAttribute::create("NumberType", numberType))
        .addAttribute(XmlAttribute::create("Precision", precision))
        .addAttribute(XmlAttribute::create("Format", format))
        .addAttribute(XmlAttribute::create("Dimensions", dimensionStream.str()));
    data->setImmediate(writer::WriteInline::createString(payload));
    return data;
  }

  [[nodiscard]] std::shared_ptr<XmlEntry>
      dataItem(const std::shared_ptr<datatype::Datatype>& datatype,
               const std::string& format,
               const std::vector<std::size_t>& dimensions,
               const std::string& payload) const {
    if (dynamic_cast<datatype::IntegerDatatype*>(datatype.get()) != nullptr) {
      return dataItem("Int", datatype->size(), format, dimensions, payload);
    } else if (dynamic_cast<datatype::F32Datatype*>(datatype.get()) != nullptr) {
      return dataItem("Float", datatype->size(), format, dimensions, payload);
    } else if (dynamic_cast<datatype::F64Datatype*>(datatype.get()) != nullptr) {
      return dataItem("Float", datatype->size(), format, dimensions, payload);
    } else {
      logError() << "Internal error while writing an Xdmf file.";
      throw;
    }
  }

  [[nodiscard]] std::shared_ptr<XmlNode>
      hyperslab(const std::shared_ptr<datatype::Datatype>& datatype,
                const std::string& format,
                const std::vector<std::size_t>& offset,
                const std::vector<std::size_t>& dimensions,
                const std::string& payload) const {
    auto hyperslab = std::make_shared<XmlNode>("DataItem");
    std::ostringstream dimensionStream;
    // offset
    dimensionStream << offset[0];
    for (std::size_t i = 1; i < offset.size(); ++i) {
      dimensionStream << " " << offset[i];
    }
    // stride (always uniform)
    for (std::size_t _ = 0; _ < dimensions.size(); ++_) {
      dimensionStream << " 1";
    }
    // dimensions
    for (std::size_t i = 0; i < dimensions.size(); ++i) {
      dimensionStream << " " << dimensions[i];
    }
    hyperslab->addNode(dataItem("Int", 8, "XML", {3, dimensions.size()}, dimensionStream.str()));
    hyperslab->addNode(dataItem(datatype, format, dimensions, payload));
    return hyperslab;
  }

  [[nodiscard]] std::shared_ptr<XmlNode> makeNode() const {
    std::shared_ptr<XmlNode> node(nullptr);

    if (type == "Topology") {
      node = std::make_shared<XmlNode>("Topology");
      node->addAttribute(XmlAttribute::create("TopologyType", name))
          .addAttribute(XmlAttribute::create("NumberOfElements", dimensions[0]));
    } else if (type == "Geometry") {
      node = std::make_shared<XmlNode>("Geometry");
      node->addAttribute(XmlAttribute::create("GeometryType", name))
          .addAttribute(XmlAttribute::create("NumberOfElements", dimensions[0]));
    } else if (type == "AttributeCell") {
      node = std::make_shared<XmlNode>("Attribute");
      node->addAttribute(XmlAttribute::create("Name", name))
          .addAttribute(XmlAttribute::create("Center", "Cell"));
    } else if (type == "AttributeNode") {
      node = std::make_shared<XmlNode>("Attribute");
      node->addAttribute(XmlAttribute::create("Name", name))
          .addAttribute(XmlAttribute::create("Center", "Node"));
    }

    if (offset.empty()) {
      node->addNode(dataItem(datatype, format, dimensions, location));
    } else {
      node->addNode(hyperslab(datatype, format, offset, dimensions, location));
    }

    return node;
  }
};
struct XdmfGrid {
  std::vector<XdmfDataset> datasets;
  std::string name;
};

struct XdmfMeta {
  XmlFile getxml() {

    auto root = std::make_shared<XmlNode>("Xdmf");
    root->addAttribute(XmlAttribute::create("Version", "2.0"));

    auto root2 = std::make_shared<XmlNode>("Domain");
    root->addNode(root2);

    auto root3 = std::make_shared<XmlNode>("Grid");
    root3->addAttribute(XmlAttribute::create("Name", "TimeSeries"))
        .addAttribute(XmlAttribute::create("GridType", "Collection"))
        .addAttribute(XmlAttribute::create("CollectionType", "Temporal"));

    for (const auto& entry : entries) {
      auto grid = std::make_shared<XmlNode>("Grid");
      grid->addAttribute(XmlAttribute::create("Name", entry.name))
          .addAttribute(XmlAttribute::create("GridType", "Uniform"));
      for (const auto& datasetEntry : entriesConst) {
        grid->addNode(datasetEntry.makeNode());
      }
      for (const auto& datasetEntry : entry.datasets) {
        grid->addNode(datasetEntry.makeNode());
      }
      root3->addNode(grid);
    }

    XmlFile file;
    file.setRoot(root);

    return file;
  }

  std::vector<XdmfGrid> entries;
  std::vector<XdmfDataset> entriesConst;
};

} // namespace

namespace seissol::io::instance::mesh {

XdmfWriter::XdmfWriter(const std::string& name,
                       std::size_t localElementCount,
                       geometry::Shape shape,
                       std::size_t targetDegree)
    : localElementCount(localElementCount), globalElementCount(localElementCount), name(name),
      pointsPerElement(geometry::numPoints(targetDegree, shape)),
      type(geometry::xdmfType(shape, targetDegree)) {
  MPI_Exscan(&localElementCount,
             &elementOffset,
             1,
             datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
             MPI_SUM,
             seissol::MPI::mpi.comm());
  MPI_Allreduce(&localElementCount,
                &globalElementCount,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::MPI::mpi.comm());
  pointOffset = elementOffset * pointsPerElement;
  localPointCount = localElementCount * pointsPerElement;
  globalPointCount = globalElementCount * pointsPerElement;

  const auto selfPointOffset = pointOffset;

  addData(
      type,
      "Topology",
      true,
      localPointCount,
      writer::GeneratedBuffer::createElementwise<int64_t>(
          localPointCount, 1, std::vector<std::size_t>(), [=](int64_t* target, std::size_t index) {
            target[0] = index + selfPointOffset;
          }));
}

void XdmfWriter::addData(const std::string& name,
                         const std::string& type,
                         bool isConst,
                         std::size_t localCount,
                         const std::shared_ptr<writer::DataSource>& data) {
  auto& instrarray = isConst ? instructionsConst : instructions;

  const auto datasetId = datasetCount;
  ++datasetCount;

  auto binary{this->binary};

  std::size_t globalCount = localCount;

  MPI_Allreduce(&localCount,
                &globalCount,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::MPI::mpi.comm());

  instrarray.emplace_back([=](const std::string& filename, std::size_t counter) {
    WriteResult result{};
    result.name = name;
    result.type = type;
    result.datatype = data->datatype();
    result.format = binary ? "Binary" : "HDF";

    // two modes:
    // * binary: append to the same file at some position
    // * HDF5: add a new dataset

    if (binary) {
      const auto trueFilename = filename + "-dataset" + std::to_string(datasetId);
      result.format = "Binary";
      result.location = trueFilename;
      result.instruction = std::make_shared<writer::instructions::BinaryWrite>(trueFilename, data);
      if (!isConst) {
        // write into a single file
        result.offset.emplace_back(counter);

        result.offset.emplace_back(0);
        for (std::size_t _ = 0; _ < data->shape().size(); ++_) {
          result.offset.emplace_back(0);
        }

        result.dimensions.emplace_back(counter + 1);
      }
    } else {
      const std::string groupName = "dataset" + std::to_string(datasetId);
      const std::string datasetName = "data" + std::to_string(counter);
      result.format = "HDF";
      result.location = filename + ":/" + groupName + "/" + datasetName;
      result.instruction = std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(filename, {groupName}),
          datasetName,
          data,
          data->datatype());
    }

    result.dimensions.emplace_back(globalCount);
    result.dimensions.insert(result.dimensions.end(), data->shape().begin(), data->shape().end());

    return result;
  });
}

void XdmfWriter::addHook(const std::function<void(std::size_t, double)>& hook) {
  hooks.push_back(hook);
}

std::function<writer::Writer(const std::string&, std::size_t, double)> XdmfWriter::makeWriter() {
  logInfo() << "Adding Xdmf writer" << name;
  const auto self = *this;
  return [self, meta = XdmfMeta()](const std::string& prefix,
                                   std::size_t counter,
                                   double time) mutable -> writer::Writer {
    for (const auto& hook : self.hooks) {
      hook(counter, time);
    }

    auto writer = writer::Writer();

    const auto filenameMeta = prefix + "-" + self.name + ".xdmf";
    const auto foldernameData = prefix + "-" + self.name + "-data";

    auto grid = XdmfGrid{};
    grid.name = "step-" + std::to_string(counter);
    if (counter == 0) {
      for (const auto& instruction : self.instructionsConst) {
        const auto invoked = instruction(foldernameData, counter);
        writer.addInstruction(invoked.instruction);
        meta.entriesConst.emplace_back(invoked);
      }
    }
    for (const auto& instruction : self.instructions) {
      const auto invoked = instruction(foldernameData, counter);
      writer.addInstruction(invoked.instruction);
      grid.datasets.emplace_back(invoked);
    }
    meta.entries.emplace_back(grid);

    writer.addInstructions(meta.getxml().instructions(filenameMeta));
    return writer;
  };
}

} // namespace seissol::io::instance::mesh
