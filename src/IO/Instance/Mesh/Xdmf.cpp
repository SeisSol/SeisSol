// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Xdmf.h"

#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"
#include "IO/Datatype/MPIType.h"
#include "IO/Instance/Geometry/Typedefs.h"
#include "IO/Instance/Metadata/Xml.h"
#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Hdf5.h"
#include "IO/Writer/Writer.h"
#include "Parallel/MPI.h"

#include <cstddef>
#include <cstdint>
#include <functional>
#include <memory>
#include <mpi.h>
#include <sstream>
#include <string>
#include <utils/logger.h>
#include <utils/stringutils.h>
#include <vector>

namespace seissol::io::instance::mesh {

namespace {
using namespace seissol::io::instance::metadata;

struct XdmfDataset {
  std::string name;
  std::string type;

  std::shared_ptr<datatype::Datatype> datatype;
  std::string format;
  std::vector<std::size_t> offset;
  std::vector<std::size_t> dimensions;
  std::string location;

  explicit XdmfDataset(const instance::mesh::XdmfWriter::WriteResult& result)
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
        .addAttribute(XmlAttribute::create("Precision", std::to_string(precision)))
        .addAttribute(XmlAttribute::create("Format", format))
        .addAttribute(XmlAttribute::create("Dimensions", dimensionStream.str()));
    data->setDataSource(writer::WriteInline::createString(payload));
    return data;
  }

  [[nodiscard]] std::shared_ptr<XmlEntry>
      dataItem(const std::shared_ptr<datatype::Datatype>& datatype,
               const std::string& format,
               const std::vector<std::size_t>& dimensions,
               const std::string& payload) const {

    std::string itemType = "invalid";
    if (dynamic_cast<datatype::IntegerDatatype*>(datatype.get()) != nullptr) {
      itemType = "Int";
    } else if (dynamic_cast<datatype::F32Datatype*>(datatype.get()) != nullptr) {
      itemType = "Float";
    } else if (dynamic_cast<datatype::F64Datatype*>(datatype.get()) != nullptr) {
      itemType = "Float";
    } else {
      logError() << "Internal error while writing an Xdmf file.";
    }

    return dataItem(itemType, datatype->size(), format, dimensions, payload);
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
    hyperslab->addAttribute(XmlAttribute::create("ItemType", "HyperSlab"))
        .addAttribute(XmlAttribute::create("Dimensions", std::to_string(dimensions[1])));
    hyperslab->addNode(dataItem("UInt", 8, "XML", {3, dimensions.size()}, dimensionStream.str()));
    hyperslab->addNode(dataItem(datatype, format, dimensions, payload));
    return hyperslab;
  }

  [[nodiscard]] std::shared_ptr<XmlNode> makeNode() const {
    std::shared_ptr<XmlNode> node(nullptr);

    if (type == "Topology") {
      node = std::make_shared<XmlNode>("Topology");
      node->addAttribute(XmlAttribute::create("TopologyType", name))
          .addAttribute(XmlAttribute::create("NumberOfElements", std::to_string(dimensions[0])));
    } else if (type == "Geometry") {
      node = std::make_shared<XmlNode>("Geometry");
      node->addAttribute(XmlAttribute::create("GeometryType", name))
          .addAttribute(XmlAttribute::create("NumberOfElements", std::to_string(dimensions[0])));
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
  double time{0};
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
    root2->addNode(root3);

    for (const auto& entry : entries) {
      auto grid = std::make_shared<XmlNode>("Grid");
      grid->addAttribute(XmlAttribute::create("Name", entry.name))
          .addAttribute(XmlAttribute::create("GridType", "Uniform"));

      // add the time node here already
      auto gridTime = std::make_shared<XmlNode>("Time");
      gridTime->addAttribute(XmlAttribute::create("Value", std::to_string(entry.time)));
      grid->addNode(gridTime);

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

XdmfWriter::XdmfWriter(const std::string& name,
                       std::size_t localElementCount,
                       geometry::Shape shape,
                       std::size_t targetDegree,
                       bool binary,
                       int32_t compress)
    : name(name), type(geometry::xdmfType(shape, targetDegree)), binary(binary), compress(compress),
      localElementCount(localElementCount), globalElementCount(localElementCount),
      pointsPerElement(
          geometry::numPoints(std::max(targetDegree, static_cast<std::size_t>(1)), shape)) {
  MPI_Exscan(&localElementCount,
             &elementOffset,
             1,
             datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
             MPI_SUM,
             seissol::Mpi::mpi.comm());
  MPI_Allreduce(&localElementCount,
                &globalElementCount,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::Mpi::mpi.comm());
  pointOffset = elementOffset * pointsPerElement;
  localPointCount = localElementCount * pointsPerElement;
  globalPointCount = globalElementCount * pointsPerElement;

  const auto selfPointOffset = pointOffset;
  const auto selfPointsPerElement = pointsPerElement;

  addData(type,
          "Topology",
          true,
          localElementCount,
          writer::GeneratedBuffer::createElementwise<int64_t>(
              localElementCount, 1, {pointsPerElement}, [=](int64_t* target, std::size_t index) {
                for (std::size_t i = 0; i < selfPointsPerElement; ++i) {
                  target[i] = selfPointsPerElement * index + i + selfPointOffset;
                }
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

  const auto binary{this->binary};
  const auto compress{this->compress};

  std::size_t globalCount = localCount;

  MPI_Allreduce(&localCount,
                &globalCount,
                1,
                datatype::convertToMPI(datatype::inferDatatype<std::size_t>()),
                MPI_SUM,
                seissol::Mpi::mpi.comm());

  instrarray.emplace_back([=](const std::string& preFilename, std::size_t counter) {
    WriteResult result{};
    result.name = name;
    result.type = type;
    result.datatype = data->datatype();
    result.format = binary ? "Binary" : "HDF";

    const auto filenameParts = utils::StringUtils::split(preFilename, '/');
    const auto& filename = filenameParts.back();

    // for both binary and HDF5: append to the same dataset (add a dummy dimension)

    if (!isConst) {
      // write into a single file
      result.offset.emplace_back(counter);

      result.offset.emplace_back(0);
      for (std::size_t _ = 0; _ < data->shape().size(); ++_) {
        result.offset.emplace_back(0);
      }

      result.dimensions.emplace_back(1);
    }

    if (binary) {
      const auto trueFilepath = preFilename + "-dataset" + std::to_string(datasetId);
      const auto trueFilename = filename + "-dataset" + std::to_string(datasetId);
      result.format = "Binary";
      result.location = trueFilename;
      result.instruction =
          std::make_shared<writer::instructions::BinaryWrite>(trueFilepath, data, 0, counter > 0);
    } else {
      const std::string datasetName = "dataset" + std::to_string(datasetId);
      result.format = "HDF";
      result.location = filename + ":/" + datasetName;
      result.instruction = std::make_shared<writer::instructions::Hdf5DataWrite>(
          writer::instructions::Hdf5Location(preFilename, {}),
          datasetName,
          data,
          data->datatype(),
          true,
          compress);
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
    grid.time = time;
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
