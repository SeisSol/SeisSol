// SPDX-FileCopyrightText: 2025 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Xdmf.h"
#include <IO/Datatype/Datatype.h>
#include <IO/Instance/Geometry/Typedefs.h>
#include <IO/Instance/Metadata/Xml.h>
#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Hdf5.h>
#include <IO/Writer/Instructions/Instruction.h>
#include <IO/Writer/Writer.h>
#include <cstddef>
#include <functional>
#include <memory>
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
    dimensionStream << offset[0];
    for (std::size_t i = 1; i < offset.size(); ++i) {
      dimensionStream << " " << offset[i];
    }
    for (std::size_t i = 0; i < dimensions.size(); ++i) {
      dimensionStream << " " << dimensions[i];
    }
    hyperslab->addNode(dataItem("Int", 8, "XML", {dimensions.size(), 2}, dimensionStream.str()));
    hyperslab->addNode(dataItem(datatype, format, dimensions, payload));
    return hyperslab;
  }

  [[nodiscard]] std::shared_ptr<XmlNode> makeNode() const {
    std::shared_ptr<XmlNode> node(nullptr);

    if (type == "Topology") {
      node = std::make_shared<XmlNode>("Topology");
      node->addAttribute(XmlAttribute::create("TopologyType", "Tetrahedron"))
          .addAttribute(XmlAttribute::create("NumberOfElements", dimensions[0]));
    } else if (type == "Geometry") {
      node = std::make_shared<XmlNode>("Geometry");
      node->addAttribute(XmlAttribute::create("GeometryType", "XYZ"))
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
      type(geometry::xdmfType(shape, targetDegree)) {}

void XdmfWriter::addData(const std::string& name,
                         const std::string& type,
                         bool isConst,
                         const std::shared_ptr<writer::DataSource>& data) {
  auto& instrarray = isConst ? instructionsConst : instructions;

  auto binary{this->binary};

  instrarray.emplace_back([=](const std::string& filename, double time) {
    const auto location = binary ? filename : (filename + ":/data");
    const std::string format = binary ? "Binary" : "HDF";
    const auto writeInstruction =
        binary ? std::dynamic_pointer_cast<writer::instructions::WriteInstruction>(
                     std::make_shared<writer::instructions::BinaryWrite>(filename, data))
               : std::dynamic_pointer_cast<writer::instructions::WriteInstruction>(
                     std::make_shared<writer::instructions::Hdf5DataWrite>(
                         writer::instructions::Hdf5Location(filename, {}),
                         "data",
                         data,
                         data->datatype()));
    return WriteResult{
        name,
        type,
        data->datatype(),
        format,
        {},
        data->shape(),
        location,
    };
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
      std::size_t ccount = 0;
      for (const auto& instruction : self.instructionsConst) {
        const auto invoked = instruction(prefix + "c" + std::to_string(ccount), time);
        writer.addInstruction(invoked.instruction);
        meta.entriesConst.emplace_back(invoked);
        ++ccount;
      }
    }
    std::size_t count = 0;
    for (const auto& instruction : self.instructions) {
      const auto invoked = instruction(prefix + std::to_string(count), time);
      writer.addInstruction(invoked.instruction);
      grid.datasets.emplace_back(invoked);
      ++count;
    }
    meta.entries.emplace_back(grid);

    writer.addInstructions(meta.getxml().instructions(filenameMeta));
    return writer;
  };
}

} // namespace seissol::io::instance::mesh
