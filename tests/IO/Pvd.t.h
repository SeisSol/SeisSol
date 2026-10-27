// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include <doctest.h>

#include "IO/Instance/Metadata/Pvd.h"
#include "IO/Instance/Metadata/Xml.h"
#include "IO/Writer/Instructions/Binary.h"

#include <memory>
#include <string>
#include <vector>

namespace seissol::unit_test {
using namespace seissol::io::instance::metadata;
using namespace seissol::io::writer;

// Helper: extract concatenated inline text from write instructions.
// BinaryWrite instructions with non-distributed DataSources carry inline data.
static std::string
    extractInlineText(const std::vector<std::shared_ptr<instructions::WriteInstruction>>& instrs) {
  std::string result;
  for (const auto& instr : instrs) {
    const auto* bw = dynamic_cast<instructions::BinaryWrite*>(instr.get());
    if ((bw != nullptr) && bw->dataSource) {
      const auto* ptr = reinterpret_cast<const char*>(bw->dataSource->getLocalPointer());
      const auto size = bw->dataSource->getLocalSize();
      if (ptr != nullptr && size > 0) {
        result.append(ptr, size);
      }
    }
  }
  return result;
}

// ---------------------------------------------------------------------------
// Times in the metadata
// ---------------------------------------------------------------------------

TEST_CASE("formatExact keeps a time as it is" * doctest::test_suite("io")) {
  CHECK(formatExact(0.1) == "0.1");
  CHECK(formatExact(2.0) == "2");
  // an output interval far below a microsecond must not collapse to zero
  for (const double value : {1e-7, 2.5e-9, 1.0 / 3.0, 12345.678901234567}) {
    CHECK(std::stod(formatExact(value)) == value);
  }
}

// ---------------------------------------------------------------------------
// XmlAttribute
// ---------------------------------------------------------------------------

TEST_CASE("XmlAttribute create and getImmediate" * doctest::test_suite("io")) {
  SUBCASE("String attribute") {
    const auto attr = XmlAttribute::create("name", "test");
    CHECK(attr.getImmediate<std::string>() == std::string("test"));
  }

  // TODO: those tests were different before; maybe make different again (if we need the Xml
  // interface more for that)

  SUBCASE("Int attribute") {
    const auto attr = XmlAttribute::create("count", "42");
    CHECK(attr.getImmediate<std::string>() == std::string("42"));
  }

  SUBCASE("Double attribute") {
    const auto attr = XmlAttribute::create("timestep", "1.5");
    CHECK(attr.getImmediate<std::string>() == std::string("1.5"));
  }
}

// ---------------------------------------------------------------------------
// XmlFile: simple node with attributes
// ---------------------------------------------------------------------------

TEST_CASE("XmlFile simple node produces XML" * doctest::test_suite("io")) {
  const auto root = std::make_shared<XmlNode>("Root");
  root->addAttribute(XmlAttribute::create("version", std::string("1.0")));

  XmlFile file;
  file.setRoot(root);

  const auto instrs = file.instructions("out.xml");
  // Should produce at least one instruction
  CHECK_FALSE(instrs.empty());

  const std::string xml = extractInlineText(instrs);
  // Must contain the XML header
  CHECK(xml.find("<?xml version=\"1.0\"?>") != std::string::npos);
  // Must contain the root element
  CHECK(xml.find("<Root ") != std::string::npos);
  // Must contain the attribute name
  CHECK(xml.find("version=") != std::string::npos);
}

// ---------------------------------------------------------------------------
// XmlFile: nested nodes
// ---------------------------------------------------------------------------

TEST_CASE("XmlFile nested nodes" * doctest::test_suite("io")) {
  const auto root = std::make_shared<XmlNode>("Parent");
  const auto child = std::make_shared<XmlNode>("Child");
  child->addAttribute(XmlAttribute::create("id", "1"));
  root->addNode(child);

  XmlFile file;
  file.setRoot(root);

  const std::string xml = extractInlineText(file.instructions("test.xml"));
  CHECK(xml.find("<Parent>") != std::string::npos);
  CHECK(xml.find("<Child ") != std::string::npos);
  CHECK(xml.find("</Parent>") != std::string::npos);
}

// ---------------------------------------------------------------------------
// makePvu: PVD file generation
// ---------------------------------------------------------------------------

TEST_CASE("makePvu empty entries" * doctest::test_suite("io")) {
  const auto file = makePvu({});
  const auto instrs = file.instructions("test.pvd");
  CHECK_FALSE(instrs.empty());

  const std::string xml = extractInlineText(instrs);
  CHECK(xml.find("<?xml version=\"1.0\"?>") != std::string::npos);
  CHECK(xml.find("<VTKFile ") != std::string::npos);
  CHECK(xml.find("type=") != std::string::npos);
  CHECK(xml.find("<Collection />") != std::string::npos);
  CHECK(xml.find("</VTKFile>") != std::string::npos);
  // No DataSet elements
  CHECK(xml.find("<DataSet ") == std::string::npos);
}

TEST_CASE("makePvu with entries" * doctest::test_suite("io")) {
  const std::vector<PvuEntry> entries = {
      {"output-0001.vtu", 0.0},
      {"output-0002.vtu", 0.5},
      {"output-0003.vtu", 1.0},
  };
  const auto file = makePvu(entries);
  const auto instrs = file.instructions("test.pvd");
  const std::string xml = extractInlineText(instrs);

  CHECK(xml.find("<VTKFile ") != std::string::npos);
  CHECK(xml.find("<Collection>") != std::string::npos);

  // Each entry should produce a DataSet element with file and timestep attributes
  // The file names appear as attribute values via inline data
  CHECK(xml.find("file=") != std::string::npos);
  CHECK(xml.find("timestep=") != std::string::npos);

  // 3 DataSet entries
  std::size_t datasetCount = 0;
  std::size_t pos = 0;
  while ((pos = xml.find("<DataSet ", pos)) != std::string::npos) {
    ++datasetCount;
    ++pos;
  }
  CHECK(datasetCount == 3);
}

TEST_CASE("makePvu single entry" * doctest::test_suite("io")) {
  const auto file = makePvu({{"step.vtu", 42.0}});
  const auto instrs = file.instructions("out.pvd");
  const std::string xml = extractInlineText(instrs);

  CHECK(xml.find("<VTKFile ") != std::string::npos);
  CHECK(xml.find("<Collection>") != std::string::npos);

  CHECK(xml.find("<DataSet ") != std::string::npos);

  CHECK(xml.find("</Collection>") != std::string::npos);
  CHECK(xml.find("</VTKFile>") != std::string::npos);

  // All instructions reference the same file
  for (const auto& instr : instrs) {
    auto* bw = dynamic_cast<instructions::BinaryWrite*>(instr.get());
    if (bw != nullptr) {
      CHECK(bw->filename == "out.pvd");
    }
  }
}

// ---------------------------------------------------------------------------
// GrowingXmlFile: the file written by step is the file written in full
// ---------------------------------------------------------------------------

// Carries out the writes of @p instrs on @p file as BinaryWriter does.
static void
    applyWrites(std::string& file,
                const std::vector<std::shared_ptr<instructions::WriteInstruction>>& instrs) {
  for (const auto& instr : instrs) {
    auto* bw = dynamic_cast<instructions::BinaryWrite*>(instr.get());
    REQUIRE(bw != nullptr);
    REQUIRE_FALSE(bw->dataSource->distributed());
    const auto* ptr = reinterpret_cast<const char*>(bw->dataSource->getLocalPointer());
    const std::string text(ptr, ptr + bw->dataSource->getLocalSize());
    if (!bw->append) {
      file.clear();
    }
    const auto at = bw->offset.value_or(file.size());
    if (file.size() < at + text.size()) {
      file.resize(at + text.size());
    }
    file.replace(at, text.size(), text);
  }
}

TEST_CASE("GrowingXmlFile writes the file a step at a time" * doctest::test_suite("io")) {
  // nested as an Xdmf file, whose time series takes a grid per step
  const auto document = [](const std::vector<std::shared_ptr<XmlEntry>>& steps) {
    const auto root = std::make_shared<XmlNode>("Outer");
    const auto series = std::make_shared<XmlNode>("Series");
    series->addAttribute(XmlAttribute::create("kind", "temporal"));
    root->addNode(series);
    for (const auto& step : steps) {
      series->addNode(step);
    }
    XmlFile file;
    file.setRoot(root);
    return file;
  };
  const auto makeStep = [](int index) {
    auto step = std::make_shared<XmlNode>("Step");
    step->addAttribute(XmlAttribute::create("index", std::to_string(index)));
    const auto value = std::make_shared<XmlData>("Value");
    value->setImmediate<std::string>("data of step " + std::to_string(index));
    step->addNode(value);
    return step;
  };

  GrowingXmlFile growing(document);
  std::vector<std::shared_ptr<XmlEntry>> steps;
  std::string file;
  for (int index = 0; index < 4; ++index) {
    steps.emplace_back(makeStep(index));
    const auto instrs = growing.add("out.xdmf", {steps.back()});
    REQUIRE(instrs.size() == 1);
    const auto* bw = dynamic_cast<instructions::BinaryWrite*>(instrs.front().get());
    REQUIRE(bw != nullptr);
    CHECK(bw->filename == "out.xdmf");
    // only the first write writes the head, and from the start
    CHECK(bw->append == (index > 0));
    CHECK(bw->offset.has_value() == (index > 0));
    CHECK((extractInlineText(instrs).find("<?xml") == std::string::npos) == (index > 0));

    applyWrites(file, instrs);
    CHECK(file == document(steps).text());
  }
}

TEST_CASE("GrowingXmlFile writes a pvd a step at a time" * doctest::test_suite("io")) {
  GrowingXmlFile growing(makePvuCollection);
  std::vector<PvuEntry> entries;
  std::string file;
  for (int index = 0; index < 3; ++index) {
    entries.push_back({"output-" + std::to_string(index) + ".vtkhdf", 0.5 * index});
    applyWrites(file, growing.add("out.pvd", {makePvuDataset(entries.back())}));
    CHECK(file == makePvu(entries).text());
  }
}

} // namespace seissol::unit_test
