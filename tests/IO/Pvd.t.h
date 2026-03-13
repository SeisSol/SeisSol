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
    auto* bw = dynamic_cast<instructions::BinaryWrite*>(instr.get());
    if (bw && bw->dataSource) {
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
// XmlAttribute
// ---------------------------------------------------------------------------

TEST_CASE("XmlAttribute create and getImmediate" * doctest::test_suite("io")) {
  SUBCASE("String attribute") {
    auto attr = XmlAttribute::create("name", "test");
    CHECK(attr.getImmediate<std::string>() == std::string("test"));
  }

  // TODO: those tests were different before; maybe make different again (if we need the Xml
  // interface more for that)

  SUBCASE("Int attribute") {
    auto attr = XmlAttribute::create("count", "42");
    CHECK(attr.getImmediate<std::string>() == std::string("42"));
  }

  SUBCASE("Double attribute") {
    auto attr = XmlAttribute::create("timestep", "1.5");
    CHECK(attr.getImmediate<std::string>() == std::string("1.5"));
  }
}

// ---------------------------------------------------------------------------
// XmlFile: simple node with attributes
// ---------------------------------------------------------------------------

TEST_CASE("XmlFile simple node produces XML" * doctest::test_suite("io")) {
  auto root = std::make_shared<XmlNode>("Root");
  root->addAttribute(XmlAttribute::create("version", std::string("1.0")));

  XmlFile file;
  file.setRoot(root);

  auto instrs = file.instructions("out.xml");
  // Should produce at least one instruction
  CHECK_FALSE(instrs.empty());

  std::string xml = extractInlineText(instrs);
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
  auto root = std::make_shared<XmlNode>("Parent");
  auto child = std::make_shared<XmlNode>("Child");
  child->addAttribute(XmlAttribute::create("id", "1"));
  root->addNode(child);

  XmlFile file;
  file.setRoot(root);

  std::string xml = extractInlineText(file.instructions("test.xml"));
  CHECK(xml.find("<Parent>") != std::string::npos);
  CHECK(xml.find("<Child ") != std::string::npos);
  CHECK(xml.find("</Parent>") != std::string::npos);
}

// ---------------------------------------------------------------------------
// makePvu: PVD file generation
// ---------------------------------------------------------------------------

TEST_CASE("makePvu empty entries" * doctest::test_suite("io")) {
  auto file = makePvu({});
  auto instrs = file.instructions("test.pvd");
  CHECK_FALSE(instrs.empty());

  std::string xml = extractInlineText(instrs);
  CHECK(xml.find("<?xml version=\"1.0\"?>") != std::string::npos);
  CHECK(xml.find("<VTKFile ") != std::string::npos);
  CHECK(xml.find("type=") != std::string::npos);
  CHECK(xml.find("<Collection />") != std::string::npos);
  CHECK(xml.find("</VTKFile>") != std::string::npos);
  // No DataSet elements
  CHECK(xml.find("<DataSet ") == std::string::npos);
}

TEST_CASE("makePvu with entries" * doctest::test_suite("io")) {
  std::vector<PvuEntry> entries = {
      {"output-0001.vtu", 0.0},
      {"output-0002.vtu", 0.5},
      {"output-0003.vtu", 1.0},
  };
  auto file = makePvu(entries);
  auto instrs = file.instructions("test.pvd");
  std::string xml = extractInlineText(instrs);

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
  auto file = makePvu({{"step.vtu", 42.0}});
  auto instrs = file.instructions("out.pvd");
  std::string xml = extractInlineText(instrs);

  CHECK(xml.find("<VTKFile ") != std::string::npos);
  CHECK(xml.find("<Collection>") != std::string::npos);

  CHECK(xml.find("<DataSet ") != std::string::npos);

  CHECK(xml.find("</Collection>") != std::string::npos);
  CHECK(xml.find("</VTKFile>") != std::string::npos);

  // All instructions reference the same file
  for (const auto& instr : instrs) {
    auto* bw = dynamic_cast<instructions::BinaryWrite*>(instr.get());
    if (bw) {
      CHECK(bw->filename == "out.pvd");
    }
  }
}

} // namespace seissol::unit_test
