// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Xml.h"

#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Instruction.h"

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
namespace seissol::io::instance::metadata {

XmlInstructor::XmlInstructor(const std::string& file) : file(file) {}

void XmlInstructor::addText(const std::string& text) { cache << text; }

void XmlInstructor::newLine() {
  cache << '\n';
  for (int i = 0; i < indent; ++i) {
    cache << ' ';
  }
}

void XmlInstructor::indentLeft() { indent -= 1; }

void XmlInstructor::indentRight() { indent += 1; }

void XmlInstructor::addBuffer(const std::shared_ptr<writer::DataSource>& dataSource) {
  if (dataSource->distributed()) {
    flush();
    instructionList.push_back(
        std::make_shared<writer::instructions::BinaryWrite>(file, dataSource, 0, written));
    written = true;
  } else {
    cache << dataSource->datatype()->toStringRaw(dataSource->getLocalPointer());
  }
}

void XmlInstructor::flush() {
  const auto data = cache.str();
  if (!data.empty()) {
    instructionList.push_back(std::make_shared<writer::instructions::BinaryWrite>(
        file, writer::WriteInline::createString(data), 0, written));
    written = true;
  }
  cache.clear();
}

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> XmlInstructor::instructions() {
  flush();
  return instructionList;
}

XmlAttribute::XmlAttribute(const std::string& name) : name(name) {}

void XmlAttribute::write(XmlInstructor& instructor) const {
  instructor.addText(" ");
  instructor.addText(name);
  instructor.addText("=\"");
  instructor.addBuffer(data);
  instructor.addText("\"");
}

XmlEntry::XmlEntry(const std::string& name) : name(name) {}

XmlEntry& XmlEntry::addAttribute(const XmlAttribute& attribute) {
  attributes.emplace_back(attribute);
  return *this;
}

void XmlEntry::write(XmlInstructor& instructor) const {
  instructor.addText("<");
  instructor.addText(name);
  for (const auto& attribute : attributes) {
    attribute.write(instructor);
  }
  if (isEmpty()) {
    instructor.addText(" />");
  } else {
    instructor.addText(">");
    innerWrite(instructor);
    instructor.addText("</");
    instructor.addText(name);
    instructor.addText(">");
  }
}

XmlNode::XmlNode(const std::string& name) : XmlEntry(name) {}

void XmlNode::addNode(const std::shared_ptr<XmlEntry>& entry) { entries.push_back(entry); }

void XmlNode::innerWrite(XmlInstructor& instructor) const {
  for (const auto& entry : entries) {
    instructor.indentRight();
    instructor.newLine();
    entry->write(instructor);
    instructor.indentLeft();
  }
  instructor.newLine();
}

bool XmlNode::isEmpty() const { return entries.empty(); }

XmlData::XmlData(const std::string& name) : XmlEntry(name) {}

void XmlData::innerWrite(XmlInstructor& instructor) const { instructor.addBuffer(data); }

bool XmlData::isEmpty() const { return data == nullptr; }

void XmlData::setDataSource(const std::shared_ptr<writer::DataSource>& dataSource) {
  data = dataSource;
}

const std::shared_ptr<writer::DataSource>& XmlData::getDataSource() const { return data; }

std::shared_ptr<XmlEntry> XmlFile::getRoot() { return root; }

void XmlFile::setRoot(std::shared_ptr<XmlEntry> entry) { root = std::move(entry); }

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
    XmlFile::instructions(const std::string& file) const {
  XmlInstructor instructor(file);
  instructor.addText("<?xml version=\"1.0\"?>");
  instructor.newLine();
  root->write(instructor);
  instructor.newLine();
  return instructor.instructions();
}

} // namespace seissol::io::instance::metadata
