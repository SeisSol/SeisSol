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

XmlInstructor::XmlInstructor(const std::string& file) : file_(file) {}

void XmlInstructor::addText(const std::string& text) { cache_ << text; }

void XmlInstructor::addBuffer(const std::shared_ptr<writer::DataSource>& dataSource) {
  if (dataSource->distributed()) {
    flush();
    instructionList_.push_back(
        std::make_shared<writer::instructions::BinaryWrite>(file_, dataSource));
  } else {
    dataSource->datatype()->toStringRaw(dataSource->getLocalPointer());
  }
}

void XmlInstructor::flush() {
  const auto data = cache_.str();
  if (!data.empty()) {
    instructionList_.push_back(std::make_shared<writer::instructions::BinaryWrite>(
        file_, writer::WriteInline::createString(data)));
  }
  cache_.clear();
}

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> XmlInstructor::instructions() {
  flush();
  return instructionList_;
}

XmlAttribute::XmlAttribute(const std::string& name) : name_(name) {}

void XmlAttribute::write(XmlInstructor& instructor) const {
  instructor.addText(name_);
  instructor.addText("=\"");
  instructor.addBuffer(data_);
  instructor.addText("\" ");
}

XmlEntry::XmlEntry(const std::string& name) : name_(name) {}

XmlEntry& XmlEntry::addAttribute(const XmlAttribute& attribute) {
  attributes_.emplace_back(attribute);
  return *this;
}

void XmlEntry::write(XmlInstructor& instructor) const {
  instructor.addText("<");
  instructor.addText(name_);
  instructor.addText(" ");
  for (const auto& attribute : attributes_) {
    attribute.write(instructor);
  }
  instructor.addText(">");
  innerWrite(instructor);
  instructor.addText("</");
  instructor.addText(name_);
  instructor.addText(">");
}

XmlNode::XmlNode(const std::string& name) : XmlEntry(name) {}

void XmlNode::addNode(const std::shared_ptr<XmlEntry>& entry) { entries_.push_back(entry); }

void XmlNode::innerWrite(XmlInstructor& instructor) const {
  for (const auto& entry : entries_) {
    entry->write(instructor);
  }
}

XmlData::XmlData(const std::string& name) : XmlEntry(name) {}

void XmlData::innerWrite(XmlInstructor& instructor) const { instructor.addBuffer(data_); }

std::shared_ptr<XmlEntry> XmlFile::getRoot() { return root_; }

void XmlFile::setRoot(std::shared_ptr<XmlEntry> entry) { root_ = std::move(entry); }

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
    XmlFile::instructions(const std::string& file) const {
  XmlInstructor instructor(file);
  instructor.addText("<?xml version=\"1.0\"?>");
  root_->write(instructor);
  return instructor.instructions();
}

} // namespace seissol::io::instance::metadata
