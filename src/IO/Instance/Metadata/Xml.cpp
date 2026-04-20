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

void XmlInstructor::newLine() {
  cache_ << '\n';
  for (int i = 0; i < indent_; ++i) {
    cache_ << ' ';
  }
}

void XmlInstructor::indentLeft() { indent_ -= 1; }

void XmlInstructor::indentRight() { indent_ += 1; }

void XmlInstructor::addBuffer(const std::shared_ptr<writer::DataSource>& dataSource) {
  if (dataSource->distributed()) {
    flush();
    instructionList_.push_back(
        std::make_shared<writer::instructions::BinaryWrite>(file_, dataSource, 0, written_));
    written_ = true;
  } else {
    cache_ << dataSource->datatype()->toStringRaw(dataSource->getLocalPointer());
  }
}

void XmlInstructor::flush() {
  const auto data = cache_.str();
  if (!data.empty()) {
    instructionList_.push_back(std::make_shared<writer::instructions::BinaryWrite>(
        file_, writer::WriteInline::createString(data), 0, written_));
    written_ = true;
  }
  cache_.clear();
}

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> XmlInstructor::instructions() {
  flush();
  return instructionList_;
}

XmlAttribute::XmlAttribute(const std::string& name) : name_(name) {}

void XmlAttribute::write(XmlInstructor& instructor) const {
  instructor.addText(" ");
  instructor.addText(name_);
  instructor.addText("=\"");
  instructor.addBuffer(data_);
  instructor.addText("\"");
}

XmlEntry::XmlEntry(const std::string& name) : name_(name) {}

XmlEntry& XmlEntry::addAttribute(const XmlAttribute& attribute) {
  attributes_.emplace_back(attribute);
  return *this;
}

void XmlEntry::write(XmlInstructor& instructor) const {
  instructor.addText("<");
  instructor.addText(name_);
  for (const auto& attribute : attributes_) {
    attribute.write(instructor);
  }
  if (isEmpty()) {
    instructor.addText(" />");
  } else {
    instructor.addText(">");
    innerWrite(instructor);
    instructor.addText("</");
    instructor.addText(name_);
    instructor.addText(">");
  }
}

XmlNode::XmlNode(const std::string& name) : XmlEntry(name) {}

void XmlNode::addNode(const std::shared_ptr<XmlEntry>& entry) { entries_.push_back(entry); }

void XmlNode::innerWrite(XmlInstructor& instructor) const {
  for (const auto& entry : entries_) {
    instructor.indentRight();
    instructor.newLine();
    entry->write(instructor);
    instructor.indentLeft();
  }
  instructor.newLine();
}

bool XmlNode::isEmpty() const { return entries_.empty(); }

XmlData::XmlData(const std::string& name) : XmlEntry(name) {}

void XmlData::innerWrite(XmlInstructor& instructor) const { instructor.addBuffer(data_); }

bool XmlData::isEmpty() const { return data_ == nullptr; }

void XmlData::setDataSource(const std::shared_ptr<writer::DataSource>& dataSource) {
  data_ = dataSource;
}

const std::shared_ptr<writer::DataSource>& XmlData::getDataSource() const { return data_; }

std::shared_ptr<XmlEntry> XmlFile::getRoot() { return root_; }

void XmlFile::setRoot(std::shared_ptr<XmlEntry> entry) { root_ = std::move(entry); }

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
    XmlFile::instructions(const std::string& file) const {
  XmlInstructor instructor(file);
  instructor.addText("<?xml version=\"1.0\"?>");
  instructor.newLine();
  root_->write(instructor);
  instructor.newLine();
  return instructor.instructions();
}

} // namespace seissol::io::instance::metadata
