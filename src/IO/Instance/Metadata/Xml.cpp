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
#include <utils/logger.h>
#include <vector>
namespace seissol::io::instance::metadata {

XmlInstructor::XmlInstructor(std::string file) : file_(std::move(file)) {}

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
  // (clear() alone would only reset the state of the stream, and keep the text)
  cache_.str(std::string());
  cache_.clear();
}

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> XmlInstructor::instructions() {
  flush();
  return instructionList_;
}

std::string XmlInstructor::text() const {
  if (!instructionList_.empty()) {
    logError() << "The XML file" << file_
               << "holds distributed data, which is not written as text.";
  }
  return cache_.str();
}

XmlAttribute::XmlAttribute(std::string name) : name_(std::move(name)) {}

void XmlAttribute::write(XmlInstructor& instructor) const {
  instructor.addText(" ");
  instructor.addText(name_);
  instructor.addText("=\"");
  instructor.addBuffer(data_);
  instructor.addText("\"");
}

XmlEntry::XmlEntry(std::string name) : name_(std::move(name)) {}

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

void XmlFile::write(XmlInstructor& instructor) const {
  instructor.addText("<?xml version=\"1.0\"?>");
  instructor.newLine();
  root_->write(instructor);
  instructor.newLine();
}

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
    XmlFile::instructions(const std::string& file) const {
  XmlInstructor instructor(file);
  write(instructor);
  return instructor.instructions();
}

std::string XmlFile::text() const {
  XmlInstructor instructor("");
  write(instructor);
  return instructor.text();
}

GrowingXmlFile::GrowingXmlFile(const Document& document) {
  // the file with a placeholder in the place of the entries, cut there
  const std::string placeholderName = "seissol-entries";
  const std::string placeholder = "<" + placeholderName + " />";
  const auto text = document({std::make_shared<XmlData>(placeholderName)}).text();
  const auto at = text.find(placeholder);
  // an entry starts on a line of its own, indented
  const auto lineStart = at == std::string::npos ? at : text.rfind('\n', at);
  if (lineStart == std::string::npos) {
    logError() << "The entries a growing XML file takes have to be children of one of its nodes.";
  }
  head_ = text.substr(0, lineStart);
  tail_ = text.substr(at + placeholder.size());
  indent_ = static_cast<int>(at - lineStart - 1);
}

std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
    GrowingXmlFile::add(const std::string& file,
                        const std::vector<std::shared_ptr<XmlEntry>>& entries) {
  // the same as the parent node writes for an entry
  XmlInstructor instructor(file);
  for (int i = 0; i < indent_; ++i) {
    instructor.indentRight();
  }
  for (const auto& entry : entries) {
    instructor.newLine();
    entry->write(instructor);
  }
  const auto added = instructor.text();

  if (!tailOffset_.has_value()) {
    tailOffset_ = head_.size() + added.size();
    return {std::make_shared<writer::instructions::BinaryWrite>(
        file, writer::WriteInline::createString(head_ + added + tail_), 0, false)};
  }
  const auto offset = tailOffset_.value();
  tailOffset_ = offset + added.size();
  return {std::make_shared<writer::instructions::BinaryWrite>(
      file, writer::WriteInline::createString(added + tail_), 0, true, offset)};
}

} // namespace seissol::io::instance::metadata
