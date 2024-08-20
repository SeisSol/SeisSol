#pragma once

#include <IO/Writer/Instructions/Binary.hpp>
#include <IO/Writer/Instructions/Data.hpp>
#include <IO/Writer/Instructions/Instruction.hpp>
#include <memory>
#include <sstream>
#include <string>
namespace seissol::io::instance::metadata {

class XmlInstructor {
  public:
  XmlInstructor(const std::string file) : file(file) {}

  void addText(const std::string& text) { cache << text; }

  void addBuffer(std::shared_ptr<writer::DataSource> dataSource) {
    if (dataSource->distributed()) {
      flush();
      instructionList.push_back(
          std::make_shared<writer::instructions::BinaryWrite>(file, dataSource));
    } else {
      dataSource->datatype()->toStringRaw(dataSource->getLocalPointer());
    }
  }

  void flush() {
    const auto data = cache.str();
    if (!data.empty()) {
      instructionList.push_back(std::make_shared<writer::instructions::BinaryWrite>(
          file, writer::WriteInline::createString(data)));
    }
    cache.clear();
  }

  std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> instructions() {
    flush();
    return instructionList;
  }

  private:
  std::string file;
  std::ostringstream cache;
  std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> instructionList;
};

class XmlAttribute {
  public:
  XmlAttribute(const std::string& name) : name(name) {}

  template <typename T>
  static XmlAttribute create(const std::string& name, const T& value) {
    auto attribute = XmlAttribute(name);
    attribute.setImmediate(value);
    return attribute;
  }

  template <typename T>
  void setImmediate(const T& data) {
    this->data = writer::WriteInline::create(data);
  }

  template <>
  void setImmediate<std::string>(const std::string& data) {
    this->data = writer::WriteInline::createString(data);
  }

  template <typename T>
  T getImmediate() const {
    const auto* data = this->data->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const T*>(data);
    return *dataConv;
  }

  template <>
  std::string getImmediate<std::string>() const {
    const auto* data = this->data->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const char*>(data);
    return std::string(dataConv, dataConv + this->data->getLocalSize());
  }

  void write(XmlInstructor& instructor) const {
    instructor.addText(name);
    instructor.addText("=\"");
    instructor.addBuffer(data);
    instructor.addText("\" ");
  }

  private:
  std::string name;
  std::shared_ptr<writer::DataSource> data;
};

class XmlEntry {
  public:
  XmlEntry(const std::string& name) : name(name) {}

  XmlEntry& addAttribute(XmlAttribute&& attribute) {
    attributes.emplace_back(attribute);
    return *this;
  }

  void write(XmlInstructor& instructor) const {
    instructor.addText("<");
    instructor.addText(name);
    instructor.addText(" ");
    for (const auto& attribute : attributes) {
      attribute.write(instructor);
    }
    instructor.addText(">");
    innerWrite(instructor);
    instructor.addText("</");
    instructor.addText(name);
    instructor.addText(">");
  }

  protected:
  virtual void innerWrite(XmlInstructor& instructor) const = 0;

  private:
  std::string name;
  std::vector<XmlAttribute> attributes;
};

class XmlNode : public XmlEntry {
  public:
  XmlNode(const std::string& name) : XmlEntry(name) {}

  void addNode(std::shared_ptr<XmlEntry> entry) { entries.push_back(entry); }

  void addNode(XmlEntry&& entry) {
    entries.emplace_back(std::make_shared<XmlEntry>(std::forward<XmlEntry>(entry)));
  }

  protected:
  void innerWrite(XmlInstructor& instructor) const override {
    for (const auto& entry : entries) {
      entry->write(instructor);
    }
  }

  private:
  std::vector<std::shared_ptr<XmlEntry>> entries;
};

class XmlData : public XmlEntry {
  public:
  XmlData(const std::string& name) : XmlEntry(name) {}

  template <typename T>
  void setImmediate(const T& data) {
    this->data = writer::WriteInline::create(data);
  }

  template <>
  void setImmediate<std::string>(const std::string& data) {
    this->data = writer::WriteInline::createString(data);
  }

  template <typename T>
  T getImmediate() const {
    const auto* data = this->data->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const T*>(data);
    return *dataConv;
  }

  template <>
  std::string getImmediate<std::string>() const {
    const auto* data = this->data->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const char*>(data);
    return std::string(dataConv, dataConv + this->data->getLocalSize());
  }

  template <typename T>
  void setBuffer(const T*) {}

  template <typename T>
  const T* getBuffer() const;

  protected:
  void innerWrite(XmlInstructor& instructor) const override { instructor.addBuffer(data); }

  private:
  std::shared_ptr<writer::DataSource> data;
};

class XmlFile {
  private:
  std::shared_ptr<XmlEntry> root;

  public:
  std::shared_ptr<XmlEntry> getRoot() { return root; }

  void setRoot(std::shared_ptr<XmlEntry> entry) { root = entry; }

  void setRoot(XmlEntry&& entry) {
    root = std::make_shared<XmlEntry>(std::forward<XmlEntry>(entry));
  }

  std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
      instructions(const std::string& file) const {
    XmlInstructor instructor(file);
    instructor.addText("<?xml version=\"1.0\"?>");
    root->write(instructor);
    return instructor.instructions();
  }
};

} // namespace seissol::io::instance::metadata
