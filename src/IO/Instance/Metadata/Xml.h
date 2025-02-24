// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_METADATA_XML_H_
#define SEISSOL_SRC_IO_INSTANCE_METADATA_XML_H_

#include <IO/Writer/Instructions/Binary.h>
#include <IO/Writer/Instructions/Data.h>
#include <IO/Writer/Instructions/Instruction.h>
#include <memory>
#include <sstream>
#include <string>
namespace seissol::io::instance::metadata {

class XmlInstructor {
  public:
  XmlInstructor(const std::string& file);

  void addText(const std::string& text);

  void addBuffer(const std::shared_ptr<writer::DataSource>& dataSource);

  void flush();

  std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> instructions();

  private:
  std::string file;
  std::ostringstream cache;
  std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> instructionList;
};

class XmlAttribute {
  public:
  XmlAttribute(const std::string& name);

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

  template <typename T>
  T getImmediate() const {
    const auto* data = this->data->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const T*>(data);
    return *dataConv;
  }

  void write(XmlInstructor& instructor) const;

  private:
  std::string name;
  std::shared_ptr<writer::DataSource> data;
};

template <>
inline void XmlAttribute::setImmediate<std::string>(const std::string& data) {
  this->data = writer::WriteInline::createString(data);
}

template <>
inline std::string XmlAttribute::getImmediate<std::string>() const {
  const auto* data = this->data->getLocalPointer();
  const auto* dataConv = reinterpret_cast<const char*>(data);
  return std::string(dataConv, dataConv + this->data->getLocalSize());
}

class XmlEntry {
  public:
  virtual ~XmlEntry() = default;
  XmlEntry(const std::string& name);

  XmlEntry& addAttribute(const XmlAttribute& attribute);

  void write(XmlInstructor& instructor) const;

  protected:
  virtual void innerWrite(XmlInstructor& instructor) const = 0;

  private:
  std::string name;
  std::vector<XmlAttribute> attributes;
};

class XmlNode : public XmlEntry {
  public:
  ~XmlNode() override = default;
  XmlNode(const std::string& name);

  void addNode(const std::shared_ptr<XmlEntry>& entry);

  protected:
  void innerWrite(XmlInstructor& instructor) const override;

  private:
  std::vector<std::shared_ptr<XmlEntry>> entries;
};

class XmlData : public XmlEntry {
  public:
  ~XmlData() override = default;
  XmlData(const std::string& name);

  template <typename T>
  void setImmediate(const T& data) {
    this->data = writer::WriteInline::create(data);
  }

  template <typename T>
  T getImmediate() const {
    const auto* data = this->data->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const T*>(data);
    return *dataConv;
  }

  template <typename T>
  void setBuffer(const T* /*unused*/) {}

  template <typename T>
  const T* getBuffer() const;

  protected:
  void innerWrite(XmlInstructor& instructor) const override;

  private:
  std::shared_ptr<writer::DataSource> data;
};

template <>
inline void XmlData::setImmediate<std::string>(const std::string& data) {
  this->data = writer::WriteInline::createString(data);
}

template <>
inline std::string XmlData::getImmediate<std::string>() const {
  const auto* data = this->data->getLocalPointer();
  const auto* dataConv = reinterpret_cast<const char*>(data);
  return std::string(dataConv, dataConv + this->data->getLocalSize());
}

class XmlFile {
  private:
  std::shared_ptr<XmlEntry> root;

  public:
  std::shared_ptr<XmlEntry> getRoot();

  void setRoot(std::shared_ptr<XmlEntry> entry);

  [[nodiscard]] std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
      instructions(const std::string& file) const;
};

} // namespace seissol::io::instance::metadata

#endif // SEISSOL_SRC_IO_INSTANCE_METADATA_XML_H_
