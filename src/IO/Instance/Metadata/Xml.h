// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_INSTANCE_METADATA_XML_H_
#define SEISSOL_SRC_IO_INSTANCE_METADATA_XML_H_

#include "IO/Writer/Instructions/Binary.h"
#include "IO/Writer/Instructions/Data.h"
#include "IO/Writer/Instructions/Instruction.h"

#include <cstddef>
#include <functional>
#include <iomanip>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <vector>

namespace seissol::io::instance::metadata {

/**
 * @brief A double as text that reads back as the same double, e.g. an output time.
 *
 * Fifteen significant digits where they suffice, which keeps 0.1 as "0.1", and seventeen where
 * they do not.
 */
inline std::string formatExact(double value) {
  std::ostringstream shortForm;
  shortForm << std::setprecision(15) << value;
  if (std::stod(shortForm.str()) == value) {
    return shortForm.str();
  }
  std::ostringstream longForm;
  longForm << std::setprecision(17) << value;
  return longForm.str();
}

class XmlInstructor {
  public:
  explicit XmlInstructor(const std::string& file);

  void addText(const std::string& text);

  void newLine();

  void indentLeft();
  void indentRight();

  void addBuffer(const std::shared_ptr<writer::DataSource>& dataSource);

  void flush();

  std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> instructions();

  /**
   * @brief What was written, as text; only for an XML file holding no distributed data (which
   * gets written by an instruction of its own).
   */
  [[nodiscard]] std::string text() const;

  private:
  bool written_{false};
  int indent_{0};
  std::string file_;
  std::ostringstream cache_;
  std::vector<std::shared_ptr<writer::instructions::WriteInstruction>> instructionList_;
};

class XmlAttribute {
  public:
  explicit XmlAttribute(const std::string& name);

  template <typename T>
  void setImmediate(const T& data) {
    this->data_ = writer::WriteInline::create(data);
  }

  template <typename T>
  T getImmediate() const {
    const auto* data = this->data_->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const T*>(data);
    return *dataConv;
  }

  static XmlAttribute create(const std::string& name, const std::string& value);

  void write(XmlInstructor& instructor) const;

  private:
  std::string name_;
  std::shared_ptr<writer::DataSource> data_;
};

template <>
inline void XmlAttribute::setImmediate<std::string>(const std::string& data) {
  this->data_ = writer::WriteInline::createString(data);
}

template <>
inline std::string XmlAttribute::getImmediate<std::string>() const {
  const auto* data = this->data_->getLocalPointer();
  const auto* dataConv = reinterpret_cast<const char*>(data);

  // createString stores the text without a null terminator
  return std::string(dataConv, dataConv + this->data_->getLocalSize());
}

inline XmlAttribute XmlAttribute::create(const std::string& name, const std::string& value) {
  auto attribute = XmlAttribute(name);
  attribute.setImmediate(value);
  return attribute;
}

class XmlEntry {
  public:
  virtual ~XmlEntry() = default;
  explicit XmlEntry(const std::string& name);

  XmlEntry& addAttribute(const XmlAttribute& attribute);

  void write(XmlInstructor& instructor) const;

  protected:
  virtual void innerWrite(XmlInstructor& instructor) const = 0;
  [[nodiscard]] virtual bool isEmpty() const = 0;

  private:
  std::string name_;
  std::vector<XmlAttribute> attributes_;
};

class XmlNode : public XmlEntry {
  public:
  ~XmlNode() override = default;
  explicit XmlNode(const std::string& name);

  void addNode(const std::shared_ptr<XmlEntry>& entry);

  protected:
  void innerWrite(XmlInstructor& instructor) const override;
  [[nodiscard]] bool isEmpty() const override;

  private:
  std::vector<std::shared_ptr<XmlEntry>> entries_;
};

class XmlData : public XmlEntry {
  public:
  ~XmlData() override = default;
  explicit XmlData(const std::string& name);

  template <typename T>
  void setImmediate(const T& data) {
    this->data_ = writer::WriteInline::create(data);
  }

  template <typename T>
  T getImmediate() const {
    const auto* data = this->data_->getLocalPointer();
    const auto* dataConv = reinterpret_cast<const T*>(data);
    return *dataConv;
  }

  void setDataSource(const std::shared_ptr<writer::DataSource>& dataSource);

  [[nodiscard]] const std::shared_ptr<writer::DataSource>& getDataSource() const;

  protected:
  void innerWrite(XmlInstructor& instructor) const override;
  [[nodiscard]] bool isEmpty() const override;

  private:
  std::shared_ptr<writer::DataSource> data_;
};

template <>
inline void XmlData::setImmediate<std::string>(const std::string& data) {
  this->data_ = writer::WriteInline::createString(data);
}

template <>
inline std::string XmlData::getImmediate<std::string>() const {
  const auto* data = this->data_->getLocalPointer();
  const auto* dataConv = reinterpret_cast<const char*>(data);

  // createString stores the text without a null terminator
  return std::string(dataConv, dataConv + this->data_->getLocalSize());
}

class XmlFile {
  private:
  std::shared_ptr<XmlEntry> root_;

  void write(XmlInstructor& instructor) const;

  public:
  std::shared_ptr<XmlEntry> getRoot();

  void setRoot(std::shared_ptr<XmlEntry> entry);

  [[nodiscard]] std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
      instructions(const std::string& file) const;

  //! The file as text; cf. XmlInstructor::text.
  [[nodiscard]] std::string text() const;
};

/**
 * @brief An XML file that grows by entries at one place, such as by a step per output.
 *
 * Writing the whole file anew every time it grows makes every write longer than the one before,
 * and all of them together take time growing with the square of their number. Here, only the
 * first write writes the file in full. A later one puts its entries where the closing tags behind
 * the earlier ones are, and writes these after them again; so it is as long as what it adds. The
 * file is the same as if written in full.
 *
 * The entries cannot hold distributed data (cf. XmlInstructor::text).
 */
class GrowingXmlFile {
  public:
  //! The file holding the given entries at the place where it grows.
  using Document = std::function<XmlFile(const std::vector<std::shared_ptr<XmlEntry>>&)>;

  explicit GrowingXmlFile(const Document& document);

  //! The instructions adding @p entries to @p file, behind those of the calls before.
  [[nodiscard]] std::vector<std::shared_ptr<writer::instructions::WriteInstruction>>
      add(const std::string& file, const std::vector<std::shared_ptr<XmlEntry>>& entries);

  private:
  //! The file before and behind the entries.
  std::string head_;
  std::string tail_;
  //! The indentation of the line an entry starts on.
  int indent_{0};
  //! Where the tail starts in the file; nothing before the first write.
  std::optional<std::size_t> tailOffset_;
};

} // namespace seissol::io::instance::metadata

#endif // SEISSOL_SRC_IO_INSTANCE_METADATA_XML_H_
