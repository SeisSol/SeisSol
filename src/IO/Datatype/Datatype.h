// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_DATATYPE_DATATYPE_H_
#define SEISSOL_SRC_IO_DATATYPE_DATATYPE_H_

#include <memory>
#include <optional>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::datatype {

class Datatype;

class Array {
  public:
  Array(std::shared_ptr<Datatype> type, const std::vector<std::size_t>& dimensions);

  private:
  std::shared_ptr<Datatype> type_;
  std::vector<std::size_t> dimensions_;
};

class Datatype : public std::enable_shared_from_this<Datatype> {
  public:
  virtual ~Datatype();
  [[nodiscard]] virtual std::size_t size() const = 0;
  [[nodiscard]] virtual YAML::Node serialize() const = 0;
  virtual Array unwrap(std::size_t maxDimensions);
  virtual std::string toStringRaw(const void* data) const = 0;
  [[nodiscard]] virtual std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const = 0;

  template <typename T>
  std::string toString(const T& data) const {
    return toStringRaw(&data);
  }
  template <typename T>
  std::optional<T> fromString(const std::string& str) const {
    const auto result = fromStringRaw(str);
    if (result.has_value()) {
      const char* dataRaw = result.value().data();
      const T* data = reinterpret_cast<const T*>(dataRaw);
      return std::make_optional<T>(*data);
    } else {
      return std::optional<T>();
    }
  }
  static std::shared_ptr<Datatype> deserialize(YAML::Node node);
};

class OpaqueDatatype : public Datatype {
  public:
  explicit OpaqueDatatype(std::size_t size);

  explicit OpaqueDatatype(YAML::Node node);

  [[nodiscard]] std::size_t size() const override;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;

  private:
  std::size_t sizeP_;
};

class StringDatatype : public Datatype {
  public:
  explicit StringDatatype(std::size_t size);

  explicit StringDatatype(YAML::Node node);

  [[nodiscard]] std::size_t size() const override;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;

  private:
  std::size_t sizeP_;
};

class F32Datatype : public Datatype {
  public:
  [[nodiscard]] std::size_t size() const override;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;
};

class F64Datatype : public Datatype {
  public:
  [[nodiscard]] std::size_t size() const override;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;
};

class F80Datatype : public Datatype {
  public:
  [[nodiscard]] std::size_t size() const override;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;
};

class IntegerDatatype : public Datatype {
  public:
  IntegerDatatype(std::size_t size, bool sign);

  explicit IntegerDatatype(YAML::Node node);

  [[nodiscard]] std::size_t size() const override;

  [[nodiscard]] bool sign() const;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;

  private:
  std::size_t sizeP_;
  bool signP_;
};

class ArrayDatatype : public Datatype {
  public:
  ArrayDatatype(std::shared_ptr<Datatype> base, const std::vector<std::size_t>& dimensions);

  explicit ArrayDatatype(YAML::Node node);

  [[nodiscard]] std::size_t size() const override;

  Array unwrap(std::size_t maxDimensions) override;

  [[nodiscard]] const std::vector<std::size_t>& dimensions() const;

  [[nodiscard]] std::shared_ptr<Datatype> base() const;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;

  private:
  std::shared_ptr<Datatype> baseP_;
  std::vector<std::size_t> dimensionsP_;
};

class StructDatatype : public Datatype {
  public:
  struct MemberInfo {
    std::string name;
    std::size_t offset{};
    std::shared_ptr<Datatype> datatype;
  };

  explicit StructDatatype(const std::vector<MemberInfo>& members);

  StructDatatype(const std::vector<MemberInfo>& members, std::size_t size);

  explicit StructDatatype(YAML::Node node);

  [[nodiscard]] std::size_t size() const override;

  [[nodiscard]] const std::vector<MemberInfo>& members() const;

  [[nodiscard]] YAML::Node serialize() const override;

  std::string toStringRaw(const void* data) const override;
  [[nodiscard]] std::optional<std::vector<char>>
      fromStringRaw(const std::string& str) const override;

  private:
  static std::size_t minSize(const std::vector<MemberInfo>& members);

  std::size_t sizeP_;
  std::vector<MemberInfo> membersP_;
};

} // namespace seissol::io::datatype

#endif // SEISSOL_SRC_IO_DATATYPE_DATATYPE_H_
