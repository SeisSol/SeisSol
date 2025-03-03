// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Datatype.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <iomanip>
#include <memory>
#include <optional>
#include <sstream>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {
std::vector<seissol::io::datatype::StructDatatype::MemberInfo>
    deserializeMemberInfo(const YAML::Node& members) {
  std::vector<seissol::io::datatype::StructDatatype::MemberInfo> memberInfo;
  for (YAML::Node member : members) {
    seissol::io::datatype::StructDatatype::MemberInfo info;
    info.name = member["name"].as<std::string>();
    info.offset = member["offset"].as<std::size_t>();
    info.datatype = seissol::io::datatype::Datatype::deserialize(member["datatype"]);
    memberInfo.push_back(info);
  }
  return memberInfo;
}

const std::string Base64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
const std::string Base64Pad = "=";
const std::array<int, 256> FromBase64 = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};

template <int Count>
void base64Convert(std::ostringstream& sstr, int dataPtr1, int dataPtr2, int dataPtr3) {
  const int data = dataPtr1 | (dataPtr2 << 8) | (dataPtr3 << 16);
  if constexpr (Count > 0) {
    const auto data1 = data & 0x3f;
    sstr << Base64[data1];
  }
  if constexpr (Count > 1) {
    const auto data2 = (data >> 6) & 0x3f;
    sstr << Base64[data2];
  }
  if constexpr (Count > 2) {
    const auto data3 = (data >> 12) & 0x3f;
    sstr << Base64[data3];
  }
  if constexpr (Count > 3) {
    const auto data4 = (data >> 18) & 0x3f;
    sstr << Base64[data4];
  }
}

template <int Count>
void base64Revert(const std::string& idata, std::size_t ipos, char* odata, std::size_t opos) {
  const auto idata1 = FromBase64[idata[ipos]];
  const auto idata2 = FromBase64[idata[ipos + 1]];
  const auto idata3 = FromBase64[idata[ipos + 2]];
  const auto idata4 = FromBase64[idata[ipos + 3]];
  const int idataNumber = idata1 | (idata2 << 6) | (idata3 << 12) | (idata4 << 18);
  if constexpr (Count > 0) {
    const auto data = idataNumber & 0xff;
    odata[opos] = data;
  }
  if constexpr (Count > 1) {
    const auto data = (idataNumber >> 8) & 0xff;
    odata[opos + 1] = data;
  }
  if constexpr (Count > 2) {
    const auto data = (idataNumber >> 16) & 0xff;
    odata[opos + 2] = data;
  }
}

template <typename T>
std::string toStringRawPrimitive(const void* data, int precision) {
  const auto value = *reinterpret_cast<const T*>(data);
  std::ostringstream sstr;
  sstr << std::setprecision(precision) << value;
  return sstr.str();
}

template <typename T>
std::optional<std::vector<char>> fromStringRawPrimitive(const std::string& str) {
  std::istringstream sstr(str);
  T value;
  sstr >> value;
  std::vector<char> data(sizeof(T));
  const char* valuePtr = reinterpret_cast<const char*>(&value);
  std::copy_n(valuePtr, sizeof(T), data.begin());
  return std::make_optional(data);
}
} // namespace

namespace seissol::io::datatype {

Datatype::~Datatype() = default;

Array::Array(std::shared_ptr<Datatype> type, const std::vector<std::size_t>& dimensions)
    : type(std::move(type)), dimensions(dimensions) {}

Array Datatype::unwrap(std::size_t maxDimensions) { return Array(shared_from_this(), {}); }

OpaqueDatatype::OpaqueDatatype(std::size_t size) : sizeP(size) {}

std::size_t OpaqueDatatype::size() const { return sizeP; }

YAML::Node OpaqueDatatype::serialize() const {
  YAML::Node node;
  node["type"] = "opaque";
  node["size"] = sizeP;
  return node;
}

OpaqueDatatype::OpaqueDatatype(YAML::Node node) : sizeP(node["size"].as<size_t>()) {}

std::string OpaqueDatatype::toStringRaw(const void* data) const {
  std::ostringstream sstr;
  const char* dataPtr = reinterpret_cast<const char*>(data);
  for (std::size_t i = 0; i < sizeP; i += 3) {
    base64Convert<4>(sstr, dataPtr[i], dataPtr[i + 1], dataPtr[i + 2]);
  }
  if (sizeP % 3 == 1) {
    base64Convert<2>(sstr, dataPtr[sizeP - 1], 0, 0);
    sstr << Base64Pad << Base64Pad;
  }
  if (sizeP % 3 == 2) {
    base64Convert<3>(sstr, dataPtr[sizeP - 2], dataPtr[sizeP - 3], 0);
    sstr << Base64Pad;
  }

  return sstr.str();
}
std::optional<std::vector<char>> OpaqueDatatype::fromStringRaw(const std::string& str) const {
  std::vector<char> data(sizeP);
  for (std::size_t i = 0; i < (sizeP + 3) / 4; ++i) {
    base64Revert<3>(str, 4 * i, data.data(), 3 * i);
  }
  return std::make_optional(data);
}

StringDatatype::StringDatatype(std::size_t size) : sizeP(size) {}

std::size_t StringDatatype::size() const { return sizeP; }

YAML::Node StringDatatype::serialize() const {
  YAML::Node node;
  node["type"] = "string";
  node["size"] = sizeP;
  return node;
}

StringDatatype::StringDatatype(YAML::Node node) : sizeP(node["size"].as<size_t>()) {}

std::string StringDatatype::toStringRaw(const void* data) const {
  const char* dataPtr = reinterpret_cast<const char*>(data);
  return std::string(dataPtr, dataPtr + sizeP);
}
std::optional<std::vector<char>> StringDatatype::fromStringRaw(const std::string& str) const {
  return std::make_optional(std::vector<char>(str.begin(), str.end()));
}

std::size_t F32Datatype::size() const { return 4; }

YAML::Node F32Datatype::serialize() const {
  YAML::Node node;
  node["type"] = "f32";
  return node;
}

std::string F32Datatype::toStringRaw(const void* data) const {
  return toStringRawPrimitive<float>(data, 8);
}
std::optional<std::vector<char>> F32Datatype::fromStringRaw(const std::string& str) const {
  return fromStringRawPrimitive<float>(str);
}

std::size_t F64Datatype::size() const { return 8; }

YAML::Node F64Datatype::serialize() const {
  YAML::Node node;
  node["type"] = "f64";
  return node;
}

std::string F64Datatype::toStringRaw(const void* data) const {
  return toStringRawPrimitive<double>(data, 16);
}
std::optional<std::vector<char>> F64Datatype::fromStringRaw(const std::string& str) const {
  return fromStringRawPrimitive<double>(str);
}

std::size_t F80Datatype::size() const { return 10; }

YAML::Node F80Datatype::serialize() const {
  YAML::Node node;
  node["type"] = "f80";
  return node;
}

std::string F80Datatype::toStringRaw(const void* data) const {
  return toStringRawPrimitive<long double>(data, 20);
}
std::optional<std::vector<char>> F80Datatype::fromStringRaw(const std::string& str) const {
  return fromStringRawPrimitive<long double>(str);
}

IntegerDatatype::IntegerDatatype(std::size_t size, bool sign) : sizeP(size), signP(sign) {
  assert(size > 0);
}

IntegerDatatype::IntegerDatatype(YAML::Node node)
    : sizeP(node["size"].as<std::size_t>()), signP(node["sign"].as<bool>()) {}

std::size_t IntegerDatatype::size() const { return sizeP; }

bool IntegerDatatype::sign() const { return signP; }

YAML::Node IntegerDatatype::serialize() const {
  YAML::Node node;
  node["type"] = "int";
  node["sign"] = signP;
  node["size"] = sizeP;
  return node;
}

std::string IntegerDatatype::toStringRaw(const void* data) const {
  // for now
  return toStringRawPrimitive<long long>(data, 0);
}
std::optional<std::vector<char>> IntegerDatatype::fromStringRaw(const std::string& str) const {
  // for now
  return fromStringRawPrimitive<long long>(str);
}

ArrayDatatype::ArrayDatatype(std::shared_ptr<Datatype> base,
                             const std::vector<std::size_t>& dimensions)
    : baseP(std::move(base)), dimensionsP(dimensions) {}

ArrayDatatype::ArrayDatatype(YAML::Node node)
    : baseP(Datatype::deserialize(node["base"])),
      dimensionsP(node["shape"].as<std::vector<std::size_t>>()) {}

std::size_t ArrayDatatype::size() const {
  std::size_t totalSize = baseP->size();
  for (const auto& dim : dimensionsP) {
    totalSize *= dim;
  }
  return totalSize;
}

Array ArrayDatatype::unwrap(std::size_t maxDimensions) { return Array(baseP, dimensionsP); }

const std::vector<std::size_t>& ArrayDatatype::dimensions() const { return dimensionsP; }

std::shared_ptr<Datatype> ArrayDatatype::base() const { return baseP; }

YAML::Node ArrayDatatype::serialize() const {
  YAML::Node node;
  node["type"] = "array";
  node["base"] = baseP->serialize();
  node["shape"] = dimensionsP;
  return node;
}

std::string ArrayDatatype::toStringRaw(const void* data) const { return ""; }
std::optional<std::vector<char>> ArrayDatatype::fromStringRaw(const std::string& str) const {
  return std::optional<std::vector<char>>();
}

struct MemberInfo {
  std::string name;
  std::size_t offset;
  std::shared_ptr<Datatype> datatype;
};

StructDatatype::StructDatatype(const std::vector<MemberInfo>& members)
    : StructDatatype(members, minSize(members)) {}

StructDatatype::StructDatatype(const std::vector<MemberInfo>& members, std::size_t size)
    : membersP(members), sizeP(size) {
  assert(size >= minSize(members));
}

StructDatatype::StructDatatype(YAML::Node node)
    : sizeP(node["size"].as<std::size_t>()), membersP(deserializeMemberInfo(node["members"])) {}

std::size_t StructDatatype::size() const { return sizeP; }

const std::vector<StructDatatype::MemberInfo>& StructDatatype::members() const { return membersP; }

YAML::Node StructDatatype::serialize() const {
  YAML::Node node;
  node["type"] = "struct";
  node["size"] = sizeP;
  std::vector<YAML::Node> members;
  for (const auto& member : membersP) {
    YAML::Node membernode;
    membernode["name"] = member.name;
    membernode["offset"] = member.offset;
    membernode["datatype"] = member.datatype->serialize();
    members.push_back(membernode);
  }
  node["members"] = members;
  return node;
}

std::string StructDatatype::toStringRaw(const void* data) const { return ""; }
std::optional<std::vector<char>> StructDatatype::fromStringRaw(const std::string& str) const {
  return std::optional<std::vector<char>>();
}

std::size_t StructDatatype::minSize(const std::vector<MemberInfo>& members) {
  std::size_t size = 0;
  for (const auto& member : members) {
    size = std::max(size, member.offset + member.datatype->size());
  }
  return size;
}

std::shared_ptr<Datatype> Datatype::deserialize(YAML::Node node) {
  auto type = node["type"].as<std::string>();
  if (type == "f32") {
    return std::make_shared<F32Datatype>();
  }
  if (type == "f64") {
    return std::make_shared<F64Datatype>();
  }
  if (type == "f80") {
    return std::make_shared<F80Datatype>();
  }
  if (type == "opaque") {
    return std::make_shared<OpaqueDatatype>(node);
  }
  if (type == "string") {
    return std::make_shared<StringDatatype>(node);
  }
  if (type == "int") {
    return std::make_shared<IntegerDatatype>(node);
  }
  if (type == "array") {
    return std::make_shared<ArrayDatatype>(node);
  }
  if (type == "struct") {
    return std::make_shared<StructDatatype>(node);
  }
  // error
  return std::make_shared<OpaqueDatatype>(0);
}

} // namespace seissol::io::datatype
