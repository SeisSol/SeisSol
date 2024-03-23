#include "Datatype.hpp"

#include <cassert>
#include <exception>
#include <memory>
#include <numeric>
#include <string>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace {
static std::vector<seissol::io::datatype::StructDatatype::MemberInfo>
    deserializeMemberInfo(YAML::Node members) {
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
} // namespace

namespace seissol::io::datatype {

Array::Array(std::shared_ptr<Datatype> type, const std::vector<std::size_t>& dimensions)
    : type(type), dimensions(dimensions) {}

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

StringDatatype::StringDatatype(std::size_t size) : sizeP(size) {}

std::size_t StringDatatype::size() const { return sizeP; }

YAML::Node StringDatatype::serialize() const {
  YAML::Node node;
  node["type"] = "string";
  node["size"] = sizeP;
  return node;
}

StringDatatype::StringDatatype(YAML::Node node) : sizeP(node["size"].as<size_t>()) {}

std::size_t F32Datatype::size() const { return 4; }

YAML::Node F32Datatype::serialize() const {
  YAML::Node node;
  node["type"] = "f32";
  return node;
}

std::size_t F64Datatype::size() const { return 8; }

YAML::Node F64Datatype::serialize() const {
  YAML::Node node;
  node["type"] = "f64";
  return node;
}

std::size_t F80Datatype::size() const { return 10; }

YAML::Node F80Datatype::serialize() const {
  YAML::Node node;
  node["type"] = "f80";
  return node;
}

IntegerDatatype::IntegerDatatype(std::size_t size, bool sign) : sizeP(size), signP(sign) {}

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

ArrayDatatype::ArrayDatatype(std::shared_ptr<Datatype> base,
                             const std::vector<std::size_t>& dimensions)
    : baseP(base), dimensionsP(dimensions) {}

ArrayDatatype::ArrayDatatype(YAML::Node node)
    : baseP(Datatype::deserialize(node["base"])),
      dimensionsP(node["shape"].as<std::vector<std::size_t>>()) {}

std::size_t ArrayDatatype::size() const {
  return baseP->size() * std::reduce(dimensionsP.begin(),
                                     dimensionsP.end(),
                                     static_cast<std::size_t>(1),
                                     [](auto a, auto b) { return a * b; });
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
  for (auto member : membersP) {
    YAML::Node membernode;
    membernode["name"] = member.name;
    membernode["offset"] = member.offset;
    membernode["datatype"] = member.datatype->serialize();
    members.push_back(membernode);
  }
  node["mebers"] = members;
  return node;
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
