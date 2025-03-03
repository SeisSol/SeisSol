// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Data.h"
#include <IO/Datatype/Datatype.h>
#include <algorithm>
#include <async/ExecInfo.h>
#include <cstddef>
#include <cstring>
#include <memory>
#include <string>
#include <utility>
#include <vector>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

DataSource::DataSource(std::shared_ptr<datatype::Datatype> datatype,
                       const std::vector<std::size_t>& shape)
    : datatypeP(std::move(datatype)), shapeP(shape) {}

DataSource::~DataSource() = default;

std::shared_ptr<seissol::io::datatype::Datatype> DataSource::datatype() const { return datatypeP; }

const std::vector<std::size_t>& DataSource::shape() const { return shapeP; }

WriteInline::WriteInline(const void* dataPtr,
                         std::size_t size,
                         std::shared_ptr<datatype::Datatype> datatype,
                         const std::vector<std::size_t>& shape)
    : DataSource(std::move(datatype), shape) {
  data.resize(size);
  std::memcpy(data.data(), dataPtr, size);
}

WriteInline::WriteInline(YAML::Node node)
    : DataSource(datatype::Datatype::deserialize(node["datatype"]),
                 node["shape"].as<std::vector<std::size_t>>()) {
  const auto rawData = node["data"].as<YAML::Binary>();
  data.resize(rawData.size());
  std::copy_n(rawData.data(), rawData.size(), data.begin());
}

YAML::Node WriteInline::serialize() {
  YAML::Node node;
  node["type"] = "inline";
  node["data"] = YAML::Binary(const_cast<const unsigned char*>(data.data()), data.size());
  node["datatype"] = datatype()->serialize();
  node["shape"] = shape();
  return node;
}

const void* WriteInline::getPointer(const async::ExecInfo& info) { return data.data(); }

const void* WriteInline::getLocalPointer() const { return data.data(); }
std::size_t WriteInline::getLocalSize() const { return data.size(); }

std::size_t WriteInline::count(const async::ExecInfo& info) {
  return data.size() / datatype()->size();
}

void WriteInline::assignId(int /*id*/) {}

bool WriteInline::distributed() { return false; }

bool WriteBufferRemote::distributed() { return true; }

WriteBufferRemote::WriteBufferRemote(YAML::Node node)
    : DataSource(datatype::Datatype::deserialize(node["datatype"]),
                 node["shape"].as<std::vector<std::size_t>>()) {
  id = node["id"].as<int>();
  datatypeP = datatype::Datatype::deserialize(node["datatype"]);
}

YAML::Node WriteBufferRemote::serialize() {
  YAML::Node node;
  node["id"] = id;
  node["datatype"] = datatype()->serialize();
  node["type"] = "buffer";
  node["shape"] = shape();
  return node;
}

const void* WriteBufferRemote::getPointer(const async::ExecInfo& info) { return info.buffer(id); }

std::size_t WriteBufferRemote::count(const async::ExecInfo& info) {
  return info.bufferSize(id) / datatype()->size();
}

void WriteBufferRemote::assignId(int /*id*/) {}

const void* WriteBufferRemote::getLocalPointer() const { return nullptr; }
std::size_t WriteBufferRemote::getLocalSize() const { return 0; }

bool WriteBuffer::distributed() { return true; }

WriteBuffer::WriteBuffer(const void* data,
                         size_t size,
                         std::shared_ptr<datatype::Datatype> datatype,
                         const std::vector<std::size_t>& shape)
    : data(data), size(size), DataSource(std::move(datatype), shape) {}

YAML::Node WriteBuffer::serialize() {
  YAML::Node node;
  node["id"] = id;
  node["datatype"] = datatype()->serialize();
  node["type"] = "buffer";
  node["shape"] = shape();
  return node;
}

const void* WriteBuffer::getLocalPointer() const { return data; }
size_t WriteBuffer::getLocalSize() const {
  std::size_t shapeprod = 1;
  for (auto dim : shape()) {
    shapeprod *= dim;
  }
  return shapeprod * size * datatype()->size();
}

const void* WriteBuffer::getPointer(const async::ExecInfo& info) { return data; }

std::size_t WriteBuffer::count(const async::ExecInfo& info) {
  return info.bufferSize(id) / datatype()->size();
}

void WriteBuffer::assignId(int givenId) { id = givenId; }

std::unique_ptr<DataSource> DataSource::deserialize(YAML::Node node) {
  const auto nodeType = node["type"].as<std::string>();
  if (nodeType == "inline") {
    return std::make_unique<WriteInline>(node);
  } else if (nodeType == "buffer") {
    return std::make_unique<WriteBufferRemote>(node);
  } else {
    throw std::exception();
  }
}

} // namespace seissol::io::writer
