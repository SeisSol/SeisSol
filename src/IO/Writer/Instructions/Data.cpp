// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#include "Data.h"

#include "IO/Datatype/Datatype.h"

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

namespace {
//! @brief What one entry holds: the dimensions that neither move between ranks nor grow.
std::vector<std::size_t> fixedSizes(const std::vector<Dimension>& dimensions) {
  std::vector<std::size_t> result;
  for (const auto& dimension : dimensions) {
    if (!dimension.isDistributed() && !dimension.isAppended()) {
      result.push_back(dimension.size);
    }
  }
  return result;
}
} // namespace

DataSource::DataSource(std::shared_ptr<datatype::Datatype> datatype,
                       const std::vector<std::size_t>& shape,
                       bool leadingDistributed)
    : DataSource(std::move(datatype), makeDimensions(shape, leadingDistributed)) {}

DataSource::DataSource(std::shared_ptr<datatype::Datatype> datatype,
                       std::vector<Dimension> dimensions)
    : datatypeP_(std::move(datatype)), shapeP_(fixedSizes(dimensions)),
      dimensionsP_(std::move(dimensions)) {
  checkDimensions(dimensionsP_, "a data source");
}

DataSource::~DataSource() = default;

std::shared_ptr<seissol::io::datatype::Datatype> DataSource::datatype() const { return datatypeP_; }

const std::vector<Dimension>& DataSource::dimensions() const { return dimensionsP_; }

const std::vector<std::size_t>& DataSource::shape() const { return shapeP_; }

bool DataSource::distributed() const {
  return std::any_of(dimensionsP_.begin(), dimensionsP_.end(), [](const Dimension& dimension) {
    return dimension.isDistributed();
  });
}

WriteInline::WriteInline(const void* dataPtr,
                         std::size_t size,
                         std::shared_ptr<datatype::Datatype> datatype,
                         const std::vector<std::size_t>& shape)
    : DataSource(std::move(datatype), shape, false) {
  data_.resize(size);
  std::memcpy(data_.data(), dataPtr, size);
}

WriteInline::WriteInline(const void* dataPtr,
                         std::size_t size,
                         std::shared_ptr<datatype::Datatype> datatype,
                         std::vector<Dimension> dimensions)
    : DataSource(std::move(datatype), std::move(dimensions)) {
  data_.resize(size);
  std::memcpy(data_.data(), dataPtr, size);
}

WriteInline::WriteInline(YAML::Node node)
    : DataSource(datatype::Datatype::deserialize(node["datatype"]),
                 deserializeDimensions(node["dimensions"])) {
  const auto rawData = node["data"].as<YAML::Binary>();
  data_.resize(rawData.size());
  std::copy_n(rawData.data(), rawData.size(), data_.begin());
}

YAML::Node WriteInline::serialize() {
  YAML::Node node;
  node["type"] = "inline";
  node["data"] = YAML::Binary(const_cast<const unsigned char*>(data_.data()), data_.size());
  node["datatype"] = datatype()->serialize();
  node["dimensions"] = serializeDimensions(dimensions());
  return node;
}

const void* WriteInline::getPointer(const async::ExecInfo& /*info*/) { return data_.data(); }

const void* WriteInline::getLocalPointer() const { return data_.data(); }
std::size_t WriteInline::getLocalSize() const { return data_.size(); }

std::size_t WriteInline::count(const async::ExecInfo& /*info*/) {
  return data_.size() / datatype()->size();
}

void WriteInline::assignId(int /*id*/) {}

bool WriteInline::managed() const { return false; }

bool WriteBufferRemote::managed() const { return true; }

WriteBufferRemote::WriteBufferRemote(YAML::Node node)
    : DataSource(datatype::Datatype::deserialize(node["datatype"]),
                 deserializeDimensions(node["dimensions"])) {
  id_ = node["id"].as<int>();
  datatypeP_ = datatype::Datatype::deserialize(node["datatype"]);
}

YAML::Node WriteBufferRemote::serialize() {
  YAML::Node node;
  node["id"] = id_;
  node["datatype"] = datatype()->serialize();
  node["type"] = "buffer";
  node["dimensions"] = serializeDimensions(dimensions());
  return node;
}

const void* WriteBufferRemote::getPointer(const async::ExecInfo& info) { return info.buffer(id_); }

std::size_t WriteBufferRemote::count(const async::ExecInfo& info) {
  return info.bufferSize(id_) / datatype()->size();
}

void WriteBufferRemote::assignId(int /*id*/) {}

const void* WriteBufferRemote::getLocalPointer() const { return nullptr; }
std::size_t WriteBufferRemote::getLocalSize() const { return 0; }

bool WriteBuffer::managed() const { return true; }

WriteBuffer::WriteBuffer(const void* data,
                         size_t size,
                         std::shared_ptr<datatype::Datatype> datatype,
                         const std::vector<std::size_t>& shape)
    : DataSource(std::move(datatype), shape, true), data_(data), size_(size) {}

WriteBuffer::WriteBuffer(const void* data,
                         size_t size,
                         std::shared_ptr<datatype::Datatype> datatype,
                         std::vector<Dimension> dimensions)
    : DataSource(std::move(datatype), std::move(dimensions)), data_(data), size_(size) {}

YAML::Node WriteBuffer::serialize() {
  YAML::Node node;
  node["id"] = id_;
  node["datatype"] = datatype()->serialize();
  node["type"] = "buffer";
  node["dimensions"] = serializeDimensions(dimensions());
  return node;
}

const void* WriteBuffer::getLocalPointer() const { return data_; }
size_t WriteBuffer::getLocalSize() const {
  std::size_t shapeprod = 1;
  for (auto dim : shape()) {
    shapeprod *= dim;
  }
  return shapeprod * size_ * datatype()->size();
}

const void* WriteBuffer::getPointer(const async::ExecInfo& /*info*/) { return data_; }

std::size_t WriteBuffer::count(const async::ExecInfo& info) {
  return info.bufferSize(id_) / datatype()->size();
}

void WriteBuffer::assignId(int givenId) { id_ = givenId; }

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
