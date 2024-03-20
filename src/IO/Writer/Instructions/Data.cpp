#include "Data.hpp"

namespace seissol::io::writer {

DataSource::DataSource(std::shared_ptr<datatype::Datatype> datatype) : datatypeP(datatype) {}

std::shared_ptr<seissol::io::datatype::Datatype> DataSource::datatype() const { return datatypeP; }

WriteInline::WriteInline(void* dataPtr,
                         std::size_t size,
                         std::shared_ptr<datatype::Datatype> datatype)
    : DataSource(datatype) {
  data.resize(size);
  std::memcpy(data.data(), dataPtr, size);
}

WriteInline::WriteInline(YAML::Node node)
    : DataSource(datatype::Datatype::deserialize(node["datatype"])) {
  YAML::Binary rawData = node["data"].as<YAML::Binary>();
  data.resize(rawData.size());
  std::copy_n(rawData.data(), rawData.size(), data.begin());
}

YAML::Node WriteInline::serialize() {
  YAML::Node node;
  node["type"] = "inline";
  node["data"] = YAML::Binary(const_cast<const unsigned char*>(data.data()), data.size());
  node["datatype"] = datatype()->serialize();
  return node;
}

const void* WriteInline::getPointer(const async::ExecInfo& info) { return data.data(); }

std::size_t WriteInline::count(const async::ExecInfo& info) { return 1; }

void WriteInline::assignId(int) {}

WriteBufferRemote::WriteBufferRemote(YAML::Node node)
    : DataSource(datatype::Datatype::deserialize(node["datatype"])) {
  id = node["id"].as<int>();
  datatypeP = datatype::Datatype::deserialize(node["datatype"]);
}

YAML::Node WriteBufferRemote::serialize() {
  YAML::Node node;
  node["id"] = id;
  node["datatype"] = datatype()->serialize();
  node["type"] = "buffer";
  return node;
}

const void* WriteBufferRemote::getPointer(const async::ExecInfo& info) { return info.buffer(id); }

std::size_t WriteBufferRemote::count(const async::ExecInfo& info) {
  return info.bufferSize(id) / datatype()->size();
}

void WriteBufferRemote::assignId(int) {}

WriteBuffer::WriteBuffer(void* data, size_t size, std::shared_ptr<datatype::Datatype> datatype)
    : id(-1), data(data), size(size), DataSource(datatype) {}

YAML::Node WriteBuffer::serialize() {
  YAML::Node node;
  node["id"] = id;
  node["datatype"] = datatype()->serialize();
  node["type"] = "buffer";
  return node;
}

void* WriteBuffer::getLocalPointer() { return data; }
size_t WriteBuffer::getLocalSize() { return size; }

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
