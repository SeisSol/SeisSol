#pragma once

#include "async/ExecInfo.h"
#include <IO/Datatype/Datatype.hpp>
#include <IO/Datatype/Inference.hpp>
#include <cstring>
#include <exception>
#include <functional>
#include <yaml-cpp/yaml.h>

namespace seissol::io::writer {

class DataSource {
  public:
  DataSource(std::shared_ptr<datatype::Datatype> datatype);

  virtual YAML::Node serialize() = 0;
  virtual const void* getPointer(const async::ExecInfo& info) = 0;
  virtual std::size_t count(const async::ExecInfo& info) = 0;
  virtual void assignId(int id) = 0;

  std::shared_ptr<seissol::io::datatype::Datatype> datatype() const;

  static std::unique_ptr<DataSource> deserialize(YAML::Node node);

  protected:
  std::shared_ptr<seissol::io::datatype::Datatype> datatypeP;
};

class WriteInline : public DataSource {
  public:
  WriteInline(void* dataPtr, std::size_t size, std::shared_ptr<datatype::Datatype> datatype);

  explicit WriteInline(YAML::Node node);

  YAML::Node serialize() override;

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  void assignId(int) override;

  template <typename T>
  static std::shared_ptr<DataSource>
      create(const T& data,
             std::shared_ptr<datatype::Datatype> datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteInline>(&data, sizeof(T), datatype);
  }

  private:
  std::vector<unsigned char> data;
};

class WriteBufferRemote : public DataSource {
  public:
  explicit WriteBufferRemote(YAML::Node node);

  YAML::Node serialize() override;

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  void assignId(int) override;

  private:
  int id;
};

class WriteBuffer : public DataSource {
  public:
  WriteBuffer(void* data, size_t size, std::shared_ptr<datatype::Datatype> datatype);

  YAML::Node serialize() override;

  void* getLocalPointer();
  size_t getLocalSize();

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  void assignId(int givenId) override;

  template <typename T>
  static std::shared_ptr<DataSource>
      create(T* data,
             size_t count,
             std::shared_ptr<datatype::Datatype> datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteBuffer>(data, sizeof(T) * count, datatype);
  }

  private:
  void* data;
  size_t size;
  int id;
};

class AdhocBuffer : public DataSource {
  public:
  virtual std::size_t getTargetSize() = 0;
  virtual void setData(void* target) = 0;

  YAML::Node serialize() override {
    YAML::Node node;
    node["id"] = id;
    node["datatype"] = datatype()->serialize();
    node["type"] = "buffer";
    return node;
  }

  const void* getPointer(const async::ExecInfo& info) override { return nullptr; }

  std::size_t count(const async::ExecInfo& info) override {
    return getTargetSize() / datatype()->size();
  }

  void assignId(int givenId) override { id = givenId; }

  private:
  int id;
};

template <std::size_t TargetCount, typename T, typename S, typename F>
class ElementTransformBuffer : public AdhocBuffer {
  public:
  ElementTransformBuffer(F&& handler, T* source, std::size_t count)
      : handler(std::forward<F>(handler)), source(source), sourceCount(count) {}

  void setData(void* targetPtr) override {
    S* target = targetPtr;
    for (std::size_t i = 0; i < sourceCount; ++i) {
      using TargetArray = S[TargetCount];
      std::invoke(handler, *reinterpret_cast<TargetArray*>(&target[i * TargetCount]), source[i]);
    }
  }

  std::size_t getTargetSize() override { return TargetCount * sizeof(S) * sourceCount; }

  private:
  F handler;
  T* source;
  size_t sourceCount;
};

} // namespace seissol::io::writer
