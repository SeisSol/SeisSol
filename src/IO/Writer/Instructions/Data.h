// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DATA_H_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DATA_H_

#include <IO/Datatype/Datatype.h>
#include <IO/Datatype/Inference.h>
#include <cstring>
#include <functional>
#include <memory>
#include <utility>
#include <yaml-cpp/yaml.h>

namespace async {
class ExecInfo;
} // namespace async

namespace seissol::io::writer {

class DataSource {
  public:
  DataSource(std::shared_ptr<datatype::Datatype> datatype, const std::vector<std::size_t>& shape);
  virtual ~DataSource();

  virtual YAML::Node serialize() = 0;
  virtual const void* getPointer(const async::ExecInfo& info) = 0;
  virtual std::size_t count(const async::ExecInfo& info) = 0;
  virtual void assignId(int id) = 0;
  virtual bool distributed() = 0;
  [[nodiscard]] virtual const void* getLocalPointer() const = 0;
  [[nodiscard]] virtual size_t getLocalSize() const = 0;

  [[nodiscard]] const std::vector<std::size_t>& shape() const;
  [[nodiscard]] std::shared_ptr<seissol::io::datatype::Datatype> datatype() const;

  static std::unique_ptr<DataSource> deserialize(YAML::Node node);

  protected:
  std::shared_ptr<seissol::io::datatype::Datatype> datatypeP;
  std::vector<std::size_t> shapeP;
};

class WriteInline : public DataSource {
  public:
  WriteInline(const void* dataPtr,
              std::size_t size,
              std::shared_ptr<datatype::Datatype> datatype,
              const std::vector<std::size_t>& shape);

  explicit WriteInline(YAML::Node node);

  YAML::Node serialize() override;

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  bool distributed() override;

  void assignId(int /*id*/) override;

  [[nodiscard]] const void* getLocalPointer() const override;
  [[nodiscard]] size_t getLocalSize() const override;

  template <typename T>
  static std::shared_ptr<DataSource>
      create(const T& data,
             std::shared_ptr<datatype::Datatype> datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteInline>(&data, sizeof(T), datatype, std::vector<std::size_t>());
  }

  static std::shared_ptr<DataSource> createString(const std::string& data) {
    return std::make_shared<WriteInline>(
        data.data(),
        (data.size() + 1) * sizeof(char),
        std::make_shared<datatype::StringDatatype>(data.size() + 1),
        std::vector<std::size_t>());
  }

  template <typename T>
  static std::shared_ptr<DataSource>
      createArray(const std::vector<std::size_t>& shape,
                  const std::vector<T>& data,
                  std::shared_ptr<datatype::Datatype> datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteInline>(data.data(), sizeof(T) * data.size(), datatype, shape);
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

  void assignId(int /*id*/) override;

  bool distributed() override;

  [[nodiscard]] const void* getLocalPointer() const override;
  [[nodiscard]] size_t getLocalSize() const override;

  private:
  int id;
};

class WriteBuffer : public DataSource {
  public:
  WriteBuffer(const void* data,
              size_t size,
              std::shared_ptr<datatype::Datatype> datatype,
              const std::vector<std::size_t>& shape);

  YAML::Node serialize() override;

  [[nodiscard]] const void* getLocalPointer() const override;
  [[nodiscard]] size_t getLocalSize() const override;

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  void assignId(int givenId) override;

  bool distributed() override;

  template <typename T>
  static std::shared_ptr<DataSource>
      create(const T* data,
             size_t count,
             const std::vector<std::size_t>& shape = {},
             std::shared_ptr<datatype::Datatype> datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteBuffer>(data, count, datatype, shape);
  }

  private:
  const void* data;
  size_t size;
  int id{-1};
};

class AdhocBuffer : public DataSource {
  public:
  [[nodiscard]] virtual std::size_t getTargetSize() const = 0;
  virtual void setData(void* target) = 0;

  AdhocBuffer(std::shared_ptr<datatype::Datatype> datatype, const std::vector<std::size_t>& shape)
      : DataSource(std::move(datatype), shape) {}

  YAML::Node serialize() override {
    YAML::Node node;
    node["id"] = id;
    node["datatype"] = datatype()->serialize();
    node["type"] = "buffer";
    node["shape"] = shape();
    return node;
  }

  const void* getPointer(const async::ExecInfo& info) override { return nullptr; }

  [[nodiscard]] const void* getLocalPointer() const override { return nullptr; }
  [[nodiscard]] size_t getLocalSize() const override { return getTargetSize(); }

  std::size_t count(const async::ExecInfo& info) override {
    return getTargetSize() / datatype()->size();
  }

  void assignId(int givenId) override { id = givenId; }

  bool distributed() override { return true; }

  private:
  int id;
};

class GeneratedBuffer : public AdhocBuffer {
  public:
  GeneratedBuffer(std::size_t sourceCount,
                  std::size_t targetCount,
                  std::function<void(void*)> generator,
                  std::shared_ptr<datatype::Datatype> datatype,
                  const std::vector<std::size_t>& shape)
      : generator(std::move(generator)), sourceCount(sourceCount), targetCount(targetCount),
        targetStride(targetCount), AdhocBuffer(std::move(datatype), shape) {

    for (auto dim : shape) {
      targetStride *= dim;
    }
  }

  [[nodiscard]] std::size_t getTargetSize() const override {
    return targetStride * datatype()->size() * sourceCount;
  }

  void setData(void* targetPtr) override { std::invoke(generator, targetPtr); }

  template <typename T, typename F>
  static std::shared_ptr<GeneratedBuffer> createElementwise(
      std::size_t sourceCount,
      std::size_t targetCount,
      const std::vector<std::size_t>& shape,
      F handler,
      std::shared_ptr<datatype::Datatype> datatype = datatype::inferDatatype<T>()) {
    std::size_t localTargetStride = targetCount;
    for (auto dim : shape) {
      localTargetStride *= dim;
    }
    return std::make_shared<GeneratedBuffer>(
        sourceCount,
        targetCount,
        [=](void* targetPtr) {
          T* target = reinterpret_cast<T*>(targetPtr);
#ifdef _OPENMP
#pragma omp parallel for schedule(static)
#endif
          for (std::size_t i = 0; i < sourceCount; ++i) {
            std::invoke(handler, &target[i * localTargetStride], i);
          }
        },
        datatype,
        shape);
  }

  private:
  std::function<void(void*)> generator;
  std::size_t sourceCount;
  std::size_t targetCount;
  std::size_t targetStride;
};

} // namespace seissol::io::writer

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DATA_H_
