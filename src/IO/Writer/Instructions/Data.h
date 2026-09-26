// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DATA_H_
#define SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DATA_H_

#include "Dimension.h"
#include "IO/Datatype/Datatype.h"
#include "IO/Datatype/Inference.h"

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
  DataSource(std::shared_ptr<datatype::Datatype> datatype,
             const std::vector<std::size_t>& shape,
             bool leadingDistributed);
  DataSource(std::shared_ptr<datatype::Datatype> datatype, std::vector<Dimension> dimensions);
  virtual ~DataSource();

  virtual YAML::Node serialize() = 0;
  virtual const void* getPointer(const async::ExecInfo& info) = 0;
  virtual std::size_t count(const async::ExecInfo& info) = 0;
  virtual void assignId(int id) = 0;
  [[nodiscard]] virtual const void* getLocalPointer() const = 0;
  [[nodiscard]] virtual size_t getLocalSize() const = 0;

  /**
   * @brief Whether the data has to travel to the executor in a buffer of its own.
   *
   * This is a property of where the memory comes from, not of how the data is laid out across the
   * ranks: data carried inside the plan needs no buffer, everything else does.
   */
  [[nodiscard]] virtual bool managed() const = 0;

  //! @brief The full shape: at most one dimension distributed, at most one appended.
  [[nodiscard]] const std::vector<Dimension>& dimensions() const;
  /**
   * @brief The sizes of the dimensions that neither move between ranks nor grow between writes.
   *
   * What one entry of the data holds, in other words. An attribute or an Xdmf payload describes
   * itself with these; the dataset writer works off dimensions() instead.
   */
  [[nodiscard]] const std::vector<std::size_t>& shape() const;
  [[nodiscard]] bool distributed() const;
  [[nodiscard]] std::shared_ptr<seissol::io::datatype::Datatype> datatype() const;

  static std::unique_ptr<DataSource> deserialize(YAML::Node node);

  protected:
  std::shared_ptr<seissol::io::datatype::Datatype> datatypeP_;
  std::vector<std::size_t> shapeP_;
  std::vector<Dimension> dimensionsP_;
};

class WriteInline : public DataSource {
  public:
  WriteInline(const void* dataPtr,
              std::size_t size,
              std::shared_ptr<datatype::Datatype> datatype,
              const std::vector<std::size_t>& shape);

  WriteInline(const void* dataPtr,
              std::size_t size,
              std::shared_ptr<datatype::Datatype> datatype,
              std::vector<Dimension> dimensions);

  explicit WriteInline(YAML::Node node);

  YAML::Node serialize() override;

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  [[nodiscard]] bool managed() const override;

  void assignId(int /*id*/) override;

  [[nodiscard]] const void* getLocalPointer() const override;
  [[nodiscard]] size_t getLocalSize() const override;

  template <typename T>
  static std::shared_ptr<DataSource>
      create(const T& data,
             const std::shared_ptr<datatype::Datatype>& datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteInline>(&data, sizeof(T), datatype, std::vector<std::size_t>());
  }

  static std::shared_ptr<DataSource> createString(const std::string& data) {
    return std::make_shared<WriteInline>(data.data(),
                                         (data.size()) * sizeof(char),
                                         std::make_shared<datatype::StringDatatype>(data.size()),
                                         std::vector<std::size_t>());
  }

  template <typename T>
  static std::shared_ptr<DataSource> createArray(
      const std::vector<std::size_t>& shape,
      const std::vector<T>& data,
      const std::shared_ptr<datatype::Datatype>& datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteInline>(data.data(), sizeof(T) * data.size(), datatype, shape);
  }

  //! @brief An array whose shape says how it joins what is already in the dataset.
  template <typename T>
  static std::shared_ptr<DataSource> createShaped(
      const std::vector<Dimension>& dimensions,
      const std::vector<T>& data,
      const std::shared_ptr<datatype::Datatype>& datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteInline>(
        data.data(), sizeof(T) * data.size(), datatype, dimensions);
  }

  private:
  std::vector<unsigned char> data_;
};

class WriteBufferRemote : public DataSource {
  public:
  explicit WriteBufferRemote(YAML::Node node);

  YAML::Node serialize() override;

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  void assignId(int /*id*/) override;

  [[nodiscard]] bool managed() const override;

  [[nodiscard]] const void* getLocalPointer() const override;
  [[nodiscard]] size_t getLocalSize() const override;

  private:
  int id_;
};

class WriteBuffer : public DataSource {
  public:
  WriteBuffer(const void* data,
              size_t size,
              std::shared_ptr<datatype::Datatype> datatype,
              const std::vector<std::size_t>& shape);

  WriteBuffer(const void* data,
              size_t size,
              std::shared_ptr<datatype::Datatype> datatype,
              std::vector<Dimension> dimensions);

  YAML::Node serialize() override;

  [[nodiscard]] const void* getLocalPointer() const override;
  [[nodiscard]] size_t getLocalSize() const override;

  const void* getPointer(const async::ExecInfo& info) override;

  std::size_t count(const async::ExecInfo& info) override;

  void assignId(int givenId) override;

  [[nodiscard]] bool managed() const override;

  template <typename T>
  static std::shared_ptr<DataSource>
      create(const T* data,
             size_t count,
             const std::vector<std::size_t>& shape = {},
             const std::shared_ptr<datatype::Datatype>& datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteBuffer>(data, count, datatype, shape);
  }

  //! @brief A buffer whose shape says how it joins what is already in the dataset.
  template <typename T>
  static std::shared_ptr<DataSource> createShaped(
      const T* data,
      size_t count,
      const std::vector<Dimension>& dimensions,
      const std::shared_ptr<datatype::Datatype>& datatype = datatype::inferDatatype<T>()) {
    return std::make_shared<WriteBuffer>(data, count, datatype, dimensions);
  }

  private:
  const void* data_;
  size_t size_;
  int id_{-1};
};

class AdhocBuffer : public DataSource {
  public:
  [[nodiscard]] virtual std::size_t getTargetSize() const = 0;
  virtual void setData(void* target) = 0;

  AdhocBuffer(std::shared_ptr<datatype::Datatype> datatype, const std::vector<std::size_t>& shape)
      : DataSource(std::move(datatype), shape, true) {}

  AdhocBuffer(std::shared_ptr<datatype::Datatype> datatype, std::vector<Dimension> dimensions)
      : DataSource(std::move(datatype), std::move(dimensions)) {}

  YAML::Node serialize() override {
    YAML::Node node;
    node["id"] = id_;
    node["datatype"] = datatype()->serialize();
    node["type"] = "buffer";
    node["dimensions"] = serializeDimensions(dimensions());
    return node;
  }

  const void* getPointer(const async::ExecInfo& /*info*/) override { return nullptr; }

  [[nodiscard]] const void* getLocalPointer() const override { return nullptr; }
  [[nodiscard]] size_t getLocalSize() const override { return getTargetSize(); }

  std::size_t count(const async::ExecInfo& /*info*/) override {
    return getTargetSize() / datatype()->size();
  }

  void assignId(int givenId) override { id_ = givenId; }

  [[nodiscard]] bool managed() const override { return true; }

  private:
  int id_{-1};
};

class GeneratedBuffer : public AdhocBuffer {
  public:
  GeneratedBuffer(std::size_t sourceCount,
                  std::size_t targetCount,
                  std::function<void(void*)> generator,
                  std::shared_ptr<datatype::Datatype> datatype,
                  const std::vector<std::size_t>& shape)
      : GeneratedBuffer(sourceCount,
                        targetCount,
                        std::move(generator),
                        std::move(datatype),
                        makeDimensions(shape, true)) {}

  GeneratedBuffer(std::size_t sourceCount,
                  std::size_t targetCount,
                  std::function<void(void*)> generator,
                  std::shared_ptr<datatype::Datatype> datatype,
                  std::vector<Dimension> dimensions)
      : AdhocBuffer(std::move(datatype), std::move(dimensions)), generator_(std::move(generator)),
        sourceCount_(sourceCount), targetStride_(targetCount) {

    for (auto dim : shape()) {
      targetStride_ *= dim;
    }
  }

  [[nodiscard]] std::size_t getTargetSize() const override {
    return targetStride_ * datatype()->size() * sourceCount_;
  }

  void setData(void* targetPtr) override { std::invoke(generator_, targetPtr); }

  template <typename T, typename F>
  static std::shared_ptr<GeneratedBuffer> createElementwise(
      std::size_t sourceCount,
      std::size_t targetCount,
      const std::vector<std::size_t>& shape,
      const F& handler,
      const std::shared_ptr<datatype::Datatype>& datatype = datatype::inferDatatype<T>()) {
    return createElementwiseShaped<T>(
        sourceCount, targetCount, makeDimensions(shape, true), handler, datatype);
  }

  //! @brief As createElementwise, with the shape stated rather than assumed.
  template <typename T, typename F>
  static std::shared_ptr<GeneratedBuffer> createElementwiseShaped(
      std::size_t sourceCount,
      std::size_t targetCount,
      const std::vector<Dimension>& dimensions,
      const F& handler,
      const std::shared_ptr<datatype::Datatype>& datatype = datatype::inferDatatype<T>()) {
    std::size_t localTargetStride = targetCount;
    for (const auto& dimension : dimensions) {
      if (!dimension.isDistributed() && !dimension.isAppended()) {
        localTargetStride *= dimension.size;
      }
    }
    return std::make_shared<GeneratedBuffer>(
        sourceCount,
        targetCount,
        [=](void* targetPtr) {
          T* target = reinterpret_cast<T*>(targetPtr);

#pragma omp parallel for schedule(static)
          for (std::size_t i = 0; i < sourceCount; ++i) {
            std::invoke(handler, &target[i * localTargetStride], i);
          }
        },
        datatype,
        dimensions);
  }

  private:
  std::function<void(void*)> generator_;
  std::size_t sourceCount_;
  std::size_t targetStride_;
};

} // namespace seissol::io::writer

#endif // SEISSOL_SRC_IO_WRITER_INSTRUCTIONS_DATA_H_
