#pragma once

#include <IO/Datatype/Datatype.hpp>
#include <IO/Datatype/Inference.hpp>
#include <hdf5.h>
#include <memory>
#include <stack>
#include <string>
#include <vector>

namespace seissol::io::reader::file {
class Hdf5Reader {
  public:
  Hdf5Reader(MPI_Comm comm);
  void openFile(const std::string& name);
  void openGroup(const std::string& name);
  template <typename T>
  std::vector<T>
      readAttribute(const std::string& name,
                    std::shared_ptr<datatype::Datatype> type = datatype::inferDatatype<T>()) {
    const auto count = attributeCount(name);
    std::vector<T> output(count);
    readAttributeRaw(output.data(), name, type);
    return output;
  }
  template <typename T>
  T readAttributeScalar(const std::string& name,
                        std::shared_ptr<datatype::Datatype> type = datatype::inferDatatype<T>()) {
    T attr;
    readAttributeRaw(&attr, name, type);
    return attr;
  }
  std::size_t attributeCount(const std::string& name);
  void readAttributeRaw(void* data,
                        const std::string& name,
                        std::shared_ptr<datatype::Datatype> type);
  template <typename T>
  std::vector<T>
      readData(const std::string& name,
               std::shared_ptr<datatype::Datatype> targetType = datatype::inferDatatype<T>()) {
    const auto count = dataCount(name);
    std::vector<T> output(count);
    readDataRaw(output.data(), name, count, targetType);
    return output;
  }
  std::size_t dataCount(const std::string& name);
  void readDataRaw(void* data,
                   const std::string& name,
                   std::size_t count,
                   std::shared_ptr<datatype::Datatype> targetType);
  void closeGroup();
  void closeFile();

  void checkExistence(const std::string& name, const std::string& type);

  private:
  std::stack<hid_t> handles;
  MPI_Comm comm;
};
} // namespace seissol::io::reader::file
