#ifndef SEISSOL_PARAMETER_READER_H
#define SEISSOL_PARAMETER_READER_H

#include <string>
#include <unordered_set>

#include <utils/logger.h>
#include <utils/stringutils.h>
#include <yaml-cpp/yaml.h>

#include "Parallel/MPI.h"

namespace seissol::initializer::parameters {

// converts a string to lower case, and trims it.
inline void sanitize(std::string& input) {
  utils::StringUtils::trim(input);
  utils::StringUtils::toLower(input);
}

// A small helper class which reads a YAML node dictionary. It keeps track of all items that have
// been read and reports all values which are not used or not used anymore.
// TODO(David): maybe make the reader more tree-like (i.e. keep a central set on which nodes have
// been visited), and output all non-understood values at the end and not between sections
class ParameterReader {
  public:
  ParameterReader(const YAML::Node& node, bool empty) : node(node), empty(empty) {}

  template <typename T>
  T readWithDefault(const std::string& field, const T& defaultValue) {
    T value = defaultValue;
    if (hasField(field)) {
      value = readUnsafe<T>(field);
    } else {
      logDebug(seissol::MPI::mpi.rank())
          << "The field" << field << "was not specified, using fallback.";
    }
    return value;
  }

  // TODO(David): long-term (if we don't switch to another format first), merge readWithDefaultEnum
  // with readWithDefaultStringEnum, i.e. allow both numerical and textual values for an enum (can
  // we maybe auto-generate a parser from an enum definition?)
  template <typename T>
  T readWithDefaultEnum(const std::string& field,
                        const T& defaultValue,
                        const std::unordered_set<T>& validValues) {
    int value = readWithDefault(field, static_cast<int>(defaultValue));
    if (validValues.find(static_cast<T>(value)) == validValues.end()) {
      logError() << "The field" << field << "had an invalid enum value:" << value;
    }
    return static_cast<T>(value);
  }

  template <typename T>
  T readWithDefaultStringEnum(const std::string& field,
                              const std::string& defaultValue,
                              const std::unordered_map<std::string, T>& validValues) {
    std::string value = readWithDefault(field, defaultValue); // TODO(David): sanitize string
    sanitize(value);
    if (validValues.find(value) == validValues.end()) {
      logError() << "The field" << field << "had an invalid enum value:" << value;
    }
    return validValues.at(value);
  }

  template <typename T>
  T readIfRequired(const std::string& field, bool required) {
    if (required) {
      const std::string failMessage =
          "The field " + field + " is required, but not found in the parameters file.";
      return readOrFail<T>(field, failMessage);
    } else {
      markUnused({field});
      return T{};
    }
  }

  template <typename T>
  T readOrFail(const std::string& field, const std::string& failMessage) {
    if (hasField(field)) {
      return readUnsafe<T>(field);
    } else {
      logError() << "The field" << field << "was not found, but it is required.";
      return T(); // unreachable. TODO(David): use compiler hint instead
    }
  }

  void warnDeprecatedSingle(const std::string& field);
  void warnDeprecated(const std::vector<std::string>& fields);
  void warnUnknown(const std::string& prefix = "") const;
  void markUnused(const std::vector<std::string>& fields);
  ParameterReader* readSubNode(const std::string& subnodeName);
  bool hasField(const std::string& field);

  private:
  template <typename T>
  T readUnsafe(const std::string& field) {
    visited.emplace(field);
    logDebug(seissol::MPI::mpi.rank()) << "The field" << field << "was read.";
    try {
      // booleans are stored as integers
      if constexpr (std::is_same<T, bool>::value) {
        return node[field].as<int>() > 0;
      } else {
        return node[field].as<T>();
      }
    } catch (std::exception& e) {
      logError() << "Error while reading field" << field << ":" << e.what();
      return T(); // unreachable. TODO(David): use compiler hint instead
    }
  }

  YAML::Node node; // apparently the YAML nodes use a reference semantic. Hence, we do it like this.
  bool empty;
  std::unordered_set<std::string> visited;
  std::unordered_map<std::string, std::shared_ptr<ParameterReader>> subreaders;
};
} // namespace seissol::initializer::parameters

#endif
