#ifndef SEISSOL_PARAMETER_READER_H
#define SEISSOL_PARAMETER_READER_H

#include <string>
#include <unordered_set>

#include <utils/stringutils.h>
#include <utils/logger.h>
#include <yaml-cpp/yaml.h>

#include "Parallel/MPI.h"

namespace seissol::initializers::parameters {

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
  T readOrFail(const std::string& field, const std::string& failMessage) {
    if (hasField(field)) {
      return readUnsafe<T>(field);
    } else {
      logError() << "The field" << field << "was not found, but it is required.";
      return T(); // unreachable. TODO(David): use compiler hint instead
    }
  }

  void warnDeprecatedSingle(const std::string& field) {
    if (hasField(field)) {
      visited.emplace(field);
      logInfo(seissol::MPI::mpi.rank())
          << "The field" << field
          << "is no longer in use. You may safely remove it from your parameters file.";
    }
  }

  void warnDeprecated(const std::vector<std::string>& fields) {
    for (const auto& field : fields) {
      warnDeprecatedSingle(field);
    }
  }

  void warnUnknown() {
    for (const auto& subnodes : node) {
      auto field = subnodes.first.as<std::string>();
      if (visited.find(field) == visited.end()) {
        logWarning(seissol::MPI::mpi.rank()) << "The field" << field << "is not known to SeisSol.";
      }
    }
  }

  template <typename... Args>
  void markUnused(const Args&... argFields) {
    for (const auto& field : {argFields...}) {
      logDebug(seissol::MPI::mpi.rank()) << "The field" << field << "is ignored (if it is found).";
      visited.emplace(field);
    }
  }

  ParameterReader readSubNode(const std::string& subnodeName) {
    visited.emplace(subnodeName);
    logDebug(seissol::MPI::mpi.rank()) << "Entering section" << subnodeName;
    if (hasField(subnodeName)) {
      return ParameterReader(node[subnodeName], false);
    } else {
      logDebug(seissol::MPI::mpi.rank())
          << "Section" << subnodeName
          << "not found in the given parameter file. Using an empty reader.";
      return ParameterReader(node[subnodeName], true);
    }
  }

  bool hasField(const std::string& field) { return !empty && node[field]; }

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

  bool empty;
  YAML::Node node; // apparently the YAML nodes use a reference semantic. Hence, we do it like this.
  std::unordered_set<std::string> visited;
};
} // namespace seissol::initializers::parameters

#endif
