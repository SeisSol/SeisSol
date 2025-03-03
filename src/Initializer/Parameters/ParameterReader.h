// SPDX-FileCopyrightText: 2023-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_PARAMETERS_PARAMETERREADER_H_
#define SEISSOL_SRC_INITIALIZER_PARAMETERS_PARAMETERREADER_H_

#include <optional>
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
  ParameterReader(const YAML::Node& node, const std::string& rootPath, bool empty);

  template <typename T>
  std::optional<T> read(const std::string& field) {
    if (hasField(field)) {
      return readUnsafe<T>(field);
    } else {
      return std::optional<T>();
    }
  }

  template <typename T>
  T readWithDefault(const std::string& field, const T& defaultValue) {
    const auto value = read<T>(field);
    if (value.has_value()) {
      return value.value();
    } else {
      logDebug() << "The field" << field << "was not specified, using fallback.";
      return defaultValue;
    }
  }

  // TODO(David): long-term (if we don't switch to another format first), merge readWithDefaultEnum
  // with readWithDefaultStringEnum, i.e. allow both numerical and textual values for an enum (can
  // we maybe auto-generate a parser from an enum definition?)
  template <typename T>
  T readWithDefaultEnum(const std::string& field,
                        const T& defaultValue,
                        const std::unordered_set<T>& validValues) {
    const int value = readWithDefault(field, static_cast<int>(defaultValue));
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
      const std::string failMessage = "The field " + field + " is required.";
      return readOrFail<T>(field, failMessage);
    } else {
      markUnused({field});
      return T{};
    }
  }

  template <typename T>
  T readIfRequiredAlternatives(const std::vector<std::string>& fields, bool required) {
    if (required) {
      bool found = false;
      T value = T{};
      for (std::size_t i = 0; i < fields.size(); ++i) {
        if (found) {
          if (i > 0) {
            warnDeprecatedSingle(fields[i]);
          } else {
            markUnused({fields[0]});
          }
        } else {
          const auto tryRead = read<T>(fields[i]);
          if (tryRead.has_value()) {
            found = true;
            value = tryRead.value();
            if (i > 0) {
              logWarning() << "The field name" << fields[i]
                           << "is deprecated; consider changing its name to" << fields[0];
            }
          }
        }
      }
      if (!found) {
        logError() << "The field" << fields[0] << "was not found, but it is required";
      }
      return value;
    } else {
      markUnused(fields);
      return T{};
    }
  }

  template <typename T>
  T readOrFail(const std::string& field, const std::string& failMessage) {
    if (hasField(field)) {
      return readUnsafe<T>(field);
    } else {
      logError() << "The field" << field
                 << "was not found, but it is required (details: " << failMessage.c_str() << ")";
      return T(); // unreachable. TODO(David): use compiler hint instead
    }
  }

  std::optional<std::string> readPath(const std::string& field);
  std::string readPathOrFail(const std::string& field, const std::string& failMessage);

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
    logDebug() << "The field" << field << "was read.";
    try {
      // booleans are stored as integers
      if constexpr (std::is_same_v<T, bool>) {
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
  std::string rootPath;
  bool empty;
  std::unordered_set<std::string> visited;
  std::unordered_map<std::string, std::shared_ptr<ParameterReader>> subreaders;
};
} // namespace seissol::initializer::parameters

#endif // SEISSOL_SRC_INITIALIZER_PARAMETERS_PARAMETERREADER_H_
