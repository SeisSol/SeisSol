#ifndef INITIALIZER_INPUTAUX_H_
#define INITIALIZER_INPUTAUX_H_

#include <type_traits>
#include <string>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Typedefs.hpp"

namespace YAML {
template <>
struct convert<seissol::dr::FrictionLawType> {
  static Node encode(const seissol::dr::FrictionLawType& rhs) {
    Node node;
    node.push_back(static_cast<unsigned int>(rhs));
    return node;
  }

  static bool decode(const Node& node, seissol::dr::FrictionLawType& rhs) {
    if (node.IsSequence() || node.size() != 0) {
      return false;
    }

    rhs = static_cast<seissol::dr::FrictionLawType>(node.as<unsigned int>());
    return true;
  }
};
}; // namespace YAML
namespace seissol::initializers {
  template <typename T>
  void updateIfExists(const YAML::Node& param, std::string&& field, T& value) {
    // if params stores a node with name field override value
    if (param[field]) {
      // booleans are stored as integers
      if constexpr(std::is_same<T, bool>::value) {
        value = param[field].as<int>() > 0;
      } else {
        value = param[field].as<T>();
      }
    }
  }
}

#endif //INITIALIZER_INPUTAUX_H_