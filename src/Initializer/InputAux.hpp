#ifndef INITIALIZER_INPUTAUX_H_
#define INITIALIZER_INPUTAUX_H_

#include <type_traits>
#include <string>
#include <yaml-cpp/yaml.h>

namespace seissol {
  namespace initializers {
    template <typename T>
    T getParamIfExists(const YAML::Node& Param, std::string&& Field, T DefaultValue) {
      if (std::is_same<T, bool>::value) {
        T Value{DefaultValue};
        if (Param[Field]) {
          Value = Param[Field].as<int>() > 0;
        }
        return Value;
      }
      else {
        return Param[Field] ? Param[Field].as<T>() : DefaultValue;
      }
    }
  }
}

#endif //INITIALIZER_INPUTAUX_H_