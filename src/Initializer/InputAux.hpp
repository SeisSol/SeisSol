#ifndef INITIALIZER_INPUTAUX_H_
#define INITIALIZER_INPUTAUX_H_

#include <type_traits>
#include <string>
#include <fstream>
#include <yaml-cpp/yaml.h>

#include "DynamicRupture/Typedefs.hpp"
#include "utils/logger.h"

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
} // namespace YAML
namespace seissol::initializers {
/*
   * If param stores a node with name field, return. Otherwise error.
   * @param param: YAML Node, which we want to read from.
   * @param field: Name of the field, we would like to read
 */
template <typename T>
T getUnsafe(const YAML::Node& param, const std::string& field) {
    try {
      // booleans are stored as integers
      if constexpr(std::is_same<T, bool>::value) {
        return param[field].as<int>() > 0;
      } else {
        return param[field].as<T>();
      }
    } catch(std::exception& e) {
      logError() << "Error while reading field " << field << ": " << e.what();
    }
}

 /*
  * If param stores a node with name field override value
  * @param param: YAML Node, which we want to read from.
  * @param field: Name of the field, we would like to read
  * @param value: Reference to the value, which we want to override
 */
template <typename T>
T getWithDefault(const YAML::Node& param, const std::string& field, T defaultValue) {
  T value = defaultValue;
  if (param[field]) {
    value = getUnsafe<T>(param, field);
  }
  return value;
}

/**
 * \brief Returns true if number elements in the input string (separated by the white space)
 *  is less or equal to the size of a container
 * */
template<typename OutputType, typename ContainerT>
bool isCapacityEnough(const std::string& inputString, ContainerT& outputMask) {
  std::istringstream inputStream(inputString);
  auto begin = std::istream_iterator<OutputType>(inputStream);
  auto end = std::istream_iterator<OutputType>();

  const size_t numInputElements = std::distance(begin, end);
  return numInputElements <= outputMask.size();
}

/**
 * \brief Initializes the given input mask with values from the input string
 *
 * \throws runtime_error if an input string contains more parameters than the capacity of a provided container
 * */
template<typename ContainerT>
void convertStringToMask(const std::string& stringMask, ContainerT& mask) {
  using T = typename std::iterator_traits<typename ContainerT::iterator>::value_type;

  if (!isCapacityEnough<T>(stringMask, mask))
    throw std::runtime_error("Number of input elements is more than the mask capacity");

  std::istringstream inputStream(stringMask);
  auto it = std::istream_iterator<T>(inputStream);
  auto end = std::istream_iterator<T>();

  for (int index = 0; it != end; ++index, ++it) {
    if (std::is_same<T, bool>::value) {
      mask[index] = (*it) > 0;
    }
    else {
      mask[index] = (*it);
    }
  }
}

/**
 * \brief Converts an input string to a vector of the datatype T. If the input string does not contain enough
 * items to convert to an array of length n, we fill with default values. If the input string contains
 * more than n items, we only consider the first n ones.
 *
 * \throws runtime_error if the input string contains parameters, which can not be converted to T.
 * */
template<typename T, size_t n>
std::array<T, n> convertStringToArray(const std::string& inputString) {
  auto result = std::array<T, n>();
  if (inputString.empty()) {
    return result;
  }

  size_t begin = 0;
  size_t wordCount = 0;
  enum class State { Word, Delimiter };
  State s = inputString.at(0) == ' ' ? State::Delimiter : State::Word;

  auto convert = [&inputString] (size_t begin, size_t end) {
    size_t count = end - begin;
    std::string word = inputString.substr(begin, count);
    if constexpr (std::is_integral<T>::value) {
      return std::stoi(word);
    } else if constexpr (std::is_floating_point<T>::value) {
      return std::stod(word);
    } else {
      return static_cast<T>(word);
    }
  };

  for (size_t i = 0; i < inputString.size(); i++) {
    if (inputString.at(i) == ' ') {
      if (s == State::Word) {
        result.at(wordCount) = convert(begin, i);
        wordCount++;
        if (wordCount >= n) {
          break;
        }
      }
      begin = i;
      s = State::Delimiter;
    } else {
      s = State::Word;
    }
  }
  if (wordCount < n) {
    result.at(wordCount) = convert(begin, inputString.size());
  }

  return result;
}

using StringsType = std::list<std::string>;
class FileProcessor {
public:
  static StringsType getFileAsStrings(const std::string & fileName) {
    StringsType content;
    std::fstream paramFile(fileName, std::ios_base::in);
    if (!paramFile.is_open()) {
      throw std::runtime_error("cannot open file: " + fileName);
    }

    std::string tmpString;
    while (std::getline(paramFile, tmpString)) {
      content.push_back(tmpString);
    }

    paramFile.close();
    return content;
  }

  static void removeEmptyLines(StringsType& content) {
    const std::string WHITESPACE = " \n\r\t\f\v";
    auto isEmptyString = [&WHITESPACE](const std::string & string) -> bool {
      size_t start = string.find_first_not_of(WHITESPACE);
      return start == std::string::npos;
    };

    std::vector<StringsType::iterator> deletees;
    for (auto itr = content.begin(); itr != content.end(); ++itr) {
      if (isEmptyString(*itr))
        deletees.push_back(itr);
    }

    for (auto &itr : deletees) {
      content.erase(itr);
    }
  }
};
} // seissol::initializers

#endif //INITIALIZER_INPUTAUX_H_
