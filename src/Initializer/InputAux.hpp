#ifndef INITIALIZER_INPUTAUX_H_
#define INITIALIZER_INPUTAUX_H_

#include <type_traits>
#include <string>
#include <fstream>
#include <yaml-cpp/yaml.h>

namespace seissol::initializers {
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
/**
 * \brief Returns true if number elements in the input string (separated by the white space)
 *  is less or equal to the size of a container
 * */
template<typename OutputType, typename ContainerT>
bool isCapacityEnough(const std::string& InputString, ContainerT& OutputMask) {
  std::istringstream InputStream(InputString);
  auto Begin = std::istream_iterator<OutputType>(InputStream);
  auto End = std::istream_iterator<OutputType>();

  const size_t NumInputElements = std::distance(Begin, End);
  return NumInputElements <= OutputMask.size();
}

/**
 * \brief Initializes the given input mask with values from the input string
 *
 * \throws runtime_error if an input string contains more parameters than the capacity of a provided container
 * */
template<typename ContainerT>
void convertStringToMask(const std::string& StringMask, ContainerT& Mask) {
  using T = typename std::iterator_traits<typename ContainerT::iterator>::value_type;

  if (!isCapacityEnough<T>(StringMask, Mask))
    throw std::runtime_error("Num. input elements is more than the Mask capacity");

  std::istringstream InputStream(StringMask);
  auto Begin = std::istream_iterator<T>(InputStream);
  auto End = std::istream_iterator<T>();

  for (int Index = 0; Begin != End; ++Index, ++Begin) {
    if (std::is_same<T, bool>::value) {
      Mask[Index] = (*Begin) > 0;
    }
    else {
      Mask[Index] = (*Begin);
    }
  }
}

using StringsT = std::list<std::string>;
class FileProcessor {
public:
  static StringsT getFileAsStrings(const std::string &FileName) {
    StringsT Content;
    std::fstream ParamFile(FileName, std::ios_base::in);
    if (!ParamFile.is_open()) {
      throw std::runtime_error("cannot open file: " + FileName);
    }

    std::string TmpString;
    while (std::getline(ParamFile, TmpString)) {
      Content.push_back(TmpString);
    }

    ParamFile.close();
    return Content;
  }

  static void removeEmptyLines(StringsT &Content) {

    const std::string WHITESPACE = " \n\r\t\f\v";
    auto isEmptyString = [&WHITESPACE](const std::string &String) -> bool {
      size_t Start = String.find_first_not_of(WHITESPACE);
      return Start == std::string::npos;
    };

    std::vector<StringsT::iterator> Deletees;
    for (auto Itr = Content.begin(); Itr != Content.end(); ++Itr) {
      if (isEmptyString(*Itr))
        Deletees.push_back(Itr);
    }

    for (auto &Itr : Deletees) {
      Content.erase(Itr);
    }
  }
};
} // seissol::initializers

#endif //INITIALIZER_INPUTAUX_H_