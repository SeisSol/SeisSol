// SPDX-FileCopyrightText: 2020-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_INITIALIZER_INPUTAUX_H_
#define SEISSOL_SRC_INITIALIZER_INPUTAUX_H_

#include <fstream>
#include <iterator>
#include <list>
#include <sstream>
#include <utils/logger.h>
#include <vector>

namespace seissol::initializer {
/**
 * \brief Returns true if number elements in the input string (separated by the white space)
 *  is less or equal to the size of a container
 * */
template <typename OutputType, typename ContainerT>
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
 * \throws runtime_error if an input string contains more parameters than the capacity of a provided
 * container
 * */
template <typename ContainerT>
void convertStringToMask(const std::string& stringMask, ContainerT& mask) {
  using T = typename std::iterator_traits<typename ContainerT::iterator>::value_type;

  if (!isCapacityEnough<T>(stringMask, mask)) {
    throw std::runtime_error("Number of input elements is more than the mask capacity");
  }

  std::istringstream inputStream(stringMask);
  auto it = std::istream_iterator<T>(inputStream);
  auto end = std::istream_iterator<T>();

  for (int index = 0; it != end; ++index, ++it) {
    if (std::is_same<T, bool>::value) {
      mask[index] = (*it) > 0;
    } else {
      mask[index] = (*it);
    }
  }
}

/**
 * \brief Converts an input string to a vector of the datatype T and length n. Unless explicitly
 * ignored, this contains: we want exactly n elements, and empty elements are not considered (i.e.
 * two delimiters following on each other are treated as a single one).
 *
 * \throws runtime_error if the input string contains parameters, which can not be converted to T.
 * Or if the length is not equal to n (unless ignored).
 * */
template <typename T, size_t N>
std::array<T, N> convertStringToArray(const std::string& inputString,
                                      bool exactLength = true,
                                      bool skipEmpty = true,
                                      char delimiter = ' ') {
  auto result = std::array<T, N>();
  if (inputString.empty()) {
    if (exactLength && N > 0) {
      throw std::runtime_error(
          std::string("Insufficient number of elements in array. Given: 0. Required: ") +
          std::to_string(N) + std::string("."));
    } else {
      return result;
    }
  }

  auto convert = [&inputString](size_t begin, size_t end) {
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

  size_t begin = 0;
  size_t wordCount = 0;
  enum class State { Word, Delimiter };
  State s = inputString.at(0) == delimiter ? State::Delimiter : State::Word;

  // iterate over all words. We need to start at zero here: suppose we had ";;;" with delimiter ';'.
  // when !skipEmpty, this denotes an array with four elements; the first three are found with this
  // loop here.
  for (size_t i = 0; i < inputString.size(); i++) {
    if (inputString.at(i) == delimiter) {
      // either we have a word, or two subsequent delimiters
      if (s == State::Word || !skipEmpty) {
        result.at(wordCount) = convert(begin, i);
        ++wordCount;
        if (wordCount >= N) {
          break;
        }
      }

      // exclude the delimiter, hence i+1
      begin = i + 1;
      s = State::Delimiter;
    } else {
      s = State::Word;
    }
  }

  // handle rest. Note that if a line ends with a delimiter, we consider the last element to be an
  // empty one again.
  if ((s == State::Word || !skipEmpty) && wordCount < N) {
    result.at(wordCount) = convert(begin, inputString.size());
    ++wordCount;
  }

  if (wordCount != N && exactLength) {
    throw std::runtime_error(std::string("Insufficient number of elements in array. Given: ") +
                             std::to_string(wordCount) + std::string(". Required: ") +
                             std::to_string(N));
  }

  return result;
}
//
using StringsType = std::list<std::string>;
class FileProcessor {
  public:
  static StringsType getFileAsStrings(const std::string& fileName, const std::string& what) {
    StringsType content;
    std::fstream paramFile(fileName, std::ios_base::in);
    if (!paramFile.is_open()) {
      logError() << "Cannot open file (" << what.c_str() << "):" << fileName;
    }

    std::string tmpString;
    while (std::getline(paramFile, tmpString)) {
      content.push_back(tmpString);
    }

    paramFile.close();
    return content;
  }

  static void removeEmptyLines(StringsType& content) {
    const std::string whitespace = " \n\r\t\f\v";
    auto isEmptyString = [&whitespace](const std::string& string) -> bool {
      size_t start = string.find_first_not_of(whitespace);
      return start == std::string::npos;
    };

    std::vector<StringsType::iterator> deletees;
    for (auto itr = content.begin(); itr != content.end(); ++itr) {
      if (isEmptyString(*itr))
        deletees.push_back(itr);
    }

    for (auto& itr : deletees) {
      content.erase(itr);
    }
  }
};
} // namespace seissol::initializer

#endif // SEISSOL_SRC_INITIALIZER_INPUTAUX_H_
