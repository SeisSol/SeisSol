// SPDX-FileCopyrightText: 2026 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_JSON_H_
#define SEISSOL_SRC_COMMON_JSON_H_

#include <cstddef>
#include <cstdint>
#include <iosfwd>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <variant>
#include <vector>

namespace seissol::json {

/**
 * @brief A JSON value, built up in memory and written out as text.
 *
 * Meant for the files SeisSol writes for other programs to read: an object keeps its keys in the
 * order they were added, so that what a person reads follows what the code wrote; numbers that
 * are not finite, which JSON has no spelling for, are written as null.
 */
class Value {
  public:
  using Array = std::vector<Value>;
  using Object = std::vector<std::pair<std::string, Value>>;

  //! @brief null
  Value() = default;
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(std::nullptr_t /*null*/) {}
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(bool value) : data_(value) {}

  template <typename T>
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(T value)
    requires(std::is_integral_v<T> && !std::is_same_v<T, bool>)
  {
    if constexpr (std::is_signed_v<T>) {
      data_ = static_cast<std::int64_t>(value);
    } else {
      data_ = static_cast<std::uint64_t>(value);
    }
  }

  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(double value) : data_(value) {}
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(float value) : data_(static_cast<double>(value)) {}
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(const char* value) : data_(std::string(value)) {}
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(std::string value) : data_(std::move(value)) {}
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(std::string_view value) : data_(std::string(value)) {}
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(Array value) : data_(std::move(value)) {}
  // NOLINTNEXTLINE(misc-explicit-constructor)
  Value(Object value) : data_(std::move(value)) {}

  //! @brief An array of @p values, one element each.
  template <typename T>
  static Value arrayOf(const std::vector<T>& values) {
    Array array;
    array.reserve(values.size());
    for (const auto& value : values) {
      array.emplace_back(value);
    }
    return {std::move(array)};
  }

  static Value array() { return {Array{}}; }
  static Value object() { return {Object{}}; }

  [[nodiscard]] bool isNull() const { return std::holds_alternative<std::nullptr_t>(data_); }
  [[nodiscard]] bool isBool() const { return std::holds_alternative<bool>(data_); }
  [[nodiscard]] bool isNumber() const;
  [[nodiscard]] bool isString() const { return std::holds_alternative<std::string>(data_); }
  [[nodiscard]] bool isArray() const { return std::holds_alternative<Array>(data_); }
  [[nodiscard]] bool isObject() const { return std::holds_alternative<Object>(data_); }

  /**
   * @brief The member @p key of an object, added as null if there is none yet.
   *
   * A null value turns into an empty object first, so that nested objects can be filled in
   * without declaring each level. Anything else that is not an object is an error.
   */
  Value& operator[](std::string_view key);

  //! @brief The member @p key of an object, or nullptr if there is none.
  [[nodiscard]] const Value* find(std::string_view key) const;

  /**
   * @brief Appends @p value to an array and returns the element.
   *
   * A null value turns into an empty array first.
   */
  Value& push(Value value);

  //! @brief The element @p index of an array; fails for anything else, or past its end.
  Value& at(std::size_t index);

  //! @brief The number of elements of an array or of members of an object; zero otherwise.
  [[nodiscard]] std::size_t size() const;

  //! @brief The value of a string; fails for anything else.
  [[nodiscard]] const std::string& asString() const;

  //! @brief The value of a number, as a double; fails for anything else.
  [[nodiscard]] double asDouble() const;

  //! @brief The value of a boolean; fails for anything else.
  [[nodiscard]] bool asBool() const;

  [[nodiscard]] const Array& asArray() const;
  [[nodiscard]] const Object& asObject() const;

  /**
   * @brief Writes the value as JSON text.
   *
   * @param indent The number of spaces a nesting level is indented by, with one member or element
   * per line; a negative number writes everything on one line, without any whitespace.
   */
  void write(std::ostream& stream, int indent = 2) const;

  //! @brief The value as JSON text, as write() produces it.
  [[nodiscard]] std::string dump(int indent = 2) const;

  private:
  void write(std::ostream& stream, int indent, int depth) const;

  std::
      variant<std::nullptr_t, bool, std::int64_t, std::uint64_t, double, std::string, Array, Object>
          data_;
};

//! @brief @p text as a JSON string literal, quotes included.
std::string quote(std::string_view text);

/**
 * @brief The shortest decimal spelling of @p value that reads back as the same double.
 *
 * Not finite values come out as null, which is what JSON allows instead.
 */
std::string formatNumber(double value);

} // namespace seissol::json

#endif // SEISSOL_SRC_COMMON_JSON_H_
