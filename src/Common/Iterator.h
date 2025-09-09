// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_ITERATOR_H_
#define SEISSOL_SRC_COMMON_ITERATOR_H_

#include <functional>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

// NOTE: the following are essentially backports of C++20 and C++23 features

namespace seissol::common {

// TODO: remove all in here once C++23 lands in SeisSol

// cf. https://stackoverflow.com/a/44661987
template <typename RangeT>
using IteratorType = decltype(std::begin(std::declval<RangeT&>()));

/**
  Runs multiple iterators at the same time in a tuple, akin to Python `zip`.

  Can most likely be simplified with C++20 view transforms.
  Will be replaced with C++23 std::views::zip .
 */
template <typename... RangeTs>
class Zip {
  public:
  // only a discount implementation for now: we only support LegacyInputIterator.

  template <typename... IteratorTs>
  class Iterator {
public:
    // NOLINTNEXTLINE
    using iterator_category = std::input_iterator_tag;
    // NOLINTNEXTLINE
    using difference_type =
        std::tuple<typename std::iterator_traits<IteratorTs>::difference_type...>;
    // NOLINTNEXTLINE
    using value_type = std::tuple<typename std::iterator_traits<IteratorTs>::value_type...>;
    // NOLINTNEXTLINE
    using pointer = std::tuple<typename std::iterator_traits<IteratorTs>::pointer...>;
    // NOLINTNEXTLINE
    using reference = std::tuple<typename std::iterator_traits<IteratorTs>::reference...>;

    constexpr static auto makeReference(std::tuple<IteratorTs...>& iterators) -> reference {
      return std::apply([](auto&... values) { return reference{*values...}; }, iterators);
    }

    Iterator(bool lenient,
             std::tuple<IteratorTs...> iterators,
             std::tuple<IteratorTs...> iteratorEnds)
        : iterators(iterators), iteratorEnds(iteratorEnds),
          ended(tupleAnyEqual(iterators, iteratorEnds)), lenient(lenient) {
      if (!lenient && ended && iterators != iteratorEnds) {
        throw std::runtime_error("Not all iterators are finished (but they should be).");
      }
    }

    constexpr auto operator*() -> reference { return makeReference(iterators); }

    constexpr auto operator*() const -> reference { return makeReference(iterators); }

    constexpr auto operator++() -> Iterator& {
      iterators = tupleTransform([](auto& value) { return ++value; }, iterators);

      ended = tupleAnyEqual(iterators, iteratorEnds);
      if (!lenient && ended && iterators != iteratorEnds) {
        throw std::runtime_error("Not all iterators are finished (but they should be).");
      }

      return *this;
    }

    constexpr auto operator==(const Iterator& other) const -> bool {
      return (ended && other.ended) || (!ended && !other.ended && iterators == other.iterators);
    }

    constexpr auto operator!=(const Iterator& other) const -> bool { return !(*this == other); }

private:
    std::tuple<IteratorTs...> iterators;
    std::tuple<IteratorTs...> iteratorEnds;
    bool ended;
    bool lenient;
  };

  Zip(bool lenient, RangeTs&&... ranges) : lenient(lenient), ranges(ranges...) {}

  constexpr auto begin() {
    return Iterator<IteratorType<RangeTs>...>(
        lenient,
        tupleTransform([](auto& value) { return std::begin(value); }, ranges),
        tupleTransform([](auto&& value) { return std::end(value); }, ranges));
  }

  constexpr auto end() {
    return Iterator<IteratorType<RangeTs>...>(
        lenient,
        tupleTransform([](auto& value) { return std::end(value); }, ranges),
        tupleTransform([](auto&& value) { return std::end(value); }, ranges));
  }

  constexpr auto begin() const {
    return Iterator<IteratorType<RangeTs>...>(
        lenient,
        tupleTransform([](const auto& value) { return std::cbegin(value); }, ranges),
        tupleTransform([](const auto& value) { return std::cend(value); }, ranges));
  }

  constexpr auto end() const {
    return Iterator<IteratorType<RangeTs>...>(
        lenient,
        tupleTransform([](const auto& value) { return std::cend(value); }, ranges),
        tupleTransform([](const auto& value) { return std::cend(value); }, ranges));
  }

  private:
  template <typename TupleT, std::size_t... Idx>
  constexpr static bool
      tupleAnyEqualImpl(const TupleT& tuple1, const TupleT& tuple2, std::index_sequence<Idx...>) {
    return ((std::get<Idx>(tuple1) == std::get<Idx>(tuple2)) || ...);
  }

  template <typename TupleT>
  constexpr static bool tupleAnyEqual(const TupleT& tuple1, const TupleT& tuple2) {
    using Indices = std::make_index_sequence<std::tuple_size_v<TupleT>>;
    return tupleAnyEqualImpl(tuple1, tuple2, Indices());
  }

  template <typename F, typename TupleA>
  constexpr static auto tupleTransform(F function, TupleA& tuple) {
    return std::apply([function](auto&... args) { return std::tuple{function(args)...}; }, tuple);
  }

  template <typename F, typename TupleA>
  constexpr static auto tupleTransform(F function, const TupleA& tuple) {
    return std::apply([function](const auto&... args) { return std::tuple{function(args)...}; },
                      tuple);
  }

  std::tuple<RangeTs...> ranges;
  bool lenient;
};

/**
  Runs multiple iterators at the same time in a tuple, akin to Python `zip`.

  Usage:

  for (const auto [first, second, third] : zip(firstIt, secondIt, thirdIt)) {
    // ...
  }

  Can most likely be simplified with C++20 view transforms.
  Will be replaced with C++23 std::views::zip .
 */
template <typename... RangeTs>
constexpr auto zip(RangeTs&&... ranges) {
  return Zip<RangeTs...>(false, std::forward<RangeTs>(ranges)...);
}

/**
  A simple incrementing range.

  TODO: replace/remove, once C++20 lands in SeisSol
 */
template <typename T>
class Range {
  static_assert(std::is_integral_v<T>, "For now, T needs to be integer");

  public:
  class Iterator {
public:
    // NOLINTNEXTLINE
    using iterator_category = std::input_iterator_tag;
    // NOLINTNEXTLINE
    using difference_type = std::make_signed_t<T>;
    // NOLINTNEXTLINE
    using value_type = T;
    // NOLINTNEXTLINE
    using pointer = const T*;
    // NOLINTNEXTLINE
    using reference =
        std::conditional_t<std::is_trivially_copyable_v<T> && sizeof(T) <= sizeof(std::size_t),
                           T,
                           const T&>;

    Iterator(T value, T step, T end) : value(value), step(step), end(end) {
      if (value > end) {
        value = end;
      }
    }

    constexpr auto operator++() {
      value += step;
      if (value > end) {
        value = end;
      }
      return *this;
    }

    constexpr auto operator*() -> reference { return value; }

    constexpr auto operator*() const -> reference { return value; }

    constexpr auto operator==(const Iterator& other) const -> bool {
      return value == other.value && step == other.step && end == other.end;
    }

    constexpr auto operator!=(const Iterator& other) const -> bool { return !(*this == other); }

private:
    T value;
    T step;
    T end;
  };

  Range(T start, T stop, T step) : startVal(start), stopVal(stop), stepVal(step) {}

  [[nodiscard]] constexpr auto begin() const { return Iterator(startVal, stepVal, stopVal); }

  [[nodiscard]] constexpr auto end() const { return Iterator(stopVal, stepVal, stopVal); }

  [[nodiscard]] constexpr auto cbegin() const { return Iterator(startVal, stepVal, stopVal); }

  [[nodiscard]] constexpr auto cend() const { return Iterator(stopVal, stepVal, stopVal); }

  private:
  T startVal;
  T stopVal;
  T stepVal;
};

/**
  A range to an iterator, running from 0 to `stop - 1`. (similar to Python)

  With C++20, replace by std::views::iota
 */
template <typename T>
constexpr auto range(T stop) {
  return Range(0, stop, 1);
}

/**
  A range to an iterator, running from `start` to `stop - 1`. (similar to Python)

  With C++20, replace by std::views::iota
 */
template <typename T>
constexpr auto range(T start, T stop) {
  return Range(start, stop, 1);
}

/**
  A range to an iterator, running from `start` to `stop - 1` in increment `step`. (similar to
  Python)

  With C++20, replace by std::views::iota
 */
template <typename T>
constexpr auto range(T start, T stop, T step) {
  return Range(start, stop, step);
}

/**
  Counts while running through an iterator.

  Usage:

  for (const auto [i, item] : zip(it)) {
    // ...
  }

  Will give:
  0, it[0]
  1, it[1]
  ...

  Can most likely be simplified with C++20 view transforms.
  Will be replaced with C++23 std::views::enumerate .
 */
template <typename RangeT>
constexpr auto enumerate(RangeT&& iterator) {
  // a tiny bit hacky to use both zip and the int range like that. But it should work.
  return Zip<Range<std::size_t>, RangeT>(
      true,
      Range<std::size_t>(0, std::numeric_limits<std::size_t>::max(), 1),
      std::forward<RangeT>(iterator));
}

/**
  Filter an iterator by a function; i.e. skip certain elements while iterating over it. E.g.

  auto beginIt = FilteredIterator(obj.begin(), obj.end(), filter);
  auto endIt = FilteredIterator(obj.end(), obj.end(), filter);

  for (auto it = beginIt; it != endIt; ++it) {
    // ...
  }

  is equivalent to

  for (auto it = obj.begin(); it != obj.end(); ++it) {
    if (filter(*it)) {
      // ...
    }
  }

  (NOTE: might benefit from a short-hand notation like for zip and enumerate in this file; however
  we didn't need that so far; and C++20 is around the corner to be adopted anyways)
 */
template <typename T>
class FilteredIterator {
  public:
  // NOLINTNEXTLINE
  using iterator_category = typename std::iterator_traits<T>::iterator_category;
  // NOLINTNEXTLINE
  using difference_type = typename std::iterator_traits<T>::difference_type;
  // NOLINTNEXTLINE
  using value_type = typename std::iterator_traits<T>::value_type;
  // NOLINTNEXTLINE
  using pointer = typename std::iterator_traits<T>::pointer;
  // NOLINTNEXTLINE
  using reference = typename std::iterator_traits<T>::reference;

  FilteredIterator(T base, T end, std::function<bool(reference)> filter)
      : base(base), end(end), filter(filter) {
    // skip initially-filtered elements
    while (this->base != end && !filter(*this->base)) {
      ++this->base;
    }
  }

  constexpr auto operator++() {
    // advance, and skip if needed
    ++base;
    while (base != end && !filter(*base)) {
      ++base;
    }
    return *this;
  }

  constexpr auto operator*() -> reference { return *base; }

  constexpr auto operator*() const -> reference { return *base; }

  constexpr auto operator==(const T& other) const -> bool { return other == base; }

  constexpr auto operator!=(const T& other) const -> bool { return !(*this == other); }

  constexpr auto operator==(const FilteredIterator<T>& other) const -> bool {
    return other.base == base;
  }

  constexpr auto operator!=(const FilteredIterator<T>& other) const -> bool {
    return !(*this == other);
  }

  private:
  T base;
  T end;
  std::function<bool(reference)> filter;
};

} // namespace seissol::common

#endif // SEISSOL_SRC_COMMON_ITERATOR_H_
