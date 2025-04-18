// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_COMMON_ITERATOR_H_
#define SEISSOL_SRC_COMMON_ITERATOR_H_

#include <iterator>
#include <limits>
#include <stdexcept>
#include <tuple>
#include <utility>

// NOTE: the following are essentially backports of C++20 and C++23 features

namespace seissol::common {

// TODO: remove once C++23 lands in SeisSol

// cf. https://stackoverflow.com/a/44661987
template <typename RangeT>
using IteratorType = decltype(std::begin(std::declval<RangeT&>()));

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

template <typename... RangeTs>
constexpr auto zip(RangeTs&&... ranges) {
  return Zip<RangeTs...>(false, std::forward<RangeTs>(ranges)...);
}

// TODO: replace/remove, once C++20 lands in SeisSol

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

template <typename T>
constexpr auto range(T stop) {
  return Range(0, stop, 1);
}

template <typename T>
constexpr auto range(T start, T stop) {
  return Range(start, stop, 1);
}

template <typename T>
constexpr auto range(T start, T stop, T step) {
  return Range(start, stop, step);
}

template <typename RangeT>
constexpr auto enumerate(RangeT&& iterator) {
  // a tiny bit hacky to use both zip and the int range like that. But it should work.
  return Zip<Range<std::size_t>, RangeT>(
      true,
      Range<std::size_t>(0, std::numeric_limits<std::size_t>::max(), 1),
      std::forward<RangeT>(iterator));
}

} // namespace seissol::common

#endif // SEISSOL_SRC_COMMON_ITERATOR_H_
