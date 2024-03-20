#pragma once

#include "Datatype.hpp"
#include <limits>
#include <memory>
#include <type_traits>

namespace seissol::io::datatype {
template <typename S>
static std::vector<std::size_t> arrayTotalExtent() {
  if constexpr (std::is_array_v<S>) {
    static_assert(std::extent_v<S> != 0,
                  "No arrays of unspecified extent (T[]) are supported at the moment.");
    auto subExtent = arrayTotalExtent<std::remove_extent_t<S>>();
    subExtent.insert(subExtent.begin(), std::extent_v<S>);
    return subExtent;
  } else {
    return {};
  }
}

template <typename T, typename = void>
struct InferDatatypeStruct {
  const static inline std::shared_ptr<Datatype> Type = std::make_shared<OpaqueDatatype>(sizeof(T));
};

template <typename T>
struct InferDatatypeStruct<T, std::enable_if_t<std::is_array_v<T>>> {
  const static inline std::shared_ptr<Datatype> Type = std::make_shared<ArrayDatatype>(
      InferDatatypeStruct<std::remove_all_extents_t<T>>::Type, arrayTotalExtent<T>());
};

template <typename T, std::size_t N>
struct InferDatatypeStruct<std::array<T, N>> {
  const static inline std::shared_ptr<Datatype> Type =
      std::make_shared<ArrayDatatype>(InferDatatypeStruct<T>::Type, {N});
};

template <typename... Ts>
struct InferDatatypeStruct<std::tuple<Ts...>> {
  private:
  template <std::size_t I, typename T>
  static StructDatatype::MemberInfo makeMemberInfo() {
    std::tuple<Ts...> instance;
    return StructDatatype::MemberInfo{std::string("_") + std::to_string(I),
                                      reinterpret_cast<std::size_t>(&std::get<I>(instance)) -
                                          reinterpret_cast<std::size_t>(&instance),
                                      InferDatatypeStruct<T>::Type};
  }

  template <std::size_t I, typename... T>
  static std::vector<StructDatatype::MemberInfo> makeMemberInfos();

  template <std::size_t I, typename T, typename... TRest>
  static std::vector<StructDatatype::MemberInfo> makeMemberInfos() {
    auto prevVector = makeMemberInfos<I + 1, TRest...>();
    prevVector.insert(prevVector.begin(), makeMemberInfo<I, T>());
    return prevVector;
  }

  template <std::size_t I>
  static std::vector<StructDatatype::MemberInfo> makeMemberInfos() {
    return {};
  }

  public:
  const static inline std::shared_ptr<Datatype> Type =
      std::make_shared<StructDatatype>(makeMemberInfos<0, Ts...>());
};

/*template<typename T>
struct InferDatatypeStruct<T, std::enable_if_t<std::is_same_v<T,
decltype(std::declval<T>().layout())>>> { const static inline std::shared_ptr<Datatype> Type =
std::make_shared<StructDatatype>(T::layout());
};*/

template <typename T>
struct InferDatatypeStruct<T, std::enable_if_t<std::is_integral_v<T>>> {
  const static inline std::shared_ptr<Datatype> Type =
      std::make_shared<IntegerDatatype>(sizeof(T), std::numeric_limits<T>::is_signed);
};

template <>
struct InferDatatypeStruct<float> {
  const static inline std::shared_ptr<Datatype> Type = std::make_shared<F32Datatype>();
};

template <>
struct InferDatatypeStruct<double> {
  const static inline std::shared_ptr<Datatype> Type = std::make_shared<F64Datatype>();
};

template <>
struct InferDatatypeStruct<long double> {
  const static inline std::shared_ptr<Datatype> Type = std::make_shared<F80Datatype>();
};

template <typename T>
std::shared_ptr<Datatype> inferDatatype() {
  return InferDatatypeStruct<T>::Type;
}
} // namespace seissol::io::datatype
