// SPDX-FileCopyrightText: 2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_IO_DATATYPE_INFERENCE_H_
#define SEISSOL_SRC_IO_DATATYPE_INFERENCE_H_

#include "Datatype.h"
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

template <typename T>
std::shared_ptr<Datatype> inferDatatype();

template <typename T, typename = void>
struct InferDatatypeStruct {
  static auto type() { return OpaqueDatatype(sizeof(T)); }
};

template <typename T>
struct InferDatatypeStruct<T, std::enable_if_t<std::is_array_v<T>>> {
  static auto type() {
    return ArrayDatatype(inferDatatype<std::remove_all_extents_t<T>>(), arrayTotalExtent<T>());
  }
};

template <typename T, std::size_t N>
struct InferDatatypeStruct<std::array<T, N>> {
  static auto type() { return ArrayDatatype(inferDatatype<T>(), {N}); }
};

template <typename... Ts>
struct InferDatatypeStruct<std::tuple<Ts...>> {
  private:
  using TupleT = std::tuple<Ts...>;
  template <std::size_t I, typename T>
  static StructDatatype::MemberInfo makeMemberInfo() {
    TupleT instance;
    auto elemname = std::string("_") + std::to_string(I);
    auto elemoffset = reinterpret_cast<std::size_t>(&std::get<I>(instance)) -
                      reinterpret_cast<std::size_t>(&instance);
    auto elemtype = inferDatatype<T>();
    return StructDatatype::MemberInfo{elemname, elemoffset, elemtype};
  }

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
  static auto type() { return StructDatatype(makeMemberInfos<0, Ts...>(), sizeof(TupleT)); }
};

template <typename T1, typename T2>
struct InferDatatypeStruct<std::pair<T1, T2>> {
  private:
  using PairT = std::pair<T1, T2>;

  public:
  static auto type() {
    return StructDatatype(
        std::vector<StructDatatype::MemberInfo>{
            StructDatatype::MemberInfo{"first", offsetof(PairT, first), inferDatatype<T1>()},
            StructDatatype::MemberInfo{"second", offsetof(PairT, second), inferDatatype<T2>()},
        },
        sizeof(PairT));
  }
};

template <typename T>
struct InferDatatypeStruct<
    T,
    std::enable_if_t<
        std::is_same_v<decltype(T::datatypeLayout()), std::vector<StructDatatype::MemberInfo>>>> {
  static auto type() { return StructDatatype(T::datatypeLayout(), sizeof(T)); }
};

template <typename T>
struct InferDatatypeStruct<T, std::enable_if_t<std::is_integral_v<T>>> {
  static auto type() { return IntegerDatatype(sizeof(T), std::numeric_limits<T>::is_signed); }
};

template <>
struct InferDatatypeStruct<float> {
  static auto type() { return F32Datatype(); }
};

template <>
struct InferDatatypeStruct<double> {
  static auto type() { return F64Datatype(); }
};

template <>
struct InferDatatypeStruct<long double> {
  static auto type() { return F80Datatype(); }
};

template <typename T>
std::shared_ptr<Datatype> inferDatatype() {
  return std::make_shared<decltype(InferDatatypeStruct<T>::type())>(InferDatatypeStruct<T>::type());
}
} // namespace seissol::io::datatype

#endif // SEISSOL_SRC_IO_DATATYPE_INFERENCE_H_
