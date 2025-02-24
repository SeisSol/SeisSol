// SPDX-FileCopyrightText: 2019-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff
// SPDX-FileContributor: Carsten Uphoff

#ifndef SEISSOL_SRC_INITIALIZER_TREE_INTERFACEHELPER_H_
#define SEISSOL_SRC_INITIALIZER_TREE_INTERFACEHELPER_H_

#include "easi/util/Magic.h"

namespace seissol {
template <typename X>
struct extract_type {
  typedef X type;
};

template <template <typename> class F, typename X>
struct extract_type<F<X>> {
  typedef X type;
};
} // namespace seissol

#define _LTSTREE_MEMBER_PTR(N, HANDLE_STRUCT, X)                                                   \
  seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type* X = nullptr;
#define _LTSTREE_MEMBER_REF(N, HANDLE_STRUCT, X)                                                   \
  seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type& X;
#define _LTSTREE_MEMBER_REF_CS(HANDLE_STRUCT, X)                                                   \
  seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type& X
#define _LTSTREE_MEMBER_PTR_CS(HANDLE_STRUCT, X)                                                   \
  seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type* X
#define _LTSTREE_MEMBER_GETTER(N, HANDLE_STRUCT, X)                                                \
  seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type& X() noexcept { return *X##_ptr; }
#define _LTSTREE_MEMBER_CONST_GETTER(N, HANDLE_STRUCT, X)                                          \
  const seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type& X() const noexcept {              \
    return *X##_ptr;                                                                               \
  }
#define _LTSTREE_MEMBER_SUFFIXED_PTR(N, HANDLE_STRUCT, X)                                          \
  seissol::extract_type<decltype(HANDLE_STRUCT::X)>::type* X##_ptr = nullptr;
#define _LTSTREE_MEMBER_GETTERED_PTR(N, HANDLE_STRUCT, X)                                          \
  _LTSTREE_MEMBER_SUFFIXED_PTR(N, HANDLE_STRUCT, X)                                                \
  _LTSTREE_MEMBER_GETTER(N, HANDLE_STRUCT, X)                                                      \
  _LTSTREE_MEMBER_CONST_GETTER(N, HANDLE_STRUCT, X)
#define _LTSTREE_LOAD(N, HANDLE_STRUCT, X) X = tree.var(handleStruct.X, place);
#define _LTSTREE_INIT_MEMBER(HANDLE_STRUCT, X) X(X)
#define _LTSTREE_INIT_MEMBER_PTR(HANDLE_STRUCT, X) X##_ptr(X)
#define _LTSTREE_INIT_MEMBER_PTRREF(HANDLE_STRUCT, X) X##_ptr(&X)
#define _LTSTREE_ACCESS(HANDLE_STRUCT, X) X[index]
#define _LTSTREE_PTR_ACCESS(HANDLE_STRUCT, X) X + index
#define _LTSTREE_LOOKUP(HANDLE_STRUCT, X) lut.lookup(handleStruct.X, meshId, place)
#define LTSTREE_GENERATE_INTERFACE(NAME, HANDLE_STRUCT, ...)                                       \
  struct NAME {                                                                                    \
    MAGIC_FOR_EACH(_LTSTREE_MEMBER_REF, HANDLE_STRUCT, __VA_ARGS__)                                \
    NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_MEMBER_REF_CS, HANDLE_STRUCT, __VA_ARGS__))       \
        : MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_INIT_MEMBER, HANDLE_STRUCT, __VA_ARGS__) {}      \
    template <typename T>                                                                          \
    static NAME lookup(HANDLE_STRUCT const& handleStruct,                                          \
                       T const& lut,                                                               \
                       unsigned meshId,                                                            \
                       seissol::initializer::AllocationPlace place =                               \
                           seissol::initializer::AllocationPlace::Host) {                          \
      return NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_LOOKUP, HANDLE_STRUCT, __VA_ARGS__));    \
    }                                                                                              \
    struct Loader {                                                                                \
      MAGIC_FOR_EACH(_LTSTREE_MEMBER_PTR, HANDLE_STRUCT, __VA_ARGS__)                              \
      template <typename T>                                                                        \
      void load(HANDLE_STRUCT const& handleStruct,                                                 \
                T& tree,                                                                           \
                seissol::initializer::AllocationPlace place =                                      \
                    seissol::initializer::AllocationPlace::Host){                                  \
          MAGIC_FOR_EACH(_LTSTREE_LOAD, HANDLE_STRUCT, __VA_ARGS__)} NAME entry(unsigned index) {  \
        return NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_ACCESS, HANDLE_STRUCT, __VA_ARGS__));  \
      }                                                                                            \
      NAME entry(unsigned index) const {                                                           \
        return NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_ACCESS, HANDLE_STRUCT, __VA_ARGS__));  \
      }                                                                                            \
    };                                                                                             \
  };

#define LTSTREE_GENERATE_INTERFACE_GETTERED(NAME, HANDLE_STRUCT, ...)                              \
  struct NAME {                                                                                    \
    MAGIC_FOR_EACH(_LTSTREE_MEMBER_GETTERED_PTR, HANDLE_STRUCT, __VA_ARGS__)                       \
    NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_MEMBER_PTR_CS, HANDLE_STRUCT, __VA_ARGS__))       \
        : MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_INIT_MEMBER_PTR, HANDLE_STRUCT, __VA_ARGS__) {}  \
    NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_MEMBER_REF_CS, HANDLE_STRUCT, __VA_ARGS__))       \
        : MAGIC_FOR_EACH_COMMA_SEPARATED(                                                          \
              _LTSTREE_INIT_MEMBER_PTRREF, HANDLE_STRUCT, __VA_ARGS__) {}                          \
    template <typename T>                                                                          \
    static NAME lookup(HANDLE_STRUCT const& handleStruct,                                          \
                       T const& lut,                                                               \
                       unsigned meshId,                                                            \
                       seissol::initializer::AllocationPlace place =                               \
                           seissol::initializer::AllocationPlace::Host) {                          \
      return NAME(MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_LOOKUP, HANDLE_STRUCT, __VA_ARGS__));    \
    }                                                                                              \
    struct Loader {                                                                                \
      MAGIC_FOR_EACH(_LTSTREE_MEMBER_PTR, HANDLE_STRUCT, __VA_ARGS__)                              \
      template <typename T>                                                                        \
      void load(HANDLE_STRUCT const& handleStruct,                                                 \
                T& tree,                                                                           \
                seissol::initializer::AllocationPlace place =                                      \
                    seissol::initializer::AllocationPlace::Host){                                  \
          MAGIC_FOR_EACH(_LTSTREE_LOAD, HANDLE_STRUCT, __VA_ARGS__)} NAME entry(unsigned index) {  \
        return NAME(                                                                               \
            MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_PTR_ACCESS, HANDLE_STRUCT, __VA_ARGS__));      \
      }                                                                                            \
      NAME entry(unsigned index) const {                                                           \
        return NAME(                                                                               \
            MAGIC_FOR_EACH_COMMA_SEPARATED(_LTSTREE_PTR_ACCESS, HANDLE_STRUCT, __VA_ARGS__));      \
      }                                                                                            \
    };                                                                                             \
  };

#endif // SEISSOL_SRC_INITIALIZER_TREE_INTERFACEHELPER_H_
