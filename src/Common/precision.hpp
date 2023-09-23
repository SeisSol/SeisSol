#pragma once

#include <string>
#ifdef USE_MPI
#include <mpi.h>
#endif
#ifdef USE_HDF
#include <hdf5.h>
#endif

namespace seissol {

enum class PrecisionType {
  // TODO: add half precision(s) here?
  Float,
  Double
};

template <PrecisionType precision>
struct PrecisionToType {};
template <>
struct PrecisionToType<PrecisionType::Float> {
  using Type = float;
};
template <>
struct PrecisionToType<PrecisionType::Double> {
  using Type = double;
};
template <typename RealType>
struct PrecisionFromType {};
template <>
struct PrecisionFromType<float> {
  static constexpr PrecisionType Precision = PrecisionType::Float;
  static inline const std::string Text = "float";
#ifdef USE_MPI
  static inline const MPI_Datatype MPIType = MPI_FLOAT;
#endif
#ifdef USE_HDF
  static inline const hid_t HDF5Type = H5T_IEEE_F32LE;
#endif
};
template <>
struct PrecisionFromType<double> {
  static constexpr PrecisionType Precision = PrecisionType::Double;
  static inline const std::string Text = "double";
#ifdef USE_MPI
  static inline const MPI_Datatype MPIType = MPI_DOUBLE;
#endif
#ifdef USE_HDF
  static inline const hid_t HDF5Type = H5T_IEEE_F64LE;
#endif
};

} // namespace seissol
