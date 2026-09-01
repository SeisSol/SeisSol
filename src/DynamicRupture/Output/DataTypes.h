// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_DATATYPES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_DATATYPES_H_

#include "Common/Iterator.h"
#include "GeneratedCode/tensor.h"
#include "Geometry.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Memory/Tree/Backmap.h"
#include "Parallel/DataCollector.h"
#include "Parallel/Runtime/Stream.h"

#include <Eigen/Dense>
#include <array>
#include <cassert>
#include <cstddef>
#include <cstring>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <unordered_map>
#include <vector>

namespace seissol::dr::output {
template <std::size_t Dim>
struct VarT {
  VarT() = default;
  [[nodiscard]] constexpr std::size_t dim() const { return Dim; }

  real* operator[](std::size_t dim) {
    assert(dim < Dim && "access is out of the Dim. bounds");
    return data[dim].data();
  }

  real& operator()(std::size_t dim, size_t level, size_t index) {
    assert(dim < Dim && "access is out of Dim. bounds");
    assert(level < maxCacheLevel && "access is out of cache bounds");
    assert(index < size && "access is out of size bounds");
    assert(!data[dim].empty() && "data has not been initialized yet");
    return data[dim][index + level * size];
  }

  real& operator()(size_t level, size_t index) {
    static_assert(Dim == 1, "access of the overload is allowed only for 1 dim variables");
    return this->operator()(0, level, index);
  }

  const real* operator[](std::size_t dim) const {
    assert(dim < Dim && "access is out of the Dim. bounds");
    return data[dim].data();
  }

  const real& operator()(std::size_t dim, size_t level, size_t index) const {
    assert(dim < Dim && "access is out of Dim. bounds");
    assert(level < maxCacheLevel && "access is out of cache bounds");
    assert(index < size && "access is out of size bounds");
    assert(!data[dim].empty() && "data has not been initialized yet");
    return data[dim][index + level * size];
  }

  const real& operator()(size_t level, size_t index) const {
    static_assert(Dim == 1, "access of the overload is allowed only for 1 dim variables");
    return this->operator()(0, level, index);
  }

  // allocates data for a var (for all dimensions and cache levels)
  // initialized to zeros if var is active.
  // Otherwise, inits with nullptr
  void allocateData(size_t dataSize) {
    size = dataSize;
    if (isActive) {
      for (std::size_t dim = 0; dim < Dim; ++dim) {
        data[dim].resize(size * maxCacheLevel);
      }
    } else {
      for (std::size_t dim = 0; dim < Dim; ++dim) {
        data[dim].resize(0);
      }
    }
  }

  void resizeCache(size_t newMaxCacheLevel) {
    maxCacheLevel = newMaxCacheLevel;
    if (isActive) {
      for (std::size_t dim = 0; dim < Dim; ++dim) {
        data[dim].resize(size * maxCacheLevel);
      }
    }
  }

  std::array<std::vector<real>, Dim> data;
  bool isActive{false};
  size_t size{0};
  size_t maxCacheLevel{1};
};

using Var1D = VarT<1>;
using Var2D = VarT<2>;
using Var3D = VarT<3>;

// Description is given in `enum VariableID`
using DrVarsT =
    std::tuple<Var2D, Var3D, Var1D, Var2D, Var3D, Var2D, Var1D, Var1D, Var1D, Var1D, Var1D, Var2D>;

enum DirectionID { Strike = 0, Dip = 1, Normal = 2 };
enum TPID { Pressure = 0, Temperature = 1 };
enum ParamID { FrictionCoefficient = 0, State = 1 };

enum VariableID {
  SlipRate = 0,
  TransientTractions,
  NormalVelocity,
  FrictionAndState,
  TotalTractions,
  Slip,
  RuptureVelocity,
  AccumulatedSlip,
  PeakSlipRate,
  RuptureTime,
  DynamicStressTime,
  ThermalPressurizationVariables,
  Size
};

const inline std::vector<std::vector<std::string>> VariableLabels = {{"SRs", "SRd"},
                                                                     {"T_s", "T_d", "P_n"},
                                                                     {"u_n"},
                                                                     {"Mud", "StV"},
                                                                     {"Ts0", "Td0", "Pn0"},
                                                                     {"Sls", "Sld"},
                                                                     {"Vr"},
                                                                     {"ASl"},
                                                                     {"PSR"},
                                                                     {"RT"},
                                                                     {"DS"},
                                                                     {"P_f", "Tmp"}};

} // namespace seissol::dr::output

namespace seissol::dr {
struct PlusMinusBasisFunctions {
  std::vector<real> plusSide;
  std::vector<real> minusSide;
};

/**
  Data of a single fault face taking part in the output.

  Everything stored here is a function of the fault face alone, so it is set up once per face and
  evaluated once per face and time step, no matter how many output points the face carries.
 */
struct OutputFace {
  std::size_t faultFaceIndex{};
  ::seissol::initializer::StoragePosition position{};

  std::size_t elementIndex{};
  std::size_t localFaceSideId{};

  FaultDirections faultDirections{};
  std::array<real, seissol::tensor::stressRotationMatrix::size()> stressGlbToDipStrikeAligned{};
  std::array<real, seissol::tensor::stressRotationMatrix::size()> stressFaceAlignedToGlb{};
  std::array<real, seissol::tensor::Tinv::size()> glbToFaceAlignedData{};
  Eigen::Matrix<real, 2, 2> jacobianT2d{Eigen::Matrix<real, 2, 2>::Zero()};

  // gather indices into ReceiverOutputData::deviceDataCollector
  std::size_t deviceDataPlus{};
  std::size_t deviceDataMinus{};
};

/**
  Data of a single output location on a fault face; i.e. everything which depends on the position
  within the face, but not on the simulation a receiver belongs to.
 */
struct OutputPoint {
  std::size_t faceId{};
  PlusMinusBasisFunctions basisFunctions;
  std::size_t nearestGpIndex{};
  std::size_t nearestInternalGpIndex{};
};

/**
  The face -> point -> receiver hierarchy of the on-fault output, in CSR form.

  A face owns a contiguous range of points, and a point owns a contiguous range of receivers. The
  receivers are renumbered while the topology is built, which is what makes the ranges contiguous
  and reduces the structure to two offset vectors -- no index arrays are needed.

    points of face f:     [pointOffset[f],     pointOffset[f + 1])
    receivers of point p: [receiverOffset[p],  receiverOffset[p + 1])
 */
struct OutputTopology {
  std::vector<OutputFace> faces;
  std::vector<OutputPoint> points;
  std::vector<std::size_t> pointOffset{0};
  std::vector<std::size_t> receiverOffset{0};

  [[nodiscard]] std::size_t faceCount() const { return faces.size(); }
  [[nodiscard]] std::size_t pointCount() const { return points.size(); }
  [[nodiscard]] std::size_t receiverCount() const { return receiverOffset.back(); }

  [[nodiscard]] auto pointsOf(std::size_t face) const {
    return common::range(pointOffset[face], pointOffset[face + 1]);
  }

  [[nodiscard]] auto receiversOf(std::size_t point) const {
    return common::range(receiverOffset[point], receiverOffset[point + 1]);
  }

  /**
    The first receiver of a point; usable wherever a representative of the point is needed (e.g.
    for its coordinates, which all its receivers share).
   */
  [[nodiscard]] std::size_t representative(std::size_t point) const {
    return receiverOffset[point];
  }

  void addFace(const OutputFace& face) {
    faces.push_back(face);
    pointOffset.push_back(points.size());
  }

  void addPoint(const OutputPoint& point, std::size_t receiverCount) {
    points.push_back(point);
    receiverOffset.push_back(receiverOffset.back() + receiverCount);
    pointOffset.back() = points.size();
  }
};

struct ReceiverOutputData {
  output::DrVarsT vars;
  std::vector<Receiver> receivers;
  OutputTopology topology;

  std::vector<double> cachedTime;
  size_t currentCacheLevel{0};
  size_t maxCacheLevel{50};
  bool isActive{false};
  std::optional<int64_t> clusterId;

  std::unique_ptr<parallel::DataCollector<real>> deviceDataCollector;
  std::size_t cellCount{0};

  std::unordered_map<std::size_t, std::unique_ptr<parallel::DataCollectorUntyped>> deviceVariables;
  std::optional<parallel::runtime::StreamRuntime> extraRuntime;
};
} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_DATATYPES_H_
