// SPDX-FileCopyrightText: 2021 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_DATATYPES_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_DATATYPES_H_

#include "GeneratedCode/tensor.h"
#include "Geometry.h"
#include "Initializer/Parameters/DRParameters.h"
#include "Kernels/Precision.h"
#include "Memory/Descriptor/DynamicRupture.h"
#include "Parallel/DataCollector.h"
#include <Common/Templating.h>
#include <Config.h>
#include <Eigen/Dense>
#include <Parallel/Runtime/Stream.h>
#include <array>
#include <cassert>
#include <cstring>
#include <limits>
#include <string>
#include <tuple>
#include <vector>

namespace seissol::dr::output {
template <std::size_t DIM, typename RealT>
struct VarT {
  VarT() = default;
  ~VarT() { releaseData(); }
  constexpr std::size_t dim() { return DIM; }

  VarT(const VarT&) = delete;
  auto operator=(const VarT&) -> VarT& = delete;

  VarT(VarT&&) = default;
  auto operator=(VarT&&) -> VarT& = default;

  RealT* operator[](std::size_t dim) {
    assert(dim < DIM && "access is out of the DIM. bounds");
    assert(data[dim] != nullptr && "data has not been initialized yet");
    return data[dim];
  }

  RealT& operator()(std::size_t dim, size_t level, size_t index) {
    assert(dim < DIM && "access is out of DIM. bounds");
    assert(level < maxCacheLevel && "access is out of cache bounds");
    assert(index < size && "access is out of size bounds");
    assert(data[dim] != nullptr && "data has not been initialized yet");
    return data[dim][index + level * size];
  }

  RealT& operator()(size_t level, size_t index) {
    static_assert(DIM == 1, "access of the overload is allowed only for 1 dim variables");
    return this->operator()(0, level, index);
  }

  // allocates data for a var (for all dimensions and cache levels)
  // initialized to zeros if var is active.
  // Otherwise, inits with nullptr
  void allocateData(size_t dataSize) {
    size = dataSize;
    if (isActive) {
      for (std::size_t dim = 0; dim < DIM; ++dim) {
        assert(data[dim] == nullptr && "double allocation is not allowed");
        data[dim] = new RealT[size * maxCacheLevel];
        std::memset(static_cast<void*>(data[dim]), 0, size * maxCacheLevel * sizeof(RealT));
      }
    } else {
      for (std::size_t dim = 0; dim < DIM; ++dim) {
        data[dim] = nullptr;
      }
    }
  }

  void releaseData() {
    if (isActive) {
      for (auto& item : data) {
        delete[] item;
        item = nullptr;
      }
    }
  }

  std::array<RealT*, DIM> data{};
  bool isActive{false};
  size_t size{};
  size_t maxCacheLevel{1};
};

using Var1D = VarT<1, double>;
using Var2D = VarT<2, double>;
using Var3D = VarT<3, double>;

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

using FaceToLtsMapType = std::vector<std::pair<DynamicRupture::Layer*, size_t>>;

} // namespace seissol::dr::output

namespace seissol::dr {
template<typename RealT>
struct PlusMinusBasisFunctions {
  std::vector<RealT> plusSide;
  std::vector<RealT> minusSide;
};

template<typename Cfg>
struct TransformData {
  PlusMinusBasisFunctions<Real<Cfg>> basisFunctions;
  std::array<Real<Cfg>, seissol::tensor::stressRotationMatrix<Cfg>::size()> stressGlbToDipStrikeAligned;
  std::array<Real<Cfg>, seissol::tensor::stressRotationMatrix<Cfg>::size()> stressFaceAlignedToGlb;
  std::array<Real<Cfg>, seissol::tensor::T<Cfg>::size()> faceAlignedToGlbData;
  std::array<Real<Cfg>, seissol::tensor::Tinv<Cfg>::size()> glbToFaceAlignedData;
  Eigen::Matrix<Real<Cfg>, 2, 2> jacobianT2d;
};

using TransformDataVariant = TransformVariadicT<TransformData, ConfigVariant>;

struct ReceiverOutputData {
  output::DrVarsT vars;

  std::vector<TransformDataVariant, Eigen::aligned_allocator<TransformDataVariant>> transformData;
  std::vector<ReceiverPoint> receiverPoints;

  std::vector<FaultDirections> faultDirections;
  std::vector<double> cachedTime;
  size_t currentCacheLevel{0};
  size_t maxCacheLevel{50};
  bool isActive{false};
  std::optional<int64_t> clusterId;

  std::unique_ptr<parallel::DataCollector<Real<Cfg>>> deviceDataCollector;
  std::vector<std::size_t> deviceDataPlus;
  std::vector<std::size_t> deviceDataMinus;
  std::size_t cellCount{0};

  std::unordered_map<std::size_t, std::unique_ptr<parallel::DataCollector<Real<Cfg>>>> deviceVariables;
  std::vector<std::size_t> deviceIndices;
  std::optional<parallel::runtime::StreamRuntime> extraRuntime;
};
} // namespace seissol::dr

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_DATATYPES_H_
