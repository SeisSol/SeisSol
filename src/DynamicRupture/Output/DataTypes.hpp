#ifndef SEISSOL_DR_OUTPUT_DATA_TYPES_HPP
#define SEISSOL_DR_OUTPUT_DATA_TYPES_HPP

#include "Geometry.hpp"
#include "Initializer/tree/Layer.hpp"
#include <array>
#include <cassert>
#include <cstring>
#include <limits>
#include <tuple>
#include <vector>

namespace seissol::dr::output {
template <int DIM>
struct VarT {
  constexpr int dim() { return DIM; }

  real* operator[](int dim) {
    assert(dim < DIM && "access is out of the DIM. bounds");
    assert(data[dim] != nullptr && "data has not been initialized yet");
    return data[dim];
  }

  real& operator()(int dim, size_t level, size_t index) {
    assert(dim < DIM && "access is out of DIM. bounds");
    assert(level < maxCacheLevel && "access is out of cache bounds");
    assert(index < size && "access is out of size bounds");
    assert(data[dim] != nullptr && "data has not been initialized yet");
    return data[dim][index + level * size];
  }

  real& operator()(size_t level, size_t index) {
    static_assert(DIM == 1, "access of the overload is allowed only for 1 dim variables");
    return this->operator()(0, level, index);
  }

  // allocates data for a var (for all dimensions and cache levels) initialized to zeros
  // if var is active. Otherwise, inits with nullptr
  void allocateData(size_t dataSize) {
    size = dataSize;
    if (isActive) {
      for (int dim = 0; dim < DIM; ++dim) {
        data[dim] = new real[size * maxCacheLevel];
        std::memset(static_cast<void*>(data[dim]), 0, size * maxCacheLevel * sizeof(real));
      }
    } else {
      for (int dim = 0; dim < DIM; ++dim)
        data[dim] = nullptr;
    }
  }

  void releaseData() {
    if (isActive) {
      for (auto item : data) {
        delete[] item;
      }
    }
  }

  std::array<real*, DIM> data{};
  bool isActive{false};
  size_t size{};
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
  TpVariables, // Thermal Pressurization
  Size
};

enum class OutputType : int {
  None = 0,
  AtPickpoint = 3,
  Elementwise = 4,
  AtPickpointAndElementwise = 5
};

struct GeneralParamsT {
  OutputType outputPointType{OutputType::None};
  int slipRateOutputType{1};
  int frictionLawType{0};
  bool isRfOutputOn{false};
  bool isDsOutputOn{false};
  bool isMagnitudeOutputOn{false};
  bool isEnergyRateOutputOn{false};
  bool isThermalPressurizationOn{false};
  int energyRatePrintTimeInterval{30};
  bool isRfTimeOn{false};
  bool faultOutputFlag{false};
  std::string outputFilePrefix{"data"};
  std::string xdmfWriterBackend{"hdf5"};
  std::string checkPointBackend{"none"};
  double endTime{15.0};
  size_t maxIteration{1000000000};
};

struct PickpointParamsT {
  std::array<bool, std::tuple_size<DrVarsT>::value> outputMask{
      true, true, true}; // the rest is false by default
  int printTimeInterval{1};
  int maxPickStore{50};
  std::string ppFileName{};
};

enum class RefinerType { Triple = 1, Quad = 2, Invalid = 3 };

struct ElementwiseFaultParamsT {
  double printTimeIntervalSec{1.0};
  std::array<bool, std::tuple_size<DrVarsT>::value> outputMask{true, true, true, true};
  RefinerType refinementStrategy{RefinerType::Quad};
  int refinement{2};
};
} // namespace seissol::dr::output

namespace seissol::dr {
struct PlusMinusBasisFunctionsT {
  std::vector<real> plusSide;
  std::vector<real> minusSide;
};

struct IntialTraction {
  real p0{0.0};
  real ts0{0.0};
  real td0{0.0};
};
using ConstantsT = std::vector<IntialTraction>;

struct OutputData {
  output::DrVarsT vars;
  std::vector<PlusMinusBasisFunctionsT> basisFunctions;
  std::vector<ReceiverPointT> receiverPoints;
  std::vector<std::vector<real>> glbToDipStrikeAligned;
  std::vector<FaultDirectionsT> faultDirections{};
  std::vector<IntialTraction> intialTractions;
  std::vector<double> cachedTime{};
  size_t currentCacheLevel{0};
  size_t maxCacheLevel{50};
  bool isActive{false};
};
} // namespace seissol::dr

#endif // SEISSOL_DR_OUTPUT_DATA_TYPES_HPP
