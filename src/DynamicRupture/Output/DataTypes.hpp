#ifndef SEISSOL_DROUTOUT_DRDATATYPES_HPP
#define SEISSOL_DROUTOUT_DRDATATYPES_HPP

#include "Geometry.hpp"
#include "Initializer/tree/Layer.hpp"
#include <vector>
#include <array>
#include <tuple>
#include <cassert>
#include <limits>
#include <cstring>

namespace seissol::dr::output {
template <int DIM>
struct VarT {
  constexpr int dim() { return DIM; }

  real* operator[](int dim) {
    assert(dim < DIM && "access is out of the DIM. bounds");
    assert(data[dim] != nullptr && "data has been initialized yet");
    return data[dim];
  }

  real& operator()(int dim, size_t level, size_t index) {
    assert(dim < DIM && "access is out of DIM. bounds");
    assert(level < maxCacheLevel && "access is out of cache bounds");
    assert(index < size && "access is out of size bounds");
    assert(data[dim] != nullptr && "data has been initialized yet");
    return data[dim][index + level * size];
  }

  real& operator()(size_t level, size_t index) {
    static_assert(DIM == 1, "access of the overload is allowed only for 1 dim variables");
    assert(level < maxCacheLevel && "access is out of cache bounds");
    assert(index < size && "access is out of size bounds");
    assert(data[0] != nullptr && "data has been initialized yet");
    return data[0][index + level * size];
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
using DrVarsT =
    std::tuple<Var2D, Var3D, Var1D, Var2D, Var3D, Var2D, Var1D, Var1D, Var1D, Var1D, Var1D, Var2D>;

enum DirectionID { STRIKE = 0, DIP = 1, NORMAL = 2 };
enum ThermoID { PRESSURE = 0, TEMPERATURE = 1 };
enum ParamID { FUNCTION = 0, STATE = 1 };

enum VariableID {
  SlipRate = 0,
  TransientShierStress,
  NormalVelocity,
  FunctionAndState,
  TotalStresses,
  Slip,
  RuptureVelocity,
  AbsoluteSlip,
  PeakSlipsRate,
  RuptureTime,
  Ds,
  Thermo,
  Size
};

enum class OutputType : int {
  None = 0,
  AtPickpoint = 3,
  Elementwise = 4,
  AtPickpointAndElementwise = 5
};

struct GeneralParamsT {
  OutputType outputPointType{OutputType::AtPickpoint};
  int slipRateOutputType{1};
  int frictionLawType{0};
  int backgroundType{1};
  bool isRfOutputOn{false};
  bool isDsOutputOn{false};
  bool isMagnitudeOutputOn{false};
  bool isEnergyRateOutputOn{false};
  bool isGpWiseOutput{false};
  bool isTermalPressureOn{false};
  int energyRatePrintTimeInterval{1};
  bool isRfTimeOn{false};
  bool faultOutputFlag{false};
  std::string outputFilePrefix{"data"};
  std::string xdmfWriterBackend{"hdf5"};
  std::string checkPointBackend{"none"};
  double endTime{15.0};
  size_t maxIteration{10000000};
};

struct PickpointParamsT {
  std::array<bool, std::tuple_size<DrVarsT>::value> outputMask{
      true, true, true}; // the rest is false by default
  int printTimeInterval{1};
  int numOutputPoints{0};
  int maxPickStore{50};
  std::string ppFileName{};
};

struct ElementwiseFaultParamsT {
  int printTimeInterval{2};
  double printTimeIntervalSec{1.0};
  int printIntervalCriterion{1};
  int maxPickStore{50};
  std::array<bool, std::tuple_size<DrVarsT>::value> outputMask{true, true, true, true};
  int refinementStrategy{2};
  int refinement{2};
};
} // namespace seissol::dr::output

namespace seissol::dr {
struct PlusMinusBasisFunctionsT {
  std::vector<real> plusSide;
  std::vector<real> minusSide;
};

struct ConstantT {
  real p0{0.0};
  real ts0{0.0};
  real td0{0.0};
};
using ConstantsT = std::vector<ConstantT>;

struct OutputData {
  output::DrVarsT vars;
  std::vector<PlusMinusBasisFunctionsT> basisFunctions;
  std::vector<ReceiverPointT> receiverPoints;
  std::vector<std::vector<real>> rotationMatrices;
  std::vector<FaultDirectionsT> faultDirections{};
  std::vector<ConstantT> constrains;
  std::vector<double> cachedTime{};
  size_t currentCacheLevel{0};
  size_t maxCacheLevel{50};
  bool isActive{false};
};
} // namespace seissol::dr

#endif // SEISSOL_DROUTOUT_DRDATATYPES_HPP
