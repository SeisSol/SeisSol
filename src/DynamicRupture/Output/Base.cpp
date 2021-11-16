#include "Solver/Interoperability.h"
#include "Initializer/tree/Layer.hpp"
#include "Initializer/DynamicRupture.h"
#include "SeisSol.h"
#include "Base.hpp"
#include "ResultWriter/common.hpp"
#include <fstream>
#include <unordered_map>
#include <filesystem>
#include <type_traits>

struct NativeFormat {};
struct WideFormat {};
template<typename T, typename U = NativeFormat>
struct FormattedBuildInType {
  T value;
};

template<typename T, typename U = NativeFormat>
auto makeFormatted(T value) {
  return FormattedBuildInType<T, U>{value};
}

template<typename T, typename U = NativeFormat>
std::ostream& operator<<(std::ostream& stream, FormattedBuildInType<T, U> obj) {
  if constexpr (std::is_floating_point_v<T>) {
    stream << std::setprecision(16) << std::scientific << obj.value;
  }
  else if constexpr (std::is_integral_v<T> && std::is_same_v<U, WideFormat>) {
    stream << std::setw(5) << std::setfill('0') << obj.value;
  }
  else {
    stream << obj.value;
  }
  return stream;
}


std::string seissol::dr::output::Base::constructPpReveiverFileName(const int receiverGlobalIndex) const {
  std::stringstream fileName;
  fileName << generalParams.outputFilePrefix
           << "-new-faultreceiver-"
           << makeFormatted<int, WideFormat>(receiverGlobalIndex);
#ifdef PARALLEL
  fileName << '-' << makeFormatted<int, WideFormat>(MPI::mpi.rank());
#endif
  fileName << ".dat";
  return fileName.str();
}


void seissol::dr::output::Base::initEwOutput(const std::unordered_map<std::string, double*>& FaultParams) {
  ewOutputBuilder->init(FaultParams);

  const auto& receiverPoints = ewOutputBuilder->outputData.receiverPoints;
  auto cellConnectivity = getCellConnectivity(receiverPoints);
  auto vertices = getAllVertices(receiverPoints);
  constexpr auto maxNumVars = std::tuple_size<DrVarsT>::value;
  auto intMask = convertMaskFromBoolToInt<maxNumVars>(ewOutputBuilder->elementwiseParams.outputMask);

  std::string newFaultFilePrefix = generalParams.outputFilePrefix + std::string("-new");
  double printTime = ewOutputBuilder->elementwiseParams.printTimeIntervalSec;
  auto backendType = seissol::writer::backendType(generalParams.xdmfWriterBackend.c_str());

  std::vector<real*> dataPointers;
  auto recordPointers = [&dataPointers](auto& var, int) {
    if (var.isActive) {
      for (int dim = 0; dim < var.dim(); ++dim)
        dataPointers.push_back(var.data[dim]);
    }
  };
  aux::forEach(ewOutputBuilder->outputData.vars, recordPointers);


  seissol::SeisSol::main.secondFaultWriter().init(cellConnectivity.data(),
                                                  vertices.data(),
                                                  static_cast<unsigned int>(receiverPoints.size()),
                                                  static_cast<unsigned int>(3 * receiverPoints.size()),
                                                  &intMask[0],
                                                  const_cast<const real**>(dataPointers.data()),
                                                  newFaultFilePrefix.data(),
                                                  printTime,
                                                  backendType);

  seissol::SeisSol::main.secondFaultWriter().setupCallbackObject(this);
}


void seissol::dr::output::Base::initPickpointOutput() {
  ppOutputBuilder->init();

  std::stringstream baseHeader;
  baseHeader << "VARIABLES = \"Time\"";
  size_t labelCounter = 0;
  auto collectVariableNames = [&baseHeader, &labelCounter](auto& var, int i) {
    if (var.isActive) {
      for (int dim = 0; dim < var.dim(); ++dim) {
        baseHeader << " ,\"" << writer::FaultWriterExecutor::getLabelName(labelCounter) << '\"';
        ++labelCounter;
      }
    }
    else {
      labelCounter += var.dim();
    }
  };
  aux::forEach(ppOutputBuilder->outputData.vars, collectVariableNames);

  auto& outputData = ppOutputBuilder->outputData;
  for (const auto& point: outputData.receiverPoints) {
    const size_t globalIndex = point.globalReceiverIndex;
    auto fileName = constructPpReveiverFileName(globalIndex);
    if (!std::filesystem::exists(fileName)) {

      std::ofstream file(fileName, std::ios_base::out);
      if (file.is_open()) {
        std::stringstream title;
        title << "TITLE = \"Temporal Signal for fault receiver number " << globalIndex << "\"";

        file << title.str() << '\n';
        file << baseHeader.str() << '\n';
        file << "# x1\t" << makeFormatted(point.global.x) << '\n';
        file << "# x2\t" << makeFormatted(point.global.y) << '\n';
        file << "# x3\t" << makeFormatted(point.global.z) << '\n';
      }
      else {
        logError() << "cannot open " << fileName;
      }
      file.close();
    }
  }
}


void seissol::dr::output::Base::init(const std::unordered_map<std::string, double*>& FaultParams) {
  if (ewOutputBuilder) {
    initEwOutput(FaultParams);
  }
  if (ppOutputBuilder) {
    initPickpointOutput();
  }
}


void seissol::dr::output::Base::initFaceToLtsMap() {
  faceToLtsMap.resize(drTree->getNumberOfCells(Ghost));
  for (auto it = drTree->beginLeaf(seissol::initializers::LayerMask(Ghost));
    it != drTree->endLeaf(); ++it) {

    DRFaceInformation *faceInformation = it->var(dynRup->faceInformation);
    for (size_t ltsFace = 0; ltsFace < it->getNumberOfCells(); ++ltsFace) {
      faceToLtsMap[faceInformation[ltsFace].meshFace] = std::make_pair(&(*it), ltsFace);
    }
  }
}


bool seissol::dr::output::Base::isAtPickpoint(double time, double dt) {
  bool isFirstStep = iterationStep == 0;

  const double abortTime = std::min(generalParams.endTime, generalParams.maxIteration * dt);
  constexpr double timeMargin = 1.005;
  const bool isCloseToTimeOut = (abortTime - time) < (dt * timeMargin);

  const int printTimeInterval = ppOutputBuilder->pickpointParams.printTimeInterval;
  const bool isOutputIteration = iterationStep % printTimeInterval == 0;

  return ppOutputBuilder && (isFirstStep || isOutputIteration || isCloseToTimeOut);
}


void seissol::dr::output::Base::writePickpointOutput(double time, double dt) {

  if (this->isAtPickpoint(time, dt)) {

    auto& outputData = ppOutputBuilder->outputData;
    calcFaultOutput(OutputType::AtPickpoint, outputData, time);

    const bool isMaxCacheLevel = outputData.currentCacheLevel >= static_cast<size_t>(ppOutputBuilder->pickpointParams.maxPickStore);
    const bool isCloseToEnd = (generalParams.endTime - time) < dt * 1.005;
    if (isMaxCacheLevel || isCloseToEnd) {
      for (size_t pointId = 0; pointId < outputData.receiverPoints.size(); ++pointId) {

        std::stringstream data;
        for (size_t level = 0; level < outputData.currentCacheLevel; ++level) {
          data << makeFormatted(outputData.cachedTime[level]) << '\t';
          auto recordResults = [pointId, level, &data](auto &var, int) {
            if (var.isActive) {
              for (int dim = 0; dim < var.dim(); ++dim) {
                data << makeFormatted(var(dim, level, pointId)) << '\t';
              }
            }
          };
          aux::forEach(outputData.vars, recordResults);
          data << '\n';
        }

        auto fileName = constructPpReveiverFileName(outputData.receiverPoints[pointId].globalReceiverIndex);
        std::ofstream file(fileName, std::ios_base::app);
        if (file.is_open()) {
          file << data.str();
        }
        else {
          logError() << "cannot open " << fileName;
        }
        file.close();
      }
      outputData.currentCacheLevel = 0;
    }
  }
  iterationStep += 1;
}


void seissol::dr::output::Base::updateElementwiseOutput() {
  calcFaultOutput(OutputType::Elementwise, ewOutputBuilder->outputData);
}


using DrPaddedArrayT = real (*)[seissol::init::QInterpolated::Stop[0]];
void seissol::dr::output::Base::calcFaultOutput(const OutputType type, OutputData& outputData, double time) {

  size_t level = (type == OutputType::AtPickpoint) ? outputData.currentCacheLevel : 0;
  for (size_t i = 0; i < outputData.receiverPoints.size(); ++i) {

    auto ltsMap = faceToLtsMap[outputData.receiverPoints[i].faultFaceIndex];
    const auto layer = ltsMap.first;
    const auto ltsId = ltsMap.second;

    auto mu = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->mu));
    auto rt = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->rupture_time));
    auto slip = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->slip));
    auto peakSR = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->peakSR));

    const auto nearestGaussPoint = outputData.receiverPoints[i].nearestGpIndex;
    const auto faceIndex = outputData.receiverPoints[i].faultFaceIndex;

    const auto normal = outputData.faultDirections[faceIndex].faceNormal;
    const auto tangent1 = outputData.faultDirections[faceIndex].tangent1;
    const auto tangent2 = outputData.faultDirections[faceIndex].tangent2;
    const auto strike = outputData.faultDirections[faceIndex].strike;
    const auto dip = outputData.faultDirections[faceIndex].dip;

    auto &functionAndState = std::get<VariableID::FunctionAndState>(outputData.vars);
    if (functionAndState.isActive) {
      functionAndState(ParamID::FUNCTION, level, i) = mu[ltsId][nearestGaussPoint];
    }


    auto &ruptureTime = std::get<VariableID::RuptureTime>(outputData.vars);
    if (ruptureTime.isActive) {
      ruptureTime(level, i) = rt[ltsId][nearestGaussPoint];
    }


    auto &absoluteSlip = std::get<VariableID::AbsoluteSlip>(outputData.vars);
    if (absoluteSlip.isActive) {
      absoluteSlip(level, i) = slip[ltsId][nearestGaussPoint];
    }


    auto &peakSlipsRate = std::get<VariableID::PeakSlipsRate>(outputData.vars);
    if (peakSlipsRate.isActive) {
      peakSlipsRate(level, i) = peakSR[ltsId][nearestGaussPoint];
    }


    auto &slipVectors = std::get<VariableID::Slip>(outputData.vars);
    if (slipVectors.isActive) {
      VrtxCoords crossProduct = {0.0, 0.0, 0.0};
      MeshTools::cross(strike, tangent1, crossProduct);

      double cos1 = MeshTools::dot(strike, tangent1);
      double scalarProd = MeshTools::dot(crossProduct, normal);

      // Note: cos1**2 can be greater than 1.0 because of rounding errors -> min
      double sin1 = std::sqrt(1.0 - std::min(1.0, cos1 * cos1));
      sin1 = (scalarProd > 0) ? sin1 : - sin1;

      auto slip1 = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->slipStrike));
      auto slip2 = reinterpret_cast<DrPaddedArrayT>(layer->var(dynRup->slipDip));

      slipVectors(DirectionID::STRIKE, level, i) = cos1 * slip1[ltsId][nearestGaussPoint] -
                                                        sin1 * slip2[ltsId][nearestGaussPoint];

      slipVectors(DirectionID::DIP, level, i) = sin1 * slip1[ltsId][nearestGaussPoint] +
                                                     cos1 * slip2[ltsId][nearestGaussPoint];
    }
  }

  if (type == OutputType::AtPickpoint) {
    outputData.cachedTime[outputData.currentCacheLevel] = time;
    outputData.currentCacheLevel += 1;
  }
}