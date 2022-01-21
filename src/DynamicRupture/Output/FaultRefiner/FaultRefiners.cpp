#ifndef SEISSOL_TRIPLEFAULTFACEREFINER_HPP
#define SEISSOL_TRIPLEFAULTFACEREFINER_HPP

#include "FaultRefiners.hpp"
#include "DynamicRupture/Output/OutputAux.hpp"
#include "utils/logger.h"
#include <memory>

namespace seissol::dr::output::refiner {
RefinerType convertToType(int strategy) {
  switch (strategy) {
  case 1:
    return RefinerType::Triple;
  case 2:
    return RefinerType::Quad;
  default:
    logError() << "Unknown refinement strategy for Fault Face Refiner";
    // return something to suppress a compiler warning
    return RefinerType::Invalid;
  }
}

std::unique_ptr<FaultRefiner> get(RefinerType strategy) {
  switch (strategy) {
  case RefinerType::Triple:
    return std::unique_ptr<FaultRefiner>(new TripleFaultFaceRefiner);
  case RefinerType::Quad:
    return std::unique_ptr<FaultRefiner>(new QuadFaultFaceRefiner);
  case RefinerType::Invalid:
  default:
    return nullptr;
  }
}

void FaultRefiner::repeat(Data data, PointsPair& point1, PointsPair& point2, PointsPair& point3) {
  ExtTriangle subGlobalFace(
      std::get<global>(point1), std::get<global>(point2), std::get<global>(point3));

  ExtTriangle subReferenceFace(
      std::get<reference>(point1), std::get<reference>(point2), std::get<reference>(point3));

  refineAndAccumulate({data.refinementLevel - 1, data.faultFaceIndex, data.localFaceSideId},
                      std::make_pair(subGlobalFace, subReferenceFace));
}

void FaultRefiner::addReceiver(Data data, TrianglePair& face) {
  ReceiverPointT receiver{};
  receiver.isInside = true;
  receiver.faultFaceIndex = data.faultFaceIndex;
  receiver.localFaceSideId = data.localFaceSideId;
  receiver.globalReceiverIndex = points.size();
  receiver.global = getMidTrianglePoint(std::get<global>(face));
  receiver.reference = getMidTrianglePoint(std::get<reference>(face));
  receiver.globalTriangle = std::get<global>(face);

  points.push_back(receiver);
}

void TripleFaultFaceRefiner::refineAndAccumulate(Data data, TrianglePair face) {

  if (data.refinementLevel == 0) {
    addReceiver(data, face);
    return;
  }

  auto& globalFace = std::get<global>(face);
  auto& referenceFace = std::get<reference>(face);

  auto midPoint =
      std::make_pair(getMidTrianglePoint(globalFace), getMidTrianglePoint(referenceFace));
  std::array<PointsPair, 3> points{};
  for (size_t i = 0; i < 3; ++i) {
    points[i] = std::make_pair(globalFace[i], referenceFace[i]);
  }

  repeat(data, points[0], points[1], midPoint);
  repeat(data, midPoint, points[1], points[2]);
  repeat(data, points[0], midPoint, points[2]);
}

void QuadFaultFaceRefiner::refineAndAccumulate(Data data, TrianglePair face) {

  if (data.refinementLevel == 0) {
    addReceiver(data, face);
    return;
  }

  auto& globalFace = std::get<global>(face);
  auto& referenceFace = std::get<reference>(face);

  auto split = [&globalFace, &referenceFace](size_t pointIndex1, size_t pointIndex2) {
    return std::make_pair(getMidPoint(globalFace[pointIndex1], globalFace[pointIndex2]),
                          getMidPoint(referenceFace[pointIndex1], referenceFace[pointIndex2]));
  };

  auto midPoint1 = split(0, 1);
  auto midPoint2 = split(1, 2);
  auto midPoint3 = split(2, 0);

  PointsPair trianglePoint = std::make_pair(globalFace.p1, referenceFace.p1);
  repeat(data, trianglePoint, midPoint1, midPoint3);

  trianglePoint = std::make_pair(globalFace.p2, referenceFace.p2);
  repeat(data, midPoint1, trianglePoint, midPoint2);

  repeat(data, midPoint1, midPoint2, midPoint3);

  trianglePoint = std::make_pair(globalFace.p3, referenceFace.p3);
  repeat(data, midPoint3, midPoint2, trianglePoint);
}
} // namespace seissol::dr::output::refiner

#endif // SEISSOL_TRIPLEFAULTFACEREFINER_HPP
