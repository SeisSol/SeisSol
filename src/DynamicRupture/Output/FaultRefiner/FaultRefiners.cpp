#include "FaultRefiners.hpp"

#include <memory>

#include "utils/logger.h"

#include "DynamicRupture/Output/OutputAux.hpp"
#include "Initializer/parameters/OutputParameters.h"

namespace seissol::dr::output::refiner {
seissol::initializers::parameters::FaultRefinement castToRefinerType(int strategy) {
  switch (strategy) {
  case 1:
    return seissol::initializers::parameters::FaultRefinement::Triple;
  case 2:
    return seissol::initializers::parameters::FaultRefinement::Quad;
  default:
    logError() << "Unknown refinement strategy for Fault Face Refiner";
    // return something to suppress a compiler warning
    return seissol::initializers::parameters::FaultRefinement::None;
  }
}

std::unique_ptr<FaultRefiner> get(seissol::initializers::parameters::FaultRefinement strategy) {
  switch (strategy) {
  case seissol::initializers::parameters::FaultRefinement::Triple:
    return std::make_unique<FaultFaceTripleRefiner>();
  case seissol::initializers::parameters::FaultRefinement::Quad:
    return std::make_unique<FaultFaceQuadRefiner>();
  case seissol::initializers::parameters::FaultRefinement::None:
  default:
    return nullptr;
  }
}

void FaultRefiner::repeatRefinement(Data data,
                                    PointsPair& point1,
                                    PointsPair& point2,
                                    PointsPair& point3) {
  const ExtTriangle subGlobalFace(
      std::get<global>(point1), std::get<global>(point2), std::get<global>(point3));

  const ExtTriangle subReferenceFace(
      std::get<reference>(point1), std::get<reference>(point2), std::get<reference>(point3));

  auto updatedData = data;
  updatedData.refinementLevel -= 1;
  refineAndAccumulate(updatedData, std::make_pair(subGlobalFace, subReferenceFace));
}

void FaultRefiner::addReceiver(Data data, TrianglePair& face) {
  ReceiverPoint receiver{};
  receiver.isInside = true;
  receiver.faultFaceIndex = data.faultFaceIndex;
  receiver.localFaceSideId = data.localFaceSideId;
  receiver.elementIndex = data.elementId;
  receiver.globalReceiverIndex = points.size();
  receiver.global = getMidPointTriangle(std::get<global>(face));
  receiver.reference = getMidPointTriangle(std::get<reference>(face));
  receiver.globalTriangle = std::get<global>(face);

  points.push_back(receiver);
}

void FaultFaceTripleRefiner::refineAndAccumulate(Data data, TrianglePair face) {

  if (data.refinementLevel == 0) {
    addReceiver(data, face);
    return;
  }

  auto& globalFace = std::get<global>(face);
  auto& referenceFace = std::get<reference>(face);

  auto midPoint =
      std::make_pair(getMidPointTriangle(globalFace), getMidPointTriangle(referenceFace));
  std::array<PointsPair, 3> points{};
  for (size_t i = 0; i < 3; ++i) {
    points[i] = std::make_pair(globalFace.point(i), referenceFace.point(i));
  }

  repeatRefinement(data, points[0], points[1], midPoint);
  repeatRefinement(data, midPoint, points[1], points[2]);
  repeatRefinement(data, points[0], midPoint, points[2]);
}

void FaultFaceQuadRefiner::refineAndAccumulate(Data data, TrianglePair face) {

  if (data.refinementLevel == 0) {
    addReceiver(data, face);
    return;
  }

  auto& globalFace = std::get<global>(face);
  auto& referenceFace = std::get<reference>(face);

  auto split = [&globalFace, &referenceFace](size_t pointIndex1, size_t pointIndex2) {
    return std::make_pair(
        getMidPoint(globalFace.point(pointIndex1), globalFace.point(pointIndex2)),
        getMidPoint(referenceFace.point(pointIndex1), referenceFace.point(pointIndex2)));
  };

  auto midPoint1 = split(0, 1);
  auto midPoint2 = split(1, 2);
  auto midPoint3 = split(2, 0);

  PointsPair trianglePoint = std::make_pair(globalFace.point(0), referenceFace.point(0));
  repeatRefinement(data, trianglePoint, midPoint1, midPoint3);

  trianglePoint = std::make_pair(globalFace.point(1), referenceFace.point(1));
  repeatRefinement(data, midPoint1, trianglePoint, midPoint2);

  repeatRefinement(data, midPoint1, midPoint2, midPoint3);

  trianglePoint = std::make_pair(globalFace.point(2), referenceFace.point(2));
  repeatRefinement(data, midPoint3, midPoint2, trianglePoint);
}
} // namespace seissol::dr::output::refiner
