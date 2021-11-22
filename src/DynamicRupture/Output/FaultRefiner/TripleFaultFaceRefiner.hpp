#ifndef SEISSOL_TRIPLEFAULTFACEREFINER_HPP
#define SEISSOL_TRIPLEFAULTFACEREFINER_HPP

#include "FaultRefinerInterface.hpp"

namespace seissol::dr::output {
class TripleFaultFaceRefiner : public FaultRefinerInterface {
  int getNumSubTriangles() override { return 3; }
  void refineAndAccumulate(int refinementLevel,
                           int faultFaceIndex,
                           int localFaceSideId,
                           ExtTriangle referenceFace,
                           ExtTriangle globalFace) override {

    if (refinementLevel == 0) {
      ReceiverPointT receiver{};
      receiver.isInside = true;
      receiver.faultFaceIndex = faultFaceIndex;
      receiver.localFaceSideId = localFaceSideId;
      receiver.globalReceiverIndex = points.size();
      receiver.global = getMidTrianglePoint(globalFace);
      receiver.referece = getMidTrianglePoint(referenceFace);
      receiver.globalSubTet = globalFace;

      points.push_back(receiver);
      return;
    }

    auto refMidTriangle = getMidTrianglePoint(referenceFace);
    auto glbMidTriangle = getMidTrianglePoint(globalFace);

    {
      // First sub-face (-triangle)
      ExtTriangle subReferenceFace(referenceFace.p1, referenceFace.p2, refMidTriangle);
      ExtTriangle subGlobalFace(globalFace.p1, globalFace.p2, glbMidTriangle);
      this->refineAndAccumulate(
          refinementLevel - 1, faultFaceIndex, localFaceSideId, subReferenceFace, subGlobalFace);
    }

    {
      // Second sub-face (-triangle)
      ExtTriangle subReferenceFace(refMidTriangle, referenceFace.p2, referenceFace.p3);
      ExtTriangle subGlobalFace(glbMidTriangle, referenceFace.p2, referenceFace.p3);
      this->refineAndAccumulate(
          refinementLevel - 1, faultFaceIndex, localFaceSideId, subReferenceFace, subGlobalFace);
    }

    {
      // Third sub-face (-triangle)
      ExtTriangle subReferenceFace(referenceFace.p1, refMidTriangle, referenceFace.p3);
      ExtTriangle subGlobalFace(globalFace.p1, glbMidTriangle, globalFace.p3);
      this->refineAndAccumulate(
          refinementLevel - 1, faultFaceIndex, localFaceSideId, subReferenceFace, subGlobalFace);
    }
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_TRIPLEFAULTFACEREFINER_HPP
