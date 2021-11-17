#ifndef SEISSOL_QUADFAULTFACEREFINER_HPP
#define SEISSOL_QUADFAULTFACEREFINER_HPP

#include "FaultRefinerInterface.hpp"
#include "Numerical_aux/Transformation.h"
#include "DynamicRupture/Output/OutputAux.hpp"

namespace seissol::dr::output {
class QuadFaultFaceRefiner : public FaultRefinerInterface {
  int getNumSubTriangles() override { return 4; }
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

    ExtVrtxCoords refMidPoint1 = getMidPoint(referenceFace.p1, referenceFace.p2);
    ExtVrtxCoords refMidPoint2 = getMidPoint(referenceFace.p2, referenceFace.p3);
    ExtVrtxCoords refMidPoint3 = getMidPoint(referenceFace.p3, referenceFace.p1);

    ExtVrtxCoords glbMidPoint1 = getMidPoint(globalFace.p1, globalFace.p2);
    ExtVrtxCoords glbMidPoint2 = getMidPoint(globalFace.p2, globalFace.p3);
    ExtVrtxCoords glbMidPoint3 = getMidPoint(globalFace.p3, globalFace.p1);

    {
      // First sub-face (-triangle)
      ExtTriangle subReferenceFace(referenceFace.p1, refMidPoint1, refMidPoint3);
      ExtTriangle subGlobalFace(globalFace.p1, glbMidPoint1, glbMidPoint3);
      this->refineAndAccumulate(
          refinementLevel - 1, faultFaceIndex, localFaceSideId, subReferenceFace, subGlobalFace);
    }

    {
      // Second sub-face (-triangle)
      ExtTriangle subReferenceFace(refMidPoint1, referenceFace.p2, refMidPoint2);
      ExtTriangle subGlobalFace(glbMidPoint1, globalFace.p2, glbMidPoint2);
      this->refineAndAccumulate(
          refinementLevel - 1, faultFaceIndex, localFaceSideId, subReferenceFace, subGlobalFace);
    }

    {
      // Third sub-face (-triangle)
      ExtTriangle subReferenceFace(refMidPoint1, refMidPoint2, refMidPoint3);
      ExtTriangle subGlobalFace(glbMidPoint1, glbMidPoint2, glbMidPoint3);
      this->refineAndAccumulate(
          refinementLevel - 1, faultFaceIndex, localFaceSideId, subReferenceFace, subGlobalFace);
    }

    {
      // Fourth sub-face (-triangle)
      ExtTriangle subReferenceFace(refMidPoint3, refMidPoint2, referenceFace.p3);
      ExtTriangle subGlobalFace(glbMidPoint3, glbMidPoint2, globalFace.p3);
      this->refineAndAccumulate(
          refinementLevel - 1, faultFaceIndex, localFaceSideId, subReferenceFace, subGlobalFace);
    }
  }
};
} // namespace seissol::dr::output

#endif // SEISSOL_QUADFAULTFACEREFINER_HPP
