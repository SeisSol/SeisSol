#ifndef SEISSOL_INTERFACE_HPP
#define SEISSOL_INTERFACE_HPP

#include "DynamicRupture/Output/DataTypes.hpp"
#include <memory>

namespace seissol::dr::output {
class FaultRefinerInterface {
  public:
  virtual int getNumSubTriangles() = 0;
  virtual void refineAndAccumulate(int refinementLevel,
                                   int faultFaceIndex,
                                   int localFaceSideId,
                                   ExtTriangle referenceFace,
                                   ExtTriangle globalFace) = 0;

  ReceiverPointsT&& moveAllReceiverPoints() { return std::move(points); }
  ReceiverPointsT getAllReceiverPoints() { return points; }

  protected:
  ReceiverPointsT points{};
};
} // namespace seissol::dr::output
#endif // SEISSOL_INTERFACE_HPP
