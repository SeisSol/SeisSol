#ifndef SEISSOL_DR_OUTPUT_REFINERS_HPP
#define SEISSOL_DR_OUTPUT_REFINERS_HPP

#include "DynamicRupture/Output/DataTypes.hpp"
#include <memory>
#include <tuple>

namespace seissol::dr::output::refiner {
RefinerType convertToType(int strategy);

class FaultRefiner;
std::unique_ptr<FaultRefiner> get(RefinerType strategy);

class FaultRefiner {
  public:
  struct Data {
    int refinementLevel{};
    int faultFaceIndex{};
    int localFaceSideId{};
  };
  using PointsPair = std::pair<ExtVrtxCoords, ExtVrtxCoords>;
  using TrianglePair = std::pair<ExtTriangle, ExtTriangle>;

  virtual int getNumSubTriangles() = 0;
  virtual void refineAndAccumulate(Data data, TrianglePair Face) = 0;

  ReceiverPointsT&& moveAllReceiverPoints() { return std::move(points); }
  ReceiverPointsT getAllReceiverPoints() { return points; }

  protected:
  ReceiverPointsT points{};

  static constexpr size_t GLOBAL = 0;
  static constexpr size_t REFERENCE = 1;

  inline void repeat(Data data, PointsPair& point1, PointsPair& point2, PointsPair& point3);
  inline void addReceiver(Data data, TrianglePair& Face);
};

class TripleFaultFaceRefiner : public FaultRefiner {
  public:
  int getNumSubTriangles() final { return 3; }
  void refineAndAccumulate(Data data, TrianglePair Face) final;
};

class QuadFaultFaceRefiner : public FaultRefiner {
  public:
  int getNumSubTriangles() final { return 4; }
  void refineAndAccumulate(Data data, TrianglePair Face) final;
};
} // namespace seissol::dr::output::refiner
#endif // SEISSOL_DR_OUTPUT_REFINERS_HPP
