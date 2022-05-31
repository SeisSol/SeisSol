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
    int elementId{-1};
  };
  using PointsPair = std::pair<ExtVrtxCoords, ExtVrtxCoords>;
  using TrianglePair = std::pair<ExtTriangle, ExtTriangle>;

  virtual int getNumSubTriangles() = 0;
  virtual void refineAndAccumulate(Data data, TrianglePair face) = 0;
  virtual ~FaultRefiner(){};

  ReceiverPointsT&& moveAllReceiverPoints() { return std::move(points); }
  ReceiverPointsT getAllReceiverPoints() { return points; }

  protected:
  ReceiverPointsT points{};

  static constexpr size_t global = 0;
  static constexpr size_t reference = 1;

  inline void repeat(Data data, PointsPair& point1, PointsPair& point2, PointsPair& point3);
  inline void addReceiver(Data data, TrianglePair& face);
};

class TripleFaultFaceRefiner : public FaultRefiner {
  public:
  int getNumSubTriangles() final { return 3; }
  void refineAndAccumulate(Data data, TrianglePair face) final;
};

class QuadFaultFaceRefiner : public FaultRefiner {
  public:
  int getNumSubTriangles() final { return 4; }
  void refineAndAccumulate(Data data, TrianglePair face) final;
};
} // namespace seissol::dr::output::refiner
#endif // SEISSOL_DR_OUTPUT_REFINERS_HPP
