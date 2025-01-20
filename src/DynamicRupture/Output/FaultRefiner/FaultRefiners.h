// SPDX-FileCopyrightText: 2021-2024 SeisSol Group
//
// SPDX-License-Identifier: BSD-3-Clause
// SPDX-LicenseComments: Full text under /LICENSE and /LICENSES/
//
// SPDX-FileContributor: Author lists in /AUTHORS and /CITATION.cff

#ifndef SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_FAULTREFINER_FAULTREFINERS_H_
#define SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_FAULTREFINER_FAULTREFINERS_H_

#include "DynamicRupture/Output/DataTypes.h"
#include "Initializer/Parameters/OutputParameters.h"
#include <memory>
#include <tuple>

namespace seissol::dr::output::refiner {
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

  [[nodiscard]] virtual int getNumSubTriangles() const = 0;
  virtual void refineAndAccumulate(Data data, TrianglePair face) = 0;
  virtual ~FaultRefiner() = default;

  ReceiverPoints&& moveAllReceiverPoints() { return std::move(points); }

  protected:
  ReceiverPoints points;

  static constexpr size_t Global = 0;
  static constexpr size_t Reference = 1;

  inline void
      repeatRefinement(Data data, PointsPair& point1, PointsPair& point2, PointsPair& point3);
  inline void addReceiver(Data data, TrianglePair& face);
};

class NoRefiner : public FaultRefiner {
  public:
  [[nodiscard]] int getNumSubTriangles() const final { return 1; }
  void refineAndAccumulate(Data data, TrianglePair face) final;
};

class FaultFaceTripleRefiner : public FaultRefiner {
  public:
  [[nodiscard]] int getNumSubTriangles() const final { return 3; }
  void refineAndAccumulate(Data data, TrianglePair face) final;
};

class FaultFaceQuadRefiner : public FaultRefiner {
  public:
  [[nodiscard]] int getNumSubTriangles() const final { return 4; }
  void refineAndAccumulate(Data data, TrianglePair face) final;
};

std::unique_ptr<FaultRefiner> get(seissol::initializer::parameters::FaultRefinement strategy);
} // namespace seissol::dr::output::refiner

#endif // SEISSOL_SRC_DYNAMICRUPTURE_OUTPUT_FAULTREFINER_FAULTREFINERS_H_
