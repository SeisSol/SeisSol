#ifndef SEISSOL_DR_RECEIVER_BASED_OUTPUT_BUILDER_HPP
#define SEISSOL_DR_RECEIVER_BASED_OUTPUT_BUILDER_HPP

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.hpp"
#include "DynamicRupture/Output/OutputAux.hpp"
#include "Geometry/MeshReader.h"
#include "Initializer/InputAux.hpp"
#include "Model/common.hpp"
#include "Numerical_aux/Transformation.h"
#include "Parallel/MPI.h"

namespace seissol::dr::output {
class ReceiverBasedOutputBuilder {
  public:
  virtual ~ReceiverBasedOutputBuilder() = default;
  virtual void build(std::shared_ptr<ReceiverOutputData> outputData) = 0;

  void setMeshReader(const MeshReader* reader);

  protected:
  virtual void initTimeCaching() = 0;

  void initBasisFunctions();
  void initFaultDirections();
  void initRotationMatrices();
  void initOutputVariables(std::array<bool, std::tuple_size<DrVarsT>::value>& outputMask);
  void initJacobian2dMatrices();
  void assignNearestInternalGaussianPoints();

  protected:
  const MeshReader* meshReader{};
  std::shared_ptr<ReceiverOutputData> outputData;
  int localRank{-1};
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_RECEIVER_BASED_OUTPUT_BUILDER_HPP
