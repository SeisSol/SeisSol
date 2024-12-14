#ifndef SEISSOL_DR_RECEIVER_BASED_OUTPUT_BUILDER_HPP
#define SEISSOL_DR_RECEIVER_BASED_OUTPUT_BUILDER_HPP

#include "DynamicRupture/Misc.h"
#include "DynamicRupture/Output/DataTypes.h"
#include "DynamicRupture/Output/OutputAux.h"
#include "Geometry/MeshReader.h"
#include "Initializer/DynamicRupture.h"
#include "Initializer/InputAux.h"
#include "Initializer/LTS.h"
#include "Initializer/Tree/LTSTree.h"
#include "Initializer/Tree/Lut.h"
#include "Kernels/Precision.h"
#include "Model/Common.h"
#include "Numerical/Transformation.h"
#include "Parallel/MPI.h"
#include <vector>

namespace seissol::dr::output {
class ReceiverBasedOutputBuilder {
  public:
  virtual ~ReceiverBasedOutputBuilder() = default;
  virtual void build(std::shared_ptr<ReceiverOutputData> outputData) = 0;

  void setMeshReader(const seissol::geometry::MeshReader* reader);
  void setLtsData(seissol::initializer::LTSTree* userWpTree,
                  seissol::initializer::LTS* userWpDescr,
                  seissol::initializer::Lut* userWpLut,
                  seissol::initializer::LTSTree* userDrTree,
                  seissol::initializer::DynamicRupture* userDrDescr);

  void setVariableList(const std::vector<std::size_t>& variables);
  void setFaceToLtsMap(std::vector<std::size_t>* faceToLtsMap);

  protected:
  virtual void initTimeCaching() = 0;

  void initBasisFunctions();
  void initFaultDirections();
  void initRotationMatrices();
  void initOutputVariables(std::array<bool, std::tuple_size<DrVarsT>::value>& outputMask);
  void initJacobian2dMatrices();
  void assignNearestInternalGaussianPoints();
  void assignFaultTags();

  protected:
  const seissol::geometry::MeshReader* meshReader{};
  seissol::initializer::LTSTree* wpTree;
  seissol::initializer::LTS* wpDescr;
  seissol::initializer::Lut* wpLut;
  seissol::initializer::LTSTree* drTree;
  seissol::initializer::DynamicRupture* drDescr;
  std::shared_ptr<ReceiverOutputData> outputData;
  std::vector<std::size_t> variables;
  std::vector<std::size_t>* faceToLtsMap{nullptr};
  int localRank{-1};
};
} // namespace seissol::dr::output

#endif // SEISSOL_DR_RECEIVER_BASED_OUTPUT_BUILDER_HPP
