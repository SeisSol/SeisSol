#ifndef SEISSOL_DR_OUTPUT_BASE_HPP
#define SEISSOL_DR_OUTPUT_BASE_HPP

#include "DynamicRupture/Output/ParametersInitializer.hpp"
#include "Initializer/tree/Lut.hpp"
#include "Initializer/LTS.h"
#include "Initializer/DynamicRupture.h"
#include "Geometry/MeshReader.h"
#include "Solver/Interoperability.h"

namespace seissol::dr::output {
class ReceiverBasedOutput {
  public:
  virtual ~ReceiverBasedOutput() = default;

  void setLtsData(seissol::initializers::LTSTree* userWpTree,
                  seissol::initializers::LTS* userWpDescr,
                  seissol::initializers::Lut* userWpLut,
                  seissol::initializers::LTSTree* userDrTree,
                  seissol::initializers::DynamicRupture* userDrDescr);

  void setMeshReader(MeshReader* userMeshReader) { meshReader = userMeshReader; }
  void setFaceToLtsMap(FaceToLtsMapT* map) { faceToLtsMap = map; }
  virtual void tiePointers(seissol::initializers::Layer& layerData,
                           seissol::initializers::DynamicRupture* description,
                           seissol::Interoperability& eInteroperability);
  void calcFaultOutput(OutputType type,
                       ReceiverBasedOutputData& state,
                       const GeneralParamsT& generalParams,
                       double time = 0.0);

  protected:
  void getDofs(real dofsPlus[tensor::Q::size()], int meshId, int side);
  void computeLocalStresses();
  virtual real computeLocalStrength() = 0;
  virtual real computePf() { return 0.0; }
  void computeLocalTraction(real strength);
  virtual void computeSlipAndRate(std::array<real, 6>&, std::array<real, 6>&);
  void computeSlipAndRate(const double* tangent1,
                          const double* tangent2,
                          const double* strike,
                          const double* dip);

  virtual void adjustRotatedTractionAndStresses(std::array<real, 6>& rotatedTraction,
                                                std::array<real, 6>& rotatedLocalStress){};

  int getClosestInternalGp(int nearestGpIndex, int nPoly);

  virtual void
      outputSpecifics(ReceiverBasedOutputData& data, size_t outputSpecifics, size_t receiverIdx) {}
  real computeRuptureVelocity(Eigen::Matrix<real, 2, 2>& jacobiT2d);

  seissol::initializers::LTS* wpDescr{nullptr};
  seissol::initializers::LTSTree* wpTree{nullptr};
  seissol::initializers::Lut* wpLut{nullptr};
  seissol::initializers::LTSTree* drTree{nullptr};
  seissol::initializers::DynamicRupture* drDescr{nullptr};
  MeshReader* meshReader{nullptr};
  FaceToLtsMapT* faceToLtsMap{nullptr};

  struct LocalInfo {
    seissol::initializers::Layer* layer{};
    size_t ltsId{};
    int nearestGpIndex{};

    real pf{};
    real mu{};
    real sXY{};
    real sXZ{};
    real p0{};

    real p{};
    real u{};
    real yyStress{};
    real zzStress{};
    real xyStress{};
    real xzStress{};
    real yzStress{};
    real tracEla{};

    real xyTraction{};
    real xzTraction{};

    real srS{};
    real srD{};

    real faceAlignedValuesPlus[tensor::QAtPoint::size()]{};
    real faceAlignedValuesMinus[tensor::QAtPoint::size()]{};

    model::IsotropicWaveSpeeds* waveSpeedsPlus{};
    model::IsotropicWaveSpeeds* waveSpeedsMinus{};
  } local{};
};
} // namespace seissol::dr::output
#endif // SEISSOL_DR_OUTPUT_BASE_HPP
