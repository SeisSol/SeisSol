#ifndef SEISSOL_TENSOR_H_
#define SEISSOL_TENSOR_H_
namespace seissol {
  namespace tensor {
    struct AminusT {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct AplusT {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct I {
      constexpr static unsigned const Shape[2] = {56, 9};
      constexpr static unsigned const Size = 504;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct INodal {
      constexpr static unsigned const Shape[2] = {21, 9};
      constexpr static unsigned const Size = 216;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct INodalUpdate {
      constexpr static unsigned const Shape[2] = {21, 9};
      constexpr static unsigned const Size = 216;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct M2inv {
      constexpr static unsigned const Shape[2] = {21, 21};
      constexpr static unsigned const Size = 441;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct M3inv {
      constexpr static unsigned const Shape[2] = {56, 56};
      constexpr static unsigned const Size = 3136;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct Q {
      constexpr static unsigned const Shape[2] = {56, 9};
      constexpr static unsigned const Size = 504;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QAtPoint {
      constexpr static unsigned const Shape[1] = {9};
      constexpr static unsigned const Size = 9;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QEtaModal {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QEtaNodal {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QFortran {
      constexpr static unsigned const Shape[2] = {56, 9};
      constexpr static unsigned const Size = 504;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QInterpolated {
      constexpr static unsigned const Shape[2] = {49, 9};
      constexpr static unsigned const Size = 468;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QInterpolatedMinus {
      constexpr static unsigned const Shape[2] = {49, 9};
      constexpr static unsigned const Size = 468;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QInterpolatedPlus {
      constexpr static unsigned const Shape[2] = {49, 9};
      constexpr static unsigned const Size = 468;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QStress {
      constexpr static unsigned const Shape[2] = {56, 6};
      constexpr static unsigned const Size = 336;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QStressNodal {
      constexpr static unsigned const Shape[2] = {56, 6};
      constexpr static unsigned const Size = 336;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QgodLocal {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct QgodNeighbor {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct T {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct Tinv {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct TinvT {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct V2nTo2JacobiQuad {
      constexpr static unsigned const Shape[2] = {49, 21};
      constexpr static unsigned const Size = 1029;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct V3mTo2n {
      constexpr static unsigned const Shape[][2] = {{49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}, {49, 56}};
      constexpr static unsigned const Size[] = {2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912, 2912};
      constexpr static unsigned index(unsigned i0, unsigned i1) {
        return 1*i0 + 4*i1;
      }
      constexpr static unsigned size(unsigned i0, unsigned i1) {
        return Size[index(i0, i1)];
      }
      template<typename T>
      struct Container {
        T data[16];
        Container() : data{} {}
        inline T& operator()(unsigned i0, unsigned i1) {
          return data[index(i0, i1)];
        }
        inline T const& operator()(unsigned i0, unsigned i1) const {
          return data[index(i0, i1)];
        }
      };
    };
    struct V3mTo2nTWDivM {
      constexpr static unsigned const Shape[][2] = {{56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}, {56, 49}};
      constexpr static unsigned const Size[] = {2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744, 2744};
      constexpr static unsigned index(unsigned i0, unsigned i1) {
        return 1*i0 + 4*i1;
      }
      constexpr static unsigned size(unsigned i0, unsigned i1) {
        return Size[index(i0, i1)];
      }
      template<typename T>
      struct Container {
        T data[16];
        Container() : data{} {}
        inline T& operator()(unsigned i0, unsigned i1) {
          return data[index(i0, i1)];
        }
        inline T const& operator()(unsigned i0, unsigned i1) const {
          return data[index(i0, i1)];
        }
      };
    };
    struct averageNormalDisplacement {
      constexpr static unsigned const Shape[1] = {21};
      constexpr static unsigned const Size = 24;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct basisFunctionsAtPoint {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct christoffel {
      constexpr static unsigned const Shape[2] = {3, 3};
      constexpr static unsigned const Size = 9;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct dQ {
      constexpr static unsigned const Shape[][2] = {{56, 9}, {56, 9}, {56, 9}, {56, 9}, {56, 9}, {56, 9}};
      constexpr static unsigned const Size[] = {504, 324, 180, 108, 36, 36};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[6];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct direction {
      constexpr static unsigned const Shape[1] = {3};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct displacementRotationMatrix {
      constexpr static unsigned const Shape[2] = {3, 3};
      constexpr static unsigned const Size = 12;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct dofsQP {
      constexpr static unsigned const Shape[2] = {343, 9};
      constexpr static unsigned const Size = 3096;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct easiBoundaryConstant {
      constexpr static unsigned const Shape[2] = {9, 21};
      constexpr static unsigned const Size = 189;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct easiBoundaryMap {
      constexpr static unsigned const Shape[3] = {9, 9, 21};
      constexpr static unsigned const Size = 1701;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct easiIdentMap {
      constexpr static unsigned const Shape[3] = {9, 9, 21};
      constexpr static unsigned const Size = 1701;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct evalAtQP {
      constexpr static unsigned const Shape[2] = {343, 56};
      constexpr static unsigned const Size = 19264;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct fMrT {
      constexpr static unsigned const Shape[][2] = {{21, 56}, {21, 56}, {21, 56}, {21, 56}};
      constexpr static unsigned const Size[] = {1344, 1344, 1344, 1344};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct fP {
      constexpr static unsigned const Shape[][2] = {{21, 21}, {21, 21}, {21, 21}};
      constexpr static unsigned const Size[] = {504, 21, 504};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct faceDisplacement {
      constexpr static unsigned const Shape[2] = {21, 3};
      constexpr static unsigned const Size = 72;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct fluxSolver {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct frictionalEnergy {
      constexpr static unsigned const Size = 1;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct identityT {
      constexpr static unsigned const Shape[2] = {9, 9};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct iniCond {
      constexpr static unsigned const Shape[2] = {343, 9};
      constexpr static unsigned const Size = 3096;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct initialLoading {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 6;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct kDivM {
      constexpr static unsigned const Shape[][2] = {{56, 56}, {56, 56}, {56, 56}};
      constexpr static unsigned const Size[] = {1960, 1960, 1960};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct kDivMT {
      constexpr static unsigned const Shape[][2] = {{56, 56}, {56, 56}, {56, 56}};
      constexpr static unsigned const Size[] = {1980, 1980, 1980};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct mInvJInvPhisAtSources {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct mNormal {
      constexpr static unsigned const Shape[1] = {3};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct mSlip {
      constexpr static unsigned const Shape[1] = {3};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct meanStress {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct momentFSRM {
      constexpr static unsigned const Shape[1] = {9};
      constexpr static unsigned const Size = 9;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct momentToNRF {
      constexpr static unsigned const Shape[3] = {9, 3, 3};
      constexpr static unsigned const Size = 54;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct project2nFaceTo3m {
      constexpr static unsigned const Shape[][2] = {{56, 21}, {56, 21}, {56, 21}, {56, 21}};
      constexpr static unsigned const Size[] = {1176, 1176, 1176, 1176};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct projectQP {
      constexpr static unsigned const Shape[2] = {56, 343};
      constexpr static unsigned const Size = 19208;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct rDivM {
      constexpr static unsigned const Shape[][2] = {{56, 21}, {56, 21}, {56, 21}, {56, 21}};
      constexpr static unsigned const Size[] = {192, 1176, 1176, 1176};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct rT {
      constexpr static unsigned const Shape[][2] = {{21, 56}, {21, 56}, {21, 56}, {21, 56}};
      constexpr static unsigned const Size[] = {1344, 1344, 1344, 1344};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct replicateInitialLoading {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct resample {
      constexpr static unsigned const Shape[2] = {49, 49};
      constexpr static unsigned const Size = 2401;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct rotatedFaceDisplacement {
      constexpr static unsigned const Shape[2] = {21, 3};
      constexpr static unsigned const Size = 72;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct rotatedFaceDisplacementAtQuadratureNodes {
      constexpr static unsigned const Shape[2] = {49, 3};
      constexpr static unsigned const Size = 156;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct samplingDirections {
      constexpr static unsigned const Shape[2] = {200, 3};
      constexpr static unsigned const Size = 600;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct secondInvariant {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct selectBulkAverage {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct selectBulkNegative {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct selectVelocity {
      constexpr static unsigned const Shape[2] = {9, 3};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct slipInterpolated {
      constexpr static unsigned const Shape[2] = {49, 3};
      constexpr static unsigned const Size = 147;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct slipRateInterpolated {
      constexpr static unsigned const Shape[2] = {49, 3};
      constexpr static unsigned const Size = 147;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct spaceTimePredictor {
      constexpr static unsigned const Shape[3] = {56, 9, 6};
      constexpr static unsigned const Size = 3024;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct spaceWeights {
      constexpr static unsigned const Shape[1] = {49};
      constexpr static unsigned const Size = 49;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct squaredNormSlipRateInterpolated {
      constexpr static unsigned const Shape[1] = {49};
      constexpr static unsigned const Size = 49;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct star {
      constexpr static unsigned const Shape[][2] = {{9, 9}, {9, 9}, {9, 9}};
      constexpr static unsigned const Size[] = {81, 81, 81};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[3];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct stiffnessTensor {
      constexpr static unsigned const Shape[4] = {3, 3, 3, 3};
      constexpr static unsigned const Size = 81;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct subTriangleDofs {
      constexpr static unsigned const Shape[][2] = {{1, 3}, {4, 3}, {16, 3}, {64, 3}};
      constexpr static unsigned const Size[] = {12, 12, 48, 192};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct subTriangleProjection {
      constexpr static unsigned const Shape[][2] = {{1, 56}, {4, 56}, {16, 56}, {64, 56}};
      constexpr static unsigned const Size[] = {224, 224, 896, 3584};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct subTriangleProjectionFromFace {
      constexpr static unsigned const Shape[][2] = {{1, 21}, {4, 21}, {16, 21}, {64, 21}};
      constexpr static unsigned const Size[] = {84, 84, 336, 1344};
      constexpr static unsigned index(unsigned i0) {
        return 1*i0;
      }
      constexpr static unsigned size(unsigned i0) {
        return Size[index(i0)];
      }
      template<typename T>
      struct Container {
        T data[4];
        Container() : data{} {}
        inline T& operator()(unsigned i0) {
          return data[index(i0)];
        }
        inline T const& operator()(unsigned i0) const {
          return data[index(i0)];
        }
      };
    };
    struct timeBasisFunctionsAtPoint {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 6;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct tractionInterpolated {
      constexpr static unsigned const Shape[2] = {49, 3};
      constexpr static unsigned const Size = 147;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct tractionMinusMatrix {
      constexpr static unsigned const Shape[2] = {9, 3};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct tractionPlusMatrix {
      constexpr static unsigned const Shape[2] = {9, 3};
      constexpr static unsigned const Size = 3;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct v {
      constexpr static unsigned const Shape[2] = {56, 56};
      constexpr static unsigned const Size = 3136;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct vInv {
      constexpr static unsigned const Shape[2] = {56, 56};
      constexpr static unsigned const Size = 3136;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct weightSecondInvariant {
      constexpr static unsigned const Shape[1] = {6};
      constexpr static unsigned const Size = 6;
      constexpr static unsigned size() {
        return Size;
      }
    };
    struct yieldFactor {
      constexpr static unsigned const Shape[1] = {56};
      constexpr static unsigned const Size = 56;
      constexpr static unsigned size() {
        return Size;
      }
    };
  } // namespace tensor
  namespace nodal {
    namespace tensor {
      struct MV2nTo2m {
        constexpr static unsigned const Shape[2] = {21, 21};
        constexpr static unsigned const Size = 504;
        constexpr static unsigned size() {
          return Size;
        }
      };
      struct V3mTo2nFace {
        constexpr static unsigned const Shape[][2] = {{21, 56}, {21, 56}, {21, 56}, {21, 56}};
        constexpr static unsigned const Size[] = {1344, 1344, 1344, 1344};
        constexpr static unsigned index(unsigned i0) {
          return 1*i0;
        }
        constexpr static unsigned size(unsigned i0) {
          return Size[index(i0)];
        }
        template<typename T>
        struct Container {
          T data[4];
          Container() : data{} {}
          inline T& operator()(unsigned i0) {
            return data[index(i0)];
          }
          inline T const& operator()(unsigned i0) const {
            return data[index(i0)];
          }
        };
      };
      struct nodes2D {
        constexpr static unsigned const Shape[2] = {21, 2};
        constexpr static unsigned const Size = 48;
        constexpr static unsigned size() {
          return Size;
        }
      };
    } // namespace tensor
  } // namespace nodal
} // namespace seissol
#endif
