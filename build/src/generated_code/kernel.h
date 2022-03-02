#ifndef SEISSOL_KERNEL_H_
#define SEISSOL_KERNEL_H_
#include <cmath>
#include <limits>
#include "yateto.h"
#include "tensor.h"
namespace seissol {
  namespace kernel {
    struct computeFluxSolverLocal {
      constexpr static unsigned long const NonZeroFlops = 3186;
      constexpr static unsigned long const HardwareFlops = 5103;
      constexpr static unsigned long const TmpMemRequiredInBytes = 1296;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 1296;

      double fluxScale = std::numeric_limits<double>::signaling_NaN();
      double* AplusT{};
      double const* QgodLocal{};
      double const* T{};
      double const* Tinv{};
      tensor::star::Container<double const*> star;


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct computeFluxSolverNeighbor {
      constexpr static unsigned long const NonZeroFlops = 3186;
      constexpr static unsigned long const HardwareFlops = 5103;
      constexpr static unsigned long const TmpMemRequiredInBytes = 1296;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 1296;

      double fluxScale = std::numeric_limits<double>::signaling_NaN();
      double* AminusT{};
      double const* QgodNeighbor{};
      double const* T{};
      double const* Tinv{};
      tensor::star::Container<double const*> star;


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct copyQToQFortran {
      constexpr static unsigned long const NonZeroFlops = 0;
      constexpr static unsigned long const HardwareFlops = 0;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* Q{};
      double* QFortran{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct computeChristoffel {
      constexpr static unsigned long const NonZeroFlops = 162;
      constexpr static unsigned long const HardwareFlops = 171;
      constexpr static unsigned long const TmpMemRequiredInBytes = 72;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 72;

      double* christoffel{};
      double const* direction{};
      double const* stiffnessTensor{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct projectIniCond {
      constexpr static unsigned long const NonZeroFlops = 345240;
      constexpr static unsigned long const HardwareFlops = 345744;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double* Q{};
      double const* iniCond{};
      double const* projectQP{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct evalAtQP {
      constexpr static unsigned long const NonZeroFlops = 342657;
      constexpr static unsigned long const HardwareFlops = 346752;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* Q{};
      double* dofsQP{};
      double const* evalAtQP{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct volume {
      constexpr static unsigned long const NonZeroFlops = 34839;
      constexpr static unsigned long const HardwareFlops = 123336;
      constexpr static unsigned long const TmpMemRequiredInBytes = 2592;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 2592;

      double const* I{};
      double* Q{};
      tensor::kDivM::Container<double const*> kDivM;
      tensor::star::Container<double const*> star;


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plConvertToNodal {
      constexpr static unsigned long const NonZeroFlops = 33720;
      constexpr static unsigned long const HardwareFlops = 38304;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* QStress{};
      double* QStressNodal{};
      double const* initialLoading{};
      double const* replicateInitialLoading{};
      double const* v{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plConvertToNodalNoLoading {
      constexpr static unsigned long const NonZeroFlops = 33048;
      constexpr static unsigned long const HardwareFlops = 37632;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* QStress{};
      double* QStressNodal{};
      double const* v{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plConvertEtaModal2Nodal {
      constexpr static unsigned long const NonZeroFlops = 5508;
      constexpr static unsigned long const HardwareFlops = 6272;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* QEtaModal{};
      double* QEtaNodal{};
      double const* v{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plConvertEtaNodal2Modal {
      constexpr static unsigned long const NonZeroFlops = 5190;
      constexpr static unsigned long const HardwareFlops = 6272;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double* QEtaModal{};
      double const* QEtaNodal{};
      double const* vInv{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plComputeMean {
      constexpr static unsigned long const NonZeroFlops = 280;
      constexpr static unsigned long const HardwareFlops = 336;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* QStressNodal{};
      double* meanStress{};
      double const* selectBulkAverage{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plSubtractMean {
      constexpr static unsigned long const NonZeroFlops = 336;
      constexpr static unsigned long const HardwareFlops = 336;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double* QStressNodal{};
      double const* meanStress{};
      double const* selectBulkNegative{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plComputeSecondInvariant {
      constexpr static unsigned long const NonZeroFlops = 952;
      constexpr static unsigned long const HardwareFlops = 1008;
      constexpr static unsigned long const TmpMemRequiredInBytes = 2688;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 2688;

      double const* QStressNodal{};
      double* secondInvariant{};
      double const* weightSecondInvariant{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct plAdjustStresses {
      constexpr static unsigned long const NonZeroFlops = 31812;
      constexpr static unsigned long const HardwareFlops = 37968;
      constexpr static unsigned long const TmpMemRequiredInBytes = 2688;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 2688;

      double* QStress{};
      double const* QStressNodal{};
      double const* vInv{};
      double const* yieldFactor{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct createEasiBoundaryGhostCells {
      constexpr static unsigned long const NonZeroFlops = 3591;
      constexpr static unsigned long const HardwareFlops = 6804;
      constexpr static unsigned long const TmpMemRequiredInBytes = 1512;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 1512;

      double* INodal{};
      double const* easiBoundaryConstant{};
      double const* easiBoundaryMap{};
      double const* easiIdentMap{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct updateINodal {
      constexpr static unsigned long const NonZeroFlops = 378;
      constexpr static unsigned long const HardwareFlops = 378;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double factor = std::numeric_limits<double>::signaling_NaN();
      double* INodal{};
      double const* INodalUpdate{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct rotateFaceDisplacement {
      constexpr static unsigned long const NonZeroFlops = 315;
      constexpr static unsigned long const HardwareFlops = 432;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* displacementRotationMatrix{};
      double const* faceDisplacement{};
      double* rotatedFaceDisplacement{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct computeMInvJInvPhisAtSources {
      constexpr static unsigned long const NonZeroFlops = 112;
      constexpr static unsigned long const HardwareFlops = 9408;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double JInv = std::numeric_limits<double>::signaling_NaN();
      double const* M3inv{};
      double const* basisFunctionsAtPoint{};
      double* mInvJInvPhisAtSources{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct sourceNRF {
      constexpr static unsigned long const NonZeroFlops = 1125;
      constexpr static unsigned long const HardwareFlops = 1287;
      constexpr static unsigned long const TmpMemRequiredInBytes = 144;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 144;

      double mArea = std::numeric_limits<double>::signaling_NaN();
      double* Q{};
      double const* mInvJInvPhisAtSources{};
      double const* mNormal{};
      double const* mSlip{};
      double const* momentToNRF{};
      double const* stiffnessTensor{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct sourceFSRM {
      constexpr static unsigned long const NonZeroFlops = 1512;
      constexpr static unsigned long const HardwareFlops = 1512;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double stfIntegral = std::numeric_limits<double>::signaling_NaN();
      double* Q{};
      double const* mInvJInvPhisAtSources{};
      double const* momentFSRM{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct evaluateDOFSAtPoint {
      constexpr static unsigned long const NonZeroFlops = 999;
      constexpr static unsigned long const HardwareFlops = 1008;
      constexpr static unsigned long const TmpMemRequiredInBytes = 0;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* Q{};
      double* QAtPoint{};
      double const* basisFunctionsAtPoint{};


      void execute();
    };
  } // namespace kernel
  namespace kernel {
    struct evaluateDOFSAtPointSTP {
      constexpr static unsigned long const NonZeroFlops = 6093;
      constexpr static unsigned long const HardwareFlops = 6156;
      constexpr static unsigned long const TmpMemRequiredInBytes = 432;
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 432;

      double* QAtPoint{};
      double const* basisFunctionsAtPoint{};
      double const* spaceTimePredictor{};
      double const* timeBasisFunctionsAtPoint{};


      void execute();
    };
  } // namespace kernel
  namespace dynamicRupture {
    namespace kernel {
      struct transposeTinv {
        constexpr static unsigned long const NonZeroFlops = 0;
        constexpr static unsigned long const HardwareFlops = 0;
        constexpr static unsigned long const TmpMemRequiredInBytes = 0;
        constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

        double const* Tinv{};
        double* TinvT{};


        void execute();
      };
    } // namespace kernel
  } // namespace dynamicRupture
  namespace dynamicRupture {
    namespace kernel {
      struct rotateFluxMatrix {
        constexpr static unsigned long const NonZeroFlops = 432;
        constexpr static unsigned long const HardwareFlops = 2187;
        constexpr static unsigned long const TmpMemRequiredInBytes = 0;
        constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

        double fluxScale = std::numeric_limits<double>::signaling_NaN();
        double const* T{};
        double* fluxSolver{};
        tensor::star::Container<double const*> star;


        void execute();
      };
    } // namespace kernel
  } // namespace dynamicRupture
  namespace kernel {
    struct localFlux {
      constexpr static unsigned long const NonZeroFlops[] = {9936, 10080, 31968, 27216};
      constexpr static unsigned long const HardwareFlops[] = {31536, 49248, 49248, 49248};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {3456, 3456, 3456, 3456};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 3456;

      double const* AplusT{};
      double const* I{};
      double* Q{};
      tensor::fMrT::Container<double const*> fMrT;
      tensor::rDivM::Container<double const*> rDivM;


      struct Prefetch {
        double const* I{};
        double const* Q{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      using member_function_ptr = void (localFlux::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&localFlux::execute0, &localFlux::execute1, &localFlux::execute2, &localFlux::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct localFluxNodal {
      constexpr static unsigned long const NonZeroFlops[] = {24381, 22941, 24381, 24345};
      constexpr static unsigned long const HardwareFlops[] = {25056, 25056, 25056, 25056};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {1728, 1728, 1728, 1728};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 1728;

      double const* AminusT{};
      double const* INodal{};
      double* Q{};
      tensor::project2nFaceTo3m::Container<double const*> project2nFaceTo3m;


      struct Prefetch {
        double const* I{};
        double const* Q{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      using member_function_ptr = void (localFluxNodal::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&localFluxNodal::execute0, &localFluxNodal::execute1, &localFluxNodal::execute2, &localFluxNodal::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct neighboringFlux {
      constexpr static unsigned long const NonZeroFlops[] = {11349, 10125, 11349, 11421, 10197, 11421, 22365, 21141, 22365, 19989, 18765, 19989, 11421, 10197, 11421, 11493, 10269, 11493, 22437, 21213, 22437, 20061, 18837, 20061, 22365, 21141, 22365, 22437, 21213, 22437, 33381, 32157, 33381, 31005, 29781, 31005, 19989, 18765, 19989, 20061, 18837, 20061, 31005, 29781, 31005, 28629, 27405, 28629};
      constexpr static unsigned long const HardwareFlops[] = {40608, 31428, 40608, 40608, 31428, 40608, 40608, 31428, 40608, 40608, 31428, 40608, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320, 58320, 49140, 58320};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456, 3456, 3240, 3456};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 3456;

      double const* AminusT{};
      double const* I{};
      double* Q{};
      tensor::fP::Container<double const*> fP;
      tensor::rDivM::Container<double const*> rDivM;
      tensor::rT::Container<double const*> rT;


      struct Prefetch {
        double const* I{};
      };
      Prefetch _prefetch;

      void execute0();
      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      void execute6();
      void execute7();
      void execute8();
      void execute9();
      void execute10();
      void execute11();
      void execute12();
      void execute13();
      void execute14();
      void execute15();
      void execute16();
      void execute17();
      void execute18();
      void execute19();
      void execute20();
      void execute21();
      void execute22();
      void execute23();
      void execute24();
      void execute25();
      void execute26();
      void execute27();
      void execute28();
      void execute29();
      void execute30();
      void execute31();
      void execute32();
      void execute33();
      void execute34();
      void execute35();
      void execute36();
      void execute37();
      void execute38();
      void execute39();
      void execute40();
      void execute41();
      void execute42();
      void execute43();
      void execute44();
      void execute45();
      void execute46();
      void execute47();
      using member_function_ptr = void (neighboringFlux::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&neighboringFlux::execute0, &neighboringFlux::execute1, &neighboringFlux::execute2, &neighboringFlux::execute3, &neighboringFlux::execute4, &neighboringFlux::execute5, &neighboringFlux::execute6, &neighboringFlux::execute7, &neighboringFlux::execute8, &neighboringFlux::execute9, &neighboringFlux::execute10, &neighboringFlux::execute11, &neighboringFlux::execute12, &neighboringFlux::execute13, &neighboringFlux::execute14, &neighboringFlux::execute15, &neighboringFlux::execute16, &neighboringFlux::execute17, &neighboringFlux::execute18, &neighboringFlux::execute19, &neighboringFlux::execute20, &neighboringFlux::execute21, &neighboringFlux::execute22, &neighboringFlux::execute23, &neighboringFlux::execute24, &neighboringFlux::execute25, &neighboringFlux::execute26, &neighboringFlux::execute27, &neighboringFlux::execute28, &neighboringFlux::execute29, &neighboringFlux::execute30, &neighboringFlux::execute31, &neighboringFlux::execute32, &neighboringFlux::execute33, &neighboringFlux::execute34, &neighboringFlux::execute35, &neighboringFlux::execute36, &neighboringFlux::execute37, &neighboringFlux::execute38, &neighboringFlux::execute39, &neighboringFlux::execute40, &neighboringFlux::execute41, &neighboringFlux::execute42, &neighboringFlux::execute43, &neighboringFlux::execute44, &neighboringFlux::execute45, &neighboringFlux::execute46, &neighboringFlux::execute47};
      constexpr static member_function_ptr findExecute(unsigned i0, unsigned i1, unsigned i2) {
        return ExecutePtrs[1*i0 + 3*i1 + 12*i2];
      }
      inline void execute(unsigned i0, unsigned i1, unsigned i2) {
        (this->*findExecute(i0, i1, i2))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0, unsigned i1, unsigned i2) {
        return NonZeroFlops[1*i0 + 3*i1 + 12*i2];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0, unsigned i1, unsigned i2) {
        return HardwareFlops[1*i0 + 3*i1 + 12*i2];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0, unsigned i1, unsigned i2) {
        return TmpMemRequiredInBytes[1*i0 + 3*i1 + 12*i2];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct derivativeTaylorExpansion {
      constexpr static unsigned long const NonZeroFlops[] = {504, 630, 360, 180, 72, 18};
      constexpr static unsigned long const HardwareFlops[] = {504, 630, 360, 180, 72, 18};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {0, 0, 0, 0, 0, 0};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double power = std::numeric_limits<double>::signaling_NaN();
      double* I{};
      tensor::dQ::Container<double const*> dQ;


      void execute0();
      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      using member_function_ptr = void (derivativeTaylorExpansion::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&derivativeTaylorExpansion::execute0, &derivativeTaylorExpansion::execute1, &derivativeTaylorExpansion::execute2, &derivativeTaylorExpansion::execute3, &derivativeTaylorExpansion::execute4, &derivativeTaylorExpansion::execute5};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct derivative {
      constexpr static unsigned long const NonZeroFlops[] = {0, 34524, 13806, 4716, 1260, 216};
      constexpr static unsigned long const HardwareFlops[] = {0, 124416, 46440, 18144, 3888, 2592};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {0, 2592, 1440, 864, 288, 288};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 2592;

      tensor::dQ::Container<double*> dQ;
      tensor::kDivMT::Container<double const*> kDivMT;
      tensor::star::Container<double const*> star;


      void execute1();
      void execute2();
      void execute3();
      void execute4();
      void execute5();
      using member_function_ptr = void (derivative::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {nullptr, &derivative::execute1, &derivative::execute2, &derivative::execute3, &derivative::execute4, &derivative::execute5};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct projectToNodalBoundary {
      constexpr static unsigned long const NonZeroFlops[] = {19557, 19287, 16965, 17019};
      constexpr static unsigned long const HardwareFlops[] = {24192, 24192, 24192, 24192};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {0, 0, 0, 0};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 0;

      double const* I{};
      double* INodal{};
      nodal::tensor::V3mTo2nFace::Container<double const*> V3mTo2nFace;


      void execute0();
      void execute1();
      void execute2();
      void execute3();
      using member_function_ptr = void (projectToNodalBoundary::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&projectToNodalBoundary::execute0, &projectToNodalBoundary::execute1, &projectToNodalBoundary::execute2, &projectToNodalBoundary::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct projectToNodalBoundaryRotated {
      constexpr static unsigned long const NonZeroFlops[] = {22770, 22500, 20178, 20232};
      constexpr static unsigned long const HardwareFlops[] = {28080, 28080, 28080, 28080};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {1728, 1728, 1728, 1728};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 1728;

      double const* I{};
      double* INodal{};
      double const* Tinv{};
      nodal::tensor::V3mTo2nFace::Container<double const*> V3mTo2nFace;


      void execute0();
      void execute1();
      void execute2();
      void execute3();
      using member_function_ptr = void (projectToNodalBoundaryRotated::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&projectToNodalBoundaryRotated::execute0, &projectToNodalBoundaryRotated::execute1, &projectToNodalBoundaryRotated::execute2, &projectToNodalBoundaryRotated::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct subTriangleDisplacement {
      constexpr static unsigned long const NonZeroFlops[] = {906, 2841, 4317, 10221};
      constexpr static unsigned long const HardwareFlops[] = {4032, 3528, 5040, 11088};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {672, 576, 576, 576};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 672;

      double const* MV2nTo2m{};
      double const* faceDisplacement{};
      tensor::subTriangleDofs::Container<double*> subTriangleDofs;
      tensor::subTriangleProjectionFromFace::Container<double const*> subTriangleProjectionFromFace;


      void execute0();
      void execute1();
      void execute2();
      void execute3();
      using member_function_ptr = void (subTriangleDisplacement::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&subTriangleDisplacement::execute0, &subTriangleDisplacement::execute1, &subTriangleDisplacement::execute2, &subTriangleDisplacement::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct subTriangleVelocity {
      constexpr static unsigned long const NonZeroFlops[] = {336, 1344, 5376, 21480};
      constexpr static unsigned long const HardwareFlops[] = {1368, 1368, 5472, 21840};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {96, 96, 384, 1344};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 1344;

      double const* Q{};
      double const* selectVelocity{};
      tensor::subTriangleDofs::Container<double*> subTriangleDofs;
      tensor::subTriangleProjection::Container<double const*> subTriangleProjection;


      void execute0();
      void execute1();
      void execute2();
      void execute3();
      using member_function_ptr = void (subTriangleVelocity::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&subTriangleVelocity::execute0, &subTriangleVelocity::execute1, &subTriangleVelocity::execute2, &subTriangleVelocity::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace kernel {
    struct addVelocity {
      constexpr static unsigned long const NonZeroFlops[] = {6645, 6555, 5781, 5799};
      constexpr static unsigned long const HardwareFlops[] = {8208, 8208, 8208, 8208};
      constexpr static unsigned long const TmpMemRequiredInBytes[] = {576, 576, 576, 576};
      constexpr static unsigned long const TmpMaxMemRequiredInBytes = 576;

      double const* I{};
      nodal::tensor::V3mTo2nFace::Container<double const*> V3mTo2nFace;
      double* faceDisplacement{};
      double const* selectVelocity{};


      void execute0();
      void execute1();
      void execute2();
      void execute3();
      using member_function_ptr = void (addVelocity::*)();
      constexpr static member_function_ptr ExecutePtrs[] = {&addVelocity::execute0, &addVelocity::execute1, &addVelocity::execute2, &addVelocity::execute3};
      constexpr static member_function_ptr findExecute(unsigned i0) {
        return ExecutePtrs[1*i0];
      }
      inline void execute(unsigned i0) {
        (this->*findExecute(i0))();
      }
      constexpr static unsigned long nonZeroFlops(unsigned i0) {
        return NonZeroFlops[1*i0];
      }
      constexpr static unsigned long hardwareFlops(unsigned i0) {
        return HardwareFlops[1*i0];
      }
      constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0) {
        return TmpMemRequiredInBytes[1*i0];
      }
    };
  } // namespace kernel
  namespace dynamicRupture {
    namespace kernel {
      struct evaluateAndRotateQAtInterpolationPoints {
        constexpr static unsigned long const NonZeroFlops[] = {53676, 56448, 56448, 56448, 56448, 53676, 56448, 56448, 56448, 56448, 56448, 56448, 53676, 56448, 56448, 56448};
        constexpr static unsigned long const HardwareFlops[] = {60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840, 60840};
        constexpr static unsigned long const TmpMemRequiredInBytes[] = {3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744};
        constexpr static unsigned long const TmpMaxMemRequiredInBytes = 3744;

        double const* Q{};
        double* QInterpolated{};
        double const* TinvT{};
        tensor::V3mTo2n::Container<double const*> V3mTo2n;


        struct Prefetch {
          double const* QInterpolated{};
        };
        Prefetch _prefetch;

        void execute0();
        void execute1();
        void execute2();
        void execute3();
        void execute4();
        void execute5();
        void execute6();
        void execute7();
        void execute8();
        void execute9();
        void execute10();
        void execute11();
        void execute12();
        void execute13();
        void execute14();
        void execute15();
        using member_function_ptr = void (evaluateAndRotateQAtInterpolationPoints::*)();
        constexpr static member_function_ptr ExecutePtrs[] = {&evaluateAndRotateQAtInterpolationPoints::execute0, &evaluateAndRotateQAtInterpolationPoints::execute1, &evaluateAndRotateQAtInterpolationPoints::execute2, &evaluateAndRotateQAtInterpolationPoints::execute3, &evaluateAndRotateQAtInterpolationPoints::execute4, &evaluateAndRotateQAtInterpolationPoints::execute5, &evaluateAndRotateQAtInterpolationPoints::execute6, &evaluateAndRotateQAtInterpolationPoints::execute7, &evaluateAndRotateQAtInterpolationPoints::execute8, &evaluateAndRotateQAtInterpolationPoints::execute9, &evaluateAndRotateQAtInterpolationPoints::execute10, &evaluateAndRotateQAtInterpolationPoints::execute11, &evaluateAndRotateQAtInterpolationPoints::execute12, &evaluateAndRotateQAtInterpolationPoints::execute13, &evaluateAndRotateQAtInterpolationPoints::execute14, &evaluateAndRotateQAtInterpolationPoints::execute15};
        constexpr static member_function_ptr findExecute(unsigned i0, unsigned i1) {
          return ExecutePtrs[1*i0 + 4*i1];
        }
        inline void execute(unsigned i0, unsigned i1) {
          (this->*findExecute(i0, i1))();
        }
        constexpr static unsigned long nonZeroFlops(unsigned i0, unsigned i1) {
          return NonZeroFlops[1*i0 + 4*i1];
        }
        constexpr static unsigned long hardwareFlops(unsigned i0, unsigned i1) {
          return HardwareFlops[1*i0 + 4*i1];
        }
        constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0, unsigned i1) {
          return TmpMemRequiredInBytes[1*i0 + 4*i1];
        }
      };
    } // namespace kernel
  } // namespace dynamicRupture
  namespace dynamicRupture {
    namespace kernel {
      struct nodalFlux {
        constexpr static unsigned long const NonZeroFlops[] = {54117, 56889, 56889, 56889, 56889, 54117, 56889, 56889, 56889, 56889, 56889, 56889, 54117, 56889, 56889, 56889};
        constexpr static unsigned long const HardwareFlops[] = {57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816, 57816};
        constexpr static unsigned long const TmpMemRequiredInBytes[] = {3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744, 3744};
        constexpr static unsigned long const TmpMaxMemRequiredInBytes = 3744;

        double* Q{};
        double const* QInterpolated{};
        tensor::V3mTo2nTWDivM::Container<double const*> V3mTo2nTWDivM;
        double const* fluxSolver{};


        struct Prefetch {
          double const* I{};
        };
        Prefetch _prefetch;

        void execute0();
        void execute1();
        void execute2();
        void execute3();
        void execute4();
        void execute5();
        void execute6();
        void execute7();
        void execute8();
        void execute9();
        void execute10();
        void execute11();
        void execute12();
        void execute13();
        void execute14();
        void execute15();
        using member_function_ptr = void (nodalFlux::*)();
        constexpr static member_function_ptr ExecutePtrs[] = {&nodalFlux::execute0, &nodalFlux::execute1, &nodalFlux::execute2, &nodalFlux::execute3, &nodalFlux::execute4, &nodalFlux::execute5, &nodalFlux::execute6, &nodalFlux::execute7, &nodalFlux::execute8, &nodalFlux::execute9, &nodalFlux::execute10, &nodalFlux::execute11, &nodalFlux::execute12, &nodalFlux::execute13, &nodalFlux::execute14, &nodalFlux::execute15};
        constexpr static member_function_ptr findExecute(unsigned i0, unsigned i1) {
          return ExecutePtrs[1*i0 + 4*i1];
        }
        inline void execute(unsigned i0, unsigned i1) {
          (this->*findExecute(i0, i1))();
        }
        constexpr static unsigned long nonZeroFlops(unsigned i0, unsigned i1) {
          return NonZeroFlops[1*i0 + 4*i1];
        }
        constexpr static unsigned long hardwareFlops(unsigned i0, unsigned i1) {
          return HardwareFlops[1*i0 + 4*i1];
        }
        constexpr static unsigned long tmpMemRequiredInBytes(unsigned i0, unsigned i1) {
          return TmpMemRequiredInBytes[1*i0 + 4*i1];
        }
      };
    } // namespace kernel
  } // namespace dynamicRupture
} // namespace seissol
#endif
