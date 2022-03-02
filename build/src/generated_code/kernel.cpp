#include <cassert>
#include <cstring>
#include <cstdlib>
#include <limits>
#include "subroutine.h"
#include "kernel.h"
namespace seissol {
  constexpr unsigned long const kernel::computeFluxSolverLocal::NonZeroFlops;
  constexpr unsigned long const kernel::computeFluxSolverLocal::HardwareFlops;
  void kernel::computeFluxSolverLocal::execute() {
    assert(!std::isnan(fluxScale));
    assert(AplusT != nullptr);
    assert(QgodLocal != nullptr);
    assert(T != nullptr);
    assert(Tinv != nullptr);
    assert(star(0) != nullptr);
    double *_tmp0, *_tmp1;
    alignas(32) double _buffer0[81] ;
    alignas(32) double _buffer1[81] ;
    _tmp0 = _buffer0;
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 9; ++m) {
        _tmp0[0 + 1*m + 9*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 9; ++m) {
          _tmp0[0 + 1*m + 9*n] += 1.0 * Tinv[0 + 9*m + 1*k] * QgodLocal[0 + 1*k + 9*n];
        }
      }
    }
    _tmp1 = _buffer1;
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 9; ++m) {
        _tmp1[0 + 1*m + 9*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 9; ++m) {
          _tmp1[0 + 1*m + 9*n] += 1.0 * star(0)[0 + 1*m + 9*k] * T[0 + 9*k + 1*n];
        }
      }
    }
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 9; ++m) {
        AplusT[0 + 1*m + 9*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 9; ++m) {
          AplusT[0 + 1*m + 9*n] += fluxScale * _tmp0[0 + 1*m + 9*k] * _tmp1[0 + 1*k + 9*n];
        }
      }
    }
  }
  constexpr unsigned long const kernel::computeFluxSolverNeighbor::NonZeroFlops;
  constexpr unsigned long const kernel::computeFluxSolverNeighbor::HardwareFlops;
  void kernel::computeFluxSolverNeighbor::execute() {
    assert(!std::isnan(fluxScale));
    assert(AminusT != nullptr);
    assert(QgodNeighbor != nullptr);
    assert(T != nullptr);
    assert(Tinv != nullptr);
    assert(star(0) != nullptr);
    double *_tmp0, *_tmp1;
    alignas(32) double _buffer0[81] ;
    alignas(32) double _buffer1[81] ;
    _tmp0 = _buffer0;
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 9; ++m) {
        _tmp0[0 + 1*m + 9*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 9; ++m) {
          _tmp0[0 + 1*m + 9*n] += 1.0 * Tinv[0 + 9*m + 1*k] * QgodNeighbor[0 + 1*k + 9*n];
        }
      }
    }
    _tmp1 = _buffer1;
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 9; ++m) {
        _tmp1[0 + 1*m + 9*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 9; ++m) {
          _tmp1[0 + 1*m + 9*n] += 1.0 * star(0)[0 + 1*m + 9*k] * T[0 + 9*k + 1*n];
        }
      }
    }
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 9; ++m) {
        AminusT[0 + 1*m + 9*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 9; ++m) {
          AminusT[0 + 1*m + 9*n] += fluxScale * _tmp0[0 + 1*m + 9*k] * _tmp1[0 + 1*k + 9*n];
        }
      }
    }
  }
  constexpr unsigned long const kernel::copyQToQFortran::NonZeroFlops;
  constexpr unsigned long const kernel::copyQToQFortran::HardwareFlops;
  void kernel::copyQToQFortran::execute() {
    assert(Q != nullptr);
    assert(QFortran != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 56; ++_a) {
        QFortran[1*_a + 56*_b] = Q[1*_a + 56*_b];
      }
    }
  }
  constexpr unsigned long const kernel::computeChristoffel::NonZeroFlops;
  constexpr unsigned long const kernel::computeChristoffel::HardwareFlops;
  void kernel::computeChristoffel::execute() {
    assert(christoffel != nullptr);
    assert(direction != nullptr);
    assert(stiffnessTensor != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[9] ;
    _tmp0 = _buffer0;
    for (int _l = 0; _l < 3; ++_l) {
      #pragma omp simd
      for (int _j = 0; _j < 3; ++_j) {
        _tmp0[1*_j + 3*_l] = direction[1*_j] * direction[1*_l];
      }
    }
    for (int _k = 0; _k < 3; ++_k) {
      double const* _A = stiffnessTensor + 9*_k;
      double const* _B = _tmp0;
      double * _C = christoffel + 3*_k;
      for (int _l = 0; _l < 1; ++_l) {
        double const* _Ain = _A + 27*_l;
        double const* _Bin = _B + 3*_l;
        double * _Cin = _C;
        libxsmm_m3_n1_k3_ldA3_ldB3_ldC3_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_Ain, _Bin, _Cin, nullptr, nullptr, nullptr);
      }
      for (int _l = 1; _l < 3; ++_l) {
        double const* _Ain = _A + 27*_l;
        double const* _Bin = _B + 3*_l;
        double * _Cin = _C;
        libxsmm_m3_n1_k3_ldA3_ldB3_ldC3_alpha1_beta1_alignedA0_alignedC0_pfsigonly(_Ain, _Bin, _Cin, nullptr, nullptr, nullptr);
      }
    }
  }
  constexpr unsigned long const kernel::projectIniCond::NonZeroFlops;
  constexpr unsigned long const kernel::projectIniCond::HardwareFlops;
  void kernel::projectIniCond::execute() {
    assert(Q != nullptr);
    assert(iniCond != nullptr);
    assert(projectQP != nullptr);
    libxsmm_m56_n9_k343_ldA56_ldB344_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(projectQP, iniCond, Q, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::evalAtQP::NonZeroFlops;
  constexpr unsigned long const kernel::evalAtQP::HardwareFlops;
  void kernel::evalAtQP::execute() {
    assert(Q != nullptr);
    assert(dofsQP != nullptr);
    assert(evalAtQP != nullptr);
    libxsmm_m344_n9_k56_ldA344_ldB56_ldC344_alpha1_beta0_alignedA1_alignedC1_pfsigonly(evalAtQP, Q, dofsQP, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::volume::NonZeroFlops;
  constexpr unsigned long const kernel::volume::HardwareFlops;
  void kernel::volume::execute() {
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(kDivM(0) != nullptr);
    assert(kDivM(1) != nullptr);
    assert(kDivM(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(star(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp4;
    alignas(32) double _buffer0[324] ;
    _tmp0 = _buffer0;
    libxsmm_m36_n9_k9_ldA56_ldB9_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(I, star(0), _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k35_ldA56_ldB36_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(kDivM(0), _tmp0, Q, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m36_n9_k9_ldA56_ldB9_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(I, star(1), _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k35_ldA56_ldB36_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(kDivM(1), _tmp2, Q, nullptr, nullptr, nullptr);
    _tmp4 = _buffer0;
    libxsmm_m36_n9_k9_ldA56_ldB9_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(I, star(2), _tmp4, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k35_ldA56_ldB36_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(kDivM(2), _tmp4, Q, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::plConvertToNodal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertToNodal::HardwareFlops;
  void kernel::plConvertToNodal::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(initialLoading != nullptr);
    assert(replicateInitialLoading != nullptr);
    assert(v != nullptr);
    libxsmm_m56_n6_k56_ldA56_ldB56_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(v, QStress, QStressNodal, nullptr, nullptr, nullptr);
    for (int _p = 0; _p < 6; ++_p) {
      #pragma omp simd
      for (int _k = 0; _k < 56; ++_k) {
        QStressNodal[1*_k + 56*_p] += replicateInitialLoading[1*_k] * initialLoading[1*_p];
      }
    }
  }
  constexpr unsigned long const kernel::plConvertToNodalNoLoading::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertToNodalNoLoading::HardwareFlops;
  void kernel::plConvertToNodalNoLoading::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(v != nullptr);
    libxsmm_m56_n6_k56_ldA56_ldB56_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(v, QStress, QStressNodal, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::plConvertEtaModal2Nodal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertEtaModal2Nodal::HardwareFlops;
  void kernel::plConvertEtaModal2Nodal::execute() {
    assert(QEtaModal != nullptr);
    assert(QEtaNodal != nullptr);
    assert(v != nullptr);
    libxsmm_m56_n1_k56_ldA56_ldB56_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(v, QEtaModal, QEtaNodal, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::plConvertEtaNodal2Modal::NonZeroFlops;
  constexpr unsigned long const kernel::plConvertEtaNodal2Modal::HardwareFlops;
  void kernel::plConvertEtaNodal2Modal::execute() {
    assert(QEtaModal != nullptr);
    assert(QEtaNodal != nullptr);
    assert(vInv != nullptr);
    libxsmm_m56_n1_k56_ldA56_ldB56_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(vInv, QEtaNodal, QEtaModal, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::plComputeMean::NonZeroFlops;
  constexpr unsigned long const kernel::plComputeMean::HardwareFlops;
  void kernel::plComputeMean::execute() {
    assert(QStressNodal != nullptr);
    assert(meanStress != nullptr);
    assert(selectBulkAverage != nullptr);
    libxsmm_m56_n1_k3_ldA56_ldB3_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QStressNodal, selectBulkAverage, meanStress, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::plSubtractMean::NonZeroFlops;
  constexpr unsigned long const kernel::plSubtractMean::HardwareFlops;
  void kernel::plSubtractMean::execute() {
    assert(QStressNodal != nullptr);
    assert(meanStress != nullptr);
    assert(selectBulkNegative != nullptr);
    for (int _p = 0; _p < 3; ++_p) {
      #pragma omp simd
      for (int _k = 0; _k < 56; ++_k) {
        QStressNodal[1*_k + 56*_p] += meanStress[1*_k] * selectBulkNegative[1*_p];
      }
    }
  }
  constexpr unsigned long const kernel::plComputeSecondInvariant::NonZeroFlops;
  constexpr unsigned long const kernel::plComputeSecondInvariant::HardwareFlops;
  void kernel::plComputeSecondInvariant::execute() {
    assert(QStressNodal != nullptr);
    assert(secondInvariant != nullptr);
    assert(weightSecondInvariant != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[336] ;
    _tmp0 = _buffer0;
    for (int _q = 0; _q < 6; ++_q) {
      #pragma omp simd
      for (int _k = 0; _k < 56; ++_k) {
        _tmp0[1*_k + 56*_q] = QStressNodal[1*_k + 56*_q] * QStressNodal[1*_k + 56*_q];
      }
    }
    libxsmm_m56_n1_k6_ldA56_ldB6_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, weightSecondInvariant, secondInvariant, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::plAdjustStresses::NonZeroFlops;
  constexpr unsigned long const kernel::plAdjustStresses::HardwareFlops;
  void kernel::plAdjustStresses::execute() {
    assert(QStress != nullptr);
    assert(QStressNodal != nullptr);
    assert(vInv != nullptr);
    assert(yieldFactor != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[336] ;
    _tmp0 = _buffer0;
    for (int _p = 0; _p < 6; ++_p) {
      #pragma omp simd
      for (int _l = 0; _l < 56; ++_l) {
        _tmp0[1*_l + 56*_p] = QStressNodal[1*_l + 56*_p] * yieldFactor[1*_l];
      }
    }
    libxsmm_m56_n6_k56_ldA56_ldB56_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(vInv, _tmp0, QStress, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::createEasiBoundaryGhostCells::NonZeroFlops;
  constexpr unsigned long const kernel::createEasiBoundaryGhostCells::HardwareFlops;
  void kernel::createEasiBoundaryGhostCells::execute() {
    assert(INodal != nullptr);
    assert(easiBoundaryConstant != nullptr);
    assert(easiBoundaryMap != nullptr);
    assert(easiIdentMap != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[189] ;
    _tmp0 = _buffer0;
    for (int _l = 0; _l < 21; ++_l) {
      double const* _A = easiBoundaryMap + 81*_l;
      double const* _B = INodal + 1*_l;
      double * _C = _tmp0 + 9*_l;
      for (int n = 0; n < 1; ++n) {
        for (int m = 0; m < 9; ++m) {
          _C[0 + 1*m + 9*n] = 0.0;
        }
        for (int k = 0; k < 9; ++k) {
          for (int m = 0; m < 9; ++m) {
            _C[0 + 1*m + 9*n] += 1.0 * _A[0 + 1*m + 9*k] * _B[0 + 24*k + 216*n];
          }
        }
      }
    }
    for (int _l = 0; _l < 21; ++_l) {
      double const* _A = easiIdentMap + 81*_l;
      double const* _B = easiBoundaryConstant + 9*_l;
      double * _C = _tmp0 + 9*_l;
      libxsmm_m9_n1_k9_ldA9_ldB9_ldC9_alpha1_beta1_alignedA0_alignedC0_pfsigonly(_A, _B, _C, nullptr, nullptr, nullptr);
    }
    memset(INodal + 21, 0, 3 * sizeof(double));
    memset(INodal + 45, 0, 3 * sizeof(double));
    memset(INodal + 69, 0, 3 * sizeof(double));
    memset(INodal + 93, 0, 3 * sizeof(double));
    memset(INodal + 117, 0, 3 * sizeof(double));
    memset(INodal + 141, 0, 3 * sizeof(double));
    memset(INodal + 165, 0, 3 * sizeof(double));
    memset(INodal + 189, 0, 3 * sizeof(double));
    memset(INodal + 213, 0, 3 * sizeof(double));
    for (int _a = 0; _a < 9; ++_a) {
      #pragma omp simd
      for (int _l = 0; _l < 21; ++_l) {
        INodal[1*_l + 24*_a] = _tmp0[1*_a + 9*_l];
      }
    }
  }
  constexpr unsigned long const kernel::updateINodal::NonZeroFlops;
  constexpr unsigned long const kernel::updateINodal::HardwareFlops;
  void kernel::updateINodal::execute() {
    assert(!std::isnan(factor));
    assert(INodal != nullptr);
    assert(INodalUpdate != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 21; ++_a) {
        INodal[1*_a + 24*_b] += factor * INodalUpdate[1*_a + 24*_b];
      }
    }
  }
  constexpr unsigned long const kernel::rotateFaceDisplacement::NonZeroFlops;
  constexpr unsigned long const kernel::rotateFaceDisplacement::HardwareFlops;
  void kernel::rotateFaceDisplacement::execute() {
    assert(displacementRotationMatrix != nullptr);
    assert(faceDisplacement != nullptr);
    assert(rotatedFaceDisplacement != nullptr);
    for (int n = 0; n < 3; ++n) {
      for (int m = 0; m < 24; ++m) {
        rotatedFaceDisplacement[0 + 1*m + 24*n] = 0.0;
      }
      for (int k = 0; k < 3; ++k) {
        for (int m = 0; m < 24; ++m) {
          rotatedFaceDisplacement[0 + 1*m + 24*n] += 1.0 * faceDisplacement[0 + 1*m + 24*k] * displacementRotationMatrix[0 + 4*k + 1*n];
        }
      }
    }
  }
  constexpr unsigned long const kernel::computeMInvJInvPhisAtSources::NonZeroFlops;
  constexpr unsigned long const kernel::computeMInvJInvPhisAtSources::HardwareFlops;
  void kernel::computeMInvJInvPhisAtSources::execute() {
    assert(!std::isnan(JInv));
    assert(M3inv != nullptr);
    assert(basisFunctionsAtPoint != nullptr);
    assert(mInvJInvPhisAtSources != nullptr);
    for (int n = 0; n < 1; ++n) {
      for (int m = 0; m < 56; ++m) {
        mInvJInvPhisAtSources[0 + 1*m + 56*n] = 0.0;
      }
      for (int k = 0; k < 56; ++k) {
        for (int m = 0; m < 56; ++m) {
          mInvJInvPhisAtSources[0 + 1*m + 56*n] += JInv * M3inv[0 + 1*m + 56*k] * basisFunctionsAtPoint[0 + 1*k + 56*n];
        }
      }
    }
  }
  constexpr unsigned long const kernel::sourceNRF::NonZeroFlops;
  constexpr unsigned long const kernel::sourceNRF::HardwareFlops;
  void kernel::sourceNRF::execute() {
    assert(!std::isnan(mArea));
    assert(Q != nullptr);
    assert(mInvJInvPhisAtSources != nullptr);
    assert(mNormal != nullptr);
    assert(mSlip != nullptr);
    assert(momentToNRF != nullptr);
    assert(stiffnessTensor != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[9] ;
    alignas(32) double _buffer1[9] ;
    _tmp0 = _buffer0;
    for (int _j = 0; _j < 3; ++_j) {
      #pragma omp simd
      for (int _i = 0; _i < 3; ++_i) {
        _tmp0[1*_i + 3*_j] = mSlip[1*_i] * mNormal[1*_j];
      }
    }
    _tmp1 = _buffer1;
    libxsmm_m9_n1_k9_ldA9_ldB9_ldC9_alpha1_beta0_alignedA0_alignedC0_pfsigonly(stiffnessTensor, _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m6_n1_k9_ldA6_ldB9_ldC6_alpha1_beta0_alignedA0_alignedC0_pfsigonly(momentToNRF, _tmp1, _tmp2, nullptr, nullptr, nullptr);
    for (int _t = 0; _t < 6; ++_t) {
      #pragma omp simd
      for (int _k = 0; _k < 56; ++_k) {
        Q[1*_k + 56*_t] += mArea * mInvJInvPhisAtSources[1*_k] * _tmp2[1*_t];
      }
    }
  }
  constexpr unsigned long const kernel::sourceFSRM::NonZeroFlops;
  constexpr unsigned long const kernel::sourceFSRM::HardwareFlops;
  void kernel::sourceFSRM::execute() {
    assert(!std::isnan(stfIntegral));
    assert(Q != nullptr);
    assert(mInvJInvPhisAtSources != nullptr);
    assert(momentFSRM != nullptr);
    for (int _p = 0; _p < 9; ++_p) {
      #pragma omp simd
      for (int _k = 0; _k < 56; ++_k) {
        Q[1*_k + 56*_p] += stfIntegral * mInvJInvPhisAtSources[1*_k] * momentFSRM[1*_p];
      }
    }
  }
  constexpr unsigned long const kernel::evaluateDOFSAtPoint::NonZeroFlops;
  constexpr unsigned long const kernel::evaluateDOFSAtPoint::HardwareFlops;
  void kernel::evaluateDOFSAtPoint::execute() {
    assert(Q != nullptr);
    assert(QAtPoint != nullptr);
    assert(basisFunctionsAtPoint != nullptr);
    for (int n = 0; n < 1; ++n) {
      for (int m = 0; m < 9; ++m) {
        QAtPoint[0 + 1*m + 9*n] = 0.0;
      }
      for (int k = 0; k < 56; ++k) {
        for (int m = 0; m < 9; ++m) {
          QAtPoint[0 + 1*m + 9*n] += 1.0 * Q[0 + 56*m + 1*k] * basisFunctionsAtPoint[0 + 1*k + 56*n];
        }
      }
    }
  }
  constexpr unsigned long const kernel::evaluateDOFSAtPointSTP::NonZeroFlops;
  constexpr unsigned long const kernel::evaluateDOFSAtPointSTP::HardwareFlops;
  void kernel::evaluateDOFSAtPointSTP::execute() {
    assert(QAtPoint != nullptr);
    assert(basisFunctionsAtPoint != nullptr);
    assert(spaceTimePredictor != nullptr);
    assert(timeBasisFunctionsAtPoint != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[54] ;
    _tmp0 = _buffer0;
    for (int n = 0; n < 1; ++n) {
      for (int m = 0; m < 54; ++m) {
        _tmp0[0 + 1*m + 54*n] = 0.0;
      }
      for (int k = 0; k < 56; ++k) {
        for (int m = 0; m < 54; ++m) {
          _tmp0[0 + 1*m + 54*n] += 1.0 * spaceTimePredictor[0 + 56*m + 1*k] * basisFunctionsAtPoint[0 + 1*k + 56*n];
        }
      }
    }
    libxsmm_m9_n1_k6_ldA9_ldB6_ldC9_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp0, timeBasisFunctionsAtPoint, QAtPoint, nullptr, nullptr, nullptr);
  }
  namespace dynamicRupture {
    constexpr unsigned long const kernel::transposeTinv::NonZeroFlops;
    constexpr unsigned long const kernel::transposeTinv::HardwareFlops;
    void kernel::transposeTinv::execute() {
      assert(Tinv != nullptr);
      assert(TinvT != nullptr);
      for (int _j = 0; _j < 9; ++_j) {
        #pragma omp simd
        for (int _i = 0; _i < 9; ++_i) {
          TinvT[1*_i + 9*_j] = Tinv[1*_j + 9*_i];
        }
      }
    }
  } // namespace dynamicRupture
  namespace dynamicRupture {
    constexpr unsigned long const kernel::rotateFluxMatrix::NonZeroFlops;
    constexpr unsigned long const kernel::rotateFluxMatrix::HardwareFlops;
    void kernel::rotateFluxMatrix::execute() {
      assert(!std::isnan(fluxScale));
      assert(T != nullptr);
      assert(fluxSolver != nullptr);
      assert(star(0) != nullptr);
      for (int n = 0; n < 9; ++n) {
        for (int m = 0; m < 9; ++m) {
          fluxSolver[0 + 1*m + 9*n] = 0.0;
        }
        for (int k = 0; k < 9; ++k) {
          for (int m = 0; m < 9; ++m) {
            fluxSolver[0 + 1*m + 9*n] += fluxScale * star(0)[0 + 1*m + 9*k] * T[0 + 9*k + 1*n];
          }
        }
      }
    }
  } // namespace dynamicRupture
  constexpr unsigned long const kernel::localFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::localFlux::HardwareFlops[];
  constexpr kernel::localFlux::member_function_ptr kernel::localFlux::ExecutePtrs[];
  void kernel::localFlux::execute0() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(0) != nullptr);
    assert(rDivM(0) != nullptr);
    double *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fMrT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, AplusT, _tmp1, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp1, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::localFlux::execute1() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(1) != nullptr);
    assert(rDivM(1) != nullptr);
    double *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fMrT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, AplusT, _tmp1, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp1, Q, nullptr, _prefetch.Q, nullptr);
  }
  void kernel::localFlux::execute2() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(2) != nullptr);
    assert(rDivM(2) != nullptr);
    double *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fMrT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, AplusT, _tmp1, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp1, Q, nullptr, nullptr, nullptr);
  }
  void kernel::localFlux::execute3() {
    assert(AplusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fMrT(3) != nullptr);
    assert(rDivM(3) != nullptr);
    double *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fMrT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, AplusT, _tmp1, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp1, Q, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::localFluxNodal::NonZeroFlops[];
  constexpr unsigned long const kernel::localFluxNodal::HardwareFlops[];
  constexpr kernel::localFluxNodal::member_function_ptr kernel::localFluxNodal::ExecutePtrs[];
  void kernel::localFluxNodal::execute0() {
    assert(AminusT != nullptr);
    assert(INodal != nullptr);
    assert(Q != nullptr);
    assert(project2nFaceTo3m(0) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(INodal, AminusT, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(project2nFaceTo3m(0), _tmp0, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::localFluxNodal::execute1() {
    assert(AminusT != nullptr);
    assert(INodal != nullptr);
    assert(Q != nullptr);
    assert(project2nFaceTo3m(1) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(INodal, AminusT, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(project2nFaceTo3m(1), _tmp0, Q, nullptr, _prefetch.Q, nullptr);
  }
  void kernel::localFluxNodal::execute2() {
    assert(AminusT != nullptr);
    assert(INodal != nullptr);
    assert(Q != nullptr);
    assert(project2nFaceTo3m(2) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(INodal, AminusT, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(project2nFaceTo3m(2), _tmp0, Q, nullptr, nullptr, nullptr);
  }
  void kernel::localFluxNodal::execute3() {
    assert(AminusT != nullptr);
    assert(INodal != nullptr);
    assert(Q != nullptr);
    assert(project2nFaceTo3m(3) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(INodal, AminusT, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(project2nFaceTo3m(3), _tmp0, Q, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::neighboringFlux::NonZeroFlops[];
  constexpr unsigned long const kernel::neighboringFlux::HardwareFlops[];
  constexpr kernel::neighboringFlux::member_function_ptr kernel::neighboringFlux::ExecutePtrs[];
  void kernel::neighboringFlux::execute0() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute1() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB21_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute2() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute3() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute4() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB21_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute5() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute6() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute7() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB21_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute8() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute9() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute10() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB21_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute11() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(0) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmmsparse_38193e2bc92c2e551f139883f9580e10_m56_n9_k21_ldA0_ldB24_ldC56_alpha1_beta1_alignedA0_alignedC1_pfsigonly(rDivM(0), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute12() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute13() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute14() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute15() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute16() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute17() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute18() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute19() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute20() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute21() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute22() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute23() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(1) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(1), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute24() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute25() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute26() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute27() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute28() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute29() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute30() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute31() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute32() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute33() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute34() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute35() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(2) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(2), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute36() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute37() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute38() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(0) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(0), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute39() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute40() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute41() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(1) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(1), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute42() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute43() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute44() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(2), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute45() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(0) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(0), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute46() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(1) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[189] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmmsparse_b8016671674e78f7e2ec1c02d0b31df8_m21_n9_k21_ldA0_ldB24_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(fP(1), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m21_n9_k9_ldA21_ldB9_ldC21_alpha1_beta0_alignedA0_alignedC0_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB21_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  void kernel::neighboringFlux::execute47() {
    assert(AminusT != nullptr);
    assert(I != nullptr);
    assert(Q != nullptr);
    assert(fP(2) != nullptr);
    assert(rDivM(3) != nullptr);
    assert(rT(3) != nullptr);
    double *_tmp2, *_tmp0, *_tmp1;
    alignas(32) double _buffer0[216] ;
    alignas(32) double _buffer1[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(rT(3), I, _tmp0, nullptr, nullptr, nullptr);
    _tmp1 = _buffer1;
    libxsmm_m24_n9_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(fP(2), _tmp0, _tmp1, nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m24_n9_k9_ldA24_ldB9_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp1, AminusT, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m56_n9_k21_ldA56_ldB24_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(rDivM(3), _tmp2, Q, nullptr, _prefetch.I, nullptr);
  }
  constexpr unsigned long const kernel::derivativeTaylorExpansion::NonZeroFlops[];
  constexpr unsigned long const kernel::derivativeTaylorExpansion::HardwareFlops[];
  constexpr kernel::derivativeTaylorExpansion::member_function_ptr kernel::derivativeTaylorExpansion::ExecutePtrs[];
  void kernel::derivativeTaylorExpansion::execute0() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(0) != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 56; ++_a) {
        I[1*_a + 56*_b] = power * dQ(0)[1*_a + 56*_b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute1() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(1) != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 35; ++_a) {
        I[1*_a + 56*_b] += power * dQ(1)[1*_a + 36*_b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute2() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(2) != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 20; ++_a) {
        I[1*_a + 56*_b] += power * dQ(2)[1*_a + 20*_b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute3() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(3) != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 10; ++_a) {
        I[1*_a + 56*_b] += power * dQ(3)[1*_a + 12*_b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute4() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(4) != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 4; ++_a) {
        I[1*_a + 56*_b] += power * dQ(4)[1*_a + 4*_b];
      }
    }
  }
  void kernel::derivativeTaylorExpansion::execute5() {
    assert(!std::isnan(power));
    assert(I != nullptr);
    assert(dQ(5) != nullptr);
    for (int _b = 0; _b < 9; ++_b) {
      #pragma omp simd
      for (int _a = 0; _a < 1; ++_a) {
        I[1*_a + 56*_b] += power * dQ(5)[1*_a + 4*_b];
      }
    }
  }
  constexpr unsigned long const kernel::derivative::NonZeroFlops[];
  constexpr unsigned long const kernel::derivative::HardwareFlops[];
  constexpr kernel::derivative::member_function_ptr kernel::derivative::ExecutePtrs[];
  void kernel::derivative::execute1() {
    assert(dQ(0) != nullptr);
    assert(dQ(1) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(star(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp4;
    alignas(32) double _buffer0[324] ;
    _tmp0 = _buffer0;
    libxsmm_m36_n9_k55_ldA36_ldB56_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(0), dQ(0) + 1, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m36_n9_k9_ldA36_ldB9_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, star(0), dQ(1), nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m36_n9_k55_ldA36_ldB56_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(1), dQ(0) + 1, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m36_n9_k9_ldA36_ldB9_ldC36_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp2, star(1), dQ(1), nullptr, nullptr, nullptr);
    _tmp4 = _buffer0;
    libxsmm_m36_n9_k55_ldA36_ldB56_ldC36_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(2), dQ(0) + 1, _tmp4, nullptr, nullptr, nullptr);
    libxsmm_m36_n9_k9_ldA36_ldB9_ldC36_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp4, star(2), dQ(1), nullptr, nullptr, nullptr);
  }
  void kernel::derivative::execute2() {
    assert(dQ(1) != nullptr);
    assert(dQ(2) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(star(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp4;
    alignas(32) double _buffer0[180] ;
    _tmp0 = _buffer0;
    libxsmm_m20_n9_k34_ldA36_ldB36_ldC20_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(0), dQ(1) + 1, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m20_n9_k9_ldA20_ldB9_ldC20_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, star(0), dQ(2), nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m20_n9_k34_ldA36_ldB36_ldC20_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(1), dQ(1) + 1, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m20_n9_k9_ldA20_ldB9_ldC20_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp2, star(1), dQ(2), nullptr, nullptr, nullptr);
    _tmp4 = _buffer0;
    libxsmm_m20_n9_k34_ldA36_ldB36_ldC20_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(2), dQ(1) + 1, _tmp4, nullptr, nullptr, nullptr);
    libxsmm_m20_n9_k9_ldA20_ldB9_ldC20_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp4, star(2), dQ(2), nullptr, nullptr, nullptr);
  }
  void kernel::derivative::execute3() {
    assert(dQ(2) != nullptr);
    assert(dQ(3) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(star(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp4;
    alignas(32) double _buffer0[108] ;
    _tmp0 = _buffer0;
    libxsmm_m12_n9_k19_ldA36_ldB20_ldC12_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(0), dQ(2) + 1, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m12_n9_k9_ldA12_ldB9_ldC12_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, star(0), dQ(3), nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m12_n9_k19_ldA36_ldB20_ldC12_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(1), dQ(2) + 1, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m12_n9_k9_ldA12_ldB9_ldC12_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp2, star(1), dQ(3), nullptr, nullptr, nullptr);
    _tmp4 = _buffer0;
    libxsmm_m12_n9_k19_ldA36_ldB20_ldC12_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(2), dQ(2) + 1, _tmp4, nullptr, nullptr, nullptr);
    libxsmm_m12_n9_k9_ldA12_ldB9_ldC12_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp4, star(2), dQ(3), nullptr, nullptr, nullptr);
  }
  void kernel::derivative::execute4() {
    assert(dQ(3) != nullptr);
    assert(dQ(4) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(star(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp4;
    alignas(32) double _buffer0[36] ;
    _tmp0 = _buffer0;
    libxsmm_m4_n9_k9_ldA36_ldB12_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(0), dQ(3) + 1, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, star(0), dQ(4), nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m4_n9_k9_ldA36_ldB12_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(1), dQ(3) + 1, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp2, star(1), dQ(4), nullptr, nullptr, nullptr);
    _tmp4 = _buffer0;
    libxsmm_m4_n9_k9_ldA36_ldB12_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(2), dQ(3) + 1, _tmp4, nullptr, nullptr, nullptr);
    libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp4, star(2), dQ(4), nullptr, nullptr, nullptr);
  }
  void kernel::derivative::execute5() {
    assert(dQ(4) != nullptr);
    assert(dQ(5) != nullptr);
    assert(kDivMT(0) != nullptr);
    assert(kDivMT(1) != nullptr);
    assert(kDivMT(2) != nullptr);
    assert(star(0) != nullptr);
    assert(star(1) != nullptr);
    assert(star(2) != nullptr);
    double *_tmp2, *_tmp0, *_tmp4;
    alignas(32) double _buffer0[36] ;
    _tmp0 = _buffer0;
    libxsmm_m4_n9_k3_ldA36_ldB4_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(0), dQ(4) + 1, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, star(0), dQ(5), nullptr, nullptr, nullptr);
    _tmp2 = _buffer0;
    libxsmm_m4_n9_k3_ldA36_ldB4_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(1), dQ(4) + 1, _tmp2, nullptr, nullptr, nullptr);
    libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp2, star(1), dQ(5), nullptr, nullptr, nullptr);
    _tmp4 = _buffer0;
    libxsmm_m4_n9_k3_ldA36_ldB4_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(kDivMT(2), dQ(4) + 1, _tmp4, nullptr, nullptr, nullptr);
    libxsmm_m4_n9_k9_ldA4_ldB9_ldC4_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp4, star(2), dQ(5), nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::projectToNodalBoundary::NonZeroFlops[];
  constexpr unsigned long const kernel::projectToNodalBoundary::HardwareFlops[];
  constexpr kernel::projectToNodalBoundary::member_function_ptr kernel::projectToNodalBoundary::ExecutePtrs[];
  void kernel::projectToNodalBoundary::execute0() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(V3mTo2nFace(0) != nullptr);
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(0), I, INodal, nullptr, nullptr, nullptr);
  }
  void kernel::projectToNodalBoundary::execute1() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(V3mTo2nFace(1) != nullptr);
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(1), I, INodal, nullptr, nullptr, nullptr);
  }
  void kernel::projectToNodalBoundary::execute2() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(V3mTo2nFace(2) != nullptr);
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(2), I, INodal, nullptr, nullptr, nullptr);
  }
  void kernel::projectToNodalBoundary::execute3() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(V3mTo2nFace(3) != nullptr);
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(3), I, INodal, nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::projectToNodalBoundaryRotated::NonZeroFlops[];
  constexpr unsigned long const kernel::projectToNodalBoundaryRotated::HardwareFlops[];
  constexpr kernel::projectToNodalBoundaryRotated::member_function_ptr kernel::projectToNodalBoundaryRotated::ExecutePtrs[];
  void kernel::projectToNodalBoundaryRotated::execute0() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(Tinv != nullptr);
    assert(V3mTo2nFace(0) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(0), I, _tmp0, nullptr, nullptr, nullptr);
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 24; ++m) {
        INodal[0 + 1*m + 24*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 24; ++m) {
          INodal[0 + 1*m + 24*n] += 1.0 * _tmp0[0 + 1*m + 24*k] * Tinv[0 + 9*k + 1*n];
        }
      }
    }
  }
  void kernel::projectToNodalBoundaryRotated::execute1() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(Tinv != nullptr);
    assert(V3mTo2nFace(1) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(1), I, _tmp0, nullptr, nullptr, nullptr);
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 24; ++m) {
        INodal[0 + 1*m + 24*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 24; ++m) {
          INodal[0 + 1*m + 24*n] += 1.0 * _tmp0[0 + 1*m + 24*k] * Tinv[0 + 9*k + 1*n];
        }
      }
    }
  }
  void kernel::projectToNodalBoundaryRotated::execute2() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(Tinv != nullptr);
    assert(V3mTo2nFace(2) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(2), I, _tmp0, nullptr, nullptr, nullptr);
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 24; ++m) {
        INodal[0 + 1*m + 24*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 24; ++m) {
          INodal[0 + 1*m + 24*n] += 1.0 * _tmp0[0 + 1*m + 24*k] * Tinv[0 + 9*k + 1*n];
        }
      }
    }
  }
  void kernel::projectToNodalBoundaryRotated::execute3() {
    assert(I != nullptr);
    assert(INodal != nullptr);
    assert(Tinv != nullptr);
    assert(V3mTo2nFace(3) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[216] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n9_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(3), I, _tmp0, nullptr, nullptr, nullptr);
    for (int n = 0; n < 9; ++n) {
      for (int m = 0; m < 24; ++m) {
        INodal[0 + 1*m + 24*n] = 0.0;
      }
      for (int k = 0; k < 9; ++k) {
        for (int m = 0; m < 24; ++m) {
          INodal[0 + 1*m + 24*n] += 1.0 * _tmp0[0 + 1*m + 24*k] * Tinv[0 + 9*k + 1*n];
        }
      }
    }
  }
  constexpr unsigned long const kernel::subTriangleDisplacement::NonZeroFlops[];
  constexpr unsigned long const kernel::subTriangleDisplacement::HardwareFlops[];
  constexpr kernel::subTriangleDisplacement::member_function_ptr kernel::subTriangleDisplacement::ExecutePtrs[];
  void kernel::subTriangleDisplacement::execute0() {
    assert(MV2nTo2m != nullptr);
    assert(faceDisplacement != nullptr);
    assert(subTriangleDofs(0) != nullptr);
    assert(subTriangleProjectionFromFace(0) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[84] ;
    _tmp0 = _buffer0;
    libxsmm_m4_n21_k21_ldA4_ldB24_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjectionFromFace(0), MV2nTo2m, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m4_n3_k21_ldA4_ldB24_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, faceDisplacement, subTriangleDofs(0), nullptr, nullptr, nullptr);
  }
  void kernel::subTriangleDisplacement::execute1() {
    assert(MV2nTo2m != nullptr);
    assert(faceDisplacement != nullptr);
    assert(subTriangleDofs(1) != nullptr);
    assert(subTriangleProjectionFromFace(1) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[72] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n3_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(MV2nTo2m, faceDisplacement, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m4_n3_k21_ldA4_ldB24_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjectionFromFace(1), _tmp0, subTriangleDofs(1), nullptr, nullptr, nullptr);
  }
  void kernel::subTriangleDisplacement::execute2() {
    assert(MV2nTo2m != nullptr);
    assert(faceDisplacement != nullptr);
    assert(subTriangleDofs(2) != nullptr);
    assert(subTriangleProjectionFromFace(2) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[72] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n3_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(MV2nTo2m, faceDisplacement, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m16_n3_k21_ldA16_ldB24_ldC16_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjectionFromFace(2), _tmp0, subTriangleDofs(2), nullptr, nullptr, nullptr);
  }
  void kernel::subTriangleDisplacement::execute3() {
    assert(MV2nTo2m != nullptr);
    assert(faceDisplacement != nullptr);
    assert(subTriangleDofs(3) != nullptr);
    assert(subTriangleProjectionFromFace(3) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[72] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n3_k21_ldA24_ldB24_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(MV2nTo2m, faceDisplacement, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m64_n3_k21_ldA64_ldB24_ldC64_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjectionFromFace(3), _tmp0, subTriangleDofs(3), nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::subTriangleVelocity::NonZeroFlops[];
  constexpr unsigned long const kernel::subTriangleVelocity::HardwareFlops[];
  constexpr kernel::subTriangleVelocity::member_function_ptr kernel::subTriangleVelocity::ExecutePtrs[];
  void kernel::subTriangleVelocity::execute0() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(0) != nullptr);
    assert(subTriangleProjection(0) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[12] ;
    _tmp0 = _buffer0;
    libxsmm_m4_n3_k56_ldA4_ldB56_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjection(0), Q + 336, _tmp0, nullptr, nullptr, nullptr);
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m4_n3_k3_ldA4_ldB0_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, selectVelocity, subTriangleDofs(0), nullptr, nullptr, nullptr);
  }
  void kernel::subTriangleVelocity::execute1() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(1) != nullptr);
    assert(subTriangleProjection(1) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[12] ;
    _tmp0 = _buffer0;
    libxsmm_m4_n3_k56_ldA4_ldB56_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjection(1), Q + 336, _tmp0, nullptr, nullptr, nullptr);
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m4_n3_k3_ldA4_ldB0_ldC4_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, selectVelocity, subTriangleDofs(1), nullptr, nullptr, nullptr);
  }
  void kernel::subTriangleVelocity::execute2() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(2) != nullptr);
    assert(subTriangleProjection(2) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[48] ;
    _tmp0 = _buffer0;
    libxsmm_m16_n3_k56_ldA16_ldB56_ldC16_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjection(2), Q + 336, _tmp0, nullptr, nullptr, nullptr);
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m16_n3_k3_ldA16_ldB0_ldC16_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, selectVelocity, subTriangleDofs(2), nullptr, nullptr, nullptr);
  }
  void kernel::subTriangleVelocity::execute3() {
    assert(Q != nullptr);
    assert(selectVelocity != nullptr);
    assert(subTriangleDofs(3) != nullptr);
    assert(subTriangleProjection(3) != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[168] ;
    _tmp0 = _buffer0;
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m56_n3_k3_ldA56_ldB0_ldC56_alpha1_beta0_alignedA1_alignedC1_pfsigonly(Q + 336, selectVelocity, _tmp0, nullptr, nullptr, nullptr);
    libxsmm_m64_n3_k56_ldA64_ldB56_ldC64_alpha1_beta0_alignedA1_alignedC1_pfsigonly(subTriangleProjection(3), _tmp0, subTriangleDofs(3), nullptr, nullptr, nullptr);
  }
  constexpr unsigned long const kernel::addVelocity::NonZeroFlops[];
  constexpr unsigned long const kernel::addVelocity::HardwareFlops[];
  constexpr kernel::addVelocity::member_function_ptr kernel::addVelocity::ExecutePtrs[];
  void kernel::addVelocity::execute0() {
    assert(I != nullptr);
    assert(V3mTo2nFace(0) != nullptr);
    assert(faceDisplacement != nullptr);
    assert(selectVelocity != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[72] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n3_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(0), I + 336, _tmp0, nullptr, nullptr, nullptr);
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m24_n3_k3_ldA24_ldB0_ldC24_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp0, selectVelocity, faceDisplacement, nullptr, nullptr, nullptr);
  }
  void kernel::addVelocity::execute1() {
    assert(I != nullptr);
    assert(V3mTo2nFace(1) != nullptr);
    assert(faceDisplacement != nullptr);
    assert(selectVelocity != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[72] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n3_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(1), I + 336, _tmp0, nullptr, nullptr, nullptr);
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m24_n3_k3_ldA24_ldB0_ldC24_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp0, selectVelocity, faceDisplacement, nullptr, nullptr, nullptr);
  }
  void kernel::addVelocity::execute2() {
    assert(I != nullptr);
    assert(V3mTo2nFace(2) != nullptr);
    assert(faceDisplacement != nullptr);
    assert(selectVelocity != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[72] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n3_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(2), I + 336, _tmp0, nullptr, nullptr, nullptr);
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m24_n3_k3_ldA24_ldB0_ldC24_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp0, selectVelocity, faceDisplacement, nullptr, nullptr, nullptr);
  }
  void kernel::addVelocity::execute3() {
    assert(I != nullptr);
    assert(V3mTo2nFace(3) != nullptr);
    assert(faceDisplacement != nullptr);
    assert(selectVelocity != nullptr);
    double *_tmp0;
    alignas(32) double _buffer0[72] ;
    _tmp0 = _buffer0;
    libxsmm_m24_n3_k56_ldA24_ldB56_ldC24_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2nFace(3), I + 336, _tmp0, nullptr, nullptr, nullptr);
    libxsmmsparse_e60658f2bd24617a41c922e81bd775a8_m24_n3_k3_ldA24_ldB0_ldC24_alpha1_beta1_alignedA1_alignedC1_pfsigonly(_tmp0, selectVelocity, faceDisplacement, nullptr, nullptr, nullptr);
  }
  namespace dynamicRupture {
    constexpr unsigned long const kernel::evaluateAndRotateQAtInterpolationPoints::NonZeroFlops[];
    constexpr unsigned long const kernel::evaluateAndRotateQAtInterpolationPoints::HardwareFlops[];
    constexpr kernel::evaluateAndRotateQAtInterpolationPoints::member_function_ptr kernel::evaluateAndRotateQAtInterpolationPoints::ExecutePtrs[];
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute0() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(0,0) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(0,0), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute1() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(1,0) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(1,0), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute2() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(2,0) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(2,0), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute3() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(3,0) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(3,0), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute4() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(0,1) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(0,1), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute5() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(1,1) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(1,1), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute6() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(2,1) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(2,1), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute7() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(3,1) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(3,1), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute8() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(0,2) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(0,2), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute9() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(1,2) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(1,2), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute10() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(2,2) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(2,2), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute11() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(3,2) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(3,2), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute12() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(0,3) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(0,3), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute13() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(1,3) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(1,3), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute14() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(2,3) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(2,3), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
    void kernel::evaluateAndRotateQAtInterpolationPoints::execute15() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(TinvT != nullptr);
      assert(V3mTo2n(3,3) != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k56_ldA52_ldB56_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(V3mTo2n(3,3), Q, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(_tmp0, TinvT, QInterpolated, nullptr, _prefetch.QInterpolated, nullptr);
    }
  } // namespace dynamicRupture
  namespace dynamicRupture {
    constexpr unsigned long const kernel::nodalFlux::NonZeroFlops[];
    constexpr unsigned long const kernel::nodalFlux::HardwareFlops[];
    constexpr kernel::nodalFlux::member_function_ptr kernel::nodalFlux::ExecutePtrs[];
    void kernel::nodalFlux::execute0() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(0,0) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(0,0), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute1() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(1,0) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(1,0), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute2() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(2,0) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(2,0), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute3() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(3,0) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(3,0), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute4() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(0,1) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(0,1), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute5() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(1,1) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(1,1), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute6() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(2,1) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(2,1), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute7() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(3,1) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(3,1), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute8() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(0,2) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(0,2), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute9() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(1,2) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(1,2), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute10() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(2,2) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(2,2), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute11() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(3,2) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(3,2), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute12() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(0,3) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(0,3), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute13() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(1,3) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(1,3), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute14() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(2,3) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(2,3), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
    void kernel::nodalFlux::execute15() {
      assert(Q != nullptr);
      assert(QInterpolated != nullptr);
      assert(V3mTo2nTWDivM(3,3) != nullptr);
      assert(fluxSolver != nullptr);
      double *_tmp0;
      alignas(32) double _buffer0[468] ;
      _tmp0 = _buffer0;
      libxsmm_m52_n9_k9_ldA52_ldB9_ldC52_alpha1_beta0_alignedA1_alignedC1_pfsigonly(QInterpolated, fluxSolver, _tmp0, nullptr, nullptr, nullptr);
      libxsmm_m56_n9_k49_ldA56_ldB52_ldC56_alpha1_beta1_alignedA1_alignedC1_pfsigonly(V3mTo2nTWDivM(3,3), _tmp0, Q, nullptr, _prefetch.I, nullptr);
    }
  } // namespace dynamicRupture
} // namespace seissol
