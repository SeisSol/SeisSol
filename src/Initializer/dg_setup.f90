!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Verena Hermann (hermann AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/hermann)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2009-2017, SeisSol Group
!! All rights reserved.
!!
!! Redistribution and use in source and binary forms, with or without
!! modification, are permitted provided that the following conditions are met:
!!
!! 1. Redistributions of source code must retain the above copyright notice,
!!    this list of conditions and the following disclaimer.
!!
!! 2. Redistributions in binary form must reproduce the above copyright notice,
!!    this list of conditions and the following disclaimer in the documentation
!!    and/or other materials provided with the distribution.
!!
!! 3. Neither the name of the copyright holder nor the names of its
!!    contributors may be used to endorse or promote products derived from this
!!    software without specific prior written permission.
!!
!! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!! POSSIBILITY OF SUCH DAMAGE.
!! @section DESCRIPTION
!=============================================================================!
!!                                                                           !!
!!  ADER Discontinuous Galerkin Module                                       !!
!!                                                                           !!
!!---------------------------------------------------------------------------!!
!!                                                                           !!
!!  Equations                :  Linear Hyperbolic Equations                  !!
!!  Number of dimensions     :  3D                                           !!
!!  Element types supported  :  tetrahedrons                                 !!
!!  Spatial accuracy         :  arbitrary                                    !!
!!  Temporal accuracy        :  arbitrary                                    !!
!!                                                                           !!
!!---------------------------------------------------------------------------!!

! define preprocessor variables
#include <Initializer/preProcessorMacros.fpp>

MODULE dg_setup_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  INTERFACE closeGalerkin3D_us
     MODULE PROCEDURE closeGalerkin3D_us
  END INTERFACE
  INTERFACE iniGalerkin3D_us_level1_new
     MODULE PROCEDURE iniGalerkin3D_us_level1_new
  END INTERFACE
  INTERFACE iniGalerkin3D_us_level2_new
     MODULE PROCEDURE iniGalerkin3D_us_level2_new
  END INTERFACE
  INTERFACE iniGalerkin3D_us_intern_new
     MODULE PROCEDURE iniGalerkin3D_us_intern_new
  END INTERFACE
  INTERFACE icGalerkin3D_us_new
     MODULE PROCEDURE icGalerkin3D_us_new
  END INTERFACE
  INTERFACE IniSparseStarMatrices3D_new
     MODULE PROCEDURE IniSparseStarMatrices3D_new
  END INTERFACE
  INTERFACE Read2dGF
     MODULE PROCEDURE Read2dGF
  END INTERFACE

!  INTERFACE NonConformingGPEvaluation3D
!     MODULE PROCEDURE NonConformingGPEvaluation3D
!  END INTERFACE
  ! Private procedures and functions
  INTERFACE BuildSpecialDGGeometry3D_new
     MODULE PROCEDURE BuildSpecialDGGeometry3D_new
  END INTERFACE
  !
  !---------------------------------------------------------------------------!
  PUBLIC  ::                               &
       closeGalerkin3D_us,                 &
       AnalyseGalerkin3D_us_new,           &
       iniGalerkin3D_us_level1_new,        &
       iniGalerkin3D_us_level2_new,        &
       iniGalerkin3D_us_intern_new,        &
       icGalerkin3D_us_new,                &
       IniSparseStarMatrices3D_new,        &
       Read2dGF
!      NonConformingGPEvaluation3D
  !---------------------------------------------------------------------------!

CONTAINS

  !===========================================================================!
  !!                                                                         !!
  !! closeGalerkin deallocates the arrays used by RKDG method                !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE closeGalerkin3D_us(EQN,MESH,DISC,IO)
    !-------------------------------------------------------------------------!

    USE COMMON_operators_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tDiscretization)    :: DISC
    TYPE(tInputOutput)       :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                  :: i,j,k,iElem,iSide
    INTEGER                  :: LocElemType
    INTEGER                  :: iLocalNeighborSide
    INTEGER                  :: iLocalNeighborVrtx
    !-------------------------------------------------------------------------!
    INTENT(IN)               :: IO
    INTENT(INOUT)            :: DISC
    !-------------------------------------------------------------------------!
    !
    IF(.NOT.DISC%Galerkin%init) THEN
       logError(*) 'closeGalerkin: SeisSol Interface not initialized!!'
       STOP
    ENDIF
    !
    logInfo(*) 'Enter closeGalerkin...'
    !
    DEALLOCATE( DISC%Galerkin%geoNormals         )
    DEALLOCATE( DISC%Galerkin%geoTangent1        )
    DEALLOCATE( DISC%Galerkin%geoTangent2        )
    DEALLOCATE( DISC%Galerkin%geoSurfaces        )
    !
    IF (ASSOCIATED( DISC%Galerkin%Faculty)) DEALLOCATE(DISC%Galerkin%Faculty)
    IF (ASSOCIATED( DISC%Galerkin%TimeGaussP)) DEALLOCATE(DISC%Galerkin%TimeGaussP)
    IF (ASSOCIATED( DISC%Galerkin%TimeGaussW)) DEALLOCATE(DISC%Galerkin%TimeGaussW)
#ifdef GENERATEDKERNELS
    ! Interoperability with C needs continuous arrays in memory.
    deallocate(disc%galerkin%dgvar)
#else
    IF (ASSOCIATED( DISC%Galerkin%dgvar))   DEALLOCATE(DISC%Galerkin%dgvar)
    IF (ASSOCIATED( DISC%Galerkin%DGwork))  DEALLOCATE(DISC%Galerkin%DGwork)
#endif

! Generated Kernels: We didn't initialize the old sparse format.
#if !defined(GENERATEDKERNELS)
    DO iElem = 1, MESH%nElem
        LocElemType = MESH%LocalElemType(iElem)
        CALL CloseSparseTensor3( DISC%Galerkin%AStar_Sp(iElem) )
        CALL CloseSparseTensor3( DISC%Galerkin%BStar_Sp(iElem) )
        CALL CloseSparseTensor3( DISC%Galerkin%CStar_Sp(iElem) )
!        IF(EQN%LocAnelastic(iElem).EQ.1) THEN
!            CALL CloseSparseTensor3( DISC%Galerkin%EStar_Sp(iElem) )
!        ENDIF
        DO iSide = 1,LocElemType
            CALL CloseSparseTensor3( DISC%Galerkin%FLStar_Sp(iElem,iSide) )
            CALL CloseSparseTensor3( DISC%Galerkin%FRStar_Sp(iElem,iSide) )
        ENDDO
    ENDDO
#endif

! Generated Kernels: Didn't initialize the old sparse format to save memory.
#if !defined(GENERATEDKERNELS)
    DEALLOCATE( DISC%Galerkin%AStar_Sp)
    DEALLOCATE( DISC%Galerkin%BStar_Sp)
    DEALLOCATE( DISC%Galerkin%CStar_Sp)
    DEALLOCATE( DISC%Galerkin%FLStar_Sp)
    DEALLOCATE( DISC%Galerkin%FRStar_Sp)
#endif


    IF (ASSOCIATED( DISC%Galerkin%cPoly3D_Hex)) DEALLOCATE(DISC%Galerkin%cPoly3D_Hex)
    IF (ASSOCIATED( DISC%Galerkin%NonZeroCPoly_Hex)) DEALLOCATE(DISC%Galerkin%NonZeroCPoly_Hex)
    IF (ASSOCIATED( DISC%Galerkin%NonZeroCPolyIndex_Hex)) DEALLOCATE(DISC%Galerkin%NonZeroCPolyIndex_Hex)
    IF (ASSOCIATED( DISC%Galerkin%intGaussP_Hex)) DEALLOCATE(DISC%Galerkin%intGaussP_Hex)
    IF (ASSOCIATED( DISC%Galerkin%intGaussW_Hex)) DEALLOCATE(DISC%Galerkin%intGaussW_Hex)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussP_Hex)) DEALLOCATE(DISC%Galerkin%bndGaussP_Hex)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussW_Hex)) DEALLOCATE(DISC%Galerkin%bndGaussW_Hex)
    IF (ASSOCIATED( DISC%Galerkin%IntGPBaseFunc_Hex)) DEALLOCATE(DISC%Galerkin%IntGPBaseFunc_Hex)
    IF (ASSOCIATED( DISC%Galerkin%IntGPBaseGrad_Hex)) DEALLOCATE(DISC%Galerkin%IntGPBaseGrad_Hex)
    IF (ASSOCIATED( DISC%Galerkin%BndGPBaseFunc_Hex)) DEALLOCATE(DISC%Galerkin%BndGPBaseFunc_Hex)
    IF (ASSOCIATED( DISC%Galerkin%BndGPBaseFunc3D_Hex)) DEALLOCATE(DISC%Galerkin%BndGPBaseFunc3D_Hex)
    IF (ASSOCIATED( DISC%Galerkin%MassMatrix_Hex)) DEALLOCATE(DISC%Galerkin%MassMatrix_Hex)
    IF (ASSOCIATED( DISC%Galerkin%iMassMatrix_Hex)) DEALLOCATE(DISC%Galerkin%iMassMatrix_Hex)
    IF (ASSOCIATED( DISC%Galerkin%Kxi_k_Hex)) DEALLOCATE(DISC%Galerkin%Kxi_k_Hex)
    IF (ASSOCIATED( DISC%Galerkin%Keta_k_Hex)) DEALLOCATE(DISC%Galerkin%Keta_k_Hex)
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_k_Hex)) DEALLOCATE(DISC%Galerkin%Kzeta_k_Hex)
    IF (ASSOCIATED( DISC%Galerkin%Kxi_m_Hex)) DEALLOCATE(DISC%Galerkin%Kxi_m_Hex)
    IF (ASSOCIATED( DISC%Galerkin%Keta_m_Hex)) DEALLOCATE(DISC%Galerkin%Keta_m_Hex)
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_m_Hex)) DEALLOCATE(DISC%Galerkin%Kzeta_m_Hex)
    IF (ASSOCIATED( DISC%Galerkin%ADGxi_Hex)) DEALLOCATE(DISC%Galerkin%ADGxi_Hex)
    IF (ASSOCIATED( DISC%Galerkin%ADGeta_Hex)) DEALLOCATE(DISC%Galerkin%ADGeta_Hex)
    IF (ASSOCIATED( DISC%Galerkin%ADGzeta_Hex)) DEALLOCATE(DISC%Galerkin%ADGzeta_Hex)
    IF (ASSOCIATED( DISC%Galerkin%FluxInt_Hex)) DEALLOCATE(DISC%Galerkin%FluxInt_Hex)
    IF (ASSOCIATED( DISC%Galerkin%Kxi_k_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kxi_k_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%Kxi_k_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Keta_k_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Keta_k_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%Keta_k_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_k_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kzeta_k_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%Kzeta_k_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Kxi_m_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kxi_m_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%Kxi_m_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Keta_m_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Keta_m_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%Keta_m_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_m_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kzeta_m_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%Kzeta_m_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%ADGxi_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%ADGxi_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%ADGxi_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%ADGeta_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%ADGeta_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%ADGeta_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%ADGzeta_Hex_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%ADGzeta_Hex_Sp)
        DEALLOCATE(DISC%Galerkin%ADGzeta_Hex_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%FluxInt_Hex_Sp)) THEN
        DO iSide = 1, MESH%nSides_Hex
            CALL CloseSparseTensor3b(DISC%Galerkin%FluxInt_Hex_Sp(0,1,iSide))
            DO iLocalNeighborSide = 1, MESH%nSides_Hex
                DO iLocalNeighborVrtx = 1, MESH%nVertices_Quad
                    CALL CloseSparseTensor3b(DISC%Galerkin%FluxInt_Hex_Sp(iLocalNeighborSide,iLocalNeighborVrtx,iSide))
                ENDDO
            ENDDO
        ENDDO
        DEALLOCATE(DISC%Galerkin%FluxInt_Hex_Sp)
    ENDIF


    IF (ASSOCIATED( DISC%Galerkin%cPoly3D_Tet)) DEALLOCATE(DISC%Galerkin%cPoly3D_Tet)
    IF (ASSOCIATED( DISC%Galerkin%NonZeroCPoly_Tet)) DEALLOCATE(DISC%Galerkin%NonZeroCPoly_Tet)
    IF (ASSOCIATED( DISC%Galerkin%NonZeroCPolyIndex_Tet)) DEALLOCATE(DISC%Galerkin%NonZeroCPolyIndex_Tet)
    IF (ASSOCIATED( DISC%Galerkin%intGaussP_Tet)) DEALLOCATE(DISC%Galerkin%intGaussP_Tet)
    IF (ASSOCIATED( DISC%Galerkin%intGaussW_Tet)) DEALLOCATE(DISC%Galerkin%intGaussW_Tet)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussP_Tet)) DEALLOCATE(DISC%Galerkin%bndGaussP_Tet)
    IF (ASSOCIATED( DISC%Galerkin%bndGaussW_Tet)) DEALLOCATE(DISC%Galerkin%bndGaussW_Tet)
    IF (ASSOCIATED( DISC%Galerkin%IntGPBaseFunc_Tet)) DEALLOCATE(DISC%Galerkin%IntGPBaseFunc_Tet)
    IF (ASSOCIATED( DISC%Galerkin%IntGPBaseGrad_Tet)) DEALLOCATE(DISC%Galerkin%IntGPBaseGrad_Tet)
    IF (ASSOCIATED( DISC%Galerkin%BndGPBaseFunc_Tet)) DEALLOCATE(DISC%Galerkin%BndGPBaseFunc_Tet)
!    IF (ASSOCIATED( DISC%Galerkin%BndGPBaseFunc3D_Tet)) DEALLOCATE(DISC%Galerkin%BndGPBaseFunc3D_Tet)
    IF (ASSOCIATED( DISC%Galerkin%MassMatrix_Tet)) DEALLOCATE(DISC%Galerkin%MassMatrix_Tet)
    IF (ASSOCIATED( DISC%Galerkin%iMassMatrix_Tet)) DEALLOCATE(DISC%Galerkin%iMassMatrix_Tet)
    IF (ASSOCIATED( DISC%Galerkin%Kxi_k_Tet)) DEALLOCATE(DISC%Galerkin%Kxi_k_Tet)
    IF (ASSOCIATED( DISC%Galerkin%Keta_k_Tet)) DEALLOCATE(DISC%Galerkin%Keta_k_Tet)
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_k_Tet)) DEALLOCATE(DISC%Galerkin%Kzeta_k_Tet)
    IF (ASSOCIATED( DISC%Galerkin%Kxi_m_Tet)) DEALLOCATE(DISC%Galerkin%Kxi_m_Tet)
    IF (ASSOCIATED( DISC%Galerkin%Keta_m_Tet)) DEALLOCATE(DISC%Galerkin%Keta_m_Tet)
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_m_Tet)) DEALLOCATE(DISC%Galerkin%Kzeta_m_Tet)
    IF (ASSOCIATED( DISC%Galerkin%ADGxi_Tet)) DEALLOCATE(DISC%Galerkin%ADGxi_Tet)
    IF (ASSOCIATED( DISC%Galerkin%ADGeta_Tet)) DEALLOCATE(DISC%Galerkin%ADGeta_Tet)
    IF (ASSOCIATED( DISC%Galerkin%ADGzeta_Tet)) DEALLOCATE(DISC%Galerkin%ADGzeta_Tet)
    IF (ASSOCIATED( DISC%Galerkin%FluxInt_Tet)) DEALLOCATE(DISC%Galerkin%FluxInt_Tet)
    IF (ASSOCIATED( DISC%Galerkin%Kxi_k_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kxi_k_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%Kxi_k_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Keta_k_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Keta_k_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%Keta_k_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_k_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kzeta_k_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%Kzeta_k_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Kxi_m_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kxi_m_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%Kxi_m_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Keta_m_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Keta_m_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%Keta_m_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%Kzeta_m_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%Kzeta_m_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%Kzeta_m_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%ADGxi_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%ADGxi_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%ADGxi_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%ADGeta_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%ADGeta_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%ADGeta_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%ADGzeta_Tet_Sp)) THEN
        CALL CloseSparseTensor3b(DISC%Galerkin%ADGzeta_Tet_Sp)
        DEALLOCATE(DISC%Galerkin%ADGzeta_Tet_Sp)
    ENDIF
    IF (ASSOCIATED( DISC%Galerkin%FluxInt_Tet_Sp)) THEN
        DO iSide = 1, MESH%nSides_Tet
            CALL CloseSparseTensor3b(DISC%Galerkin%FluxInt_Tet_Sp(0,1,iSide))
            DO iLocalNeighborSide = 1, MESH%nSides_Tet
                DO iLocalNeighborVrtx = 1, MESH%nVertices_Tri
                    CALL CloseSparseTensor3b(DISC%Galerkin%FluxInt_Tet_Sp(iLocalNeighborSide,iLocalNeighborVrtx,iSide))
                ENDDO
            ENDDO
        ENDDO
        DEALLOCATE(DISC%Galerkin%FluxInt_Tet_Sp)
    ENDIF
    CONTINUE
    !
    DISC%Galerkin%init = .FALSE.
    !
    logInfo(*) 'closeGalerkin successful '
    !
  END SUBROUTINE closeGalerkin3D_us

!
!   unified solver uses routine_new
!

  !===========================================================================!
  !!                                                                         !!
  !! iniGalerkin3D_us_level1 initializes all essential DG variables that     !!
  !! do NOT depend on the mesh                                               !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE iniGalerkin3D_us_level1_new(OptionalFields,EQN,DISC,MESH,BND,IC,SOURCE,MPI,IO)

    USE DGBasis_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tUnstructOptionalFields)   :: OptionalFields
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tDiscretization)           :: DISC
    TYPE(tBoundary)                 :: BND
    TYPE(tInitialCondition)         :: IC
    TYPE(tSource)                   :: SOURCE
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i, j, iElem, iDirac, iRicker
    REAL                            :: x,y,z
    REAL                            :: pi, earth_rad
    REAL                            :: e1(3), e2(3), e3(3)
    ! ------------------------------------------------------------------------!

    DISC%Galerkin%init = .TRUE.

    DISC%Galerkin%init     = .TRUE.
    DISC%Galerkin%nDegFr   = (DISC%Galerkin%nPoly+1)*(DISC%Galerkin%nPoly+2)*(DISC%Galerkin%nPoly+3)/6   !**3

    DISC%Galerkin%nPolyRec  = DISC%Galerkin%nPoly
    DISC%Galerkin%nDegFrRec = DISC%Galerkin%nDegFr


    DISC%Galerkin%nDegFrST = DISC%Galerkin%nDegFrRec*(DISC%Galerkin%nPolyRec+1)

    logInfo(*) 'Interface SEISSOL successful '

  END SUBROUTINE iniGalerkin3D_us_level1_new

  !===========================================================================!
  !!                                                                         !!
  !! iniGalerkin3D_us_level2 initializes all essential DG variables that     !!
  !! DO depend on the mesh                                                   !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE iniGalerkin3D_us_level2_new(OptionalFields,EQN,DISC,MESH,BND,IC,SOURCE,MPI,IO)

    USE DGBasis_mod
!    USE DGSponge_mod ! not yet done for hybrids/unified version
#ifdef PARALLEL
    USE MPIExchangeValues_mod
#endif
#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc, c_null_char, c_bool
    use f_ftoc_bind_interoperability
    use calc_deltaT_mod
#endif
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tUnstructOptionalFields)   :: OptionalFields
    TYPE(tEquations)                :: EQN
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tDiscretization)           :: DISC
    TYPE(tBoundary)                 :: BND
    TYPE(tInitialCondition)         :: IC
    TYPE(tSource)                   :: SOURCE
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i, j, k, l, iElem, iDirac, iRicker
    INTEGER                         :: iDRFace
    INTEGER                         :: iCurElem
    INTEGER                         :: nDGWorkVar
    INTEGER                         :: allocstat
    INTEGER                         :: iIntGP,iTimeGP,NestSize,BlockSize,iFace
    INTEGER                         :: iDegFr_xi,iDegFr_tau,iDegFr_tau2
    REAL                            :: x,y,z
    REAL                            :: pi_const, earth_rad
    REAL                            :: e1(3), e2(3), e3(3)
    REAL                            :: xiGP,etaGP,zetaGP,tauGP
    REAL                            :: phi_i,psi_i,phi_j,psi_j
    REAL                            :: phi_xi_i(3), phi_xi_j(3)
    REAL                            :: psi_tau_i, psi_tau_j,psi_tau_grad_i2
    REAL                            :: TimeStiff(DISC%Galerkin%nPoly+1,DISC%Galerkin%nPoly+1)
    REAL                            :: TimeMass(DISC%Galerkin%nPoly+1,DISC%Galerkin%nPoly+1)
    REAL                            :: TimeF0(DISC%Galerkin%nPoly+1)
    INTEGER, POINTER                :: MPI_Dirac_Element(:,:)
    INTEGER                         :: iError,iCPU
#ifdef GENERATEDKERNELS
    real                            :: l_timeStepWidth
    real                            :: l_loads(3), l_scalings(3), l_cuts(2), l_timeScalings(2), l_gts
    integer                         :: iObject, iSide, iNeighbor, MPIIndex
    real, target                    :: materialVal(EQN%nBackgroundVar)
    logical(kind=c_bool)                        :: enableFreeSurfaceIntegration
    
#endif
    ! ------------------------------------------------------------------------!
    !
    CALL iniGalerkin3D_us_intern_new(EQN, DISC, MESH, BND, IC, SOURCE, OptionalFields, IO)
    !
    CALL BuildSpecialDGGeometry3D_new(OptionalFields%BackgroundValue,EQN,MESH,DISC,BND,MPI,IO)

#ifdef GENERATEDKERNELS
    ! we have the material parameters and insphere diameters, let's get some time step widths out of this

    ! malloc fortran arrays
    call ini_calc_deltaT( eqn%eqType,     &
                          optionalFields, &
                          eqn,            &
                          mesh,           &
                          io                )

    ! get the time step width for every tet
    call cfl_step( optionalFields, &
                   eqn,            &
                   mesh,           &
                   disc,           &
                   io                )

    ! get gts time step width
    l_gts = minval( optionalFields%dt_convectiv(:) )
    if (l_gts .le. 0.0) then
      logError(*) 'Invalid timestep width'
      stop
    endif

#ifdef PERIODIC_LTS_SCALING
    ! compute total load per half-sapce
    !         _____________________
    !        /         /     /    /|
    !       /         /     /    / |
    !      /         /     /    /  |
    !     /         /     /    /   |
    !    /         /     /    /    |
    !   /_________/_____/____/     |
    !  |         |     |     |     |
    !  |    T    |  S  |  F  |     |
    !  |    h    |  e  |  i  |     |
    !  |    i    |  c  |  r  |     |_ _ _ _ _ _ _ _ _ _
    !  |    r    |  o  |  s  |    /                    /
    !  |    d    |  n  |  t  |   /                    /
    !  |         |  d  |     |  /       Symmetry     /
    !  |         |     |     | /                    /
    !  |_________|_____|_____|/ _ _ _ _ _ _ _ _ _ _/
    !  |s^2/2 + s|  s  |s^2/2|     |     |         |
    !  |         |     |     |     |     |         |
    ! -50     -cut2  -cut1   0    cut1  cut2      50

    ! derive scaling of time step widths
    l_timeScalings(1) = PERIODIC_LTS_SCALING * PERIODIC_LTS_SCALING
    l_timeScalings(2) = PERIODIC_LTS_SCALING

    ! derive scalings
    l_scalings(1) = l_timeScalings(2) / 2
    l_scalings(2) = 1
    l_scalings(3) = l_timeScalings(1) + l_timeScalings(2)

    ! derive cuts
    l_cuts(1) = 50.0d0 * ( l_scalings(1)                 ) / ( l_scalings(1) + l_scalings(2) + l_scalings(3) )
    l_cuts(2) = 50.0d0 * ( l_scalings(1) + l_scalings(2) ) / ( l_scalings(1) + l_scalings(2) + l_scalings(3) )
#endif

    ! propagate the time step width to time manager
    do iElem = 1, mesh%nElem
      l_timeStepWidth = optionalFields%dt_convectiv(iElem)

#ifdef PERIODIC_LTS_SCALING
      ! perform the scaling of the time step width
      if(     mesh%elem%xybary(1,iElem) > -l_cuts(2) .and.  mesh%elem%xybary(1,iElem) < -l_cuts(1) .or. \
              mesh%elem%xybary(1,iElem) >  l_cuts(1) .and.  mesh%elem%xybary(1,iElem) <  l_cuts(2) )    \
      then
        l_timeStepWidth = l_gts / ( l_timeScalings(2) - 1.0E-2 )
      elseif( mesh%elem%xybary(1,iElem) > -l_cuts(1) .and.  mesh%elem%xybary(1,iElem) <  0.d0 .or.      \
              mesh%elem%xybary(1,iElem) >  0.d0      .and.  mesh%elem%xybary(1,iElem) <  l_cuts(1) )    \
      then
        l_timeStepWidth = l_gts / ( l_timeScalings(1) + 1.0E-2 )
      else
        l_timeStepWidth = l_gts
      endif
#endif

      call c_interoperability_setTimeStepWidth( i_meshId        = iElem,          &
                                                i_timeStepWidth = l_timeStepWidth &
                                              )
    enddo

    enableFreeSurfaceIntegration = (io%surfaceOutput > 0)
    ! put the clusters under control of the time manager
    call c_interoperability_initializeClusteredLts( i_clustering = disc%galerkin%clusteredLts, i_enableFreeSurfaceIntegration = enableFreeSurfaceIntegration )
#endif

    !
    SELECT CASE(DISC%Galerkin%DGMethod)
    CASE(3)
        nDGWorkVar = EQN%nVar+EQN%nAneFuncperMech
    CASE DEFAULT
        nDGWorkVar = EQN%nVarTotal
    END SELECT
    !
    ! Allocation of arrays
    ALLOCATE(                                                                                   &
         DISC%Galerkin%dgvar( DISC%Galerkin%nDegFr,EQN%nVarTotal,MESH%nElem,DISC%Galerkin%nRK), &

#ifndef GENERATEDKERNELS
         DISC%Galerkin%DGwork(DISC%Galerkin%nDegFrRec,nDGWorkVar,MESH%nElem),                   &
#endif
         STAT = allocstat                                                                       )
    IF(allocStat .NE. 0) THEN
       logError(*) 'could not allocate all variables!'
       STOP
    END IF
    !
    IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ALLOCATE( DISC%Galerkin%DGTaylor(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly,MESH%nElem), &
                  STAT = allocstat )
        IF(allocStat .NE. 0) THEN
           logError(*) 'could not allocate DISC%Galerkin%DGTayl.'
           STOP
        END IF
    ENDIF

! GENERATEDKERNELS: We eliminate double source locations in the C++ part
#ifndef GENERATEDKERNELS
    !
    ! Initialize Dirac sources in space
    !
    SELECT CASE(SOURCE%Type)
    CASE(15)    ! Dirac in space and time
       !
       SOURCE%Dirac%Element(:) = -1
       DO iElem = 1, MESH%nElem
          DO iDirac = 1, SOURCE%Dirac%nDirac
             x = SOURCE%Dirac%SpacePosition(1,iDirac)
             y = SOURCE%Dirac%SpacePosition(2,iDirac)
             z = SOURCE%Dirac%SpacePosition(3,iDirac)
             IF(XYZInElement(x,y,z,iElem,1e-5,MESH,DISC)) THEN
                SOURCE%Dirac%Element(iDirac) = iElem
             ENDIF
          ENDDO
       ENDDO
       !
#ifdef PARALLEL
       IF(MPI%nCPU.GT.1) THEN
           logInfo(*) 'Cleaning possible double source locations for MPI... '
           ALLOCATE(  MPI_Dirac_Element(SOURCE%Dirac%nDirac,0:MPI%nCPU-1) )
             CALL MPI_ALLGATHER(SOURCE%Dirac%Element,SOURCE%Dirac%nDirac,MPI_INTEGER, &
                                MPI_Dirac_Element,   SOURCE%Dirac%nDirac,MPI_INTEGER, &
                                MPI%commWorld, iError                                )
           !
           DO iDirac = 1, SOURCE%Dirac%nDirac
              IF(SOURCE%Dirac%Element(iDirac).NE.-1) THEN
                DO iCPU = 0, MPI%myrank-1
                    IF(MPI_Dirac_Element(iDirac,iCPU).NE.-1) THEN
                        logInfo(*) 'Source number ', iDirac, ' cleaned in element ',  &
                                                    SOURCE%Dirac%Element(iDirac), ' in CPU ', MPI%myrank
                        logInfo(*) '  It was already found in element ', MPI_Dirac_Element(iDirac,iCPU), &
                                                      ' in CPU ', iCPU
                        SOURCE%Dirac%Element(iDirac) = -1
                    ENDIF
                ENDDO
              ENDIF
           ENDDO
           !
           DEALLOCATE(MPI_Dirac_Element )
           logInfo(*) 'MPI source cleaning done.  '
       ENDIF
#endif ! PARALLEL
       !
    CASE(16,18)   ! Dirac in space, Ricker or Gaussian wavelet in time
       !
       SOURCE%Ricker%Element(:) = -1
       DO iElem = 1, MESH%nElem
          DO iRicker = 1, SOURCE%Ricker%nRicker
             x = SOURCE%Ricker%SpacePosition(1,iRicker)
             y = SOURCE%Ricker%SpacePosition(2,iRicker)
             z = SOURCE%Ricker%SpacePosition(3,iRicker)
             IF(XYZInElement(x,y,z,iElem,1e-5,MESH,DISC)) THEN
                SOURCE%Ricker%Element(iRicker) = iElem
             ENDIF
          ENDDO
       ENDDO
       !
#ifdef PARALLEL
       IF(MPI%nCPU.GT.1) THEN
           logInfo(*) 'Cleaning possible double source locations for MPI... '
           ALLOCATE(  MPI_Dirac_Element(SOURCE%Ricker%nRicker,0:MPI%nCPU-1) )
             CALL MPI_ALLGATHER(SOURCE%Ricker%Element,SOURCE%Ricker%nRicker,MPI_INTEGER, &
                                MPI_Dirac_Element,    SOURCE%Ricker%nRicker,MPI_INTEGER, &
                                MPI%commWorld, iError                                          )
           !
           DO iDirac = 1, SOURCE%Ricker%nRicker
              IF(SOURCE%Ricker%Element(iDirac).NE.-1) THEN
                DO iCPU = 0, MPI%myrank-1
                    IF(MPI_Dirac_Element(iDirac,iCPU).NE.-1) THEN
                        logInfo(*) 'Source number ', iDirac, ' cleaned in element ',  &
                                                    SOURCE%Ricker%Element(iDirac), ' in CPU ', MPI%myrank
                        logInfo(*) '  It was already found in element ', MPI_Dirac_Element(iDirac,iCPU), &
                                                      ' in CPU ', iCPU
                        SOURCE%Ricker%Element(iDirac) = -1
                    ENDIF
                ENDDO
              ENDIF
           ENDDO
           !
           DEALLOCATE(MPI_Dirac_Element )
           logInfo(*) 'MPI source cleaning done.  '
       ENDIF
#endif ! PARALLEL
       !
    CASE(20)   ! Dirac in space, Ricker wavelet in time
       !
       SOURCE%RP%Element(:) = -1
       DO iElem = 1, MESH%nElem
          DO iRicker = 1, SOURCE%Ricker%nRicker
             x = SOURCE%RP%SpacePosition(1,iRicker)
             y = SOURCE%RP%SpacePosition(2,iRicker)
             z = SOURCE%RP%SpacePosition(3,iRicker)
             IF(XYZInElement(x,y,z,iElem,1e-5,MESH,DISC)) THEN
                SOURCE%RP%Element(iRicker) = iElem
             ENDIF
          ENDDO
       ENDDO
       !
#ifdef PARALLEL
       IF(MPI%nCPU.GT.1) THEN
           logInfo(*) 'Cleaning possible double source locations for MPI... '
           ALLOCATE(  MPI_Dirac_Element(SOURCE%Ricker%nRicker,0:MPI%nCPU-1) )
             CALL MPI_ALLGATHER(SOURCE%RP%Element,SOURCE%Ricker%nRicker,MPI_INTEGER, &
                                MPI_Dirac_Element,    SOURCE%Ricker%nRicker,MPI_INTEGER, &
                                MPI%commWorld, iError                                          )
           !
           DO iDirac = 1, SOURCE%Ricker%nRicker
              IF(SOURCE%RP%Element(iDirac).NE.-1) THEN
                DO iCPU = 0, MPI%myrank-1
                    IF(MPI_Dirac_Element(iDirac,iCPU).NE.-1) THEN
                        logInfo(*) 'Source number ', iDirac, ' cleaned in element ',  &
                                                    SOURCE%RP%Element(iDirac), ' in CPU ', MPI%myrank
                        logInfo(*) '  It was already found in element ', MPI_Dirac_Element(iDirac,iCPU), &
                                                      ' in CPU ', iCPU
                        SOURCE%RP%Element(iDirac) = -1
                    ENDIF
                ENDDO
              ENDIF
           ENDDO
           !
           DEALLOCATE(MPI_Dirac_Element )
           logInfo(*) 'MPI source cleaning done.  '
       ENDIF
#endif ! PARALLEL
    CASE(30)   ! Finite Source Rupture Model
       !
       ! Transformation of local Rupture Plane coordinates from input file
       ! into epicentral cartesian xyz-coordinates
       !
       ALLOCATE( SOURCE%RP%RS(3,3), SOURCE%RP%RD(3,3),       &
                 SOURCE%RP%RY(3,3), SOURCE%RP%RZ(3,3),       &
                 SOURCE%RP%Rglob(3,3),SOURCE%RP%RglobT(3,3), &
                 SOURCE%RP%n_strk(3), SOURCE%RP%n_dip(3),    &
                 SOURCE%RP%corner(3) )
       ALLOCATE( MESH%IncludesFSRP(MESH%nElem) )
       ALLOCATE( SOURCE%RP%TWindowStart(SOURCE%RP%nTWindow), &
                 SOURCE%RP%TWindowEnd(SOURCE%RP%nTWindow) )

       logInfo(*) 'Computing epicentral coordinates of the rupture plane '

       SOURCE%RP%MomentTensor(:,:) = 0.

       pi_const        = ACOS(-1.0)

       SELECT CASE(SOURCE%RP%Type)


       !CASE(1)      !Only one Rupture Plane Segment (model in earth-centered geographical coordinates)
                     !############################DOES NOT WORK NOW #####################

       CASE(2,3,4)   ! Only one Rupture Plane Segment (model in epicentral coordinates)

           SOURCE%RP%strk          =  SOURCE%RP%strk   /180.*pi_const
           SOURCE%RP%dip           =  SOURCE%RP%dip    /180.*pi_const
           SOURCE%RP%rake(:,:)     =  SOURCE%RP%rake(:,:)/180.*pi_const

           ! scaling of distances into meters
           SOURCE%Hypocenter%depth = SOURCE%Hypocenter%depth * 1000
           SOURCE%RP%length        = SOURCE%RP%length * 1000
           SOURCE%RP%width         = SOURCE%RP%width  * 1000
           SOURCE%RP%HypX          = SOURCE%RP%HypX   * 1000
           SOURCE%RP%HypZ          = SOURCE%RP%HypZ   * 1000
           SOURCE%RP%dxRP          = SOURCE%RP%dxRP   * 1000
           SOURCE%RP%dzRP          = SOURCE%RP%dzRP   * 1000

           ! Matrix for rotation around z-axis by strike-angle
           ! to bring strike-direction and x-axis in coincidence

           SOURCE%RP%RZ(:,:) =  0.
           SOURCE%RP%RZ(1,1) = -COS(SOURCE%RP%strk)
           SOURCE%RP%RZ(1,2) =  SIN(SOURCE%RP%strk)
           SOURCE%RP%RZ(2,1) = -SIN(SOURCE%RP%strk)
           SOURCE%RP%RZ(2,2) = -COS(SOURCE%RP%strk)
           SOURCE%RP%RZ(3,3) =  1.

           ! Matrix for rotation around x-axis by dip-angle
           ! to bring rupture plane and xy-plane in coincidence

           SOURCE%RP%RS(:,:) =  0.
           SOURCE%RP%RS(1,1) =  1.
           SOURCE%RP%RS(2,2) =  COS(SOURCE%RP%dip)
           SOURCE%RP%RS(2,3) =  SIN(SOURCE%RP%dip)
           SOURCE%RP%RS(3,2) = -SIN(SOURCE%RP%dip)
           SOURCE%RP%RS(3,3) =  COS(SOURCE%RP%dip)

           e1 = (/1.,0.,0./)    ! strike   vector in the rupture plane system
           e2 = (/0.,-1.,0./)   ! down-dip vector in the rupture plane system

           ! unit vectors along strike and dip direction in epicentral reference xyz-coordinate system
           SOURCE%RP%RD = MATMUL(SOURCE%RP%RS,SOURCE%RP%RZ)
           ! Compute transpose of RD
           SOURCE%RP%RS      = SOURCE%RP%RD
           SOURCE%RP%RS(1,2) = SOURCE%RP%RD(2,1)
           SOURCE%RP%RS(1,3) = SOURCE%RP%RD(3,1)
           SOURCE%RP%RS(2,1) = SOURCE%RP%RD(1,2)
           SOURCE%RP%RS(2,3) = SOURCE%RP%RD(3,2)
           SOURCE%RP%RS(3,1) = SOURCE%RP%RD(1,3)
           SOURCE%RP%RS(3,2) = SOURCE%RP%RD(2,3)

           SOURCE%RP%n_strk = MATMUL(SOURCE%RP%RS,e1)
           SOURCE%RP%n_dip  = MATMUL(SOURCE%RP%RS,e2)

           ! Local position of epicenter in epicentral reference xyz-coordinate system
           x = 0.
           y = 0.
           z = -SOURCE%Hypocenter%depth

           ! Position of top left corner (origin of the local rupture-plane-coordinate system)
           SOURCE%RP%corner(1) = x - SOURCE%RP%HypX * SOURCE%RP%n_strk(1) - SOURCE%RP%HypZ * SOURCE%RP%n_dip(1)
           SOURCE%RP%corner(2) = y - SOURCE%RP%HypX * SOURCE%RP%n_strk(2) - SOURCE%RP%HypZ * SOURCE%RP%n_dip(2)
           SOURCE%RP%corner(3) = z - SOURCE%RP%HypX * SOURCE%RP%n_strk(3) - SOURCE%RP%HypZ * SOURCE%RP%n_dip(3)

           ! Computation of subfault locations for Dirac sources in epicentral reference xyz-coordinate system
           iDirac = 0
           DO j = 1, SOURCE%RP%nzRP
              DO i = 1,SOURCE%RP%nxRP
                 iDirac = iDirac + 1
                 SOURCE%RP%SpacePosition(1,iDirac) = SOURCE%RP%corner(1) + (i-0.5)*SOURCE%RP%dxRP*SOURCE%RP%n_strk(1)   &
                                                                         + (j-0.5)*SOURCE%RP%dzRP*SOURCE%RP%n_dip(1)
                 SOURCE%RP%SpacePosition(2,iDirac) = SOURCE%RP%corner(2) + (i-0.5)*SOURCE%RP%dxRP*SOURCE%RP%n_strk(2)   &
                                                                         + (j-0.5)*SOURCE%RP%dzRP*SOURCE%RP%n_dip(2)
                 SOURCE%RP%SpacePosition(3,iDirac) = SOURCE%RP%corner(3) + (i-0.5)*SOURCE%RP%dxRP*SOURCE%RP%n_strk(3)   &
                                                                         + (j-0.5)*SOURCE%RP%dzRP*SOURCE%RP%n_dip(3)
              ENDDO
           ENDDO

       CASE DEFAULT
          logError(*)  'The format type of the Finite Source Rupture Model is unknown! '
          STOP                                     ! STOP

       END SELECT

       SELECT CASE(SOURCE%RP%Type)
       CASE(2)
           ! Source time function based on triangular approximation for each time window
           ! Pre-computation of start and stop of each triangular in time
           DO i = 1,SOURCE%RP%nTWindow
               SOURCE%RP%TWindowStart(i) = (i-1)*(SOURCE%RP%TWindowLen - SOURCE%RP%WindowShift)
               SOURCE%RP%TWindowEnd(i)   = i*SOURCE%RP%TWindowLen - (i-1)*SOURCE%RP%WindowShift
           ENDDO
           ! Computation of according slip rates
           IF(SOURCE%RP%TWindowLen.GT.0.)THEN
               SOURCE%RP%Sliprate(:,:)   = SOURCE%RP%Slip(:,:)/SOURCE%RP%TWindowLen
           ENDIF

       CASE(3,4)
           ! Computation of slip rates
           SOURCE%RP%Sliprate(:,1)       = SQRT( SOURCE%RP%SSlip(:,1)*SOURCE%RP%SSlip(:,1) +   &
                                                 SOURCE%RP%DSlip(:,1)*SOURCE%RP%DSlip(:,1)   ) /SOURCE%RP%TRise(:)

       CASE DEFAULT
           logError(*)  'The format type of the Finite Source Rupture Model is unknown! '
           STOP                                     ! STOP

       END SELECT

       ! Search including elements
       SOURCE%RP%Element(:) = -1
       MESH%IncludesFSRP(:) = .FALSE.
       DO iDirac = 1, SOURCE%RP%nSbfs(1)
           x = SOURCE%RP%SpacePosition(1,iDirac)
           y = SOURCE%RP%SpacePosition(2,iDirac)
           z = SOURCE%RP%SpacePosition(3,iDirac)
           DO iElem = 1, MESH%nElem
             IF(XYZInElement(x,y,z,iElem,1e-5,MESH,DISC)) THEN
                SOURCE%RP%Element(iDirac) = iElem
                MESH%IncludesFSRP(iElem)  = .TRUE.
                !PRINT *, '## Ricker nr. ', iDirac, ' found in element ', iElem, ' in CPU ', MPI%myrank
                !PRINT *, '   Element barycenter : ', MESH%ELEM%xyBary(:,iElem)
                EXIT
             ENDIF
           ENDDO
       ENDDO
       !
       DEALLOCATE ( SOURCE%RP%RS     )
       DEALLOCATE ( SOURCE%RP%RD     )
       DEALLOCATE ( SOURCE%RP%RY     )
       DEALLOCATE ( SOURCE%RP%RZ     )
       DEALLOCATE ( SOURCE%RP%n_strk )
       DEALLOCATE ( SOURCE%RP%n_dip  )
       DEALLOCATE ( SOURCE%RP%corner )

#ifdef PARALLEL
       IF(MPI%nCPU.GT.1) THEN
           logInfo(*) 'Cleaning possible double source locations for MPI... '
           ALLOCATE(  MPI_Dirac_Element(SOURCE%RP%nSbfs(1),0:MPI%nCPU-1) )
           MPI_Dirac_Element = -1
             CALL MPI_ALLGATHER(SOURCE%RP%Element, SOURCE%RP%nSbfs(1),MPI_INTEGER, &
                                MPI_Dirac_Element, SOURCE%RP%nSbfs(1),MPI_INTEGER, &
                                MPI%commWorld, iError                                             )
           !
           DO iDirac = 1, SOURCE%RP%nSbfs(1)
              IF(SOURCE%RP%Element(iDirac).NE.-1) THEN
                DO iCPU = 0, MPI%myrank-1
                    IF(MPI_Dirac_Element(iDirac,iCPU).NE.-1) THEN
                        logInfo(*) 'Source number ', iDirac, ' cleaned in element ',  &
                                                    SOURCE%RP%Element(iDirac), ' in CPU ', MPI%myrank
                        logInfo(*) '   It was already found in element ', MPI_Dirac_Element(iDirac,iCPU), &
                                                       ' in CPU ', iCPU
                        SOURCE%RP%Element(iDirac) = -1
                        MESH%IncludesFSRP(iElem)  = .FALSE.
                    ENDIF
                ENDDO
              ENDIF
           ENDDO
           !
           DEALLOCATE(MPI_Dirac_Element )
           logInfo(*) 'MPI source cleaning done.  '
       ENDIF
#endif ! PARALLEL
    CASE(50)

       ALLOCATE( MESH%IncludesFSRP(MESH%nElem) )

       logInfo(*) 'Check which elements contain FSRM subfaults ... '

       ! Search including elements
       SOURCE%RP%Element(:) = -1
       MESH%IncludesFSRP(:) = .FALSE.
       DO iDirac = 1, SOURCE%RP%nSbfs(1)
           x = SOURCE%RP%SpacePosition(1,iDirac)
           y = SOURCE%RP%SpacePosition(2,iDirac)
           z = SOURCE%RP%SpacePosition(3,iDirac)
           DO iElem = 1, MESH%nElem
             IF(XYZInElement(x,y,z,iElem,1e-5,MESH,DISC)) THEN
                SOURCE%RP%Element(iDirac) = iElem
                MESH%IncludesFSRP(iElem)  = .TRUE.
                EXIT
             ENDIF
           ENDDO
       ENDDO
       !
       logInfo(*) '... finished! '
       !
#ifdef PARALLEL
       IF(MPI%nCPU.GT.1) THEN
           logInfo(*) 'Cleaning possible double source locations for MPI... '
           ALLOCATE(  MPI_Dirac_Element(SOURCE%RP%nSbfs(1),0:MPI%nCPU-1) )
           MPI_Dirac_Element = -1
             CALL MPI_ALLGATHER(SOURCE%RP%Element, SOURCE%RP%nSbfs(1),MPI_INTEGER, &
                                MPI_Dirac_Element, SOURCE%RP%nSbfs(1),MPI_INTEGER, &
                                MPI%commWorld, iError                                             )
           !
           DO iDirac = 1, SOURCE%RP%nSbfs(1)
              IF(SOURCE%RP%Element(iDirac).NE.-1) THEN
                DO iCPU = 0, MPI%myrank-1
                    IF(MPI_Dirac_Element(iDirac,iCPU).NE.-1) THEN
                        logInfo(*) 'Source number ', iDirac, ' cleaned in element ',  &
                                                    SOURCE%RP%Element(iDirac), ' in CPU ', MPI%myrank
                        logInfo(*) '   It was already found in element ', MPI_Dirac_Element(iDirac,iCPU), &
                                                     ' in CPU ', iCPU
                        SOURCE%RP%Element(iDirac) = -1
                        MESH%IncludesFSRP(iElem)  = .FALSE.
                    ENDIF
                ENDDO
              ENDIF
           ENDDO
           !
           DEALLOCATE(MPI_Dirac_Element )
           logInfo(*) 'MPI source cleaning done.  '
       ENDIF
#endif ! PARALLEL
       !
    END SELECT
#endif ! GENERATEDKERNELS

#ifndef GENERATEDKERNELS
    logInfo0(*) 'Initializing DR parallelization'
    DISC%DynRup%nDRElems = 0

    IF(EQN%DR.EQ.1) THEN
      DO iDRFace=1,mesh%fault%nSide
        IF(MESH%Fault%Face(iDRFace,1,1).GT.0) THEN
          ! increase face counter
          DISC%DynRup%nDRElems = DISC%DynRup%nDRElems + 1
        ENDIF
        IF(MESH%Fault%Face(iDRFace,1,2).GT.0) THEN
          ! increase face counter
          DISC%DynRup%nDRElems = DISC%DynRup%nDRElems + 1
        ENDIF
      ENDDO
    ENDIF

    logInfo(*) 'Elementfaces found for DR: ', DISC%DynRup%nDRElems

    ! Setup faced-based DR update datascrutured needed for OMP Parallelization in fricition.f90
    IF (DISC%DynRup%nDRElems.NE.0) THEN
      ! allocate the memory for the offload datastructures
      allocate( DISC%DynRup%DRupdates(DISC%Galerkin%nDegFrRec,EQN%nVar, DISC%DynRup%nDRElems) )
      allocate( DISC%DynRup%indicesOfDRElems( DISC%DynRup%nDRElems ) )
      allocate( DISC%DynRup%DRupdatesPosition( mesh%fault%nSide, 2 ) )

      DISC%DynRup%nDRElems = 0
      DO iDRFace=1,mesh%fault%nSide
        ! reset translation array for DR updates for current face
        DISC%DynRup%DRupdatesPosition( iDRFace, : ) = 0
        IF(MESH%Fault%Face(iDRFace,1,1).GT.0) THEN
          DISC%DynRup%nDRElems = DISC%DynRup%nDRElems + 1
          ! add element index to gather datastructure
          DISC%DynRup%indicesOfDRElems( DISC%DynRup%nDRElems ) = MESH%Fault%Face(iDRFace,1,1)
          ! update translation array for DR updates
          DISC%DynRup%DRupdatesPosition( iDRFace, 1 ) = DISC%DynRup%nDRElems
        ENDIF
        IF(MESH%Fault%Face(iDRFace,1,2).GT.0) THEN
          DISC%DynRup%nDRElems = DISC%DynRup%nDRElems + 1
          ! add element index to gather datastructure
          DISC%DynRup%indicesOfDRElems( DISC%DynRup%nDRElems ) = MESH%Fault%Face(iDRFace,1,2)
          ! update translation array for DR updates
          DISC%DynRup%DRupdatesPosition( iDRFace, 2 ) = DISC%DynRup%nDRElems
        ENDIF
      ENDDO
    ELSE
      ! allocate one element to have a valid pointer for offloading
      allocate( DISC%DynRup%DRupdates(DISC%Galerkin%nDegFrRec,EQN%nVarTotal, 1) )
      allocate( DISC%DynRup%indicesOfDRElems( 1 ) )
    ENDIF

    logInfo0(*) 'Initializing DR parallelization. Done.'
#endif

    !
    !CALL ini_DGSponge(EQN,DISC,MESH,IO)  ! not yet done for hybrids
    !
    !
#ifdef PARALLEL
    IF(MPI%nCPU.GT.1) THEN                                          !
                                                                    !
    ! We comunicate background fields                               !
    ! Needed for discontinuous Riemann problems                     !
    !                                                               !
    CALL MPIExchangeBackground(DISC           = DISC,             & !
                               EQN            = EQN,              & !
                               BND            = BND,              & !
                               MESH           = MESH,             & !
                               IO             = IO,               & !
                               OptionalFields = OptionalFields,   & !
                               MPI            = MPI               ) !
                                                                    !
    ENDIF                                                           !
#endif
    !
    ! Initialize sparse star matrices
    !
#if defined(GENERATEDKERNELS)
  do iElem = 1, MESH%nElem
    iSide = 0

    materialVal = OptionalFields%BackgroundValue(iElem,:)
    call c_interoperability_setMaterial( i_elem = iElem,                                         \
                                         i_side = iSide,                                         \
                                         i_materialVal = materialVal,\
                                         i_numMaterialVals = EQN%nBackgroundVar                  )

    do iSide = 1,4
      IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
          iObject         = MESH%ELEM%BoundaryToObject(iSide,iElem)
          MPIIndex        = MESH%ELEM%MPINumber(iSide,iElem)
          materialVal = BND%ObjMPI(iObject)%NeighborBackground(1:EQN%nBackgroundVar,MPIIndex) ! rho,mu,lambda
      ELSE
          SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
          CASE(0)
              iNeighbor       = MESH%ELEM%SideNeighbor(iSide,iElem)
              materialVal = OptionalFields%BackgroundValue(iNeighbor,:)
          CASE DEFAULT ! For boundary conditions take inside material
              materialVal = OptionalFields%BackgroundValue(iElem,:)
          END SELECT
      ENDIF
      call c_interoperability_setMaterial( i_elem = iElem,                        \
                                           i_side = iSide,                        \
                                           i_materialVal = materialVal,       \
                                           i_numMaterialVals = EQN%nBackgroundVar )

    enddo
  enddo

#ifdef USE_MPI
  ! synchronize redundant cell data
  logInfo0(*) 'Synchronizing copy cell material data.';
  call c_interoperability_synchronizeCellLocalData;
#endif

  ! Initialize source terms
  select case(SOURCE%Type)
    case(0)
      ! No source terms
      ! Do nothing
    case(42)
      call c_interoperability_setupNRFPointSources(trim(SOURCE%NRFFileName) // c_null_char)
    case(50)
      call c_interoperability_setupFSRMPointSources( momentTensor     = SOURCE%RP%MomentTensor,  &
                                                     numberOfSources  = SOURCE%RP%nSbfs(1),      &
                                                     centres          = SOURCE%RP%SpacePosition, &
                                                     strikes          = SOURCE%RP%Strks,         &
                                                     dips             = SOURCE%RP%Dips,          &
                                                     rakes            = SOURCE%RP%Rake,          &
                                                     onsets           = SOURCE%RP%Tonset,        &
                                                     areas            = SOURCE%RP%Area,          &
                                                     timestep         = SOURCE%RP%t_samp,        &
                                                     numberOfSamples  = SOURCE%RP%nsteps,        &
                                                     timeHistories    = SOURCE%RP%TimeHist       )
    case default
      logError(*) 'Generated Kernels: Unsupported source type: ', SOURCE%Type
      stop
  end select

  if (DISC%Galerkin%FluxMethod .ne. 0) then
    logError(*) 'Generated kernels currently supports Godunov fluxes only.'
    stop
  endif

  logInfo0(*) 'Initializing element local matrices.'
  call c_interoperability_initializeCellLocalMatrices;
#else
  CALL IniSparseStarMatrices3D_new(EQN, DISC, MESH, BND, SOURCE, OptionalFields%BackgroundValue, IO)
#endif

    IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ALLOCATE( DISC%LocalIteration(MESH%nElem) )
        ALLOCATE( DISC%LocalTime(MESH%nElem)      )
        ALLOCATE( DISC%LocalDt(MESH%nElem)        )
        DISC%LocalIteration(:)  = 0.
        DISC%LocalTime(:)       = 0.
#ifndef GENERATEDKERNELS
        DISC%Galerkin%dgwork    = 0.
#endif
    ENDIF
    !
    IF(EQN%DR.EQ.1) THEN

      ALLOCATE(DISC%DynRup%SlipRate1(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%SlipRate2(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Slip(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Slip1(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Slip2(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%TracXY(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%TracXZ(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%Mu(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%PeakSR(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%rupture_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      ALLOCATE(DISC%DynRup%dynStress_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      
      ! TODO: Transpose StateVar
      ALLOCATE(DISC%DynRup%StateVar(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      !
      DISC%DynRup%SlipRate1     = EQN%IniSlipRate1
      DISC%DynRup%SlipRate2     = EQN%IniSlipRate2
      DISC%DynRup%Slip          = 0.0D0
      DISC%DynRup%Slip1         = 0.0D0
      DISC%DynRup%Slip2         = 0.0D0
      DISC%DynRup%TracXY        = 0.0D0
      DISC%DynRup%TracXZ        = 0.0D0
      DISC%DynRup%Mu(:,:)       = EQN%IniMu(:,:)
      do iFace=1,MESH%Fault%nSide
         DISC%DynRup%StateVar(:,iFace) = EQN%IniStateVar(iFace,:)
      enddo
      DISC%DynRup%PeakSR        = 0.0D0
      DISC%DynRup%rupture_time  = 0.0D0
      DISC%DynRup%dynStress_time = 0.0D0

      allocate(disc%DynRup%output_Mu(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Strength(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Slip(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Slip1(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_Slip2(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_rupture_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_PeakSR(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      allocate(disc%DynRup%output_dynStress_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
      
      allocate(disc%DynRup%output_StateVar(MESH%Fault%nSide,DISC%Galerkin%nBndGP))

    else
        ! Allocate dummy arrays to avoid debug errors
        allocate(DISC%DynRup%SlipRate1(0,0), &
            DISC%DynRup%SlipRate2(0,0),      &
            DISC%DynRup%Slip(0,0),           &
            DISC%DynRup%Slip1(0,0),          &
            DISC%DynRup%Slip2(0,0),          &
            DISC%DynRup%Mu(0,0),             &
            DISC%DynRup%StateVar(0,0),       &
            DISC%DynRup%PeakSR(0,0),         &
            DISC%DynRup%Strength(0,0),       &
            DISC%DynRup%rupture_time(0,0),   &
            DISC%DynRup%dynStress_time(0,0)  )
        allocate(DISC%DynRup%output_Mu(0,0),      &
            DISC%DynRup%output_StateVar(0,0),     &
            DISC%DynRup%output_Strength(0,0),     &
            DISC%DynRup%output_Slip(0,0),         &
            DISC%DynRup%output_Slip1(0,0),        &
            DISC%DynRup%output_Slip2(0,0),        &
            DISC%DynRup%output_rupture_time(0,0), &
            DISC%DynRup%output_PeakSR(0,0),       &
            DISC%DynRup%output_dynStress_time(0,0))
    ENDIF
    !
    IF(DISC%Galerkin%CKMethod.EQ.1) THEN ! not yet done for hybrids
        print*,' ERROR in SUBROUTINE iniGalerkin3D_us_level2_new'
        PRINT*,' DISC%Galerkin%CKMethod.EQ.1 not implemented'
        STOP
        !
    ENDIF
  END SUBROUTINE iniGalerkin3D_us_level2_new


 !===========================================================================!
  !!                                                                         !!
  !!  iniGalerkin3D_us_intern initializes the private data                   !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE iniGalerkin3D_us_intern_new(EQN, DISC, MESH, BND, IC, SOURCE, OptionalFields, IO)
    !-------------------------------------------------------------------------!

    USE DGBasis_mod
    USE ini_QuadFreeDG_mod
    USE COMMON_operators_mod
    USE QuadPoints_mod
    USE JacobiNormal_mod

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tBoundary)          :: BND
    TYPE(tInitialCondition)  :: IC
    TYPE(tSource)            :: SOURCE
    TYPE(tUnstructOptionalFields)   :: OptionalFields
    TYPE(tInputOutput)       :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER :: allocstat                                  ! Allocation status !
    INTEGER :: stat                                       ! IO status         !
    INTEGER :: iDegFr                                     ! Loop counter      !
    INTEGER :: iElem, iSide, iBndGP
    INTEGER :: iNeighbor, iLocalNeighborSide, iNeighborSide, iNeighborVertex
    INTEGER :: iLocalNeighborVrtx
    INTEGER :: nDegFr, MaxDegFr, nDGWorkVar
    INTEGER :: i,j,k,l,m,r,r1,r2,r3                       ! Loop counter      !
    INTEGER :: iPoly, iXi, iEta, iZeta, iIntGP            ! Loop counter      !
    REAL    :: xi, eta, zeta, tau, chi, tau1, chi1
    REAL    :: xiS, etaS, zetaS, xP, yP, zP
    REAL    :: phi, gradphixieta(3)
    REAL                            :: rho,mu,lambda,c0,ce(6,6)
    REAL                            :: x(MESH%nVertexMax)
    REAL                            :: y(MESH%nVertexMax)
    REAL                            :: z(MESH%nVertexMax)
    REAL                            :: JacobiT(EQN%Dimension,EQN%Dimension)
    REAL                            :: auxMatrix(EQN%nVar,EQN%nVar)
    REAL                            :: A_Star(EQN%nVar,EQN%nVar)
    REAL                            :: B_Star(EQN%nVar,EQN%nVar)
    REAL                            :: C_Star(EQN%nVar,EQN%nVar)
    REAL                            :: A(EQN%nVar,EQN%nVar)
    REAL                            :: B(EQN%nVar,EQN%nVar)
    REAL                            :: C(EQN%nVar,EQN%nVar)
    REAL                            :: nx, ny, nz
    REAL                            :: sx, sy, sz
    REAL                            :: tx, ty, tz
    REAL                            :: locA(EQN%nVar,EQN%nVar), locabsA(EQN%nVar,EQN%nVar)
    REAL                            :: T(EQN%nVar,EQN%nVar)
    REAL                            :: iT(EQN%nVar,EQN%nVar)
    !
    INTEGER                         :: nTens3GP, indx
    REAL,POINTER                    :: Tens3BaseFunc(:,:)
    REAL,POINTER                    :: Tens3BaseGrad(:,:,:)
    REAL,POINTER                    :: Tens3GaussP(:,:)
    REAL,POINTER                    :: Tens3GaussW(:)
    INTEGER                         :: nTFMGaussP, ngll
    REAL,POINTER                    :: TFMGaussP(:,:)
    REAL,POINTER                    :: TFMGaussW(:)
    REAL                            :: phi_m, phi_l, phi_k
    REAL                            :: phigrad(EQN%Dimension)
    REAL                            :: chiGP, tauGP
    REAL                            :: VAND(DISC%Galerkin%nPoly+1,DISC%Galerkin%nPoly+1),Temp(DISC%Galerkin%nPoly+1,1)
    REAL                            :: grad(EQN%Dimension,EQN%Dimension)
    !
    LOGICAL                         :: configexist
    CHARACTER(LEN=600)              :: FileName_Tet, FileName_Hex, FileName_Time
    CHARACTER(LEN=200)              :: DGPATH
    ! CHARACTER(LEN=600)              :: TimeFile
    !
    ! Dynamic Rupture variables
    INTEGER                         :: iFace
    INTEGER                         :: iDegFr2
    INTEGER                         :: MPIIndex, iObject
    REAL                            :: xGP,yGP,zGP
    REAL                            :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
    REAL                            :: x_host(MESH%GlobalVrtxType), y_host(MESH%GlobalVrtxType), z_host(MESH%GlobalVrtxType)
    REAL                            :: phi_iDegFr,phi_iDegFr2
    !-------------------------------------------------------------------------!

    IF(.NOT.DISC%Galerkin%init) THEN
       logError(*) 'iniGalerkin: SeisSol Interface not initialized!!'
       STOP
    ENDIF

    IF(MESH%nElem_Tet.EQ.0 .AND. MESH%nElem_Hex.EQ.0) THEN
       logError(*) 'Quadraturefree ADER-DG is only implemented for tetrahedral and hexahedral.'
       STOP
    ENDIF

    ! Reading polynomial coefficients and mass matrices

    INQUIRE(                                            & !
         FILE= 'DGPATH'                               , & !
         EXIST=configexist                              ) !
    !                                                     !
    IF (configexist) THEN                                 !
       !                                                  !
       OPEN(                                            & !
            UNIT= IO%UNIT%FileIn                      , & !
            FILE= 'DGPATH'                            , & !
            IOSTAT = STAT                               ) !
       !                                                  !
       IF (stat.NE.0) THEN                                !
          logError(*) 'cannot open DGPATH'         !
          STOP
       END IF                                             !
       !                                                  !
       READ(IO%UNIT%FileIn,'(A)') DGPATH                  !
       !                                                  !
       CLOSE(IO%UNIT%FileIn)                              !
       !                                                  !
       !                                                  !
       logInfo(*) 'Path to the DG directory is: ',TRIM(DGPATH)

       IF (MESH%nElem_Tet.GT.0) WRITE(FileName_Tet,'(a,a20)') TRIM(DGPATH), 'BasisFunctions3D.tri'
       IF (MESH%nElem_Hex.GT.0) WRITE(FileName_Hex,'(a,a20)') TRIM(DGPATH), 'BasisFunctions3D.qua'

    ELSE                                                  !
       !                                                  !
       logWarning(*) 'Configuration file DGPATH missing!'
       logWarning(*) 'Use . as default path'    !
       !                                                  !
       IF (MESH%nElem_Tet.GT.0) WRITE(FileName_Tet,'(a,a20)') TRIM(DGPATH), 'BasisFunctions3D.tri'
       IF (MESH%nElem_Hex.GT.0) WRITE(FileName_Hex,'(a,a20)') TRIM(DGPATH), 'BasisFunctions3D.qua'
       !
    ENDIF

!===================================================================================!
!--------------------------------Tetrahedral Elements-------------------------------!
!===================================================================================!

    IF(MESH%nElem_Tet .GT. 0)THEN
        logInfo0(*) 'Tetrahedral Elements '
        CALL OpenFile(IO%UNIT%FileIn,FileName_Tet,.FALSE.)

        logInfo0(*) 'Reading basis functions and mass matrices for DG method '
        logInfo0(*) '  from file ', TRIM(FileName_Tet)

        ! Read comment
        READ(IO%UNIT%FileIn,*)
        READ(IO%UNIT%FileIn,*)
        ! Read maximal degree of basis polynomials stored in the file.
        READ(IO%UNIT%FileIn,*) DISC%Galerkin%nMaxPoly
        IF(DISC%Galerkin%nPoly.GT.DISC%Galerkin%nMaxPoly) THEN
            logError(*) 'Required polynomial for DG method is higher than the ones stored in file ', FileName_Tet
            STOP
        ENDIF
        MaxDegFr = (DISC%Galerkin%nPolyRec+1)*(DISC%Galerkin%nPolyRec+2)*(DISC%Galerkin%nPolyRec+3)/6
        ALLOCATE(DISC%Galerkin%cPoly3D_Tet(0:DISC%Galerkin%nPolyRec,0:DISC%Galerkin%nPolyRec,0:DISC%Galerkin%nPolyRec, &
                                           0:MaxDegFr-1, 0:DISC%Galerkin%nPolyRec),                                    &
                 DISC%Galerkin%NonZeroCPoly_Tet(0:MaxDegFr-1,0:DISC%Galerkin%nPolyRec),                                &
                 DISC%Galerkin%NonZeroCPolyIndex_Tet(EQN%Dimension,1:(DISC%Galerkin%nPolyRec+1)**3,0:MaxDegFr,         &
                                                     0:DISC%Galerkin%nPolyRec),                                        &
                 DISC%Galerkin%MassMatrix_Tet(MaxDegFr,MaxDegFr, 0:DISC%Galerkin%nPolyRec),                            &
                 DISC%Galerkin%iMassMatrix_Tet(MaxDegFr,MaxDegFr, 0:DISC%Galerkin%nPolyRec),                           &
                 DISC%Galerkin%IntGaussP_Tet(EQN%Dimension,DISC%Galerkin%nIntGP),                                      &
                 DISC%Galerkin%IntGaussW_Tet(DISC%Galerkin%nIntGP),                                                    &
                 DISC%Galerkin%BndGaussP_Tet(EQN%Dimension-1,DISC%Galerkin%nBndGP),                                      &
                 DISC%Galerkin%BndGaussW_Tet(Disc%Galerkin%nBndGP),                                                    &
                 Tens3GaussP(EQN%Dimension,(DISC%Galerkin%nPoly+2)**3),Tens3GaussW((DISC%Galerkin%nPoly+2)**3),        &
                 DISC%Galerkin%IntGPBaseGrad_Tet(EQN%Dimension,0:MaxDegFr,(DISC%Galerkin%nPolyRec+2)**3,               &
                                                 0:DISC%Galerkin%nPolyRec),                                            &
                 DISC%Galerkin%IntGPBaseFunc_Tet(0:MaxDegFr,(DISC%Galerkin%nPolyRec+2)**3,0:DISC%Galerkin%nPolyRec),   &
!                 DISC%Galerkin%BndGPBaseFunc3D_Tet(0:MaxDegFr,(DISC%Galerkin%nPolyRec+2)**2,MESH%nSides_Tet),          &
#ifndef GENERATEDKERNELS
                 Mesh%ELEM%BndBF_GP_Tet(0:MaxDegFr,DISC%Galerkin%nBndGP,MESH%nSides_Tet),                              &
#endif
                 STAT = allocstat                                                                                      )
        IF(allocStat .NE. 0) THEN
           logError(*) 'could not allocate all variables!'
           STOP
        END IF

        DISC%Galerkin%MassMatrix_Tet(:,:,:)  = 0.0d0
        DISC%Galerkin%iMassMatrix_Tet(:,:,:) = 0.0d0
        DISC%Galerkin%cPoly3D_Tet(:,:,:,:,:) = 0.0d0

        ! Read coefficients of basis functions and mass matrices up to degree DISC%Galerkin%nMaxPoly, stored in the file
        ! BasisFunctions3D.tri
        DO iPoly = 0, DISC%Galerkin%nPolyRec
            ! Read comment in front of the basis functions' coefficients
            logInfo0(*) 'Reading basis functions of order ', iPoly
            READ(IO%UNIT%FileIn,*)
            nDegFr = (iPoly + 1)*(iPoly + 2)*(iPoly + 3)/6
            ! Read polynomial coefficients
            ! where the index of the degrees of freedom starts at zero
            DO iDegFr = 0, nDegFr-1
                DO iZeta = 0, iPoly
                    DO iEta = 0, iPoly
                        DO iXi = 0, iPoly
                            READ(IO%UNIT%FileIn,*) DISC%Galerkin%cPoly3D_Tet(iXi,iEta,iZeta,iDegFr,iPoly)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            ! Read comment in front of the entries of the mass matrix
            READ(IO%UNIT%FileIn,*)
            logInfo0(*) 'Reading mass matrices   of order ', iPoly
            ! Read entries of the mass matrix
            DO k = 1, nDegFr
                DO l = 1, nDegFr
                    READ(IO%UNIT%FileIn,*) DISC%Galerkin%MassMatrix_Tet(k,l,iPoly)
                    IF(DISC%Galerkin%MassMatrix_Tet(k,l,iPoly).NE.0) THEN
                       DISC%Galerkin%iMassMatrix_Tet(k,l,iPoly) = 1.0d0/DISC%Galerkin%MassMatrix_Tet(k,l,iPoly)
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        CLOSE(IO%UNIT%FileIn)
        ! Detecting non-zero coefficients in basis polynomials
        DISC%Galerkin%NonZeroCPoly_Tet(:,:)          = 0
        DISC%Galerkin%NonZeroCPolyIndex_Tet(:,:,:,:) = -1
        DO iPoly = 0, DISC%Galerkin%nPolyRec
            nDegFr = (iPoly + 1)*(iPoly + 2)*(iPoly + 3)/6
            DO iDegFr = 0, nDegFr-1
                DO iZeta = 0, iPoly
                    DO iEta = 0, iPoly
                        DO iXi = 0, iPoly
                            IF(ABS(DISC%Galerkin%cPoly3D_Tet(iXi,iEta,iZeta,iDegFr,iPoly)).GE.1e-6) THEN
                                DISC%Galerkin%NonZeroCPoly_Tet(iDegFr,iPoly) = DISC%Galerkin%NonZeroCPoly_Tet(iDegFr,iPoly) + 1
                                DISC%Galerkin%NonZeroCPolyIndex_Tet(1,DISC%Galerkin%NonZeroCPoly_Tet(iDegFr,iPoly),iDegFr,iPoly) = iXi
                                DISC%Galerkin%NonZeroCPolyIndex_Tet(2,DISC%Galerkin%NonZeroCPoly_Tet(iDegFr,iPoly),iDegFr,iPoly) = iEta
                                DISC%Galerkin%NonZeroCPolyIndex_Tet(3,DISC%Galerkin%NonZeroCPoly_Tet(iDegFr,iPoly),iDegFr,iPoly) = iZeta
                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDIF ! Tets

!===================================================================================!
!--------------------------------Hexahedral Elements--------------------------------!
!===================================================================================!

    IF(MESH%nElem_Hex .GT. 0)THEN
        logInfo0(*) 'Hexahedral Elements '
        CALL OpenFile(IO%UNIT%FileIn,FileName_Hex,.FALSE.)

        logInfo0(*) 'Reading basis functions and mass matrices for DG method '
        logInfo0(*) '  from file ', TRIM(FileName_Hex)

        ! Read comment
        READ(IO%UNIT%FileIn,*)
        READ(IO%UNIT%FileIn,*)
        ! Read maximal degree of basis polynomials stored in the file.
        READ(IO%UNIT%FileIn,*) DISC%Galerkin%nMaxPoly
        IF(DISC%Galerkin%nPoly.GT.DISC%Galerkin%nMaxPoly) THEN
            logError(*) 'Required polynomial for DG method is higher than the ones stored in file ', FileName_Hex
            STOP
        ENDIF
        MaxDegFr = (DISC%Galerkin%nPolyRec+1)*(DISC%Galerkin%nPolyRec+2)*(DISC%Galerkin%nPolyRec+3)/6 !**3

        ALLOCATE(DISC%Galerkin%cPoly3D_Hex(0:DISC%Galerkin%nPolyRec,0:DISC%Galerkin%nPolyRec,0:DISC%Galerkin%nPolyRec,      &
                                           0:MaxDegFr-1, 0:DISC%Galerkin%nPolyRec),                                         &
                 DISC%Galerkin%NonZeroCPoly_Hex(0:MaxDegFr-1,0:DISC%Galerkin%nPolyRec),                                     &
                 DISC%Galerkin%NonZeroCPolyIndex_Hex(EQN%Dimension,1:(DISC%Galerkin%nPolyRec+1)**3,0:MaxDegFr,              &
                                                     0:DISC%Galerkin%nPolyRec),                                             &
                 DISC%Galerkin%MassMatrix_Hex(MaxDegFr,MaxDegFr, 0:DISC%Galerkin%nPolyRec),                                 &
                 DISC%Galerkin%iMassMatrix_Hex(MaxDegFr,MaxDegFr, 0:DISC%Galerkin%nPolyRec),                                &
                 Tens3GaussP(EQN%Dimension,(DISC%Galerkin%nPoly+2)**3),Tens3GaussW((DISC%Galerkin%nPoly+2)**3),             &
                 DISC%Galerkin%IntGaussP_Hex(EQN%Dimension,DISC%Galerkin%nIntGP),                                           &
                 DISC%Galerkin%IntGaussW_Hex(DISC%Galerkin%nIntGP),                                                         &
                 DISC%Galerkin%BndGaussP_Hex(EQN%Dimension,DISC%Galerkin%nBndGP),                                           &
                 DISC%Galerkin%BndGaussW_Hex(DISC%Galerkin%nBndGP),                                                         &
                 DISC%Galerkin%IntGPBaseGrad_Hex(EQN%Dimension,0:MaxDegFr,(DISC%Galerkin%nPolyRec+2)**3,                    &
                                                 0:DISC%Galerkin%nPolyRec),                                                 &
                 DISC%Galerkin%IntGPBaseFunc_Hex(0:MaxDegFr,(DISC%Galerkin%nPolyRec+2)**3,  0:DISC%Galerkin%nPolyRec),      &
                 DISC%Galerkin%BndGPBaseFunc3D_Hex(0:MaxDegFr,(DISC%Galerkin%nPolyRec+2)**2,MESH%nSides_Hex),               &
                 MESH%ELEM%BndBF_GP_Hex(0:MaxDegFr,DISC%Galerkin%nBndGP,MESH%nSides_Hex),                                   &
                 STAT = allocstat                                                                                           )
        IF(allocStat .NE. 0) THEN
            logError(*) 'could not allocate all variables!'
            STOP
        END IF

        DISC%Galerkin%MassMatrix_Hex(:,:,:)  = 0.0d0
        DISC%Galerkin%iMassMatrix_Hex(:,:,:) = 0.0d0
        DISC%Galerkin%cPoly3D_Hex(:,:,:,:,:) = 0.0d0

        ! Read coefficients of basis functions and mass matrices up to degree DISC%Galerkin%nMaxPoly, stored in the file
        ! BasisFunctions3D.qua
        DO iPoly = 0, DISC%Galerkin%nPolyRec
            ! Read comment in front of the basis functions' coefficients
            logInfo(*) 'Reading basis functions of order ', iPoly
            READ(IO%UNIT%FileIn,*)
            nDegFr = (iPoly + 1)*(iPoly + 2)*(iPoly + 3)/6
            ! Read polynomial coefficients
            ! where the index of the degrees of freedom starts at zero
            DO iDegFr = 0, nDegFr-1
                DO iZeta = 0, iPoly
                    DO iEta = 0, iPoly
                        DO iXi = 0, iPoly
                            READ(IO%UNIT%FileIn,*) DISC%Galerkin%cPoly3D_Hex(iXi,iEta,iZeta,iDegFr,iPoly)
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            ! Read comment in front of the entries of the mass matrix
            READ(IO%UNIT%FileIn,*)
            logInfo(*) 'Reading mass matrices   of order ', iPoly
            ! Read entries of the mass matrix
            DO k = 1, nDegFr
                DO l = 1, nDegFr
                    READ(IO%UNIT%FileIn,*) DISC%Galerkin%MassMatrix_Hex(k,l,iPoly)
                    IF(DISC%Galerkin%MassMatrix_Hex(k,l,iPoly).NE.0) THEN
                       DISC%Galerkin%iMassMatrix_Hex(k,l,iPoly) = 1.0d0/DISC%Galerkin%MassMatrix_Hex(k,l,iPoly)
                    ENDIF
                ENDDO
            ENDDO
        ENDDO
        CLOSE(IO%UNIT%FileIn)
        ! Detecting non-zero coefficients in basis polynomials
        DISC%Galerkin%NonZeroCPoly_Hex(:,:)          = 0
        DISC%Galerkin%NonZeroCPolyIndex_Hex(:,:,:,:) = -1
        DO iPoly = 0, DISC%Galerkin%nPolyRec
            nDegFr = (iPoly + 1)*(iPoly + 2)*(iPoly + 3)/6
            DO iDegFr = 0, nDegFr-1
                DO iZeta = 0, iPoly
                    DO iEta = 0, iPoly
                        DO iXi = 0, iPoly
                            IF(ABS(DISC%Galerkin%cPoly3D_Hex(iXi,iEta,iZeta,iDegFr,iPoly)).GE.1e-6) THEN
                                DISC%Galerkin%NonZeroCPoly_Hex(iDegFr,iPoly) = DISC%Galerkin%NonZeroCPoly_Hex(iDegFr,iPoly) + 1
                                DISC%Galerkin%NonZeroCPolyIndex_Hex(1,DISC%Galerkin%NonZeroCPoly_Hex(iDegFr,iPoly),iDegFr,iPoly) = iXi
                                DISC%Galerkin%NonZeroCPolyIndex_Hex(2,DISC%Galerkin%NonZeroCPoly_Hex(iDegFr,iPoly),iDegFr,iPoly) = iEta
                                DISC%Galerkin%NonZeroCPolyIndex_Hex(3,DISC%Galerkin%NonZeroCPoly_Hex(iDegFr,iPoly),iDegFr,iPoly) = iZeta
                            ENDIF
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
        ENDDO
    ENDIF ! Hexas

    logInfo(*) 'Basis functions and mass matrices read. '

  !==============================================================================!
  !------------------------Source Terms------------------------------------------!
  !==============================================================================!

   !
   ! Prepare Time Integration
   !
   ALLOCATE(DISC%Galerkin%cTimePoly(0:DISC%Galerkin%nPolyRec, 0:DISC%Galerkin%nPolyRec+1, 0:DISC%Galerkin%nPolyRec),      &
            DISC%Galerkin%TimeMassMatrix(1:DISC%Galerkin%nPolyRec+1,1:DISC%Galerkin%nPolyRec+1,0:DISC%Galerkin%nPolyRec), &
            DISC%Galerkin%iTimeMassMatrix(1:DISC%Galerkin%nPolyRec+1,1:DISC%Galerkin%nPolyRec+1,0:DISC%Galerkin%nPolyRec),&
            STAT = allocstat)
   IF(allocStat .NE. 0) THEN
           logError(*) 'could not allocate all variables!'
           STOP
   END IF
   !
   ! Allocate Gausspoints for ADER-DG Time Integration of source terms
   !
#ifndef NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
   DISC%Galerkin%nTimeGP = DISC%Galerkin%nPoly+1
#else
   disc%galerkin%nTimeGp = NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
#endif
   ALLOCATE(DISC%Galerkin%TimeGaussP(DISC%Galerkin%nTimeGP),                          &
            DISC%Galerkin%TimeGaussW(DISC%Galerkin%nTimeGP),                          &
            DISC%Galerkin%dtPowerFactor(0:DISC%Galerkin%nPoly,DISC%Galerkin%nTimeGP), &
            DISC%Galerkin%Faculty(0:DISC%Galerkin%nMaxPoly+5)  )
   !
   ! Precalculate the Faculty of i because otherwise this calculation takes a loooong time...
   !
   DISC%Galerkin%Faculty(0) = 1.
   DO i = 1, DISC%Galerkin%nMaxPoly+5
       DISC%Galerkin%Faculty(i) = DISC%Galerkin%Faculty(i-1)*REAL(i)
   ENDDO

   WRITE(FileName_Time,'(a,a13)') TRIM(DGPATH), 'TimeBasis.dat'
   CALL OpenFile(IO%UNIT%FileIn,FileName_Time,.FALSE.)
   !
   logInfo(*) 'Reading time basis functions and mass matrices for ADER-DG sources'
   logInfo(*) '  from file ', TRIM(FileName_Time)
   !
   ! Read comment
   READ(IO%UNIT%FileIn,*)
   READ(IO%UNIT%FileIn,*)
   ! Read maximal degree of basis polynomials stored in the file.
   READ(IO%UNIT%FileIn,*) i
   IF(DISC%Galerkin%nPoly.GT.i) THEN
      logError(*) 'Required polynomial for DG method is higher than the ones stored in file ', FileName_Time
      STOP
   ENDIF
   !
   ! Read coefficients of basis functions and mass matrices up to degree DISC%Galerkin%nMaxPoly
   !
   DISC%Galerkin%cTimePoly(:,:,:)       = 0.0d0
   DISC%Galerkin%TimeMassMatrix(:,:,:)  = 0.0d0
   DISC%Galerkin%iTimeMassMatrix(:,:,:) = 0.0d0
   !
   DO iPoly = 0, DISC%Galerkin%nPolyRec
      ! Read comment in front of the basis functions' coefficients
      READ(IO%UNIT%FileIn,*)
      nDegFr = (iPoly + 1)
      ! Read polynomial coefficients
      ! where the index of the degrees of freedom starts at zero
      DO iDegFr = 0, nDegFr-1
         DO iXi = 0, iPoly
            READ(IO%UNIT%FileIn,*) DISC%Galerkin%cTimePoly(iXi,iDegFr,iPoly)
         ENDDO
      ENDDO
      ! Read comment in front of the entries of the time mass matrix
      READ(IO%UNIT%FileIn,*)
      ! Read entries of the time mass matrix
      DO k = 1, nDegFr
         DO l = 1, nDegFr
            READ(IO%UNIT%FileIn,*) DISC%Galerkin%TimeMassMatrix(k,l,iPoly)
         ENDDO
      ENDDO
      ! Read comment in front of the entries of the inverse time mass matrix
      READ(IO%UNIT%FileIn,*)
      ! Read entries of the inverse time mass matrix
      DO k = 1, nDegFr
         DO l = 1, nDegFr
            READ(IO%UNIT%FileIn,*) DISC%Galerkin%iTimeMassMatrix(k,l,iPoly)
         ENDDO
      ENDDO
   ENDDO
   CLOSE(IO%UNIT%FileIn)
   logInfo(*) 'Time basis functions and mass matrices read. '

!===================================================================================!
!---------------------------Quadrature-free ADER DG---------------------------------!
!===================================================================================!

    DISC%Galerkin%nRK = 1

    ! Attention: Don't change Nr of GP here since some routine depend on these numbers
    DISC%Galerkin%nIntGP = (DISC%Galerkin%nPoly + 2)**3
#ifndef GENERATEDKERNELS
    DISC%Galerkin%nBndGP = (DISC%Galerkin%nPoly + 2)**2
#endif

    SELECT CASE(DISC%Galerkin%DGMethod)
    CASE(3)
        nDGWorkVar = EQN%nVar+EQN%nAneFuncperMech
    CASE DEFAULT
        nDGWorkVar = EQN%nVarTotal
    END SELECT

    IF(MESH%nElem_Tet.GT.0) THEN

         ! Allocation of arrays
         ALLOCATE(                                                                                        &
         DISC%Galerkin%Kxi_k_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),   &
         DISC%Galerkin%Keta_k_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),  &
         DISC%Galerkin%Kzeta_k_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat), &
         DISC%Galerkin%Kxi_m_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),   &
         DISC%Galerkin%Keta_m_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),  &
         DISC%Galerkin%Kzeta_m_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat), &
         DISC%Galerkin%ADGxi_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),   &
         DISC%Galerkin%ADGeta_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),  &
         DISC%Galerkin%ADGzeta_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat), &
         DISC%Galerkin%ADGklm_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),  &
         DISC%Galerkin%FluxInt_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat,  &
                                     0:MESH%nSides_Tet,1:MESH%nVertices_Tri,1:MESH%nSides_Tet),           &
         STAT=allocstat                                                                                   )
         IF(allocStat .NE. 0) THEN
            logError(*) 'could not allocate all variables!'
            STOP
         END IF

        ! Compute and store volume gaussian integration points
         CALL TetrahedronQuadraturePoints(                     &
                 nIntGP     = DISC%Galerkin%nIntGP,            &
                 IntGaussP  = DISC%Galerkin%IntGaussP_Tet,     &
                 IntGaussW  = DISC%Galerkin%IntGaussW_Tet,     &
                 M          = DISC%Galerkin%nPoly+2,           &
                 IO         = IO,                              &
                 quiet      = .TRUE.                           )

#ifdef USE_DR_CELLAVERAGE
        call CellCentresOfSubdivision(DISC%Galerkin%nPoly + 1, DISC%Galerkin%BndGaussP_Tet)
        DISC%Galerkin%BndGaussW_Tet = 1.e99 ! blow up solution if used
#else
        ! Compute and store surface gaussian integration points
        CALL TriangleQuadraturePoints(                         &
                 nIntGP     = DISC%Galerkin%nBndGP,            &
                 IntGaussP  = DISC%Galerkin%BndGaussP_Tet,     &
                 IntGaussW  = DISC%Galerkin%BndGaussW_Tet,     &
                 M          = DISC%Galerkin%nPoly+2,           &
                 IO         = IO,                              &
                 quiet      = .TRUE.                           )
#endif

        NULLIFY( Tens3GaussP )
        NULLIFY( Tens3GaussW )
        ! Compute integration points for volume integrals
        ! Local to this subroutine
         CALL TetrahedronQuadraturePoints(                                  &
                 nIntGP     = nTens3GP,                                     &
                 IntGaussP  = Tens3GaussP,                                  &
                 IntGaussW  = Tens3GaussW,                                  &
                 M          = DISC%Galerkin%nPoly+DISC%Galerkin%nPolyMat+2, &
                 IO         = IO,                                           &
                 quiet      = .TRUE.                                        )

        NULLIFY( TFMGaussP )
        NULLIFY( TFMGaussW )
        ! Compute integration points for surface integrals
        ! Local to this subroutine
        CALL TriangleQuadraturePoints(                                      &
                 nIntGP     = nTFMGaussP,                                   &
                 IntGaussP  = TFMGaussP,                                    &
                 IntGaussW  = TFMGaussW,                                    &
                 M          = DISC%Galerkin%nPoly+DISC%Galerkin%nPolyMat+2, &
                 IO         = IO,                                           &
                 quiet      = .TRUE.                                        )

        ! Compute the value and gradient of basis functions
        ! at volume integration points Tens3GaussP
        ALLOCATE( Tens3BaseFunc(  DISC%Galerkin%nDegFr,nTens3GP) )
        ALLOCATE( Tens3BaseGrad(3,DISC%Galerkin%nDegFr,nTens3GP) )

        DO iIntGP = 1, nTens3GP
            xi   = Tens3GaussP(1,iIntGP)
            eta  = Tens3GaussP(2,iIntGP)
            zeta = Tens3GaussP(3,iIntGP)
            DO l = 1, DISC%Galerkin%nDegFr
                CALL BaseFunc3D(phi_l,l,xi,eta,zeta,DISC%Galerkin%nPoly,   &
                               DISC%Galerkin%cPoly3D_Tet,                  &
                               DISC%Galerkin%NonZeroCPoly_Tet,             &
                               DISC%Galerkin%NonZeroCPolyIndex_Tet         )
                CALL BaseGrad3D(phigrad,l,xi,eta,zeta,DISC%Galerkin%nPoly, &
                               DISC%Galerkin%cPoly3D_Tet,                  &
                               DISC%Galerkin%NonZeroCPoly_Tet,             &
                               DISC%Galerkin%NonZeroCPolyIndex_Tet         )
                Tens3BaseFunc(l,iIntGP)   = phi_l
                Tens3BaseGrad(:,l,iIntGP) = phigrad
            ENDDO
        ENDDO

        ! Compute stiffness tensors for the Tetrahedrons
        DISC%Galerkin%Kxi_k_Tet   = 0.0d0
        DISC%Galerkin%Keta_k_Tet  = 0.0d0
        DISC%Galerkin%Kzeta_k_Tet = 0.0d0
        DO iIntGP = 1, nTens3GP
            xi   = Tens3GaussP(1,iIntGP)
            eta  = Tens3GaussP(2,iIntGP)
            zeta = Tens3GaussP(3,iIntGP)
            DO m = 1, DISC%Galerkin%nDegFrMat
                phi_m = Tens3BaseFunc(m,iIntGP)
                DO l = 1, DISC%Galerkin%nDegFr
                    phi_l = Tens3BaseFunc(l,iIntGP)
                    DO k = 1, DISC%Galerkin%nDegFr
                        phigrad = Tens3BaseGrad(:,k,iIntGP)
                        DISC%Galerkin%Kxi_k_Tet(   k,l,m ) = DISC%Galerkin%Kxi_k_Tet(   k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phigrad(1)*phi_l*phi_m
                        DISC%Galerkin%Keta_k_Tet(  k,l,m ) = DISC%Galerkin%Keta_k_Tet(  k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phigrad(2)*phi_l*phi_m
                        DISC%Galerkin%Kzeta_k_Tet( k,l,m ) = DISC%Galerkin%Kzeta_k_Tet( k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phigrad(3)*phi_l*phi_m
                    ENDDO
                ENDDO
            ENDDO
            !
        ENDDO
        !
        DISC%Galerkin%Kxi_m_Tet   = 0.0d0
        DISC%Galerkin%Keta_m_Tet  = 0.0d0
        DISC%Galerkin%Kzeta_m_Tet = 0.0d0
        DO iIntGP = 1, nTens3GP
            xi   = Tens3GaussP(1,iIntGP)
            eta  = Tens3GaussP(2,iIntGP)
            zeta = Tens3GaussP(3,iIntGP)
            DO m = 1, DISC%Galerkin%nDegFrMat
                phigrad = Tens3BaseGrad(:,m,iIntGP)
                DO l = 1, DISC%Galerkin%nDegFr
                    phi_l = Tens3BaseFunc(l,iIntGP)
                    DO k = 1, DISC%Galerkin%nDegFr
                        phi_k = Tens3BaseFunc(k,iIntGP)
                        DISC%Galerkin%Kxi_m_Tet(   k,l,m ) = DISC%Galerkin%Kxi_m_Tet(   k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phi_k*phi_l*phigrad(1)
                        DISC%Galerkin%Keta_m_Tet(  k,l,m ) = DISC%Galerkin%Keta_m_Tet(  k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phi_k*phi_l*phigrad(2)
                        DISC%Galerkin%Kzeta_m_Tet( k,l,m ) = DISC%Galerkin%Kzeta_m_Tet( k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phi_k*phi_l*phigrad(3)
                    ENDDO
                ENDDO
            ENDDO
            !
        ENDDO
        !
        WHERE(ABS(DISC%Galerkin%Kxi_k_Tet).LE.1e-10)
            DISC%Galerkin%Kxi_k_Tet = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Keta_k_Tet).LE.1e-10)
            DISC%Galerkin%Keta_k_Tet = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Kzeta_k_Tet).LE.1e-10)
            DISC%Galerkin%Kzeta_k_Tet = 0.0d0
        ENDWHERE
        !
        WHERE(ABS(DISC%Galerkin%Kxi_m_Tet).LE.1e-10)
            DISC%Galerkin%Kxi_m_Tet = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Keta_m_Tet).LE.1e-10)
            DISC%Galerkin%Keta_m_Tet = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Kzeta_m_Tet).LE.1e-10)
            DISC%Galerkin%Kzeta_m_Tet = 0.0d0
        ENDWHERE
        !
        DO m = 1, DISC%Galerkin%nDegFrMat
        DO l = 1, DISC%Galerkin%nDegFr
        DO k = 1, DISC%Galerkin%nDegFr
          DISC%Galerkin%ADGxi_Tet(k,l,m)   = DISC%Galerkin%Kxi_k_Tet(l,k,m)   / DISC%Galerkin%MassMatrix_Tet(k,k,DISC%Galerkin%nPoly)
          DISC%Galerkin%ADGeta_Tet(k,l,m)  = DISC%Galerkin%Keta_k_Tet(l,k,m)  / DISC%Galerkin%MassMatrix_Tet(k,k,DISC%Galerkin%nPoly)
          DISC%Galerkin%ADGzeta_Tet(k,l,m) = DISC%Galerkin%Kzeta_k_Tet(l,k,m) / DISC%Galerkin%MassMatrix_Tet(k,k,DISC%Galerkin%nPoly)
        ENDDO
        ENDDO
        ENDDO
        !
        ! ADER-DG tensor for space dependent reaction term (= identical to identity matrix -> to be optimized!)
        DISC%Galerkin%ADGklm_Tet(:,:,:) = 0.0D0
        DO m = 1, DISC%Galerkin%nDegFrMat
        DO k = 1, DISC%Galerkin%nDegFr
          DISC%Galerkin%ADGklm_Tet(k,k,m) = 1.0D0
        ENDDO
        ENDDO
        !
        ALLOCATE(DISC%Galerkin%Kxi_k_Tet_Sp, DISC%Galerkin%Keta_k_Tet_Sp, DISC%Galerkin%Kzeta_k_Tet_Sp)
        ALLOCATE(DISC%Galerkin%Kxi_m_Tet_Sp, DISC%Galerkin%Keta_m_Tet_Sp, DISC%Galerkin%Kzeta_m_Tet_Sp)
        CALL IniSparseTensor3b(DISC%Galerkin%Kxi_k_Tet_Sp,   DISC%Galerkin%Kxi_k_Tet,   DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Keta_k_Tet_Sp,  DISC%Galerkin%Keta_k_Tet,  DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Kzeta_k_Tet_Sp, DISC%Galerkin%Kzeta_k_Tet, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Kxi_m_Tet_Sp,   DISC%Galerkin%Kxi_m_Tet,   DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Keta_m_Tet_Sp,  DISC%Galerkin%Keta_m_Tet,  DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Kzeta_m_Tet_Sp, DISC%Galerkin%Kzeta_m_Tet, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        !
        ALLOCATE(DISC%Galerkin%ADGxi_Tet_Sp, DISC%Galerkin%ADGeta_Tet_Sp, DISC%Galerkin%ADGzeta_Tet_Sp, DISC%Galerkin%ADGklm_Tet_Sp)
        CALL IniSparseTensor3b(DISC%Galerkin%ADGxi_Tet_Sp,   DISC%Galerkin%ADGxi_Tet,   DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%ADGeta_Tet_Sp,  DISC%Galerkin%ADGeta_Tet,  DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%ADGzeta_Tet_Sp, DISC%Galerkin%ADGzeta_Tet, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%ADGklm_Tet_Sp,  DISC%Galerkin%ADGklm_Tet,  DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)

        !
        DISC%Galerkin%FluxInt_Tet = 0.0d0
        !
        DO iSide = 1, MESH%nSides_Tet
            !
            ! Compute for the element itself
            !
            DO iBndGP = 1, nTFMGaussP
                chiGP  = TFMGaussP(1,iBndGP)
                tauGP  = TFMGaussP(2,iBndGP)
                CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chiGP,tauGP,iSide,0)
                DO m = 1, DISC%Galerkin%nDegFrMat
                    CALL BaseFunc3D(phi_m,m,xi,eta,zeta,DISC%Galerkin%nPolyMat,         &
                                    DISC%Galerkin%cPoly3D_Tet,                          &
                                    DISC%Galerkin%NonZeroCPoly_Tet,                     &
                                    DISC%Galerkin%NonZeroCPolyIndex_Tet                 )
                    DO l = 1, DISC%Galerkin%nDegFr
                        CALL BaseFunc3D(phi_l,l,xi,eta,zeta,DISC%Galerkin%nPoly,        &
                                        DISC%Galerkin%cPoly3D_Tet,                      &
                                        DISC%Galerkin%NonZeroCPoly_Tet,                 &
                                        DISC%Galerkin%NonZeroCPolyIndex_Tet             )
                        DO k = 1, DISC%Galerkin%nDegFr
                            CALL BaseFunc3D(phi_k,k,xi,eta,zeta,DISC%Galerkin%nPoly,    &
                                            DISC%Galerkin%cPoly3D_Tet,                  &
                                            DISC%Galerkin%NonZeroCPoly_Tet,             &
                                            DISC%Galerkin%NonZeroCPolyIndex_Tet         )
                            DISC%Galerkin%FluxInt_Tet(k,l,m,0,1,iSide) =  &
                            DISC%Galerkin%FluxInt_Tet(k,l,m,0,1,iSide) + TFMGaussW(iBndGP)*phi_k*phi_l*phi_m
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            !
            ! Compute for the neighbors
            !
            DO iLocalNeighborSide = 1, MESH%nSides_Tet
            DO iLocalNeighborVrtx = 1, MESH%nVertices_Tri
            DO iBndGP = 1, nTFMGaussP
               chiGP  = TFMGaussP(1,iBndGP)
               tauGP  = TFMGaussP(2,iBndGP)
               CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chiGP,tauGP,iLocalNeighborSide,iLocalNeighborVrtx)
               CALL TrafoChiTau2XiEtaZeta(xiS,etaS,zetaS,chiGP,tauGP,iSide,0)
               DO m = 1, DISC%Galerkin%nDegFrMat
                 CALL BaseFunc3D(phi_m,m,xiS,etaS,zetaS,DISC%Galerkin%nPolyMat,         &
                                 DISC%Galerkin%cPoly3D_Tet,                             &
                                 DISC%Galerkin%NonZeroCPoly_Tet,                        &
                                 DISC%Galerkin%NonZeroCPolyIndex_Tet                    )
                 DO l = 1, DISC%Galerkin%nDegFr
                   CALL BaseFunc3D(phi_l,l,xi,eta,zeta,DISC%Galerkin%nPoly,             &
                                   DISC%Galerkin%cPoly3D_Tet,                           &
                                   DISC%Galerkin%NonZeroCPoly_Tet,                      &
                                   DISC%Galerkin%NonZeroCPolyIndex_Tet                  )
                   DO k = 1, DISC%Galerkin%nDegFr
                      CALL BaseFunc3D(phi_k,k,xiS,etaS,zetaS,DISC%Galerkin%nPoly,       &
                                      DISC%Galerkin%cPoly3D_Tet,                        &
                                      DISC%Galerkin%NonZeroCPoly_Tet,                   &
                                      DISC%Galerkin%NonZeroCPolyIndex_Tet               )
                      DISC%Galerkin%FluxInt_Tet(k,l,m,iLocalNeighborSide,iLocalNeighborVrtx,iSide) =  &
                      DISC%Galerkin%FluxInt_Tet(k,l,m,iLocalNeighborSide,iLocalNeighborVrtx,iSide) +  &
                      TFMGaussW(iBndGP)*phi_k*phi_l*phi_m
                   ENDDO
                 ENDDO
               ENDDO
               !
            ENDDO
            ENDDO
            ENDDO
        !
        ENDDO
        WHERE(ABS(DISC%Galerkin%FluxInt_Tet).LT.1e-10)
            DISC%Galerkin%FluxInt_Tet = 0.0d0
        ENDWHERE
        !
        ALLOCATE(DISC%Galerkin%FluxInt_Tet_Sp(0:MESH%nSides_Tet,1:MESH%nVertices_Tri,1:MESH%nSides_Tet))
        DO iSide = 1, MESH%nSides_Tet
            CALL IniSparseTensor3b(DISC%Galerkin%FluxInt_Tet_Sp(0,1,iSide), DISC%Galerkin%FluxInt_Tet(:,:,:,0,1,iSide), DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
            DO iLocalNeighborSide = 1, MESH%nSides_Tet
            DO iLocalNeighborVrtx = 1, MESH%nVertices_Tri
             CALL IniSparseTensor3b(DISC%Galerkin%FluxInt_Tet_Sp(iLocalNeighborSide,iLocalNeighborVrtx,iSide),    &
                                    DISC%Galerkin%FluxInt_Tet(:,:,:,iLocalNeighborSide,iLocalNeighborVrtx,iSide), &
                                    DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat       )
            ENDDO
            ENDDO
        ENDDO

#ifdef GENERATEDKERNELS
#ifndef NDEBUG
        ! assert contant material parameters per element
        if ( disc%galerkin%nDegFrMat .ne. 1 ) then
          logError(*) 'iniGalerkin3D_us_intern_new, disc%galerkin%nDegFrMat not equal 1.', disc%galerkin%nDegFrMat
          stop
        endif

        ! assert 4 sides for tetrahedrons
        if ( mesh%nSides_tet .ne. 4 ) then
          logError(*) 'iniGalerkin3D_us_intern_new, mesh%nSides_tet not equal 4.', mesh%nSides_tet
          stop
        endif

        ! assert 3 vertices for triangles
         if ( mesh%nVertices_tri .ne. 3 ) then
          logError(*) 'iniGalerkin3D_us_intern_new, mesh%nVertices_tri not equal 3.', mesh%nVertices_tri
          stop
        endif
#endif
#endif

        DEALLOCATE( Tens3BaseFunc )
        DEALLOCATE( Tens3BaseGrad )

        DO iDegFr = 1, DISC%Galerkin%nDegFr
            DISC%Galerkin%iMassMatrix_Tet(iDegFr,iDegFr,DISC%Galerkin%nPoly) = &
                1./DISC%Galerkin%MassMatrix_Tet(iDegFr,iDegFr,DISC%Galerkin%nPoly)
        ENDDO

        DISC%Galerkin%IntGPBaseFunc_Tet(:,:,:)          = 0.0d0
        DISC%Galerkin%IntGPBaseGrad_Tet(:,:,:,:)        = 0.0d0
!        DISC%Galerkin%BndGPBaseFunc3D_Tet(:,:,:)        = 0.0d0

        iPoly = DISC%Galerkin%nPoly

        DO iIntGP = 1, DISC%Galerkin%nIntGP
            DO iDegFr = 1, DISC%Galerkin%nDegFr
                !
                CALL BaseFunc3D(phi,iDegFr,                             &
                              DISC%Galerkin%intGaussP_Tet(1,iIntGP),    &
                              DISC%Galerkin%intGaussP_Tet(2,iIntGP),    &
                              DISC%Galerkin%intGaussP_Tet(3,iIntGP),    &
                              DISC%Galerkin%nPoly,                      &
                              DISC%Galerkin%cPoly3D_Tet,                &
                              DISC%Galerkin%NonZeroCPoly_Tet,           &
                              DISC%Galerkin%NonZeroCPolyIndex_Tet       )
                !
                DISC%Galerkin%IntGPBaseFunc_Tet(iDegFr,iIntGP,iPoly)   = phi
                !
                CALL BaseGrad3D( gradphixieta,iDegFr,                   &
                              DISC%Galerkin%intGaussP_Tet(1,iIntGP),    &
                              DISC%Galerkin%intGaussP_Tet(2,iIntGP),    &
                              DISC%Galerkin%intGaussP_Tet(3,iIntGP),    &
                              DISC%Galerkin%nPoly,                      &
                              DISC%Galerkin%cPoly3D_Tet,                &
                              DISC%Galerkin%NonZeroCPoly_Tet,           &
                              DISC%Galerkin%NonZeroCPolyIndex_Tet       )
                !
                DISC%Galerkin%IntGPBaseGrad_Tet(:,iDegFr,iIntGP,iPoly) = gradphixieta
                !
            ENDDO
        ENDDO
        !
#ifndef GENERATEDKERNELS
        MESH%ELEM%BndBF_GP_Tet(:,:,:) = 0.0

        DO iSide = 1, MESH%nSides_Tet
            !
            ! Compute for the element itself
            !
            DO iBndGP = 1, DISC%Galerkin%nBndGP
               DO iDegFr = 1, DISC%Galerkin%nDegFr
                  !
                  chi  = DISC%Galerkin%bndGaussP_Tet(1,iBndGP)
                  tau  = DISC%Galerkin%bndGaussP_Tet(2,iBndGP)
                  !
                  CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
                  !           !
                  CALL BaseFunc3D( phi,iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly, &
                                   DISC%Galerkin%cPoly3D_Tet,                  &
                                   DISC%Galerkin%NonZeroCPoly_Tet,             &
                                   DISC%Galerkin%NonZeroCPolyIndex_Tet         )
                  !
!                   DISC%Galerkin%BndGPBaseFunc3D_Tet(iDegFr,iBndGP,iSide)  = phi
                  MESH%ELEM%BndBF_GP_Tet(iDegFr,iBndGP,iSide) = phi
                  !
               ENDDO
            ENDDO
        ENDDO


        ! ------------------------------------------------------------
        !  Dynamic Rupture setup - force GP matching at fault surface
        ! ------------------------------------------------------------
        !
        ! The GP at '+' and '-' side at a fault surface GPs do not match necessarily.
        ! Thus, we have to force them to in order to evaluate the friction law correct
        !
        IF (EQN%DR.EQ.1) THEN
        !
        ALLOCATE(DISC%DynRup%BndBF_GP_Tet(DISC%Galerkin%nDegFr,DISC%Galerkin%nBndGP,MESH%Fault%nSide))
        ALLOCATE(DISC%DynRup%FluxInt(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFr,MESH%Fault%nSide))
        DISC%DynRup%BndBF_GP_Tet = 0.0D0
        DISC%DynRup%FluxInt = 0.0D0
        !
#if __INTEL_COMPILER == 1600
!#if __INTEL_COMPILER_UPDATE > 1
!#ifdef OMP
!        !$omp parallel private(iFace, iElem, iSide, iNeighbor, iLocalNeighborSide, iObject, MPIIndex, xV, yV, zV, x_host, y_host, z_host, iBndGP, chi, tau, xi, eta, zeta, xGP, yGP, zGP, iDegFr, iDegFr2, phi_iDegFr, phi_iDegFr2) shared( mesh, disc, eqn, bnd) default( none )
!        !$omp do schedule(static)
!#else
        logWarning(*) 'Intel Compiler 16.x were used to compile SeisSol, DR setup is not parallized'
!#endif
!#endif
#else
#ifdef OMP
        !$omp parallel private(iFace, iElem, iSide, iNeighbor, iLocalNeighborSide, iObject, MPIIndex, xV, yV, zV, x_host, y_host, z_host, iBndGP, chi, tau, xi, eta, zeta, xGP, yGP, zGP, iDegFr, iDegFr2, phi_iDegFr, phi_iDegFr2) shared( mesh, disc, eqn, bnd) default( none )
        !$omp do schedule(static)
#endif
#endif
        DO iFace=1,MESH%Fault%nSide
            !
            iElem               = MESH%Fault%Face(iFace,1,1)          ! Remark:
            iSide               = MESH%Fault%Face(iFace,2,1)          ! iElem denotes "+" side
            !
            iNeighbor           = MESH%Fault%Face(iFace,1,2)          ! iNeighbor denotes "-" side
            iLocalNeighborSide  = MESH%Fault%Face(iFace,2,2)
            !
            IF (iElem == 0) THEN
                ! iElem is in the neighbor domain
                ! obtain parameters from "-" element
                ! The neighbor element belongs to a different MPI domain
                iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
                MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
                !
                xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
                yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
                zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
                !
                ! coordinates of "-" side which is present in local domain
                x_host(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iNeighbor))
                y_host(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iNeighbor))
                z_host(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iNeighbor))
            ELSE
                !
                ! get vertices
                xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
                yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
                zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
                !
                ! get the location in the reference element of Neighbor
                IF (iNeighbor.NE.0) THEN
                    x_host(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iNeighbor))
                    y_host(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iNeighbor))
                    z_host(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iNeighbor))
                ELSE
                    ! The neighbor element belongs to a different MPI domain
                    iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
                    MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
                    !
                    x_host(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
                    y_host(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
                    z_host(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
                ENDIF
            ENDIF
            !
            DO iBndGP = 1, DISC%Galerkin%nBndGP
                chi  = DISC%Galerkin%bndGaussP_Tet(1,iBndGP)
                tau  = DISC%Galerkin%bndGaussP_Tet(2,iBndGP)
                !
                ! get xi, eta, zeta of neighbor element
                CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
                CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
                CALL TrafoXYZ2XiEtaZeta(xi,eta,zeta,xGP,yGP,zGP,x_host(1:4),y_host(1:4),z_host(1:4),4)
                !
                ! determine new basis functions at mirrored GPs
                ! variable DISC%DynRup%BndBF_GP_Tet stores basis functions at boundary for '-' element
                DO iDegFr = 1, DISC%Galerkin%nDegFr
                    CALL BaseFunc3D(                                                       &
                      value             = phi_iDegFr,   &
                      iDegFr            = iDegFr,                                          &
                      xi                = xi,                                              &
                      eta               = eta,                                             &
                      zeta              = zeta,                                            &
                      nPoly             = DISC%Galerkin%nPoly,                             &
                      cPoly3D           = DISC%Galerkin%CPoly3D_Tet,                       &
                      NonZeroCPoly      = DISC%Galerkin%NonZeroCPoly_Tet,                  &
                      NonZeroCPolyIndex = DISC%Galerkin%NonZeroCPolyIndex_Tet              )
                    DISC%DynRup%BndBF_GP_Tet(iDegFr,iBndGP,iFace) = phi_iDegFr
                    ! Compute corresponding new Flux integration matrices for '-' element
                    DO iDegFr2 = 1, DISC%Galerkin%nDegFr
                        CALL BaseFunc3D(                                                       &
                          value             = phi_iDegFr2,                                     &
                          iDegFr            = iDegFr2,                                         &
                          xi                = xi,                                              &
                          eta               = eta,                                             &
                          zeta              = zeta,                                            &
                          nPoly             = DISC%Galerkin%nPoly,                             &
                          cPoly3D           = DISC%Galerkin%CPoly3D_Tet,                       &
                          NonZeroCPoly      = DISC%Galerkin%NonZeroCPoly_Tet,                  &
                          NonZeroCPolyIndex = DISC%Galerkin%NonZeroCPolyIndex_Tet              )
                        DISC%DynRup%FluxInt(iDegFr2,iDegFr,iFace) = DISC%DynRup%FluxInt(iDegFr2,iDegFr,iFace) &
                                                                  + MESH%ELEM%BndGW_Tri(iBndGP)*phi_iDegFr*phi_iDegFr2
                    ENDDO
                ENDDO
            ENDDO
            !
        ENDDO ! iFace
#if __INTEL_COMPILER == 1600
!#if __INTEL_COMPILER_UPDATE > 1
!#ifdef OMP
!        !$omp end parallel
!#endif
!#endif
#else
#ifdef OMP
        !$omp end parallel
#endif
#endif
        !
        ! accordingly to normal FluxInt matrix all values small than 1e-10 are set to zero
        WHERE(ABS(DISC%DynRup%FluxInt).LT.1e-10)
            DISC%DynRup%FluxInt = 0.0d0
        ENDWHERE
        !
        ENDIF ! EQN%DR.EQ.1 - force GP matching at fault surface
#endif ! GENERATEDKERNELS

!        ! ------------------------------------------------------------
!        ! ATTENUATION
!        ! ------------------------------------------------------------
!        ! Preliminary generation of sparse variables NonZero_KX, CNonZero_KX and IndexNonZero_KX (with X = xi,eta,zeta)
!        IF(EQN%Anelasticity.EQ.1) THEN
!
!            ALLOCATE(                                                                          &
!            DISC%Galerkin%FMatrix3D(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,0:4,1:3,1:4), &
!            DISC%Galerkin%Kxi(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec),                   &
!            DISC%Galerkin%Keta(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec),                  &
!            DISC%Galerkin%Kzeta(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec)                  )
!
!            CALL IniQuadFreeRKDG3D(                                          &
!                                   FMatrix3D = DISC%Galerkin%FMatrix3D,      &
!                                   Kxi       = DISC%Galerkin%Kxi,            &
!                                   Keta      = DISC%Galerkin%Keta,           &
!                                   Kzeta     = DISC%Galerkin%Kzeta,          &
!                                   nDegFr    = DISC%Galerkin%nDegFr,         &
!                                   nPoly     = DISC%Galerkin%nPoly,          &
!                                   nDegFrRec = DISC%Galerkin%nDegFrRec,      &
!                                   nPolyRec  = DISC%Galerkin%nPolyRec,       &
!                                   IO        = IO                            )
!
!
!            DISC%Galerkin%NonZero_Kxi           = 0
!            DISC%Galerkin%NonZero_Keta          = 0
!            DISC%Galerkin%NonZero_Kzeta         = 0
!
!            DO l = 1, DISC%Galerkin%nDegFr
!              DO m = 1, DISC%Galerkin%nDegFrRec
!                IF(ABS(DISC%Galerkin%Kxi(l,m)).GT.1e-6) THEN
!                  DISC%Galerkin%NonZero_Kxi   = DISC%Galerkin%NonZero_Kxi + 1
!                ENDIF
!                IF(ABS(DISC%Galerkin%Keta(l,m)).GT.1e-6) THEN
!                  DISC%Galerkin%NonZero_Keta  = DISC%Galerkin%NonZero_Keta + 1
!                ENDIF
!                IF(ABS(DISC%Galerkin%Kzeta(l,m)).GT.1e-6) THEN
!                  DISC%Galerkin%NonZero_Kzeta = DISC%Galerkin%NonZero_Kzeta + 1
!                ENDIF
!              ENDDO
!            ENDDO
!
!            ALLOCATE( DISC%Galerkin%CNonZero_Kxi(  DISC%Galerkin%NonZero_Kxi  ) )
!            ALLOCATE( DISC%Galerkin%CNonZero_Keta( DISC%Galerkin%NonZero_Keta ) )
!            ALLOCATE( DISC%Galerkin%CNonZero_Kzeta(DISC%Galerkin%NonZero_Kzeta) )
!
!            ALLOCATE( DISC%Galerkin%IndexNonZero_Kxi(  2,DISC%Galerkin%NonZero_Kxi    ) )
!            ALLOCATE( DISC%Galerkin%IndexNonZero_Keta( 2,DISC%Galerkin%NonZero_Keta   ) )
!            ALLOCATE( DISC%Galerkin%IndexNonZero_Kzeta(2,DISC%Galerkin%NonZero_Kzeta  ) )
!
!            DISC%Galerkin%NonZero_Kxi           = 0
!            DISC%Galerkin%NonZero_Keta          = 0
!            DISC%Galerkin%NonZero_Kzeta         = 0
!            !
!            DISC%Galerkin%CNonZero_Kxi          = 0
!            DISC%Galerkin%CNonZero_Keta         = 0
!            DISC%Galerkin%CNonZero_Kzeta        = 0
!            !
!            DISC%Galerkin%IndexNonZero_Kxi      = 0
!            DISC%Galerkin%IndexNonZero_Keta     = 0
!            DISC%Galerkin%IndexNonZero_Kzeta    = 0
!
!            DO l = 1, DISC%Galerkin%nDegFr
!              DO m = 1, DISC%Galerkin%nDegFrRec
!                IF(ABS(DISC%Galerkin%Kxi(l,m)).GT.1e-6) THEN
!                  DISC%Galerkin%NonZero_Kxi   = DISC%Galerkin%NonZero_Kxi + 1
!                  DISC%Galerkin%CNonZero_Kxi(     DISC%Galerkin%NonZero_Kxi )        = DISC%Galerkin%Kxi(l,m)
!                  DISC%Galerkin%IndexNonZero_Kxi( :, DISC%Galerkin%NonZero_Kxi )     = (/ l , m /)
!                ENDIF
!                IF(ABS(DISC%Galerkin%Keta(l,m)).GT.1e-6) THEN
!                  DISC%Galerkin%NonZero_Keta  = DISC%Galerkin%NonZero_Keta + 1
!                  DISC%Galerkin%CNonZero_Keta(     DISC%Galerkin%NonZero_Keta )      = DISC%Galerkin%Keta(l,m)
!                  DISC%Galerkin%IndexNonZero_Keta( :, DISC%Galerkin%NonZero_Keta )   = (/ l , m /)
!                ENDIF
!                IF(ABS(DISC%Galerkin%Kzeta(l,m)).GT.1e-6) THEN
!                  DISC%Galerkin%NonZero_Kzeta = DISC%Galerkin%NonZero_Kzeta + 1
!                  DISC%Galerkin%CNonZero_Kzeta(     DISC%Galerkin%NonZero_Kzeta )    = DISC%Galerkin%Kzeta(l,m)
!                  DISC%Galerkin%IndexNonZero_Kzeta( :, DISC%Galerkin%NonZero_Kzeta ) = (/ l , m /)
!                ENDIF
!              ENDDO
!            ENDDO
!
!        ENDIF ! IF(EQN%Anelasticity.EQ.1)


        NULLIFY(DISC%Galerkin%NonZeroCPoly,DISC%Galerkin%NonZeroCPolyIndex,DISC%Galerkin%CPoly3D)

        CONTINUE

    ENDIF ! Tetras

    IF(MESH%nElem_Hex.GT.0) THEN

         ! Allocation of arrays
         ALLOCATE(                                                                                        &
         DISC%Galerkin%Kxi_k_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),   &
         DISC%Galerkin%Keta_k_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),  &
         DISC%Galerkin%Kzeta_k_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat), &
         DISC%Galerkin%Kxi_m_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),   &
         DISC%Galerkin%Keta_m_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),  &
         DISC%Galerkin%Kzeta_m_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat), &
         DISC%Galerkin%ADGxi_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),   &
         DISC%Galerkin%ADGeta_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat),  &
         DISC%Galerkin%ADGzeta_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat), &
         DISC%Galerkin%FluxInt_Hex(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFrRec,DISC%Galerkin%nDegFrMat,  &
                                     0:MESH%nSides_Hex,1:MESH%nVertices_Quad,1:MESH%nSides_Hex),          &
         STAT=allocstat                                                                                   )
         IF(allocStat .NE. 0) THEN
            logError(*) 'could not allocate all variables!'
            STOP
         END IF

        ! Compute and store volume gaussian integration points
        CALL HypercubeQuadraturePoints(                        &
                 nIntGP     = DISC%Galerkin%nIntGP,            &
                 IntGaussP  = DISC%Galerkin%IntGaussP_Hex,     &
                 IntGaussW  = DISC%Galerkin%IntGaussW_Hex,     &
                 M          = DISC%Galerkin%nPoly+2,           &
                 nDim       = 3,                               &
                 X1         = (/ 0., 0., 0. /),                &
                 X2         = (/ 1., 1., 1. /),                &
                 IO         = IO,                              &
                 quiet      = .TRUE.                           )

        ! Compute and store surface gaussian integration points
!~         CALL HypercubeQuadraturePoints(                        &
!~                  nIntGP     = DISC%Galerkin%nBndGP,            &
!~                  IntGaussP  = DISC%Galerkin%BndGaussP_Hex,     &
!~                  IntGaussW  = DISC%Galerkin%BndGaussW_Hex,     &
!~                  M          = DISC%Galerkin%nPoly+2,           &
!~                  nDim       = 2,                               &
!~                  X1         = (/ 0., 0. /),                    &
!~                  X2         = (/ 1., 1. /),                    &
!~                  IO         = IO,                              &
!~                  quiet      = .TRUE.                           )

 !       NULLIFY( Tens3GaussP )
 !       NULLIFY( Tens3GaussW )
        ! Compute integration points for volume integrals
        ! Local to this subroutine
        CALL HypercubeQuadraturePoints(                                     &
                 nIntGP     = nTens3GP,                                     &
                 IntGaussP  = Tens3GaussP,                                  &
                 IntGaussW  = Tens3GaussW,                                  &
                 M          = DISC%Galerkin%nPoly+DISC%Galerkin%nPolyMat+2, &
                 nDim       = 3,                                            &
                 X1         = (/ 0., 0., 0. /),                             &
                 X2         = (/ 1., 1., 1. /),                             &
                 IO         = IO,                                           &
                 quiet      = .TRUE.                                        )

        NULLIFY( TFMGaussP )
        NULLIFY( TFMGaussW )
        ! Compute integration points for surface integrals
        ! Local to this subroutine
        CALL HypercubeQuadraturePoints(                                     &
                 nIntGP     = nTFMGaussP,                                   &
                 IntGaussP  = TFMGaussP,                                    &
                 IntGaussW  = TFMGaussW,                                    &
                 M          = DISC%Galerkin%nPoly+DISC%Galerkin%nPolyMat+2, &
                 nDim       = 2,                                            &
                 X1         = (/ 0., 0. /),                                 &
                 X2         = (/ 1., 1. /),                                 &
                 IO         = IO,                                           &
                 quiet      = .TRUE.                                        )

        ! Compute the value and gradient of basis functions
        ! at volume integration points Tens3GaussP
        ALLOCATE( Tens3BaseFunc(  DISC%Galerkin%nDegFr,nTens3GP) )
        ALLOCATE( Tens3BaseGrad(3,DISC%Galerkin%nDegFr,nTens3GP) )

        DO iIntGP = 1, nTens3GP
            xi   = Tens3GaussP(1,iIntGP)
            eta  = Tens3GaussP(2,iIntGP)
            zeta = Tens3GaussP(3,iIntGP)
            DO l = 1, DISC%Galerkin%nDegFr
                CALL BaseFunc3D(phi_l,l,xi,eta,zeta,DISC%Galerkin%nPoly,    &
                                DISC%Galerkin%cPoly3D_Hex,                  &
                                DISC%Galerkin%NonZeroCPoly_Hex,             &
                                DISC%Galerkin%NonZeroCPolyIndex_Hex         )
                CALL BaseGrad3D(phigrad,l,xi,eta,zeta,DISC%Galerkin%nPoly,  &
                                DISC%Galerkin%cPoly3D_Hex,                  &
                                DISC%Galerkin%NonZeroCPoly_Hex,             &
                                DISC%Galerkin%NonZeroCPolyIndex_Hex         )
                Tens3BaseFunc(l,iIntGP)   = phi_l
                Tens3BaseGrad(:,l,iIntGP) = phigrad
            ENDDO
        ENDDO

        ! Compute stiffness tensors for the hexahedrons
        DISC%Galerkin%Kxi_k_Hex   = 0.0d0
        DISC%Galerkin%Keta_k_Hex  = 0.0d0
        DISC%Galerkin%Kzeta_k_Hex = 0.0d0
        DO iIntGP = 1, nTens3GP
            xi   = Tens3GaussP(1,iIntGP)
            eta  = Tens3GaussP(2,iIntGP)
            zeta = Tens3GaussP(3,iIntGP)
            DO m = 1, DISC%Galerkin%nDegFrMat
                phi_m = Tens3BaseFunc(m,iIntGP)
                DO l = 1, DISC%Galerkin%nDegFr
                    phi_l = Tens3BaseFunc(l,iIntGP)
                    DO k = 1, DISC%Galerkin%nDegFr
                        phigrad = Tens3BaseGrad(:,k,iIntGP)
                        DISC%Galerkin%Kxi_k_Hex(   k,l,m ) = DISC%Galerkin%Kxi_k_Hex(   k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phigrad(1)*phi_l*phi_m
                        DISC%Galerkin%Keta_k_Hex(  k,l,m ) = DISC%Galerkin%Keta_k_Hex(  k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phigrad(2)*phi_l*phi_m
                        DISC%Galerkin%Kzeta_k_Hex( k,l,m ) = DISC%Galerkin%Kzeta_k_Hex( k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phigrad(3)*phi_l*phi_m
                    ENDDO
                ENDDO
            ENDDO
            !
        ENDDO
        !
        DISC%Galerkin%Kxi_m_Hex   = 0.0d0
        DISC%Galerkin%Keta_m_Hex  = 0.0d0
        DISC%Galerkin%Kzeta_m_Hex = 0.0d0
        DO iIntGP = 1, nTens3GP
            xi   = Tens3GaussP(1,iIntGP)
            eta  = Tens3GaussP(2,iIntGP)
            zeta = Tens3GaussP(3,iIntGP)
            DO m = 1, DISC%Galerkin%nDegFrMat
                phigrad = Tens3BaseGrad(:,m,iIntGP)
                DO l = 1, DISC%Galerkin%nDegFr
                    phi_l = Tens3BaseFunc(l,iIntGP)
                    DO k = 1, DISC%Galerkin%nDegFr
                        phi_k = Tens3BaseFunc(k,iIntGP)
                        DISC%Galerkin%Kxi_m_Hex(   k,l,m ) = DISC%Galerkin%Kxi_m_Hex(   k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phi_k*phi_l*phigrad(1)
                        DISC%Galerkin%Keta_m_Hex(  k,l,m ) = DISC%Galerkin%Keta_m_Hex(  k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phi_k*phi_l*phigrad(2)
                        DISC%Galerkin%Kzeta_m_Hex( k,l,m ) = DISC%Galerkin%Kzeta_m_Hex( k,l,m ) +         &
                                                           Tens3GaussW(iIntGP)*phi_k*phi_l*phigrad(3)
                    ENDDO
                ENDDO
            ENDDO
            !
        ENDDO
        !
        WHERE(ABS(DISC%Galerkin%Kxi_k_Hex).LE.1e-10)
            DISC%Galerkin%Kxi_k_Hex = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Keta_k_Hex).LE.1e-10)
            DISC%Galerkin%Keta_k_Hex = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Kzeta_k_Hex).LE.1e-10)
            DISC%Galerkin%Kzeta_k_Hex = 0.0d0
        ENDWHERE
        !
        WHERE(ABS(DISC%Galerkin%Kxi_m_Hex).LE.1e-10)
            DISC%Galerkin%Kxi_m_Hex = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Keta_m_Hex).LE.1e-10)
            DISC%Galerkin%Keta_m_Hex = 0.0d0
        ENDWHERE
        WHERE(ABS(DISC%Galerkin%Kzeta_m_Hex).LE.1e-10)
            DISC%Galerkin%Kzeta_m_Hex = 0.0d0
        ENDWHERE
        !
        DO m = 1, DISC%Galerkin%nDegFrMat
        DO l = 1, DISC%Galerkin%nDegFr
        DO k = 1, DISC%Galerkin%nDegFr
          DISC%Galerkin%ADGxi_Hex(k,l,m)   = DISC%Galerkin%Kxi_k_Hex(l,k,m)   / DISC%Galerkin%MassMatrix_Hex(k,k,DISC%Galerkin%nPoly)
          DISC%Galerkin%ADGeta_Hex(k,l,m)  = DISC%Galerkin%Keta_k_Hex(l,k,m)  / DISC%Galerkin%MassMatrix_Hex(k,k,DISC%Galerkin%nPoly)
          DISC%Galerkin%ADGzeta_Hex(k,l,m) = DISC%Galerkin%Kzeta_k_Hex(l,k,m) / DISC%Galerkin%MassMatrix_Hex(k,k,DISC%Galerkin%nPoly)
        ENDDO
        ENDDO
        ENDDO
        !
        ALLOCATE(DISC%Galerkin%Kxi_k_Hex_Sp, DISC%Galerkin%Keta_k_Hex_Sp, DISC%Galerkin%Kzeta_k_Hex_Sp)
        ALLOCATE(DISC%Galerkin%Kxi_m_Hex_Sp, DISC%Galerkin%Keta_m_Hex_Sp, DISC%Galerkin%Kzeta_m_Hex_Sp)
        CALL IniSparseTensor3b(DISC%Galerkin%Kxi_k_Hex_Sp,   DISC%Galerkin%Kxi_k_Hex,   DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Keta_k_Hex_Sp,  DISC%Galerkin%Keta_k_Hex,  DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Kzeta_k_Hex_Sp, DISC%Galerkin%Kzeta_k_Hex, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Kxi_m_Hex_Sp,   DISC%Galerkin%Kxi_m_Hex,   DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Keta_m_Hex_Sp,  DISC%Galerkin%Keta_m_Hex,  DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%Kzeta_m_Hex_Sp, DISC%Galerkin%Kzeta_m_Hex, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        !
        ALLOCATE(DISC%Galerkin%ADGxi_Hex_Sp, DISC%Galerkin%ADGeta_Hex_Sp, DISC%Galerkin%ADGzeta_Hex_Sp)
        CALL IniSparseTensor3b(DISC%Galerkin%ADGxi_Hex_Sp,   DISC%Galerkin%ADGxi_Hex,   DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%ADGeta_Hex_Sp,  DISC%Galerkin%ADGeta_Hex,  DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3b(DISC%Galerkin%ADGzeta_Hex_Sp, DISC%Galerkin%ADGzeta_Hex, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
        !
        DISC%Galerkin%FluxInt_Hex = 0.0d0
        !
        DO iSide = 1, MESH%nSides_Hex
            !
            ! Compute for the element itself
            !
            DO iBndGP = 1, nTFMGaussP
                chiGP  = TFMGaussP(1,iBndGP)
                tauGP  = TFMGaussP(2,iBndGP)
                CALL HexaTrafoChiTau2XiEtaZeta(xi,eta,zeta,chiGP,tauGP,iSide,0)
                DO m = 1, DISC%Galerkin%nDegFrMat
                    CALL BaseFunc3D(phi_m,m,xi,eta,zeta,DISC%Galerkin%nPolyMat,      &
                                DISC%Galerkin%cPoly3D_Hex,                           &
                                DISC%Galerkin%NonZeroCPoly_Hex,                      &
                                DISC%Galerkin%NonZeroCPolyIndex_Hex                  )
                    DO l = 1, DISC%Galerkin%nDegFr
                    CALL BaseFunc3D(phi_l,l,xi,eta,zeta,DISC%Galerkin%nPoly,         &
                                DISC%Galerkin%cPoly3D_Hex,                           &
                                DISC%Galerkin%NonZeroCPoly_Hex,                      &
                                DISC%Galerkin%NonZeroCPolyIndex_Hex                  )
                        DO k = 1, DISC%Galerkin%nDegFr
                            CALL BaseFunc3D(phi_k,k,xi,eta,zeta,DISC%Galerkin%nPoly, &
                                            DISC%Galerkin%cPoly3D_Hex,               &
                                            DISC%Galerkin%NonZeroCPoly_Hex,          &
                                            DISC%Galerkin%NonZeroCPolyIndex_Hex      )
                            DISC%Galerkin%FluxInt_Hex(k,l,m,0,1,iSide) =  &
                            DISC%Galerkin%FluxInt_Hex(k,l,m,0,1,iSide) + TFMGaussW(iBndGP)*phi_k*phi_l*phi_m
                        ENDDO
                    ENDDO
                ENDDO
            ENDDO
            !
            ! Compute for the neighbors
            !
            DO iLocalNeighborSide = 1, MESH%nSides_Hex
            DO iLocalNeighborVrtx = 1, MESH%nVertices_Quad
            DO iBndGP = 1, nTFMGaussP
               chiGP  = TFMGaussP(1,iBndGP)
               tauGP  = TFMGaussP(2,iBndGP)
               CALL HexaTrafoChiTau2XiEtaZeta(xi,eta,zeta,chiGP,tauGP,iLocalNeighborSide,iLocalNeighborVrtx)
               CALL HexaTrafoChiTau2XiEtaZeta(xiS,etaS,zetaS,chiGP,tauGP,iSide,0)
               DO m = 1, DISC%Galerkin%nDegFrMat
                 CALL BaseFunc3D(phi_m,m,xiS,etaS,zetaS,DISC%Galerkin%nPolyMat,   &
                                 DISC%Galerkin%cPoly3D_Hex,                       &
                                 DISC%Galerkin%NonZeroCPoly_Hex,                  &
                                 DISC%Galerkin%NonZeroCPolyIndex_Hex              )
                 DO l = 1, DISC%Galerkin%nDegFr
                   CALL BaseFunc3D(phi_l,l,xi,eta,zeta,DISC%Galerkin%nPoly,       &
                                   DISC%Galerkin%cPoly3D_Hex,                     &
                                   DISC%Galerkin%NonZeroCPoly_Hex,                &
                                   DISC%Galerkin%NonZeroCPolyIndex_Hex            )
                   DO k = 1, DISC%Galerkin%nDegFr
                      CALL BaseFunc3D(phi_k,k,xiS,etaS,zetaS,DISC%Galerkin%nPoly, &
                                      DISC%Galerkin%cPoly3D_Hex,                  &
                                      DISC%Galerkin%NonZeroCPoly_Hex,             &
                                      DISC%Galerkin%NonZeroCPolyIndex_Hex         )
                      DISC%Galerkin%FluxInt_Hex(k,l,m,iLocalNeighborSide,iLocalNeighborVrtx,iSide) =  &
                      DISC%Galerkin%FluxInt_Hex(k,l,m,iLocalNeighborSide,iLocalNeighborVrtx,iSide) +  &
                      TFMGaussW(iBndGP)*phi_k*phi_l*phi_m
                   ENDDO
                 ENDDO
               ENDDO
               !
            ENDDO
            ENDDO
            ENDDO
        !
        ENDDO
        WHERE(ABS(DISC%Galerkin%FluxInt_Hex).LT.1e-10)
            DISC%Galerkin%FluxInt_Hex = 0.0d0
        ENDWHERE
        !
        ALLOCATE(DISC%Galerkin%FluxInt_Hex_Sp(0:MESH%nSides_Hex,1:MESH%nVertices_Quad,1:MESH%nSides_Hex))
        DO iSide = 1, MESH%nSides_Hex
            CALL IniSparseTensor3b(DISC%Galerkin%FluxInt_Hex_Sp(0,1,iSide), DISC%Galerkin%FluxInt_Hex(:,:,:,0,1,iSide), DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat)
            DO iLocalNeighborSide = 1, MESH%nSides_Hex
            DO iLocalNeighborVrtx = 1, MESH%nVertices_Quad
             CALL IniSparseTensor3b(DISC%Galerkin%FluxInt_Hex_Sp(iLocalNeighborSide,iLocalNeighborVrtx,iSide),    &
                                    DISC%Galerkin%FluxInt_Hex(:,:,:,iLocalNeighborSide,iLocalNeighborVrtx,iSide), &
                                    DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFrMat       )
            ENDDO
            ENDDO
        ENDDO

        DEALLOCATE( Tens3BaseFunc )
        DEALLOCATE( Tens3BaseGrad )

        DO iDegFr = 1, DISC%Galerkin%nDegFr
            DISC%Galerkin%iMassMatrix_Hex(iDegFr,iDegFr,DISC%Galerkin%nPoly) = &
                1./DISC%Galerkin%MassMatrix_Hex(iDegFr,iDegFr,DISC%Galerkin%nPoly)
        ENDDO

        DISC%Galerkin%IntGPBaseFunc_Hex(:,:,:)          = 0.0d0
        DISC%Galerkin%IntGPBaseGrad_Hex(:,:,:,:)        = 0.0d0
        DISC%Galerkin%BndGPBaseFunc3D_Hex(:,:,:)        = 0.0d0

        iPoly = DISC%Galerkin%nPoly

        DO iIntGP = 1, DISC%Galerkin%nIntGP
            DO iDegFr = 1, DISC%Galerkin%nDegFr
                !
                CALL BaseFunc3D(phi,iDegFr,                             &
                              DISC%Galerkin%intGaussP_Hex(1,iIntGP),    &
                              DISC%Galerkin%intGaussP_Hex(2,iIntGP),    &
                              DISC%Galerkin%intGaussP_Hex(3,iIntGP),    &
                              DISC%Galerkin%nPoly,                      &
                              DISC%Galerkin%cPoly3D_Hex,                &
                              DISC%Galerkin%NonZeroCPoly_Hex,           &
                              DISC%Galerkin%NonZeroCPolyIndex_Hex       )
                !
                DISC%Galerkin%IntGPBaseFunc_Hex(iDegFr,iIntGP,iPoly)   = phi
                !
                CALL BaseGrad3D( gradphixieta,iDegFr,                   &
                              DISC%Galerkin%intGaussP_Hex(1,iIntGP),    &
                              DISC%Galerkin%intGaussP_Hex(2,iIntGP),    &
                              DISC%Galerkin%intGaussP_Hex(3,iIntGP),    &
                              DISC%Galerkin%nPoly,                      &
                              DISC%Galerkin%cPoly3D_Hex,                &
                              DISC%Galerkin%NonZeroCPoly_Hex,           &
                              DISC%Galerkin%NonZeroCPolyIndex_Hex       )
                !
                DISC%Galerkin%IntGPBaseGrad_Hex(:,iDegFr,iIntGP,iPoly) = gradphixieta
                !
            ENDDO
        ENDDO
        !
        MESH%ELEM%BndBF_GP_Hex(:,:,:) = 0.0

        DO iSide = 1, MESH%nSides_Hex
            !
            ! Compute for the element itself
            !
            DO iBndGP = 1, DISC%Galerkin%nBndGP
               DO iDegFr = 1, DISC%Galerkin%nDegFr
                  !
                  chi  = DISC%Galerkin%bndGaussP_Hex(1,iBndGP)
                  tau  = DISC%Galerkin%bndGaussP_Hex(2,iBndGP)
                  !
                  CALL HexaTrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
                  !           !
                  CALL BaseFunc3D( phi,iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly, &
                                   DISC%Galerkin%cPoly3D_Hex,                  &
                                   DISC%Galerkin%NonZeroCPoly_Hex,             &
                                   DISC%Galerkin%NonZeroCPolyIndex_Hex         )
                  !
                  DISC%Galerkin%BndGPBaseFunc3D_Hex(iDegFr,iBndGP,iSide) = phi
                  MESH%ELEM%BndBF_GP_Hex(iDegFr,iBndGP,iSide) = phi
                  !
               ENDDO
            ENDDO
        ENDDO

        NULLIFY(DISC%Galerkin%NonZeroCPoly,DISC%Galerkin%NonZeroCPolyIndex,DISC%Galerkin%CPoly3D)

        CONTINUE

    ENDIF ! Hexas
    !
    logInfo(*) 'iniGalerkin successful '
    !
  END SUBROUTINE iniGalerkin3D_us_intern_new


  !===========================================================================!
  !!                                                                         !!
  !!  icGalerkin3D_us initializes the degrees of freedom at time t=0.0       !!
  !!  by L2 projection                                                       !!
  !!                                                                         !!
  !===========================================================================!


  SUBROUTINE icGalerkin3D_us_new(EQN, DISC, MESH, IC, SOURCE, IO)
    !-------------------------------------------------------------------------!

    USE DGBasis_mod
    USE COMMON_InitialField_mod
#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc
    use f_ftoc_bind_interoperability
#endif
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tInitialCondition)  :: IC
    TYPE(tSource)            :: SOURCE
    TYPE(tInputOutput)       :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER :: iElem                                                          ! Element number
    INTEGER :: iIntGP                                                         ! Index of internal Gausspoint
    INTEGER :: iDegFr                                                         ! Degree of freedom
    INTEGER :: iVar                                                           ! Variable number
    INTEGER :: iVert                                                          ! Vertex counter
    INTEGER :: iPoly, nIntGP, nDegFr, eType
    INTEGER :: LocPoly, LocDegFr                                              ! Variables for p-adaptivity
    REAL    :: xi, eta, zeta                                                  ! Reference coordinates
    REAL    :: xGP, yGP, zGP                                                  ! Physical coordinates
    REAL    :: x(MESH%nVertexMax)                                             ! Element vertices in physical coordinates system
    REAL    :: y(MESH%nVertexMax)                                             ! Element vertices in physical coordinates system
    REAL    :: z(MESH%nVertexMax)                                             ! Element vertices in physical coordinates system
    REAL    :: iniGP(EQN%nVar)                                                ! Initial state vector at Gausspoint
    REAL    :: iniGP_Ane(EQN%nAneFuncperMech*EQN%nMechanisms)                 ! Initial anelastic state vector at Gausspoint
    REAL    :: iniGP_Plast(6)                                                 ! Initial stress loading for plastic calculations
    REAL    :: phi                                                            ! Value of the base function at GP      !
    REAL    :: Einv, v                                                        ! Inverse of Young's modulus, Poisson ratio v
    !
    REAL, POINTER :: IntGaussP(:,:)     =>NULL()
    REAL, POINTER :: IntGaussW(:)       =>NULL()
    REAL, POINTER :: IntGPBaseFunc(:,:) =>NULL()
    REAL, POINTER :: MassMatrix(:,:)    =>NULL()
#ifdef GENERATEDKERNELS
    ! temporary degrees of freedom
    real    :: l_dofsUpdate(disc%galerkin%nDegFr, eqn%nVarTotal)
    real    :: l_initialLoading( NUMBER_OF_BASIS_FUNCTIONS, 6 )
    real    :: l_plasticParameters(4)
#endif
    !-------------------------------------------------------------------------!
    !
    IF(.NOT.DISC%Galerkin%init) THEN
       logError(*) 'icGalerkin: SeisSol Interface not initialized!!'
       STOP
    ENDIF
    !
#ifndef GENERATEDKERNELS
    DISC%Galerkin%dgvar  = 0.
#endif
    !
    ALLOCATE(EQN%Energy(3,1:MESH%nElem))
    EQN%Energy = 0.


    IF(EQN%Plasticity.EQ.1) THEN
      ALLOCATE(DISC%Galerkin%DOFStress(DISC%Galerkin%nDegFr,6,MESH%nElem), DISC%Galerkin%pstrain(7, MESH%nElem),&
               DISC%Galerkin%PlasticParameters(4,1:MESH%nElem), DISC%Galerkin%Strain_Matrix(6,6))
        !Initialization
        DISC%Galerkin%DOFStress = 0.
        DISC%Galerkin%pstrain = 0.
        DISC%Galerkin%PlasticParameters = 0.
        DISC%Galerkin%Strain_Matrix = 0.

        !Initialize the stress-strain relation matrix (mu and lambda should be element dependent)
        Einv = (EQN%lambda+EQN%mu)/(EQN%mu*(3*EQN%lambda+2*EQN%mu))!Inv of the Young's modulus
        v = EQN%lambda/(2*(EQN%lambda+EQN%mu)) !Poisson's ratio

        DISC%Galerkin%Strain_Matrix(1,1) = Einv
        DISC%Galerkin%Strain_Matrix(2,2) = Einv
        DISC%Galerkin%Strain_Matrix(3,3) = Einv
        DISC%Galerkin%Strain_Matrix(4,4) = 1/(2*EQN%mu)
        DISC%Galerkin%Strain_Matrix(5,5) = 1/(2*EQN%mu)
        DISC%Galerkin%Strain_Matrix(6,6) = 1/(2*EQN%mu)
        DISC%Galerkin%Strain_Matrix(1,2) = -v*Einv
        DISC%Galerkin%Strain_Matrix(1,3) = -v*Einv
        DISC%Galerkin%Strain_Matrix(2,1) = -v*Einv
        DISC%Galerkin%Strain_Matrix(2,3) = -v*Einv
        DISC%Galerkin%Strain_Matrix(3,1) = -v*Einv
        DISC%Galerkin%Strain_Matrix(3,2) = -v*Einv
    ENDIF


    logInfo0(*) 'DG initial condition projection... '
    !
    iPoly  = DISC%Galerkin%nPoly
    nIntGP = DISC%Galerkin%nIntGP
    nDegFr = DISC%Galerkin%nDegFr

#ifdef GENERATEDKERNELS
    !$omp parallel do schedule(static) shared(eqn, disc, mesh, ic, source, io, iPoly, nIntGp, nDegFr) private(iElem, iIntGP, iDegFr, iVar, iVert, eType, locPoly, locDegFr, xi, eta, zeta, xGp, yGp, zGp, x, y, z, iniGp, iniGp_ane, phi, intGaussP, intGaussW, intGPBaseFunc, massMatrix, l_dofsUpdate,l_initialLoading,l_plasticParameters, iniGP_plast)
#endif
    DO iElem = 1,MESH%nElem
#ifdef GENERATEDKERNELS
        l_dofsUpdate = 0
        l_initialLoading=0
        l_plasticParameters=0
#endif

        x = 0.; y = 0.; z = 0.;
        eType = MESH%LocalElemType(iElem)
        SELECT CASE(eType)
        CASE(4)
            DO iVert=1,MESH%nVertices_Tet
                x(iVert) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(iVert,iElem))
                y(iVert) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(iVert,iElem))
                z(iVert) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(iVert,iElem))
            ENDDO
            ! Point to the corresponding:
            ! Integration points
            intGaussP     => DISC%Galerkin%intGaussP_Tet
            intGaussW     => DISC%Galerkin%intGaussW_Tet
            ! Basis func values
            IntGPBaseFunc => DISC%Galerkin%IntGPBaseFunc_Tet(1:nDegFr,1:nIntGp,iPoly)
            ! Mass matrix
            MassMatrix    => DISC%Galerkin%MassMatrix_Tet(1:nDegFr,1:nDegFr,iPoly)
        CASE(6)
            DO iVert=1,MESH%nVertices_Hex
                x(iVert) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(iVert,iElem))
                y(iVert) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(iVert,iElem))
                z(iVert) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(iVert,iElem))
            ENDDO
            ! Point to the corresponding:
            ! Integration points
            intGaussP     => DISC%Galerkin%intGaussP_Hex
            intGaussW     => DISC%Galerkin%intGaussW_Hex
            ! Basis func values
            IntGPBaseFunc => DISC%Galerkin%IntGPBaseFunc_Hex(1:nDegFr,1:nIntGp,iPoly)
            ! Mass matrix
            MassMatrix    => DISC%Galerkin%MassMatrix_Hex(1:nDegFr,1:nDegFr,iPoly)
        END SELECT

        DO iIntGP = 1,nIntGP
            xi   = intGaussP(1,iIntGP)
            eta  = intGaussP(2,iIntGP)
            zeta = intGaussP(3,iIntGP)
            SELECT CASE(eType)
            CASE(4) ! Tetras
                CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,x,y,z)
            CASE(6) ! Hexas
                CALL HexaTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,x,y,z)
            END SELECT
            !
            ! Initial condition
            !
            CALL InitialField(iniGP, iniGP_ane, 0., xGP, yGP, zGP, iElem, EQN, IC, SOURCE, IO)
            !
            ! Projection onto the basis functions
            !
            DO iDegFr = 1, nDegFr
                phi = IntGPBaseFunc(iDegFr,iIntGP)
#ifdef GENERATEDKERNELS
                l_dofsUpdate(iDegFr, 1:EQN%nVar) = l_dofsUpdate(iDegFr, 1:EQN%nVar) + intGaussW(iIntGP)*iniGP(:)*phi
#else
                DISC%Galerkin%dgvar(iDegFr,1:EQN%nVar,iElem,1) =                                                    &
                DISC%Galerkin%dgvar(iDegFr,1:EQN%nVar,iElem,1) + IntGaussW(iIntGP)*iniGP(:)*phi
#endif
            ENDDO

            IF(EQN%Anelasticity.EQ.1) THEN
                ! Projection of anelastic functions
                DO iDegFr = 1, nDegFr
                  phi = IntGPBaseFunc(iDegFr,iIntGP)
#ifdef GENERATEDKERNELS
                  l_dofsUpdate(iDegFr, EQN%nVar+1:EQN%nVarTotal) = l_dofsUpdate(iDegFr, EQN%nVar+1:EQN%nVarTotal) + IntGaussW(iIntGP)*iniGP_Ane(:)*phi
#else
                  DISC%Galerkin%dgvar(iDegFr,EQN%nVar+1:EQN%nVarTotal,iElem,1) =                                  &
                  DISC%Galerkin%dgvar(iDegFr,EQN%nVar+1:EQN%nVarTotal,iElem,1) + IntGaussW(iIntGP)*iniGP_Ane(:)*phi
#endif
                ENDDO
            ENDIF

            IF(EQN%Plasticity.EQ.1 .AND. EQN%PlastMethod .EQ. 2) THEN !average approach for plasticity
            ! L2 projection of initial stress loading for the plastic calculations onto the DOFs
              iniGP_Plast(:) = EQN%IniStress(1:6,iElem)
              DO iDegFr = 1, nDegFr
                 phi = IntGPBaseFunc(iDegFr,iIntGP)
#ifdef GENERATEDKERNELS
                 l_initialLoading(iDegFr,1:6) = l_initialLoading(iDegFr,1:6) + IntGaussW(iIntGP)*iniGP_plast(:)*phi
#else
                DISC%Galerkin%DOFStress(iDegFr,1:6,iElem) = &
                DISC%Galerkin%DOFStress(iDegFr,1:6,iElem) + IntGaussW(iIntGP)*iniGP_plast(:)*phi
#endif
              ENDDO
           ENDIF

    ENDDO !iIntGP

            DO iDegFr = 1, nDegFr
#ifdef GENERATEDKERNELS
               l_dofsUpdate(iDegFr, :) = l_dofsUpdate( iDegFr, : ) / massMatrix(iDegFr,iDegFr)
               IF(EQN%Plasticity.EQ.1 .AND. EQN%PlastMethod .EQ. 2) THEN
                  l_initialLoading(iDegFr, :) = l_initialLoading( iDegFr, : ) / massMatrix(iDegFr,iDegFr)
               ENDIF
#else
               DISC%Galerkin%dgvar(iDegFr,:,iElem,1) = DISC%Galerkin%dgvar(iDegFr,:,iElem,1) / MassMatrix(iDegFr,iDegFr)
               IF(EQN%Plasticity.EQ.1 .AND. EQN%PlastMethod .EQ. 2) THEN
                  DISC%Galerkin%DOFStress(iDegFr,:,iElem) = DISC%Galerkin%DOFStress(iDegFr,:,iElem)/ MassMatrix(iDegFr, iDegFr)
               ENDIF
#endif
            ENDDO

            IF(EQN%Plasticity.EQ.1 .AND. EQN%PlastMethod .EQ. 0) THEN !high-order points approach
            !elementwise assignement of the initial loading
#ifdef GENERATEDKERNELS
               l_initialLoading(1,1:6) = EQN%IniStress(1:6,iElem)
#else
               DISC%Galerkin%DOFStress(1,1:6,iElem) = EQN%IniStress(1:6,iElem)
#endif
            ENDIF

#ifdef GENERATEDKERNELS
        ! write the update back
        call c_interoperability_addToDofs(  i_meshId           = iElem, \
                                            i_update           = l_dofsUpdate, \
                                            numberOfQuantities = eqn%nVarTotal )

#ifdef USE_PLASTICITY
        ! initialize the element dependent plastic parameters
        l_plasticParameters(1) = MESH%Elem%Volume(iElem)
        l_plasticParameters(2) = EQN%PlastCo(iElem) !element-dependent plastic cohesion
        l_plasticParameters(3) = EQN%Rho0    !density
        l_plasticParameters(4) = EQN%BulkFriction(iElem) !element-dependent bulk friction

        ! initialize loading in C
        call c_interoperability_setInitialLoading( i_meshId         = c_loc( iElem), \
                                                   i_initialLoading = c_loc( l_initialLoading ) )

        !initialize parameters in C
        call c_interoperability_setPlasticParameters( i_meshId         = c_loc( iElem), \
                                                   i_plasticParameters = c_loc( l_plasticParameters ) )
#endif

#else
        IF(EQN%Plasticity.EQ.1) THEN
          ! initialize plastic parameters in classic version
          DISC%Galerkin%plasticParameters(1,iElem) = MESH%Elem%Volume(iElem)
          DISC%Galerkin%plasticParameters(2,iElem) = EQN%PlastCo(iElem) !element-dependent plastic cohesion
          DISC%Galerkin%plasticParameters(3,iElem) = EQN%Rho0 !currently not needed inside the plasticity routine
          DISC%Galerkin%plasticParameters(4,iElem) = EQN%BulkFriction(iElem) !element-dependent bulk friction
        ENDIF
#endif


        NULLIFY(intGaussP)
        NULLIFY(intGaussW)
        NULLIFY(IntGPBaseFunc)
        NULLIFY(MassMatrix)


    ENDDO ! iElem

#ifdef USE_PLASTICITY
    call c_interoperability_setTv( tv = EQN%Tv )
#endif

#ifdef GENERATEDKERNELS
#ifdef USE_PLASTICITY
    ! TODO: redundant (see iniGalerkin3D_us_level2_new) call to ensure correct intitial loading in copy layers.
    call c_interoperability_synchronizeCellLocalData();
#endif
#endif

#ifndef GENERATEDKERNELS
    ! For p-adaptivity, kill all degrees of freedom that exceed the local
    ! polynomial order
    IF(DISC%Galerkin%pAdaptivity.GE.1) THEN
      DO iElem = 1, MESH%nElem
         LocPoly  = INT(DISC%Galerkin%LocPoly(iElem))
         LocDegFr = (LocPoly+1)*(LocPoly+2)*(LocPoly+3)/6
         DO iDegFr = LocDegFr+1, DISC%Galerkin%nDegFr
            DISC%Galerkin%dgvar(iDegFr,:,iElem,1) = 0.
         ENDDO
      ENDDO
    ENDIF
#endif
    !
    logInfo0(*) 'DG initial condition projection done. '
    !
  END SUBROUTINE icGalerkin3D_us_new


  !===========================================================================!
  !!                                                                         !!
  !! IniSparseStarMatrices precomputes the sparse matrices A,B,C star        !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE IniSparseStarMatrices3D_new(EQN, DISC, MESH, BND, SOURCE, MaterialVal, IO)
    !-------------------------------------------------------------------------!
    USE COMMON_operators_mod
    USE DGBasis_mod
    USE QuadPoints_mod
    USE ini_MODEL_mod
    USE JacobiNormal_mod
#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc
    use f_ftoc_bind_interoperability
#endif
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)                :: EQN
    TYPE(tDiscretization)           :: DISC
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tBoundary)                 :: BND
    TYPE(tSource)                   :: SOURCE
    TYPE(tInputOutput)              :: IO
    REAL                            :: MaterialVal(MESH%nElem,EQN%nBackgroundVar)
    REAL                            :: dt
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    REAL            :: grad(3,3)                                              ! dx/dxi    !
    REAL            :: Jacobit(3,3), JacobiDet                                ! dxi/dx    !
    REAL            :: xi,eta,zeta                                            ! Point in reference coord
    REAL            :: xGP,yGP,zGP                                            ! Point
    REAL            :: phi
    REAL            :: x(MESH%nVertexMax)                                     ! Vertex coords.                   !
    REAL            :: y(MESH%nVertexMax)                                     ! Vertex coords.                   !
    REAL            :: z(MESH%nVertexMax)                                     ! Vertex coords.                   !
    REAL            :: A(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: B(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: C(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: E(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: A_Star(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) ! Jacobians in the ref. element
    REAL            :: B_Star(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) ! Jacobians in the ref. element
    REAL            :: C_Star(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) ! Jacobians in the ref. element
    REAL            :: E_Star(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)
    REAL            :: FL(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: FR(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: FL_Star(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat,MESH%nSideMax)
    REAL            :: FR_Star(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat,MESH%nSideMax)
    REAL            :: NormalVect_n(3)
    REAL            :: NormalVect_s(3)
    REAL            :: NormalVect_t(3)
    REAL            :: T(EQN%nVarTotal,EQN%nVarTotal,MESH%nSideMax), iT(EQN%nVarTotal,EQN%nVarTotal,MESH%nSideMax)            ! Rotation matrix
    REAL            :: LocMaterialVal(EQN%nBackgroundVar)
    REAL            :: NeigMaterialVal(EQN%nBackgroundVar)
    INTEGER         :: iElem
    INTEGER         :: iSide
    INTEGER         :: iVert
    INTEGER         :: iIntGP
    INTEGER         :: iDegFr
    INTEGER         :: iNeighbor
    INTEGER         :: iObject
    INTEGER         :: MPIIndex
    INTEGER         :: LocElemType
    INTEGER         :: CleanUp
    !
    INTEGER         :: nTens3GP, JacIni, JacEnd
    REAL,POINTER    :: Tens3GaussP_Hex(:,:)    => NULL()
    REAL,POINTER    :: Tens3GaussW_Hex(:)      => NULL()
    REAL,POINTER    :: Tens3BaseFunc_Hex(:,:)  => NULL()
    REAL,POINTER    :: Tens3GaussP_Tet(:,:)    => NULL()
    REAL,POINTER    :: Tens3GaussW_Tet(:)      => NULL()
    REAL,POINTER    :: Tens3BaseFunc_Tet(:,:)  => NULL()
    !
    REAL,POINTER    :: Tens3GaussP(:,:)        => NULL()
    REAL,POINTER    :: Tens3GaussW(:)          => NULL()
    REAL,POINTER    :: Tens3BaseFunc(:,:)      => NULL()
    REAL,POINTER    :: MassMatrix(:,:,:)       => NULL()
#ifdef GENERATEDKERNELS
   real             :: l_sDivJ !< scalar \f$ |S_k| / |J_k|\f$, where \f$ | S_k | \f$ is the surface volume of a face and $\f |J_k| $\f the determinant of the Jacobian (transformation to ref. element).
   integer :: l_element, l_face, l_mpiDomain, l_mpiElement

   real, dimension(NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES, 3) :: l_starMatrices
   real, dimension(NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES, 4) :: l_nApNm1, l_nAmNm1
#ifndef NDEBUG
   integer          :: l_col, l_row !< loop variables
   real             :: l_debugVariable ! variable used for summation in debug routines
#endif
#endif

    ! register epik/scorep function IniSparseStarMatrices3D_new
    EPIK_FUNC_REG("IniSparseStarMatrices3D_new")
    SCOREP_USER_FUNC_DEFINE()

    !
    ! ------------------------------------------------------------------------!

    ! start the epik/scorep function IniSparseStarMatrices3D_new
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN( "IniSparseStarMatrices3D_new" )

    logInfo0(*) 'Initializing element local matrices.'

! Generated Kernels: Don't initialize the old sparse format to save memory.
#if !defined(GENERATEDKERNELS)
    ALLOCATE( DISC%Galerkin%AStar_Sp(MESH%nElem) )
    ALLOCATE( DISC%Galerkin%BStar_Sp(MESH%nElem) )
    ALLOCATE( DISC%Galerkin%CStar_Sp(MESH%nElem) )
    ALLOCATE( DISC%Galerkin%EStar_Sp(MESH%nElem) )
    ALLOCATE( DISC%Galerkin%FLStar_Sp(MESH%nElem,MESH%nSideMax) )
    ALLOCATE( DISC%Galerkin%FRStar_Sp(MESH%nElem,MESH%nSideMax) )
#endif

    SELECT CASE(DISC%Galerkin%CKMethod)
    CASE(1)
        ALLOCATE( DISC%Galerkin%InvSystemMatrix(MESH%nElem) )
    CASE DEFAULT
        CONTINUE
    ENDSELECT

    IF(MESH%nElem_Hex .GT. 0)THEN
        !
        NULLIFY( Tens3GaussP_Hex )
        NULLIFY( Tens3GaussW_Hex )
        CALL HypercubeQuadraturePoints(                                     &
                 nIntGP     = nTens3GP,                                     &
                 IntGaussP  = Tens3GaussP_Hex,                              &
                 IntGaussW  = Tens3GaussW_Hex,                              &
                 M          = DISC%Galerkin%nPoly+DISC%Galerkin%nPolyMat+1, &
                 nDim       = 3,                                            &
                 X1         = (/ 0., 0., 0. /),                             &
                 X2         = (/ 1., 1., 1. /),                             &
                 IO         = IO,                                           &
                 quiet      = .TRUE.                                        )
        ALLOCATE( Tens3BaseFunc_Hex(DISC%Galerkin%nDegFrMat,nTens3GP) )
        DO iIntGp = 1, nTens3GP
           xi   = Tens3GaussP_Hex(1,iIntGp)
           eta  = Tens3GaussP_Hex(2,iIntGp)
           zeta = Tens3GaussP_Hex(3,iIntGp)
           DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                CALL BaseFunc3D(Tens3BaseFunc_Hex(iDegFr,iIntGp),iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly, &
                                DISC%Galerkin%cPoly3D_Hex,                                               &
                                DISC%Galerkin%NonZeroCPoly_Hex,                                          &
                                DISC%Galerkin%NonZeroCPolyIndex_Hex                                      )
           ENDDO ! iDegFr
        ENDDO ! iIntGp
        NULLIFY(DISC%Galerkin%NonZeroCPoly,DISC%Galerkin%NonZeroCPolyIndex,DISC%Galerkin%CPoly3D)
    ENDIF
    IF(MESH%nElem_Tet .GT. 0)THEN
        !
        NULLIFY( Tens3GaussP_Tet )
        NULLIFY( Tens3GaussW_Tet )
        CALL TetrahedronQuadraturePoints(                                   &
                 nIntGP     = nTens3GP,                                     &
                 IntGaussP  = Tens3GaussP_Tet,                              &
                 IntGaussW  = Tens3GaussW_Tet,                              &
                 M          = DISC%Galerkin%nPoly+DISC%Galerkin%nPolyMat+1, &
                 IO         = IO,                                           &
                 quiet      = .TRUE.                                        )
        ALLOCATE( Tens3BaseFunc_Tet(DISC%Galerkin%nDegFrMat,nTens3GP) )
        DO iIntGp = 1, nTens3GP
           xi   = Tens3GaussP_Tet(1,iIntGp)
           eta  = Tens3GaussP_Tet(2,iIntGp)
           zeta = Tens3GaussP_Tet(3,iIntGp)
           DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                CALL BaseFunc3D(Tens3BaseFunc_Tet(iDegFr,iIntGp),iDegFr,xi,eta,zeta,DISC%Galerkin%nPoly, &
                                DISC%Galerkin%cPoly3D_Tet,                                               &
                                DISC%Galerkin%NonZeroCPoly_Tet,                                          &
                                DISC%Galerkin%NonZeroCPolyIndex_Tet                                      )
           ENDDO ! iDegFr
        ENDDO ! iIntGp
        NULLIFY(DISC%Galerkin%NonZeroCPoly,DISC%Galerkin%NonZeroCPolyIndex,DISC%Galerkin%CPoly3D)
    ENDIF

#ifdef GENERATEDKERNELS
#ifdef NDEBUG
    !$omp parallel do schedule(static) shared(eqn, disc, mesh, bnd, source, materialVal, io, Tens3GaussP_Hex, Tens3GaussW_Hex, Tens3BaseFunc_Hex,  Tens3GaussP_Tet, Tens3GaussW_Tet, Tens3BaseFunc_Tet, nTens3GP ) private(grad, Jacobit, iElem, iSide, iVert, iIntGP, iDegFr, iNeighbor, iObject, A,B,C,E,  A_Star,B_Star,C_Star,E_Star,  FL,FR,  FL_Star, FR_Star,  x,y,z, Tens3GaussP, Tens3GaussW, Tens3BaseFunc, MassMatrix, JacobiDet,  xi,eta,zeta,  xGP,yGP,zGP,  phi, NormalVect_n, NormalVect_s, NormalVect_t,  T, iT,  LocMaterialVal,  NeigMaterialVal,  MPIIndex, LocElemType, CleanUp, JacIni, JacEnd, l_sDivJ, l_starMatrices, l_nApNm1, l_nAmNm1 )
#else
    !$omp parallel do schedule(static) shared(eqn, disc, mesh, bnd, source, materialVal, io, Tens3GaussP_Hex, Tens3GaussW_Hex, Tens3BaseFunc_Hex,  Tens3GaussP_Tet, Tens3GaussW_Tet, Tens3BaseFunc_Tet, nTens3GP ) private(grad, Jacobit, iElem, iSide, iVert, iIntGP, iDegFr, iNeighbor, iObject, A,B,C,E,  A_Star,B_Star,C_Star,E_Star,  FL,FR,  FL_Star, FR_Star,  x,y,z, Tens3GaussP, Tens3GaussW, Tens3BaseFunc, MassMatrix, JacobiDet,  xi,eta,zeta,  xGP,yGP,zGP,  phi, NormalVect_n, NormalVect_s, NormalVect_t,  T, iT,  LocMaterialVal,  NeigMaterialVal,  MPIIndex, LocElemType, CleanUp, JacIni, JacEnd, l_sDivJ, l_col, l_row, l_debugVariable, l_starMatrices, l_nApNm1, l_nAmNm1 )
#endif
#endif
    DO iElem = 1, MESH%nElem
        !
        LocElemType = MESH%LocalElemType(iElem)
        !
        A(:,:)        = 0.
        B(:,:)        = 0.
        C(:,:)        = 0.
        E(:,:)        = 0.
        FL(:,:)       = 0.
        FR(:,:)       = 0.
        !
        A_Star(:,:,:)    = 0.
        B_Star(:,:,:)    = 0.
        C_Star(:,:,:)    = 0.
        E_Star(:,:,:)    = 0.
        FL_Star(:,:,:,:) = 0.
        FR_Star(:,:,:,:) = 0.
        !
        x(:) = 0.
        y(:) = 0.
        z(:) = 0.
        !
        grad(:,:)    = 0.
        JacobiT(:,:) = 0.
        SELECT CASE(LocElemType)
        CASE(4)
            DO iVert=1,MESH%nVertices_Tet
                x(iVert) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(iVert,iElem))
                y(iVert) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(iVert,iElem))
                z(iVert) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(iVert,iElem))
            ENDDO
            ! Point to the corresponding:
            ! Integration points
            Tens3GaussP   => Tens3GaussP_Tet
            Tens3GaussW   => Tens3GaussW_Tet
            ! Basis func values
            Tens3BaseFunc => Tens3BaseFunc_Tet
            ! Mass matrix
            MassMatrix    => DISC%Galerkin%MassMatrix_Tet

            ! Compute the Jacobi of the mapping x=x(xi,eta,zeta)
            ! grad = dx/dxi
            ! For tetras is constant inside the element
            CALL TetraTrafoXiEtaZetaGrad(grad,x,y,z)
            CALL MatrixInverse3x3(JacobiT,grad)
        CASE(6)
            DO iVert=1,MESH%nVertices_Hex
                x(iVert) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(iVert,iElem))
                y(iVert) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(iVert,iElem))
                z(iVert) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(iVert,iElem))
            ENDDO
            ! Point to the corresponding:
            ! Integration points
            Tens3GaussP   => Tens3GaussP_Hex
            Tens3GaussW   => Tens3GaussW_Hex
            ! Basis func values
            Tens3BaseFunc => Tens3BaseFunc_Hex
            ! Mass matrix
            MassMatrix    => DISC%Galerkin%MassMatrix_Hex

            ! Compute the Jacobi of the mapping x=x(xi,eta,zeta)
            ! grad = dx/dxi
            ! For hexas (only parallelograms) is constant inside the element
            xi = 0.5d0; eta = 0.5d0; zeta = 0.5d0; ! Baricenter
            CALL HexaTrafoXiEtaZetaGrad(grad,xi,eta,zeta,x,y,z)
            ! JacobiT = grad^-1 = dxi/dx
            CALL MatrixInverse3x3(JacobiT,grad)
        END SELECT ! LocElemType
        WHERE(ABS(grad).LT.1e-10)
         grad = 0.0d0
        ENDWHERE

        !
        DO iSide = 1,LocElemType

            ! Local side's normal and tangential vectors
            NormalVect_n = DISC%Galerkin%geoNormals(1:3,iSide,iElem)
            NormalVect_s = DISC%Galerkin%geoTangent1(1:3,iSide,iElem)
            NormalVect_t = DISC%Galerkin%geoTangent2(1:3,iSide,iElem)

            CALL RotationMatrix3D(NormalVect_n,NormalVect_s,NormalVect_t,T(:,:,iSide),iT(:,:,iSide),EQN)

        ENDDO ! iSide
#ifdef GENERATEDKERNELS
#ifndef NDEBUG
            ! set debug variable to zero
            l_debugVariable = 0.0d0
#endif
#endif
        DO iIntGP = 1, nTens3GP
            xi   = Tens3GaussP(1,iIntGP)
            eta  = Tens3GaussP(2,iIntGP)
            zeta = Tens3GaussP(3,iIntGP)
            SELECT CASE(LocElemType)
            CASE(4)
                CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,x,y,z)
            CASE(6)
                CALL HexaTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,x,y,z)
            END SELECT ! LocElemType
            LocMaterialVal(:) = MaterialVal(iElem,:) ! rho,mu,lambda
            !
            CALL ModelDefinition_new(A,B,C,E,xGP,yGP,zGP,LocMaterialVal,             &
                                     EQN%LocAnelastic(iElem),                        &
                                     EQN%LocPoroelastic(iElem),                      &
                                     EQN,MESH,DISC,SOURCE,IO                         )
            !
            DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                phi = Tens3BaseFunc(iDegFr,iIntGP)

                A_Star(:,:,iDegFr) = A_Star(:,:,iDegFr) + Tens3GaussW(iIntGP)*phi*                     &
                                    ( A(:,:)*JacobiT(1,1) + B(:,:)*JacobiT(1,2) + C(:,:)*JacobiT(1,3) )
                B_Star(:,:,iDegFr) = B_Star(:,:,iDegFr) + Tens3GaussW(iIntGP)*phi*                     &
                                    ( A(:,:)*JacobiT(2,1) + B(:,:)*JacobiT(2,2) + C(:,:)*JacobiT(2,3) )
                C_Star(:,:,iDegFr) = C_Star(:,:,iDegFr) + Tens3GaussW(iIntGP)*phi*                     &
                                    ( A(:,:)*JacobiT(3,1) + B(:,:)*JacobiT(3,2) + C(:,:)*JacobiT(3,3) )
                E_Star(:,:,iDegFr) = E_Star(:,:,iDegFr) + Tens3GaussW(iIntGP)*phi* E(:,:)
#ifdef GENERATEDKERNELS
#ifndef NDEBUG
                ! sum over gaussian weights
                ! TODO: Find literature about this procedure? Seems useless, at least in the case of constant material parameters.. (breuera)
                l_debugVariable = l_debugVariable + Tens3GaussW(iIntGP)*phi
#endif
#endif
            ENDDO ! iDegFr

            DO iSide = 1,LocElemType
                !
                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
                    iObject         = MESH%ELEM%BoundaryToObject(iSide,iElem)
                    MPIIndex        = MESH%ELEM%MPINumber(iSide,iElem)
                    NeigMaterialVal(1:3) = BND%ObjMPI(iObject)%NeighborBackground(1:3,MPIIndex) ! rho,mu,lambda
                ELSE
                    SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
                    CASE(0)
                        iNeighbor       = MESH%ELEM%SideNeighbor(iSide,iElem)
                        NeigMaterialVal(:) = MaterialVal(iNeighbor,:) ! rho,mu,lambda
                    CASE DEFAULT ! For boundary conditions take inside material
                        NeigMaterialVal(:) = LocMaterialVal(:)           ! rho,mu,lambda
                    END SELECT
                ENDIF
                !
                CALL SolveRiemannProblem( FL               = FL,                               &
                                          FR               = FR,                               &
                                          A                = A,                                &
                                          B                = B,                                &
                                          C                = C,                                &
                                          T                = T(:,:,iSide),                     &
                                          iT               = iT(:,:,iSide),                    &
                                          MaterialValLeft  = LocMaterialVal,                   &
                                          MaterialValRight = NeigMaterialVal,                  &
                                          Solver           = DISC%Galerkin%FluxMethod,         &
                                          Boundary         = MESH%ELEM%Reference(iSide,iElem)  )
                !
                DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                    phi = Tens3BaseFunc(iDegFr,iIntGP)
                    FL_Star(:,:,iDegFr,iSide) = FL_Star(:,:,iDegFr,iSide) + Tens3GaussW(iIntGP)*phi*FL
                    FR_Star(:,:,iDegFr,iSide) = FR_Star(:,:,iDegFr,iSide) + Tens3GaussW(iIntGP)*phi*FR
                ENDDO ! iDegFr
                !
            ENDDO ! iSide
            !
            CONTINUE
            !
        ENDDO ! iIntGP
        !
        DO iDegFr = 1, DISC%Galerkin%nDegFrMat
#ifdef GENERATEDKERNELS
#ifndef NDEBUG
            ! assert that the sum over the gaussian weights cancels out with the first entry of the mass matrix (volume of the reference tetrahedron)
            if( abs( l_debugVariable - MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly) ) .ge. ZEROTOLERANCE ) then
              logError(*) 'IniSparseStarMatrices3D_new, MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly) not equal to sum over gaussian weights.', MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly), l_debugVariable
              stop
            endif
#endif
#endif
            A_Star(:,:,iDegFr) = A_Star(:,:,iDegFr)/MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly)
            B_Star(:,:,iDegFr) = B_Star(:,:,iDegFr)/MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly)
            C_Star(:,:,iDegFr) = C_Star(:,:,iDegFr)/MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly)
        ENDDO ! iDegFr
        DO iSide = 1,LocElemType
            DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                FL_Star(:,:,iDegFr,iSide) = FL_Star(:,:,iDegFr,iSide)/MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly)
                FR_Star(:,:,iDegFr,iSide) = FR_Star(:,:,iDegFr,iSide)/MassMatrix(iDegFr,iDegFr,DISC%Galerkin%nPoly)
            ENDDO ! iDegFr
        ENDDO ! iSide
        !
! Generated Kernels: Don't initialize the old sparse format to save memory.
#if !defined(GENERATEDKERNELS)
        ! Build sparse tensors
        CALL IniSparseTensor3(DISC%Galerkin%AStar_Sp(iElem),A_Star,EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3(DISC%Galerkin%BStar_Sp(iElem),B_Star,EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3(DISC%Galerkin%CStar_Sp(iElem),C_Star,EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)
        CALL IniSparseTensor3(DISC%Galerkin%EStar_Sp(iElem),E_Star,EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)

        DO iSide = 1,LocElemType
            CALL IniSparseTensor3(DISC%Galerkin%FLStar_Sp(iElem,iSide),FL_Star(:,:,:,iSide),EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)
            CALL IniSparseTensor3(DISC%Galerkin%FRStar_Sp(iElem,iSide),FR_Star(:,:,:,iSide),EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)
        ENDDO ! iSide
#endif

#ifdef GENERATEDKERNELS
#ifndef NDEBUG
        ! assert contant material parameters per element
        if ( disc%galerkin%nDegFrMat .ne. 1 ) then
          logError(*) 'IniSparseStarMatrices3D_new, disc%galerkin%nDegFrMat not equal 1.', disc%galerkin%nDegFrMat
          stop
        endif

        ! assert source matrix is not set
        ! \todo @breuera What used to be the purpose of this?
!~         do l_col = lbound(E_Star, 2), ubound(E_Star, 2)
!~           do l_row = lbound(E_Star, 1), ubound(E_Star, 1)
!~             if ( abs(E_Star(l_col, l_row, 1)) .ge. ZEROTOLERANCE ) then
!~               logError(*) 'IniSparseStarMatrices3D_new, E_Star not equal zero.', E_Star(l_col, l_row, 1), l_col, l_row
!~               stop
!~             endif
!~           enddo
!~         enddo

        ! assert 4 faces (tetrahedrons)
        if ( locElemType .ne. 4 ) then
          logError(*) 'IniSparseStarMatrices3D_new, locElemType not equal 4.', locElemType
          stop
        endif
#endif
        l_starMatrices(:,1) = reshape( transpose(A_Star(:,:,1)), (/ NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES /) )
        l_starMatrices(:,2) = reshape( transpose(B_Star(:,:,1)), (/ NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES /) )
        l_starMatrices(:,3) = reshape( transpose(C_Star(:,:,1)), (/ NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES /) )

        ! iterate over faces of local element
        do iSide = 1, 4
          ! initialize "flux solvers"  \f$ N_{k,i} A_k^+ N_{k,i}^{-1} \f$ and \f$  N_{k,i} A_k^+ N_{k,i}^{-1} \f$,
          l_nApNm1(:, iSide) = reshape( transpose(FL_Star(:,:,1,iSide)), (/ NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES /) )
          l_nAmNm1(:, iSide) = reshape( transpose(FR_Star(:,:,1,iSide)), (/ NUMBER_OF_QUANTITIES*NUMBER_OF_QUANTITIES /) )

          ! scale with the scalar \f$ |S_k|/|J_k| \f$ and take care about the sign (we need to subtract the flux terms)
          l_sDivJ = -disc%galerkin%geoSurfaces(iSide,iElem) / ( 3.0d0 * mesh%elem%volume(iElem) )
          l_nApNm1(:, iSide) = l_nApNm1(:, iSide) * l_sDivJ
          l_nAmNm1(:, iSide) = l_nAmNm1(:, iSide) * l_sDivJ
        enddo

        ! initialize the constant data for this element in the overlapping structures
!~         call c_interoperability_setCellLocalData( i_meshId       = c_loc(iElem),          \
!~                                                   i_starMatrices = c_loc(l_starMatrices), \
!~                                                   i_nApNm1       = c_loc(l_nApNm1),       \
!~                                                   i_nAmNm1       = c_loc(l_nAmNm1)        )
#endif
        !
        NULLIFY(Tens3GaussP)
        NULLIFY(Tens3GaussW)
        NULLIFY(Tens3BaseFunc)
        NULLIFY(MassMatrix)
        !
        CONTINUE
        !
    ENDDO ! iElem
    !
    logInfo0(*) 'Initializing element local matrices done.'
    !

    ! end the epik/scorep function IniSparseStarMatrices3D_new
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()

  END SUBROUTINE IniSparseStarMatrices3D_new


 !===========================================================================!
  !!                                                                         !!
  !! NonConformingGPEvaluation2D precomputes the values of the            !!
  !!  basis functions at the Gauss Integration points such that these        !!
  !!  values have to be muliplied and summed only with the DOFs of the       !!
  !!  underlying elements in every time step                                 !!
  !!                                                                         !!
  !===========================================================================!

!  SUBROUTINE NonConformingGPEvaluation3D(DISC, MESH, IO, BND)
!    !-------------------------------------------------------------------------!
!
!    USE DGBasis_mod
!    !-------------------------------------------------------------------------!
!    IMPLICIT NONE
!    !-------------------------------------------------------------------------!
!    ! Argument list declaration                                               !
!    TYPE(tDiscretization)           :: DISC
!    TYPE(tUnstructMesh)             :: MESH
!    TYPE(tInputOutput)              :: IO
!    TYPE(tBoundary)                 :: BND
!    !-------------------------------------------------------------------------!
!    ! Local variable declaration                                              !
!    INTEGER                         :: i
!    INTEGER                         :: nIntGP, iIntGP
!    INTEGER                         :: iNeighbor, incNeighbor
!    INTEGER                         :: iDegFr
!    INTEGER                         :: MPIIndex
!    INTEGER                         :: iObject
!     INTEGER                         :: LocElemType, LocVrtxType
!    REAL                            :: xGP,yGP,zGP
!    REAL                            :: xi,eta,zeta
!    REAL                            :: x_host(MESH%nVertexMax), y_host(MESH%nVertexMax), z_host(MESH%nVertexMax)
!    ! ------------------------------------------------------------------------!
!
!    logInfo(*) 'Pre-computing basis functions '
!    logInfo(*) '  at Gauss points at non-conforming boundary...'
!
!    nIntGP = DISC%Galerkin%nBndGP
!    ALLOCATE(DISC%Galerkin%NC_BndBF_GP(DISC%Galerkin%nDegFr,nIntGP,MESH%nNonConformingEdges) )
!
!    DO i = 1, MESH%nNonConformingEdges
!      DO iIntGP = 1, nIntGP
!
!        xGP = MESH%ELEM%ncBndGaussP(i,iIntGP,1)
!        yGP = MESH%ELEM%ncBndGaussP(i,iIntGP,2)
!        zGP = MESH%ELEM%ncBndGaussP(i,iIntGP,3)
!
!        incNeighbor = MESH%ELEM%ncBndNeighbor(i,iIntGP)
!
!        ! Maybe it is possible to do it in a more elegant way...
!        ! ncNeighbor = -1 means MPI boundary
!        IF(incNeighbor .EQ. -1)THEN
!            iObject = MESH%ELEM%NC_BoundaryToObject(i,iIntGP)
!            MPIIndex = MESH%ELEM%MPI_NCNumber(i,iIntGP)
!            IF(MPIIndex.EQ.-1) THEN
!            PRINT*, 'NCboundary ', i, ' and GP ', iIntGP, ' have nc neighbor ', incneighbor
!                logError(*) 'Galerkin_initialize MPIIndex = -1 !'
!                STOP
!            ENDIF
!            LocElemType = BND%ObjMPI(iObject)%NC_LocalElemType(MPIIndex)
!            SELECT CASE(LocElemType)
!            CASE(4)
!                LocVrtxType = 4
!                x_host(1:LocVrtxType) = BND%ObjMPI(iObject)%NC_NeighborCoords(1,1:LocVrtxType,MPIIndex)
!                y_host(1:LocVrtxType) = BND%ObjMPI(iObject)%NC_NeighborCoords(2,1:LocVrtxType,MPIIndex)
!                z_host(1:LocVrtxType) = BND%ObjMPI(iObject)%NC_NeighborCoords(3,1:LocVrtxType,MPIIndex)
!                !
!                CALL TrafoXYZ2XiEtaZeta(xi,eta,zeta,xGP,yGP,zGP,x_host,y_host,z_host,LocVrtxType)
!                ! compute the values of the basis functions in a tetrahedron
!                DO iDegFr = 1, DISC%Galerkin%nDegFr
!                    CALL BaseFunc3D(                                                       &
!                      value             = DISC%Galerkin%NC_BndBF_GP(iDegFr,iIntGP,i),      &
!                      iDegFr            = iDegFr,                                          &
!                      xi                = xi,                                              &
!                      eta               = eta,                                             &
!                      zeta              = zeta,                                            &
!                      nPoly             = DISC%Galerkin%nPoly,                             &
!                      cPoly3D           = DISC%Galerkin%CPoly3D_Tet,                       &
!                      NonZeroCPoly      = DISC%Galerkin%NonZeroCPoly_Tet,                  &
!                      NonZeroCPolyIndex = DISC%Galerkin%NonZeroCPolyIndex_Tet              )
!                ENDDO
!            CASE(6)
!                LocVrtxType = 8
!                x_host(1:LocVrtxType) = BND%ObjMPI(iObject)%NC_NeighborCoords(1,1:LocVrtxType,MPIIndex)
!                y_host(1:LocVrtxType) = BND%ObjMPI(iObject)%NC_NeighborCoords(2,1:LocVrtxType,MPIIndex)
!                z_host(1:LocVrtxType) = BND%ObjMPI(iObject)%NC_NeighborCoords(3,1:LocVrtxType,MPIIndex)
!                CALL TrafoXYZ2XiEtaZeta(xi,eta,zeta,xGP,yGP,zGP,x_host,y_host,z_host,LocVrtxType)
!                ! compute the values of the basis functions in a tetrahedron
!                DO iDegFr = 1, DISC%Galerkin%nDegFr
!                    CALL BaseFunc3D(                                                       &
!                      value             = DISC%Galerkin%NC_BndBF_GP(iDegFr,iIntGP,i),      &
!                      iDegFr            = iDegFr,                                          &
!                      xi                = xi,                                              &
!                      eta               = eta,                                             &
!                      zeta              = zeta,                                            &
!                      nPoly             = DISC%Galerkin%nPoly,                             &
!                      cPoly3D           = DISC%Galerkin%CPoly3D_Hex,                       &
!                      NonZeroCPoly      = DISC%Galerkin%NonZeroCPoly_Hex,                  &
!                      NonZeroCPolyIndex = DISC%Galerkin%NonZeroCPolyIndex_Hex              )
!                ENDDO
!            END SELECT
!        ELSE
!            iNeighbor   = MESH%ELEM%NCB_IndexList(incNeighbor,1)
!            LocElemType = MESH%LocalElemType(iNeighbor)
!            SELECT CASE(LocElemType)
!            CASE(4)
!                LocVrtxType = 4
!                ! get the location in the reference element of Neighbor
!                x_host(1:LocVrtxType) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:LocVrtxType,iNeighbor))
!                y_host(1:LocVrtxType) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:LocVrtxType,iNeighbor))
!                z_host(1:LocVrtxType) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:LocVrtxType,iNeighbor))
!                CALL TrafoXYZ2XiEtaZeta(xi,eta,zeta,xGP,yGP,zGP,x_host,y_host,z_host,LocVrtxType)
!                ! compute the values of the basis functions in a tetrahedron
!                DO iDegFr = 1, DISC%Galerkin%nDegFr
!                    CALL BaseFunc3D(                                                       &
!                      value             = DISC%Galerkin%NC_BndBF_GP(iDegFr,iIntGP,i),      &
!                      iDegFr            = iDegFr,                                          &
!                      xi                = xi,                                              &
!                      eta               = eta,                                             &
!                      zeta              = zeta,                                            &
!                      nPoly             = DISC%Galerkin%nPoly,                             &
!                      cPoly3D           = DISC%Galerkin%CPoly3D_Tet,                       &
!                      NonZeroCPoly      = DISC%Galerkin%NonZeroCPoly_Tet,                  &
!                      NonZeroCPolyIndex = DISC%Galerkin%NonZeroCPolyIndex_Tet              )
!                ENDDO
!            CASE(6)
!                LocVrtxType = 8
!                ! get the location in the reference element of Neighbor
!                x_host(1:LocVrtxType) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:LocVrtxType,iNeighbor))
!                y_host(1:LocVrtxType) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:LocVrtxType,iNeighbor))
!                z_host(1:LocVrtxType) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:LocVrtxType,iNeighbor))
!                CALL TrafoXYZ2XiEtaZeta(xi,eta,zeta,xGP,yGP,zGP,x_host,y_host,z_host,LocVrtxType)
!                ! compute the values of the basis functions in a tetrahedron
!                DO iDegFr = 1, DISC%Galerkin%nDegFr
!                    CALL BaseFunc3D(                                                       &
!                      value             = DISC%Galerkin%NC_BndBF_GP(iDegFr,iIntGP,i),      &
!                      iDegFr            = iDegFr,                                          &
!                      xi                = xi,                                              &
!                      eta               = eta,                                             &
!                      zeta              = zeta,                                            &
!                      nPoly             = DISC%Galerkin%nPoly,                             &
!                      cPoly3D           = DISC%Galerkin%CPoly3D_Hex,                       &
!                      NonZeroCPoly      = DISC%Galerkin%NonZeroCPoly_Hex,                  &
!                      NonZeroCPolyIndex = DISC%Galerkin%NonZeroCPolyIndex_Hex              )
!                ENDDO
!            END SELECT
!        ENDIF ! MPI
!      ENDDO ! iIntGP
!    ENDDO   ! nNonConformingEdges
!    !
!    !
!    CONTINUE
!    !
!  END SUBROUTINE NonConformingGPEvaluation3D


  SUBROUTINE BuildSpecialDGGeometry3D_new(MaterialVal,EQN,MESH,DISC,BND,MPI,IO)

    USE common_operators_mod
    USE DGbasis_mod
    USE ini_faultoutput_mod
#ifdef HDF
    USE hdf_faultoutput_mod
#endif
    use, intrinsic :: iso_c_binding

    !-------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------!
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    ! Argument list declaration
    TYPE(tEquations)         :: EQN
    TYPE(tUnstructMesh)      :: MESH
    TYPE(tDiscretization)    :: DISC
    TYPE(tBoundary)          :: BND
    TYPE(tMPI)               :: MPI
    TYPE(tInputOutput)       :: IO
    REAL                     :: MaterialVal(MESH%nElem,EQN%nBackgroundVar)          !
    !-------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER :: iElem, iSide, iDomain            ! Loop counter                !
    REAL    :: sidevec(3,2)                     ! Tmp vertex connection vector!
    INTEGER :: VertexSide_Tet(MESH%nSides_Tet,MESH%nVertices_Tri)
    INTEGER :: VertexSide_Hex(MESH%nSides_Hex,MESH%nVertices_Quad)              ! # sides = 6, # vertices per side = 4
    INTEGER :: allocstat                        ! Status of allocation        !
    REAL    :: Length                           ! Length of tangent vector    !
    REAL    :: BaryVec(EQN%Dimension,MESH%nSideMax), Dist(MESH%nSideMax)
    REAL    :: minv
    INTEGER :: minl(1)
    INTEGER :: i, j, iLayer, iZone
    REAL    :: nx, ny, nz, sx, sy, sz, tx, ty, tz, rho, amax, a
    REAL    :: c(6,6), Voigt_rot(6,6), TT(6,6), T(6,6)
    REAL    :: coefficients(4), coefficients2(5)
    REAL    :: Re_solution(3), Im_solution(3)
    REAL    :: K_F, K_S, K_Mean, MM, Poro, Alpha(6)  !Porous parameters to obtain undrained c_ij parameters
    REAL    :: rho_F, rho_S, nu, Kappa(3), Tor(3)
    REAL    :: Rho1, Rho2, Beta1, Beta2, solution2(4)
    REAL, POINTER :: zone_minh(:), zone_maxh(:), zone_deltah(:), zone_deltap(:)
    COMPLEX :: solution(3)
    INTEGER :: nDOF,TotDOF, PoroFlux
    !
    INTEGER :: iErr,iPoly,iVrtx
    INTEGER :: nLocPolyElem(0:100), TempInt(MESH%nSideMax)
    REAL    :: min_h, max_h, deltah, deltap
    INTEGER :: LocPoly
    REAL    :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType), Tmp(3,3),D(2,2,2,3)
    REAL    :: xi, eta, zeta, xP, yP, zP
    INTEGER :: ngll, k, Fix1(2), Fix2(2)
    CHARACTER(LEN=200) :: Filename
    !
    ! variables for fault output
    REAL, POINTER :: S_inc(:)
    REAL, ALLOCATABLE :: chi_vector(:), tau_vector(:)
    REAL    :: S_tmp, chi, tau, phi1, phi2
    integer :: hasDR
    INTEGER :: iNeighbor, iNeighborSide, NeigBndGP, l, iFault, iDegFr, iP, iBndGP, iPlusElem
    !
    !-------------------------------------------------------------------------!
    INTENT(IN)                :: MaterialVal, EQN, IO
    INTENT(INOUT)             :: DISC, BND, MESH
    !-------------------------------------------------------------------------!
    !                                                                         !
    !
    ! The unit tetrahedron has the following 4 local vertices:                !
    !                                                                         !
    ! 1 = (0,0,0)                                                             !
    ! 2 = (1,0,0)                                                             !
    ! 3 = (0,1,0)                                                             !
    ! 4 = (0,0,1)                                                             !
    !                                                                         !
    VertexSide_Tet(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of side I        !
    VertexSide_Tet(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of side II       !
    VertexSide_Tet(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of side III      !
    VertexSide_Tet(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of side IV       !
    ! The unit hexahedron has the following 6 local vertices:                 !
    !                                                                         !
    ! 1 = (0,0,0)                                                             !
    ! 2 = (1,0,0)                                                             !
    ! 3 = (0,1,0)                                                             !
    ! 4 = (1,1,0)                                                             !
    ! 5 = (0,0,1)                                                             !
    ! 6 = (1,0,1)                                                             !
    ! 7 = (0,1,1)                                                             !
    ! 8 = (1,1,1)                                                             !
    !                                                                         !
    !         7           8         1          1,2,6,5
    !         *-----------*         2          2,4,8,6
    !        /|          /|         3          4,3,7,8
    !       / |         / |         4          3,1,5,7
    !     5*-----------*6 |         5          2,1,3,4
    !      |  |        |  |         6          5,6,8,7
    !      |  |        |  |
    !      | 3*--------|--*4        Each face is characterized by the sum of it local vertex numbers,
    !      | /         | /          i.e. side1=14, side2=20, side3=22, side4=16, side5=10, side6=26
    !      |/          |/
    !     1*-----------*2
    !
    VertexSide_Hex(1,:) =  (/ 1,2,6,5 /)   ! Local hex. vertices of side I        !
    VertexSide_Hex(2,:) =  (/ 2,4,8,6 /)   ! Local hex. vertices of side II       !
    VertexSide_Hex(3,:) =  (/ 4,3,7,8 /)   ! Local hex. vertices of side III      !
    VertexSide_Hex(4,:) =  (/ 3,1,5,7 /)   ! Local hex. vertices of side IV       !
    VertexSide_Hex(5,:) =  (/ 2,1,3,4 /)   ! Local hex. vertices of side V        !
    VertexSide_Hex(6,:) =  (/ 5,6,8,7 /)   ! Local hex. vertices of side VI       !
    !
    ALLOCATE( DISC%Galerkin%geoNormals(EQN%Dimension,MESH%nSideMax,MESH%nElem),  &
              DISC%Galerkin%geoTangent1(EQN%Dimension,MESH%nSideMax,MESH%nElem),  &
              DISC%Galerkin%geoTangent2(EQN%Dimension,MESH%nSideMax,MESH%nElem),  &
              DISC%Galerkin%geoSurfaces(MESH%nSideMax,MESH%nElem),  &
              STAT=allocstat )
    IF (allocStat .NE. 0) THEN
       logError(*) 'Interface SeisSol: could not allocate all variables!'
       STOP
    END IF

    ! Calculating boundary surfaces (3D)

    DO iElem=1,MESH%nElem
       DO iSide=1,MESH%LocalElemType(iElem)
          SELECT CASE(MESH%LocalElemType(iElem))
          CASE(4)
              ! Boundary side vector pointing in chi-direction
              sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,2),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
              ! Boundary side vector pointing in tau-direction
              sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,3),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
              ! Normal vector computed by cross product
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  sidevec(:,1).x.sidevec(:,2)
              ! Triangle surface = 0.5 * cross_product
              DISC%Galerkin%geoSurfaces(iSide,iElem)  =  0.5*SQRT(                                     &
                                                         DISC%Galerkin%geoNormals(1,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(2,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(3,iSide,iElem)**2    )
              ! Normalize normal vector to length 1
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  0.5*DISC%Galerkin%geoNormals(:,iSide,iElem) / &
                                                         DISC%Galerkin%geoSurfaces(iSide,iElem)
              !
              ! Compute MinDistBarySide :
              ! 1. Compute vector connecting barycenter of tetrahedron and local point 1 of the side
              BaryVec(:,iSide) = MESH%ELEM%xyBary(:,iElem) - MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Tet(iSide,1),iElem))
              ! 2. Compute scalar product of previous vector and normal vector (already normalized to 1)
              !    and take the absolute value.
              Dist(iSide) = ABS(DOT_PRODUCT(BaryVec(:,iSide),DISC%Galerkin%geoNormals(:,iSide,iElem)))
          CASE(6)
              ! Boundary side vector pointing in chi-direction
              sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,2),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,1),iElem))
              ! Boundary side vector pointing in tau-direction
              sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,3),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,1),iElem))
              ! Normal vector computed by cross product
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  sidevec(:,1).x.sidevec(:,2)
              ! Triangle's surface = 0.5 * cross_product
              DISC%Galerkin%geoSurfaces(iSide,iElem)  =  0.5*SQRT(                                     &
                                                         DISC%Galerkin%geoNormals(1,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(2,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(3,iSide,iElem)**2    )

              ! Boundary side vector pointing in chi-direction
              sidevec(:,1) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,4),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,3),iElem))
              ! Boundary side vector pointing in tau-direction
              sidevec(:,2) = MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,2),iElem)) - &
                             MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,3),iElem))
              ! Normal vector computed by cross product
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  sidevec(:,1).x.sidevec(:,2)
              ! Second triangle's surface = 0.5 * cross_product
              DISC%Galerkin%geoSurfaces(iSide,iElem)  =  DISC%Galerkin%geoSurfaces(iSide,iElem) +      &
                                                         0.5*SQRT(                                     &
                                                         DISC%Galerkin%geoNormals(1,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(2,iSide,iElem)**2 +  &
                                                         DISC%Galerkin%geoNormals(3,iSide,iElem)**2    )

              ! Normalize normal vector to length 1
              DISC%Galerkin%geoNormals(:,iSide,iElem) =  DISC%Galerkin%geoNormals(:,iSide,iElem) / &
                                                         SQRT(SUM(DISC%Galerkin%geoNormals(:,iSide,iElem)**2))
              !
              ! Compute MinDistBarySide :
              ! 1. Compute vector connecting barycenter of tetrahedron and local point 1 of the side
              BaryVec(:,iSide) = MESH%ELEM%xyBary(:,iElem) - MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(VertexSide_Hex(iSide,1),iElem))
              ! 2. Compute scalar product of previous vector and normal vector (already normalized to 1)
              !    and take the absolute value.
              Dist(iSide) = ABS(DOT_PRODUCT(BaryVec(:,iSide),DISC%Galerkin%geoNormals(:,iSide,iElem)))

              !
              CONTINUE
              !
          END SELECT
          ! Compute vector inside the triangle's plane for the rotation matrix
          DISC%Galerkin%geoTangent1(:,iSide,iElem) = sidevec(:,1)
          ! Normalize to 1
          Length = SQRT(DISC%Galerkin%geoTangent1(1,iSide,iElem)**2 + &
                        DISC%Galerkin%geoTangent1(2,iSide,iElem)**2 + &
                        DISC%Galerkin%geoTangent1(3,iSide,iElem)**2   )
          DISC%Galerkin%geoTangent1(:,iSide,iElem) = DISC%Galerkin%geoTangent1(:,iSide,iElem)/Length
          ! Compute second vector in the plane, orthogonal to the normal and tangent 1 vectors
          ! using the crossproduct
          DISC%Galerkin%geoTangent2(:,iSide,iElem) = DISC%Galerkin%geoNormals( :,iSide,iElem) .x. &
                                                     DISC%Galerkin%geoTangent1(:,iSide,iElem)
          ! Check: Length must be 1.
          Length = DISC%Galerkin%geoTangent2(1,iSide,iElem)**2 + &
                   DISC%Galerkin%geoTangent2(2,iSide,iElem)**2 + &
                   DISC%Galerkin%geoTangent2(3,iSide,iElem)**2
      ENDDO
      !
      SELECT CASE(MESH%LocalElemType(iElem))
      CASE(6)
          ! Minimum distance of Barycenter to a face
          MESH%ELEM%MinDistBarySide(iElem) = MINVAL(Dist(:))
      CASE(4)
          ! Insphere radius of tetrahedron (is larger the the previous minimum distance !!!)
          MESH%ELEM%MinDistBarySide(iElem) = 3*MESH%ELEM%Volume(iElem)/SUM(DISC%Galerkin%geoSurfaces(1:4,iElem))
      END SELECT
      !
    ENDDO
    !
    ! Report mesh quality
    !
    minl = MINLOC(MESH%ELEM%Volume(:))
    minv = MINVAL(MESH%ELEM%Volume(:))

    logInfo(*) 'Smallest volume found in tetraedron number : ', minl(1)
    logInfo(*) 'Smallest volume is                         : ', minv

    minl = MINLOC(MESH%ELEM%MinDistBarySide(:))
    minv = MINVAL(MESH%ELEM%MinDistBarySide(:))

    logInfo(*) 'Smallest insphere found in tetraedron number : ', minl(1)
    logInfo(*) 'Smallest insphere radius is                  : ', minv

    IF(minv.LE.1e-15) THEN
        logError(*) 'Mesh contains a singular tetrahedron with radius ', minv
        logError(*) 'Element number and position : ', minl(1), MESH%ELEM%xyBary(:,minl(1))
        STOP
    ENDIF
    DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .FALSE.
    DISC%DynRup%DynRup_out_elementwise%nDR_pick       = 0
    !
    !
    !
    ! Initialize fault rupture output
    ! only in case Dynamic rupture is turned on, and for + elements assigned to the fault
    IF(EQN%DR.EQ.1 .AND. DISC%DynRup%DR_output) THEN
        ! Case 3
        ! output at certain positions specified in the *.dyn file
        IF(DISC%DynRup%OutputPointType.EQ.3) THEN
            !
            DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_atPickpoint%nDR_pick       = DISC%DynRup%DynRup_out_atPickpoint%nOutPoints
            !
            ! test if fault pickpoints are on the fault (within a tolerance) and find corresponding "+"-element (iElem)

#ifdef HDF
            CALL ini_fault_receiver_hdf(EQN, MESH, DISC, IO, MPI)
!#else
#endif
            CALL ini_fault_receiver(EQN,MESH,BND,DISC,IO,MPI)

        ! Case 4
        ! for full fault output without pickpoints
        ELSEIF(DISC%DynRup%OutputPointType.EQ.4) THEN
            !
            DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_elementwise%nDR_pick       = 0
            CALL ini_fault_subsampled(EQN,MESH,BND,DISC,IO,MPI)
        ! Case 5
        ! for full fault output and pickpoints
        ELSEIF(DISC%DynRup%OutputPointType.EQ.5) THEN
            !
            DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_atPickpoint%nDR_pick       = DISC%DynRup%DynRup_out_atPickpoint%nOutPoints
            !
            ! test if fault pickpoints are on the fault (within a tolerance) and find corresponding "+"-element (iElem)
            CALL ini_fault_receiver(EQN,MESH,BND,DISC,IO,MPI)
            !
            !
            DISC%DynRup%DynRup_out_elementwise%DR_pick_output = .TRUE.
            DISC%DynRup%DynRup_out_elementwise%nDR_pick       = 0
            CALL ini_fault_subsampled(EQN,MESH,BND,DISC,IO,MPI)
        ENDIF ! DISC%DynRup%OutputPointType
    ENDIF ! end initialize fault output
    !
    !
    !
    !
    ! Allocate rest of MPI communication structure
    logInfo(*) 'Allocation of remaining MPI communication structure '
    logInfo(*) '  General info: ', BND%NoMPIDomains,DISC%Galerkin%nDegFr,EQN%nVar
    DO iDomain = 1, BND%NoMPIDomains
        logInfo(*) 'Bnd elements for domain ', iDomain, ' : ',  BND%ObjMPI(iDomain)%nElem
        IF(DISC%Galerkin%DGMethod.EQ.3) THEN
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDOF(DISC%Galerkin%nDegFrST,EQN%nVarTotal,BND%ObjMPI(iDomain)%nElem) )
        ELSE
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal,BND%ObjMPI(iDomain)%nElem) )
            IF (EQN%DR.EQ.1) THEN
                ALLOCATE(BND%ObjMPI(iDomain)%MPI_DR_dgvar(DISC%Galerkin%nDegFrRec,EQN%nVarTotal,BND%ObjMPI(iDomain)%nFault_MPI))
            ENDIF
        ENDIF
        ALLOCATE( BND%ObjMPI(iDomain)%NeighborBackground(EQN%nBackgroundVar,BND%ObjMPI(iDomain)%nElem)     )
        BND%ObjMPI(iDomain)%Init = .FALSE.
        ALLOCATE( BND%ObjMPI(iDomain)%AStar_Sp( BND%ObjMPI(iDomain)%nElem ) )
        ALLOCATE( BND%ObjMPI(iDomain)%BStar_Sp( BND%ObjMPI(iDomain)%nElem ) )
        ALLOCATE( BND%ObjMPI(iDomain)%CStar_Sp( BND%ObjMPI(iDomain)%nElem ) )
        ALLOCATE( BND%ObjMPI(iDomain)%EStar_Sp( BND%ObjMPI(iDomain)%nElem ) )
        ALLOCATE( BND%ObjMPI(iDomain)%FLStar_Sp(BND%ObjMPI(iDomain)%nElem,MESH%nSideMax))
        ALLOCATE( BND%ObjMPI(iDomain)%FRStar_Sp(BND%ObjMPI(iDomain)%nElem,MESH%nSideMax))
        IF(DISC%Galerkin%CKMethod.EQ.1) THEN
            ALLOCATE( BND%ObjMPI(iDomain)%InvSystemMatrix( BND%ObjMPI(iDomain)%nElem ) )
        ENDIF
        IF(DISC%Galerkin%DGMethod.EQ.3) THEN
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDuDt(DISC%Galerkin%nDegFr,EQN%nVar+EQN%nAneFuncperMech,BND%ObjMPI(iDomain)%nElem) )
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborTime( BND%ObjMPI(iDomain)%nElem) )
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborDt(   BND%ObjMPI(iDomain)%nElem) )
            ALLOCATE( BND%ObjMPI(iDomain)%NeighborUpdate(BND%ObjMPI(iDomain)%nElem))
            BND%ObjMPI(iDomain)%NeighborDuDt(:,:,:) = 0.
            BND%ObjMPI(iDomain)%NeighborTime(:)     = -1e10
            BND%ObjMPI(iDomain)%NeighborDt(:)       = -2e10
            BND%ObjMPI(iDomain)%NeighborUpdate(:)   = -1
        ENDIF
    ENDDO ! iDomain
    !
    IF(DISC%Galerkin%ZoneOrderFlag.EQ.1) THEN
        ALLOCATE( zone_minh(MESH%nZones), zone_maxh(MESH%nZones), zone_deltah(MESH%nZones), zone_deltap(MESH%nZones) )
        DO iZone = 1, MESH%nZones
            min_h = MINVAL( MESH%ELEM%Volume(:), MASK = MESH%ELEM%Reference(0,:).EQ.iZone )
            max_h = MAXVAL( MESH%ELEM%Volume(:), MASK = MESH%ELEM%Reference(0,:).EQ.iZone )
#ifdef PARALLEL
            CALL MPI_ALLREDUCE(min_h,zone_minh(iZone),1,MPI%MPI_AUTO_REAL,MPI_MIN,MPI%commWorld,iErr)
            CALL MPI_ALLREDUCE(max_h,zone_maxh(iZone),1,MPI%MPI_AUTO_REAL,MPI_MAX,MPI%commWorld,iErr)
#else
            zone_minh(iZone) = min_h
            zone_maxh(iZone) = max_h
#endif
            zone_deltah(iZone) = zone_maxh(iZone) - zone_minh(iZone)
            zone_deltap(iZone) = DISC%Galerkin%ZoneMaxPoly(iZone) - DISC%Galerkin%ZoneMinPoly(iZone)
        ENDDO
    ELSE
        min_h = MINVAL( MESH%ELEM%MinDistBarySide(:) )
        max_h = MAXVAL( MESH%ELEM%MinDistBarySide(:) )
#ifdef PARALLEL
        CALL MPI_ALLREDUCE(min_h,MESH%min_h,1,MPI%MPI_AUTO_REAL,MPI_MIN,MPI%commWorld,iErr)
        CALL MPI_ALLREDUCE(max_h,MESH%max_h,1,MPI%MPI_AUTO_REAL,MPI_MAX,MPI%commWorld,iErr)
#else
        MESH%min_h = min_h
        MESH%max_h = max_h
#endif
        deltah = MESH%max_h - MESH%min_h
        deltap = DISC%Galerkin%nPoly - DISC%Galerkin%nMinPoly
    ENDIF


    ! Allocate objects for local p-adaptivity (not yet done for hybrid meshes)
    IF(DISC%Galerkin%pAdaptivity.GE.1) THEN
        ALLOCATE( DISC%Galerkin%LocPoly(MESH%nElem) )
        !
        nLocPolyElem(:) = 0
        ! Distribute local polynomial degree
        !
        DO iElem = 1, MESH%nElem
            IF(DISC%Galerkin%ZoneOrderFlag.EQ.1)THEN
                !
                !  Distribute local polynomial degree locally by zone
                !  ----------------------------------------------------
                iLayer = MESH%ELEM%Reference(0,iElem)        ! Zone number is given by reference 0
                IF(DISC%Galerkin%ZoneMinPoly(iLayer).GT.-1)THEN
                     SELECT CASE(DISC%Galerkin%ZoneOrder(iLayer))
                     CASE(1)
                         LocPoly = NINT( DISC%Galerkin%ZoneMinPoly(iLayer) + zone_deltap(iLayer)* &
                                         ((MESH%ELEM%Volume(iElem)-zone_minh(iLayer))/zone_deltah(iLayer))**DISC%Galerkin%ZonePower(iLayer) )
                     CASE(-1)
                         LocPoly = NINT( DISC%Galerkin%ZoneMaxPoly(iLayer) - zone_deltap(iLayer)* &
                                         ((MESH%ELEM%Volume(iElem)-zone_minh(iLayer))/zone_deltah(iLayer))**DISC%Galerkin%ZonePower(iLayer) )
                     END SELECT
                ELSE
                   PRINT *, ' ERROR: local order must not be less or equal to zero! ', iLayer
                   STOP
                ENDIF
            ELSE
                !
                !  Distribute local polynomial degree globally by size
                !  ----------------------------------------------------
                LocPoly = NINT( DISC%Galerkin%nMinPoly + deltap*((MESH%ELEM%MinDistBarySide(iElem)-MESH%min_h)/deltah)**(1) )
            ENDIF
            !
            DISC%Galerkin%LocPoly(iElem) = LocPoly
            nLocPolyElem( LocPoly ) = nLocPolyElem( LocPoly ) + 1
            !
        ENDDO

        TotDOF = 0
        DO iPoly = DISC%Galerkin%nMinPoly, DISC%Galerkin%nPoly
            logInfo('(a,i3,a,i6)') 'Nr of elements with degree ', iPoly, ' : ', nLocPolyElem(iPoly)
            nDOF = nLocPolyElem(iPoly)*(iPoly+1)*(iPoly+2)*(iPoly+3)/6
            logInfo(*) 'Nr of DOF for this p-zone:   ', nDOF
            TotDOF = TotDOF + nDOF
            IF(nLocPolyElem(iPoly).GT.0) THEN
                WRITE(FileName,'(a10,i2.2,a4)') 'MeshZone-p',iPoly,'.dat'
                OPEN(UNIT=999,FILE=TRIM(FileName),RECL=300)
                WRITE(999,*) ' TITLE = "Mesh Subzone for p-Adaptivity" '
                WRITE(999,*) ' VARIABLES = "x" "y" "z" "N" '
                WRITE(999,*) ' ZONE N=  ', nLocPolyElem(iPoly)*4, &
                                  ' E=  ', nLocPolyElem(iPoly),   &
                                  ' F=FEPOINT  ET=TETRAHEDRON '
                DO iElem = 1, MESH%nElem
                   IF(DISC%Galerkin%LocPoly(iElem).EQ.iPoly) THEN
                      DO iVrtx = 1, MESH%GlobalElemType
                          WRITE(999,*) MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(iVrtx,iElem)), REAL(iPoly)
                      ENDDO
                   ENDIF
                ENDDO
                i = 0
                DO iElem = 1, MESH%nElem
                   IF(DISC%Galerkin%LocPoly(iElem).EQ.iPoly) THEN
                      DO iVrtx = 1, MESH%GlobalElemType
                        i = i + 1
                        TempInt(iVrtx) = i
                      ENDDO
                      WRITE(999,*) TempInt(:)
                   ENDIF
                ENDDO
                CLOSE(999)
            ENDIF
        ENDDO
        logInfo(*) '======================================================= '
        logInfo(*) 'Total number of degrees of freedom in domain:   ', TotDOF
        logInfo(*) '======================================================= '
    ELSE
        ALLOCATE( DISC%Galerkin%LocPoly(1) )
    ENDIF

    IF(DISC%Galerkin%ZoneOrderFlag.EQ.1) THEN
        DEALLOCATE( zone_minh, zone_maxh, zone_deltah, zone_deltap )
    ENDIF

    ALLOCATE( DISC%Galerkin%WaveSpeed(MESH%nElem,MESH%nSideMax,EQN%nNonZeroEV) )
    !
    DISC%Galerkin%WaveSpeed(:,:,:) = 0.
    !
    IF(EQN%Anisotropy.EQ.0.AND.EQN%Poroelasticity.EQ.0) THEN
      ALLOCATE( DISC%Galerkin%MaxWaveSpeed(MESH%nElem,MESH%nSideMax) )
      DO j=1,MESH%nSideMax
        DISC%Galerkin%WaveSpeed(:,j,1)=SQRT((MaterialVal(:,3)+2.*MaterialVal(:,2))/(MaterialVal(:,1)))
        DISC%Galerkin%WaveSpeed(:,j,2)=SQRT((MaterialVal(:,2))/(MaterialVal(:,1)))
        DISC%Galerkin%WaveSpeed(:,j,3)=SQRT((MaterialVal(:,2))/(MaterialVal(:,1)))
        DISC%Galerkin%MaxWaveSpeed(:,j)=SQRT((MaterialVal(:,3)+2.*MaterialVal(:,2))/(MaterialVal(:,1)))
      ENDDO
      !
    ELSE

      ALLOCATE( DISC%Galerkin%MaxWaveSpeed(MESH%nElem,MESH%nSideMax) )

      DISC%Galerkin%WaveSpeed(:,:,:) = 0.

      DO iElem = 1, MESH%nElem

       SELECT CASE(EQN%LocPoroelastic(iElem))
       CASE(0)                  ! Triclinic system
         !
         DO iSide = 1, MESH%LocalElemType(iElem)
            !
            T(:,:)          = 0.
            TT(:,:)         = 0.
            !
            nx = DISC%Galerkin%geoNormals(1,iSide,iElem)
            ny = DISC%Galerkin%geoNormals(2,iSide,iElem)
            nz = DISC%Galerkin%geoNormals(3,iSide,iElem)
            sx = DISC%Galerkin%geoTangent1(1,iSide,iElem)
            sy = DISC%Galerkin%geoTangent1(2,iSide,iElem)
            sz = DISC%Galerkin%geoTangent1(3,iSide,iElem)
            tx = DISC%Galerkin%geoTangent2(1,iSide,iElem)
            ty = DISC%Galerkin%geoTangent2(2,iSide,iElem)
            tz = DISC%Galerkin%geoTangent2(3,iSide,iElem)
            !
            ! Rotation matrices for the transformation of the Voigt matrix
            ! C' = T*C*T^T
            !
            ! Transformation matrix TO the rotated system
            T(1,1) = nx**2
            T(1,2) = ny**2
            T(1,3) = nz**2
            T(1,4) = 2*nz*ny
            T(1,5) = 2*nz*nx
            T(1,6) = 2*ny*nx
            T(2,1) = sx**2
            T(2,2) = sy**2
            T(2,3) = sz**2
            T(2,4) = 2*sz*sy
            T(2,5) = 2*sz*sx
            T(2,6) = 2*sy*sx
            T(3,1) = tx**2
            T(3,2) = ty**2
            T(3,3) = tz**2
            T(3,4) = 2*tz*ty
            T(3,5) = 2*tz*tx
            T(3,6) = 2*ty*tx
            T(4,1) = sx*tx
            T(4,2) = sy*ty
            T(4,3) = sz*tz
            T(4,4) = sz*ty+sy*tz
            T(4,5) = sz*tx+sx*tz
            T(4,6) = sy*tx+sx*ty
            T(5,1) = nx*tx
            T(5,2) = ny*ty
            T(5,3) = nz*tz
            T(5,4) = nz*ty+ny*tz
            T(5,5) = nz*tx+nx*tz
            T(5,6) = ny*tx+nx*ty
            T(6,1) = nx*sx
            T(6,2) = ny*sy
            T(6,3) = nz*sz
            T(6,4) = nz*sy+ny*sz
            T(6,5) = nz*sx+nx*sz
            T(6,6) = ny*sx+nx*sy
            !
            ! Transpose of Transformation matrix
            !
            DO i = 1, 6
              DO j = 1, 6
                TT(i,j) = T(j,i)
              ENDDO
            ENDDO
            !
            ! Material constants in the global xyz system
            !
            rho    = MaterialVal(iElem, 1)
            c(1,1) = MaterialVal(iElem, 2)
            c(1,2) = MaterialVal(iElem, 3)
            c(1,3) = MaterialVal(iElem, 4)
            c(1,4) = MaterialVal(iElem, 5)
            c(1,5) = MaterialVal(iElem, 6)
            c(1,6) = MaterialVal(iElem, 7)
            c(2,2) = MaterialVal(iElem, 8)
            c(2,3) = MaterialVal(iElem, 9)
            c(2,4) = MaterialVal(iElem,10)
            c(2,5) = MaterialVal(iElem,11)
            c(2,6) = MaterialVal(iElem,12)
            c(3,3) = MaterialVal(iElem,13)
            c(3,4) = MaterialVal(iElem,14)
            c(3,5) = MaterialVal(iElem,15)
            c(3,6) = MaterialVal(iElem,16)
            c(4,4) = MaterialVal(iElem,17)
            c(4,5) = MaterialVal(iElem,18)
            c(4,6) = MaterialVal(iElem,19)
            c(5,5) = MaterialVal(iElem,20)
            c(5,6) = MaterialVal(iElem,21)
            c(6,6) = MaterialVal(iElem,22)
            !
            ! Initialize lower half of Voigt matrix using the symmetry property
            !
            DO i = 1, 6
              DO j = 1, i-1
                c(i,j) = c(j,i)
              ENDDO
            ENDDO
            !
            ! Initialize and rotate Voigt matrix to get local material properties
            !
            Voigt_rot(:,:) = MATMUL( T(:,:), MATMUL(c(:,:),TT(:,:)) )
            !
            c(:,:) = Voigt_rot(:,:)
            !
            ! Eigenvalues of Jacobian in normal direction
            !
            coefficients(1) = rho**3
            coefficients(2) = (-rho**2*c(1,1)-c(5,5)*rho**2-rho**2*c(6,6))
            coefficients(3) = (-c(1,5)**2*rho+rho*c(6,6)*c(1,1)+c(5,5)*c(6,6)*rho-c(5,6)**2*rho- &
                              c(1,6)**2*rho+c(5,5)*rho*c(1,1))
            coefficients(4) = -c(5,5)*c(6,6)*c(1,1)+c(1,6)**2*c(5,5)-2*c(5,6)*c(1,5)*c(1,6)+         &
                               c(5,6)**2*c(1,1)+c(1,5)**2*c(6,6)
            !
            ! Normalize all coeffs w.r.t. c_11^2 (is always different from 0) to get more accurate root values (verify)
            coefficients(:) = coefficients(:)/(c(1,1)**2)
            !
           CALL ZerosPolyO3(solution,coefficients)
            !
            DO i = 1, 3
              Re_solution(i) = REAL(solution(i))
              Im_solution(i) = AIMAG(solution(i))
              IF( (Re_solution(i).LT.0.) ) THEN
                logError(*) 'Imaginary wave speed. System is not hyperbolic. '
                logError(*) ' Probably your material does not exist.'
                STOP
              ENDIF
              IF( ABS(Im_solution(i))/Re_solution(i).GT.1.e-6 ) THEN   !Roundoff error can produce larger imaginary parts than expected!
                logWarning(*) ' Imaginary wave speed. System is not hyperbolic. Wave i=', i
                logWarning(*) ' Probably your material does not exist. Im(lambda) = ', Im_solution(i)
              ENDIF
            ENDDO
            !
            !Re_Solution(:) = Re_solution(:) + Im_solution(:) !Adding the complex residues (verify)
            !
            amax         = MAXVAL( SQRT(Re_solution(:)) )
            !
            DISC%Galerkin%maxWaveSpeed(iElem,iSide) = amax
            !
            ! Descending sorting of the wavespeeds
            DO j=2,3
              a=Re_solution(j)
              DO i=j-1,1,-1
                IF (Re_solution(i) >= a) EXIT
                Re_solution(i+1)=Re_solution(i)
              ENDDO
              Re_solution(i+1)=a
            ENDDO
            !
            DISC%Galerkin%WaveSpeed(iElem,iSide,1:3) = SQRT(Re_solution(1:3))
            !
       ENDDO


      CASE(1,2,3)    !Porous case
       !
       DO iSide = 1, MESH%LocalElemType(iElem)
            !
            T(:,:)          = 0.
            TT(:,:)         = 0.
            !
            nx = DISC%Galerkin%geoNormals(1,iSide,iElem)
            ny = DISC%Galerkin%geoNormals(2,iSide,iElem)
            nz = DISC%Galerkin%geoNormals(3,iSide,iElem)
            sx = DISC%Galerkin%geoTangent1(1,iSide,iElem)
            sy = DISC%Galerkin%geoTangent1(2,iSide,iElem)
            sz = DISC%Galerkin%geoTangent1(3,iSide,iElem)
            tx = DISC%Galerkin%geoTangent2(1,iSide,iElem)
            ty = DISC%Galerkin%geoTangent2(2,iSide,iElem)
            tz = DISC%Galerkin%geoTangent2(3,iSide,iElem)
            !
            ! Rotation matrices for the transformation of the Voigt matrix
            ! C' = T*C*T^T
            !
            ! Transformation matrix TO the rotated system
            T(1,1) = nx**2
            T(1,2) = ny**2
            T(1,3) = nz**2
            T(1,4) = 2*nz*ny
            T(1,5) = 2*nz*nx
            T(1,6) = 2*ny*nx
            T(2,1) = sx**2
            T(2,2) = sy**2
            T(2,3) = sz**2
            T(2,4) = 2*sz*sy
            T(2,5) = 2*sz*sx
            T(2,6) = 2*sy*sx
            T(3,1) = tx**2
            T(3,2) = ty**2
            T(3,3) = tz**2
            T(3,4) = 2*tz*ty
            T(3,5) = 2*tz*tx
            T(3,6) = 2*ty*tx
            T(4,1) = sx*tx
            T(4,2) = sy*ty
            T(4,3) = sz*tz
            T(4,4) = sz*ty+sy*tz
            T(4,5) = sz*tx+sx*tz
            T(4,6) = sy*tx+sx*ty
            T(5,1) = nx*tx
            T(5,2) = ny*ty
            T(5,3) = nz*tz
            T(5,4) = nz*ty+ny*tz
            T(5,5) = nz*tx+nx*tz
            T(5,6) = ny*tx+nx*ty
            T(6,1) = nx*sx
            T(6,2) = ny*sy
            T(6,3) = nz*sz
            T(6,4) = nz*sy+ny*sz
            T(6,5) = nz*sx+nx*sz
            T(6,6) = ny*sx+nx*sy
            !
            ! Transpose of Transformation matrix
            !
            DO i = 1, 6
              DO j = 1, 6
                TT(i,j) = T(j,i)
              ENDDO
            ENDDO
            !
            ! Material constants in the global xyz system
            !
            rho_S  = MaterialVal(iElem, 1)
            c(1,1) = MaterialVal(iElem, 2)
            c(1,2) = MaterialVal(iElem, 3)
            c(1,3) = MaterialVal(iElem, 4)
            c(1,4) = MaterialVal(iElem, 5)
            c(1,5) = MaterialVal(iElem, 6)
            c(1,6) = MaterialVal(iElem, 7)
            c(2,2) = MaterialVal(iElem, 8)
            c(2,3) = MaterialVal(iElem, 9)
            c(2,4) = MaterialVal(iElem,10)
            c(2,5) = MaterialVal(iElem,11)
            c(2,6) = MaterialVal(iElem,12)
            c(3,3) = MaterialVal(iElem,13)
            c(3,4) = MaterialVal(iElem,14)
            c(3,5) = MaterialVal(iElem,15)
            c(3,6) = MaterialVal(iElem,16)
            c(4,4) = MaterialVal(iElem,17)
            c(4,5) = MaterialVal(iElem,18)
            c(4,6) = MaterialVal(iElem,19)
            c(5,5) = MaterialVal(iElem,20)
            c(5,6) = MaterialVal(iElem,21)
            c(6,6) = MaterialVal(iElem,22)
            rho_F  = MaterialVal(iElem,23)
            K_F    = MaterialVal(iElem,24)
            nu     = MaterialVal(iElem,25)
            K_S    = MaterialVal(iElem,26)
            Poro   = MaterialVal(iElem,27)
            ! Permeabilities and tortuosities are rotated assuming that they represent an ellipsoid in the reference system!
            ! The ellipsoid is aligned with the mesh reference system (no rotation in .par file)
            Kappa(1) = SQRT((MaterialVal(iElem,28)*nx)**2+(MaterialVal(iElem,29)*ny)**2+(MaterialVal(iElem,30)*nz)**2)
            Tor(1)   = SQRT((MaterialVal(iElem,31)*nx)**2+(MaterialVal(iElem,32)*ny)**2+(MaterialVal(iElem,33)*nz)**2)
            !
            ! Initialize lower half of Voigt matrix using the symmetry property
            !
            DO i = 1, 6
              DO j = 1, i-1
                c(i,j) = c(j,i)
              ENDDO
            ENDDO
            !
            ! Initialize and rotate Voigt matrix to get local material properties
            !
            Voigt_rot(:,:) = MATMUL( T(:,:), MATMUL(c(:,:),TT(:,:)) )
            !
            c(:,:) = Voigt_rot(:,:)
            !
            ! Derived constants for the jacobians
            rho    = rho_S * (1 - Poro) + Poro * rho_F
            K_Mean = 1./9.*(c(1,1)+c(2,2)+c(3,3)+2*(c(1,2)+c(1,3)+c(2,3)))
            MM     = K_S / ((1 - K_Mean/K_S) - Poro*(1 - K_S/K_F))
            Alpha(1) = 1 - (c(1,1)+c(1,2)+c(1,3)) / (3.*K_S)
            Alpha(2) = 1 - (c(1,2)+c(2,2)+c(2,3)) / (3.*K_S)
            Alpha(3) = 1 - (c(1,3)+c(2,3)+c(3,3)) / (3.*K_S)
            Alpha(4) = - (c(1,4)+c(2,4)+c(3,4)) / (3.*K_S)
            Alpha(5) = - (c(1,5)+c(2,5)+c(3,5)) / (3.*K_S)
            Alpha(6) = - (c(1,6)+c(2,6)+c(3,6)) / (3.*K_S)
            Rho1     = rho - (rho_F**2 / (rho_F * Tor(1) / Poro))
            Rho2     = rho_F - (rho_F * Tor(1) / Poro) * rho /rho_F
            Beta1    = rho_F / (rho_F * Tor(1) / Poro)
            Beta2    = rho / rho_F
            !
            ! Computation of undrained c(i,j) coeffs
            DO i=1,6
              DO j=1,6
                c(i,j) = c(i,j) + MM * Alpha(i)*Alpha(j)
              ENDDO
            ENDDO
            !
            ! Eigenvalues of Jacobian in normal direction
            !
            IF((c(6,6)/c(1,1)).GE.1e-10) THEN

              coefficients2(1) = Rho1**3*Rho2
              coefficients2(2) = Rho1**2*(-c(6,6)*Rho2+Alpha(1)*MM*Rho2*Beta1+Beta2*Rho1*MM-c(5,5)*Rho2-c(1,1)*Rho2-Rho1*Alpha(1)*MM)
              coefficients2(3) = Rho1*(Beta2*Rho1*Alpha(5)**2*MM**2-Beta2*Rho1*c(5,5)*MM+Beta2*Rho1*Alpha(6)**2*MM**2- &
                                 Rho1*c(1,5)*Alpha(5)*MM-Beta2*c(1,1)*Rho1*MM+c(1,1)*Rho1*MM*Beta1+c(1,1)*c(5,5)*Rho2+c(1,1)*c(6,6)*Rho2 &
                                 +Beta2*Alpha(1)**2*MM**2*Rho1-Alpha(1)**2*MM**2*Rho1*Beta1+Rho1*c(6,6)*Alpha(1)*MM-Rho1* &
                                 c(1,6)*MM*Alpha(6)-c(1,5)**2*Rho2+c(1,5)*Alpha(5)*MM*Rho2*Beta1+c(6,6)*c(5,5)*Rho2-c(1,6)**2*Rho2- &
                                 Alpha(1)*MM*c(5,5)*Rho2*Beta1+Rho1*c(5,5)*Alpha(1)*MM-c(5,6)**2*Rho2-Beta2*Rho1*c(6,6)*MM+ &
                                 c(1,6)*Alpha(6)*MM*Rho2*Beta1-Alpha(1)*MM*c(6,6)*Rho2*Beta1)
              coefficients2(4) = c(1,1)*c(5,6)**2*Rho2+c(1,6)**2*c(5,5)*Rho2+c(6,6)*c(1,5)**2*Rho2+2*Beta2*c(1,5)*Rho1*Alpha(5)* &
                                 MM**2*Alpha(1)+2*Beta2*Rho1*c(5,6)*Alpha(5)*MM**2*Alpha(6)-2*c(1,5)*Rho1*Alpha(5)*MM**2  &
                                 *Alpha(1)*Beta1+Alpha(1)*MM*c(6,6)*c(5,5)*Rho2*Beta1+Beta2*Rho1*c(6,6)*c(5,5)*MM+c(1,6)*c(5,6)* &
                                 Alpha(5)*MM*Rho2*Beta1-c(1,1)*c(6,6)*Rho1*MM*Beta1-Beta2*Alpha(1)**2*MM**2*Rho1*c(5,5)+ &
                                 Alpha(1)**2*MM**2*Rho1*c(5,5)*Beta1+Rho1*c(1,6)*MM*Alpha(6)*c(5,5)-c(1,6)*Alpha(6)*MM*c(5,5)* &
                                 Rho2*Beta1+Beta2*c(1,1)*c(6,6)*Rho1*MM+Alpha(1)**2*MM**2*c(6,6)*Rho1*Beta1-2*c(1,6)*Alpha(6) &
                                 *MM**2*Rho1*Alpha(1)*Beta1+c(1,1)*Rho1*Alpha(5)**2*MM**2*Beta1-c(1,5)*c(6,6)*Alpha(5)*MM*Rho2* &
                                 Beta1+2*Beta2*c(1,6)*Alpha(6)*MM**2*Rho1*Alpha(1)-c(1,1)*Rho1*c(5,5)*MM*Beta1-Beta2*Rho1* &
                                 c(6,6)*Alpha(5)**2*MM**2+c(1,1)*Alpha(6)**2*MM**2*Rho1*Beta1-Beta2*Alpha(1)**2*MM**2*c(6,6)*Rho1 &
                                 +c(1,5)*Alpha(6)*MM*c(5,6)*Rho2*Beta1-Rho1*c(5,6)*c(1,5)*Alpha(6)*MM+Rho1*c(6,6)*c(1,5)*Alpha(5)*MM- &
                                 Rho1*c(6,6)*c(5,5)*Alpha(1)*MM-Rho1*c(1,6)*MM*c(5,6)*Alpha(5)-Beta2*Rho1*Alpha(6)**2*MM**2*c(5,5)- &
                                 Alpha(1)*MM*c(5,6)**2*Rho2*Beta1-Beta2*c(1,1)*Rho1*Alpha(5)**2*MM**2+Beta2*c(1,1)*Rho1*c(5,5)*MM- &
                                 Beta2*c(1,1)*Alpha(6)**2*MM**2*Rho1-Beta2*c(1,6)**2*Rho1*MM-Beta2*Rho1*c(1,5)**2*MM+Rho1 &
                                 *c(1,5)**2*MM*Beta1-Beta2*Rho1*c(5,6)**2*MM+c(1,6)**2*Rho1*MM*Beta1-c(1,1)*c(6,6)*c(5,5)*Rho2+ &
                                 Rho1*c(5,6)**2*Alpha(1)*MM-2*c(1,6)*c(5,6)*c(1,5)*Rho2
              coefficients2(5) = (Beta2-Beta1)*(-Alpha(1)**2*MM*c(5,6)**2-c(1,6)**2*Alpha(5)**2*MM+c(6,6)*c(1,5)**2+c(1,1)*c(5,6)**2+ &
                                 c(1,6)**2*c(5,5)-Alpha(6)**2*MM*c(1,5)**2-2*c(1,1)*c(5,6)*Alpha(5)*MM*Alpha(6)-2*c(1,5)*c(6,6)*Alpha(5)*MM* &
                                 Alpha(1)+2*c(1,5)*Alpha(6)*MM*c(5,6)*Alpha(1)-2*c(1,6)*c(5,6)*c(1,5)+c(1,1)*Alpha(6)**2*MM*c(5,5)+2*c(1,6)*c(5,6)* &
                                 Alpha(5)*MM*Alpha(1)+2*c(1,6)*Alpha(6)*MM*c(1,5)*Alpha(5)-2*c(1,6)*Alpha(6)*MM*c(5,5)*Alpha(1)+c(1,1)*c(6,6)* &
                                 Alpha(5)**2*MM+Alpha(1)**2*MM*c(6,6)*c(5,5)-c(1,1)*c(6,6)*c(5,5))*MM
              !
              CALL ZerosPolyO4(solution2,coefficients2)
              !
              DO i = 1, 4
                IF( (solution2(i).LT.0.AND.ABS(solution2(i)).GE.MAXVAL(solution2)*1.0e-6) ) THEN
                  logError(*) 'Imaginary wave speed. System is not hyperbolic. '
                  logError(*) 'Probably your material does not exist.'
                  STOP
                ENDIF
              ENDDO
              !
              solution2 = SQRT(ABS(solution2))
              !
              amax         = MAXVAL( solution2(:) )
              !
              DISC%Galerkin%maxWaveSpeed(iElem,iSide) = amax
              !
              ! Descending sorting of the wavespeeds
              DO j=2,4
                a=solution2(j)
                DO i=j-1,1,-1
                  IF (solution2(i) >= a) EXIT
                  solution2(i+1)=solution2(i)
                ENDDO
                solution2(i+1)=a
              ENDDO
              !
              DISC%Galerkin%WaveSpeed(iElem,iSide,:) = solution2(:)
              !
            ELSE
              !Isotropic-poroacoustic: We assume K = c(1,1) and Alpha(1)=Alpha(2)=Alpha(3)
              !Solving quadratic polynomial
              coefficients2(1)=Rho2*Rho1
              coefficients2(2)=(-Rho2*c(1,1)-MM*Rho1*Alpha(1)+MM*Rho1*Beta2+Rho2*Alpha(1)*MM*Beta1)
              coefficients2(3)=-MM*Beta2*c(1,1)+Beta2*Alpha(1)**2*MM**2+MM*c(1,1)*Beta1-Alpha(1)**2*MM**2*Beta1
              !
              Re_solution(1)=SQRT(ABS((-coefficients2(2)+SQRT(coefficients2(2)**2-4*coefficients2(1)*coefficients2(3)))/(2*coefficients2(1))))
              Re_solution(2)=SQRT(ABS((-coefficients2(2)-SQRT(coefficients2(2)**2-4*coefficients2(1)*coefficients2(3)))/(2*coefficients2(1))))
              !
              DISC%Galerkin%maxWaveSpeed(iElem,iSide)=MAXVAL( Re_solution(1:2) )
              !
              IF (Re_solution(2) >= Re_solution(1)) THEN
                  a = Re_solution(2)
                  Re_solution(2) = Re_solution(1)
                  Re_solution(1) = a
              ENDIF
              !
              DISC%Galerkin%WaveSpeed(iElem,iSide,1:2) = Re_solution(1:2)
              !
           ENDIF
           !
       ENDDO
       !
     END SELECT

     ENDDO

    ENDIF !Ends the anisotropy IF
    !
    CONTINUE
    !
  END SUBROUTINE BuildSpecialDGGeometry3D_new

  SUBROUTINE AnalyseGalerkin3D_us_new(time,ANALYSE,EQN,MESH,DISC,BND,SOURCE,IC,   &
                                  IO,OptionalFields,MPI                       )
    !-------------------------------------------------------------------------!

    USE DGBasis_mod
    USE COMMON_InitialField_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
#ifdef PARALLEL
  INCLUDE 'mpif.h'
#endif
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tAnalyse)          :: ANALYSE               ! Analyse structure      !
    TYPE(tEquations)        :: EQN                   ! Equation structure     !
    TYPE(tUnstructMesh)     :: MESH                  ! Mesh structure         !
    TYPE(tDiscretization)   :: DISC                  ! Discretization struct. !
    TYPE(tBoundary)         :: BND                   ! Boundary structure     !
    TYPE(tSource)           :: SOURCE                ! Source structure       !
    TYPE(tInitialCondition) :: IC                    ! Initialcondition str.  !
    TYPE(tInputOutput)      :: IO                    ! IO structure           !
    TYPE(tUnstructOptionalFields) :: OptionalFields  ! Opt. Fields            !
    TYPE(tMPI)              :: MPI                   ! MPI structure          !
    REAL                    :: time                  ! time for analysis      !
    ! Local variable declaration
    INTEGER                 :: i,j                   ! Loop variables         !
    INTEGER                 :: iVec                  ! Variable to analyse    !
    INTEGER                 :: iVar, iDegFr, iElem   ! Loop variables         !
    INTEGER                 :: iIntGP                ! Index of internal GP   !
    INTEGER                 :: iErr                  !
    INTEGER                 :: iVert
    INTEGER                 :: LInfElement(EQN%nVar)
    INTEGER                 :: LocElemType
    REAL, POINTER           :: intGaussW(:) => NULL()
    REAL                    :: xGP(DISC%Galerkin%nIntGP)
    REAL                    :: yGP(DISC%Galerkin%nIntGP)
    REAL                    :: zGP(DISC%Galerkin%nIntGP)
    REAL                    :: xi,eta,zeta
    REAL                    :: iniGP(EQN%nVar)       ! Referencesolution      !
    REAL                    :: iniGP_ane(EQN%nAneFuncperMech*EQN%nMechanisms)  ! Initial anelastic state vector at Gausspoint
    REAL                    :: stateGP(EQN%nVar)     ! State in GP            !
    REAL                    :: gradGP(EQN%nVar,EQN%Dimension) ! State in GP   !
    REAL                    :: pstateGP(EQN%nVar)    ! Primitive state in GP  !
    REAL                    :: locError(EQN%nVar)    ! local error in GP      !
    REAL                    :: L1norm(EQN%nVar)      ! L1norm                 !
    REAL                    :: L2norm(EQN%nVar)      ! L2norm                 !
    REAL                    :: Linfnorm(EQN%nVar)    ! Max. error=L_inf norm  !
    REAL                    :: MPIL1norm(EQN%nVar)   ! Total L1norm (MPI)     !
    REAL                    :: MPIL2norm(EQN%nVar)   ! Total L2norm (MPI)     !
    REAL                    :: MPILinfnorm(EQN%nVar) ! Total L_inf norm  (MPI)!
    REAL                    :: MPI_H1_divB           ! L1 norm for MPI        !
    REAL                    :: h1, MPI_h1            ! Mesh size, MPI         !
    REAL                    :: x(MESH%GlobalVrtxType)! Vertices of triangle   !
    REAL                    :: y(MESH%GlobalVrtxType)! Vertices of triangle   !
    REAL                    :: z(MESH%GlobalVrtxType)! Vertices of triangle   !
    REAL                    :: JacobiDet             ! Mapping Jacobian det.  !
    REAL                    :: H1_divB               ! Divergence of B        !
    REAL                    :: kx, ky, kz            ! wavenumbers            !
    COMPLEX                 :: IU                    ! imaginary unit         !
    COMPLEX                 :: iniGPcomplex(EQN%nVar)! complex reference      !
    CHARACTER(LEN=20)       :: varName               ! Name of analyse var.   !
    !-------------------------------------------------------------------------!
    !INTENT(IN)              :: pvar, ANALYSE, DISC, IO
    !-------------------------------------------------------------------------!

    L1norm(:)       = 0.
    L2norm(:)       = 0.
    Linfnorm(:)     = 0.
    MPIL1norm(:)    = 0.
    MPIL2norm(:)    = 0.
    MPILinfnorm(:)  = 0.

    ! Imaginary Unit
    IU = (0.,1.)

    SELECT CASE(ANALYSE%typ)
    CASE(0)
       logInfo0(*) 'No analysis of data.'
       logInfo0(*) 'CPU-Time: ',  DISC%LoopCPUTime
       RETURN
    CASE(1)
       logInfo0(*) 'AnalyseGalerkin3D: Analyse of data. Reference is the initial condition.'
       logInfo0(*) 'CPU-Time: ',  DISC%LoopCPUTime
    CASE(3)
       logInfo0(*) 'AnalyseGalerkin3D: Reference is elastic planarwave.'
    CASE(14)                                                     !
       logInfo0(*) 'AnalyseGalerkin3D: Analyse of data. '
       logInfo0(*) 'Reference: u(x,y,t)=u0*exp[ I ( w*t - kx*x - ky*y - kz*z )]'
       logInfo0(*) 'for Anelastic Seismic Waves  '
    CASE(15)                                                     !
       logInfo0(*) 'AnalyseGalerkin3D: Analyse of data. '
       logInfo0(*) 'Reference: u(x,y,t)=u0*exp[ I ( w*t - kx*x - ky*y - kz*z )]'
       logInfo0(*) 'for Anisotropic Seismic Waves  '
    CASE DEFAULT
       logError(*) 'Analysis of Type ',ANALYSE%typ,' unknown!'
       STOP
    END SELECT

    LInfElement = -1

    DO iElem = 1, MESH%nElem
        !
        LocElemType = MESH%LocalElemType(iElem)
        !
        SELECT CASE(LocElemType)
        CASE(4) ! Tetra
            DO iVert=1,MESH%nVertices_Tet
                x(iVert) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(iVert,iElem))
                y(iVert) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(iVert,iElem))
                z(iVert) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(iVert,iElem))
            ENDDO
            DO iIntGP = 1,DISC%Galerkin%nIntGP
                xi   = DISC%Galerkin%intGaussP_Tet(1,iIntGP)
                eta  = DISC%Galerkin%intGaussP_Tet(2,iIntGP)
                zeta = DISC%Galerkin%intGaussP_Tet(3,iIntGP)
                CALL TetraTrafoXiEtaZeta2XYZ(xGP(iIntGP),yGP(iIntGP),zGP(iIntGP),xi,eta,zeta,x,y,z)
            ENDDO
            JacobiDet = 6.*MESH%ELEM%Volume(iElem)
            intGaussW => DISC%Galerkin%intGaussW_Tet
        CASE(6) ! Hexa
            DO iVert=1,MESH%nVertices_Hex
                x(iVert) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(iVert,iElem))
                y(iVert) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(iVert,iElem))
                z(iVert) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(iVert,iElem))
            ENDDO
            DO iIntGP = 1, DISC%Galerkin%nIntGP
                xi   = DISC%Galerkin%intGaussP_Hex(1,iIntGP)
                eta  = DISC%Galerkin%intGaussP_Hex(2,iIntGP)
                zeta = DISC%Galerkin%intGaussP_Hex(3,iIntGP)
                CALL HexaTrafoXiEtaZeta2XYZ(xGP(iIntGP),yGP(iIntGP),zGP(iIntGP),xi,eta,zeta,x,y,z)
            ENDDO
            JacobiDet = MESH%ELEM%Volume(iElem)
            intGaussW => DISC%Galerkin%intGaussW_Hex
        END SELECT
        !
        DO iIntGP = 1, DISC%Galerkin%nIntGP
          !
          SELECT CASE(ANALYSE%typ)
          CASE(1) ! Initialcondition
             CALL InitialField(iniGP, iniGP_ane, 0., xGP(iIntGP), yGP(iIntGP), zGP(iIntGP), iElem, EQN, IC, SOURCE, IO)
          CASE(3) ! Initialcondition
             CALL InitialField(iniGP, iniGP_ane, time, xGP(iIntGP), yGP(iIntGP), zGP(iIntGP), iElem, EQN, IC, SOURCE, IO)
          CASE(14)
             kx    = ANALYSE%PWAN%Wavenumbers(1)
             ky    = ANALYSE%PWAN%Wavenumbers(2)
             kz    = ANALYSE%PWAN%Wavenumbers(3)
             ! Analytic solution
             iniGPcomplex(:) = (0.,0.)
             iniGP(:)        = 0.
             DO iVec = 1, ANALYSE%PW%SetVar
                 iniGPcomplex(:) = iniGPcomplex(:) +                ANALYSE%PW%ampfield(iVec)  * &
                                   ANALYSE%PWAN%EigenVec(1:EQN%nVar,ANALYSE%PW%varfield(iVec)) * &
                                   exp(IU*   (ANALYSE%PWAN%EigenVal(ANALYSE%PW%varfield(iVec)) * &
                                              time - kx*xGP(iIntGP) - ky*yGP(iIntGP) - kz*zGP(iIntGP)) )
             ENDDO
             iniGP(:) = real(iniGPcomplex(:))
          CASE(15)
             ! Analytic solution
             iniGPcomplex(:) = (0.,0.)
             iniGP(:)        = 0.
             DO j = 1,3
                kx    = IC%PWANISO(j)%Wavenumbers(1)
                ky    = IC%PWANISO(j)%Wavenumbers(2)
                kz    = IC%PWANISO(j)%Wavenumbers(3)
                DO iVec = 1, IC%PWANISO(j)%SetVar
                    iniGPcomplex(:) = iniGPcomplex(:) +                 IC%PWANISO(j)%ampfield(iVec)  * &
                                      IC%PWANISO(j)%EigenVec(1:EQN%nVar,IC%PWANISO(j)%varfield(iVec)) * &
                                      exp(IU*   (IC%PWANISO(j)%EigenVal(IC%PWANISO(j)%varfield(iVec)) * &
                                              time - kx*xGP(iIntGP) - ky*yGP(iIntGP) - kz*zGP(iIntGP)) )
                ENDDO
             ENDDO
             iniGP(:) = real(iniGPcomplex(:))
          END SELECT
          !
          CALL GetStateGP_new(stateGP,iElem,iIntGP,LocElemType,EQN,DISC)
          !
          pstateGP(:) = stateGP(:)
          !
          locError(:) = ABS(iniGP(:)-pstateGP(:))
          !
          L1norm(:)   = L1norm(:) + locError(:)    * intGaussW(iIntGp)*JacobiDet
          L2norm(:)   = L2norm(:) + locError(:)**2 * intGaussW(iIntGp)*JacobiDet
          DO iVar = 1, EQN%nVar
            IF(locError(iVar).GT.Linfnorm(iVar)) THEN
                Linfnorm(iVar) = locError(iVar)
                LInfElement(iVar) = iElem
            ENDIF
          ENDDO
          !
       ENDDO
       !
    ENDDO
    !
    ! Take square-root to obtain the L2-norm

    L2norm(:) = SQRT(L2norm(:))
    !
    DO iVar = 1,EQN%nVar                                                         ! Ueberpruefe welche Variablen analysiert werden
       IF (ANALYSE%variables(iVar) ) THEN                                        ! sollen
          !
          varName = IO%TitleMask(3+iVar)
          !
          logInfo(*) 'Error analysis of variable ', TRIM(varName)
          logInfo(*) '======================================'
          logInfo(*) '  L1_norm   : ', L1norm(iVar)
          logInfo(*) '  L2_norm   : ', L2norm(iVar)
          logInfo(*) '  Linf_norm : ', Linfnorm(iVar), ' Elem/Bary: ', &
                                       LInfElement(iVar), MESH%ELEM%xyBary(:,LInfElement(iVar))
          logInfo(*) '======================================'
#ifdef PARALLEL
          CALL MPI_REDUCE(L1norm(iVar),  MPIL1norm(iVar),  1,MPI%MPI_AUTO_REAL,MPI_SUM,0,MPI%commWorld,iErr)
          L2norm(iVar) = L2norm(iVar)**2
          CALL MPI_REDUCE(L2norm(iVar),  MPIL2norm(iVar),  1,MPI%MPI_AUTO_REAL,MPI_SUM,0,MPI%commWorld,iErr)
          MPIL2norm(iVar) = SQRT(MPIL2norm(iVar))
          CALL MPI_REDUCE(Linfnorm(iVar),MPILinfnorm(iVar),1,MPI%MPI_AUTO_REAL,MPI_MAX,0,MPI%commWorld,iErr)
          logInfo0(*) 'MPI Error analysis of variable ', varName
          logInfo0(*) '======================================'
          logInfo0(*) 'MPI L1_norm   : ', MPIL1norm(iVar)
          logInfo0(*) 'MPI L2_norm   : ', MPIL2norm(iVar)
          logInfo0(*) 'MPI Linf_norm : ', MPILinfnorm(iVar)
          logInfo0(*) '======================================'
#endif
       END IF
    END DO

    h1 = MAXVAL(MESH%ELEM%Volume(:)**(1./3.))
    logInfo(*) 'h1 (3D) =  ', h1
#ifdef PARALLEL
          CALL MPI_REDUCE(h1,       MPI_h1,         1,MPI%MPI_AUTO_REAL,MPI_MAX,0,MPI%commWorld,iErr)
          logInfo(*) 'h1 (3D,MPI) =  ', MPI_h1
#endif
    !                                                                            !
    logInfo(*) '  CPU-Time: ', DISC%LoopCPUTime                   !
    logInfo(*) '--------------------------------------'           !
    logInfo(*) 'AnalyseGalerkin erfolgreich. '                    !
    logInfo(*) '--------------------------------------'           !
    !                                                                            !
  END SUBROUTINE AnalyseGalerkin3D_us_new

  ! Read the 2d Green function file
  ! Used for computing the rupture velocity

  SUBROUTINE Read2dGF(DISC,IO)
    USE QuadPoints_mod
    TYPE(tInputOutput)       :: IO
    TYPE(tDiscretization)           :: DISC
    INTEGER :: allocstat                                  ! Allocation status !
    INTEGER :: stat                                       ! IO status         !
    INTEGER nMaxPoly,MaxDegFr,iPoly,iDegFr,DegFr,iEta,iXi,k,l
    LOGICAL                         :: configexist
    CHARACTER(LEN=200)   :: DGPATH
    CHARACTER(LEN=200)   :: FileName_Tri


    INQUIRE(                                            & !
     FILE= 'DGPATH'                                   , & !
     EXIST=configexist                                  ) !
    !                                                     !
    IF (configexist) THEN                                 !
       !                                                  !
       OPEN(                                            & !
        UNIT= IO%UNIT%FileIn                          , & !
        FILE= 'DGPATH'                                , & !
        IOSTAT = STAT                                   ) !
       !                                                  !
       IF (stat.NE.0) THEN                                !
      logError(*) 'cannot open DGPATH'                    !
      STOP
       END IF                                             !
       !                                                  !
       READ(IO%UNIT%FileIn,'(A)') DGPATH                  !
       !                                                  !
       CLOSE(IO%UNIT%FileIn)                              !
       !                                                  !
       !                                                  !
       logInfo0(*) 'Path to the DG directory is: ',TRIM(DGPATH)

    WRITE(FileName_Tri,'(a,a20)') TRIM(DGPATH), 'BasisFunctions2D.tri'
    OPEN( UNIT = IO%UNIT%FileIn, FILE = TRIM(FileName_Tri), IOSTAT = STAT, STATUS='OLD' )
    !
    IF(stat.NE.0) THEN
        logError(*) ' ERROR! File ', TRIM(FileName_Tri), ' could not be opened. '
        STOP
    ENDIF
    logInfo0(*) 'Reading basis functions and mass matrices for DG method '
    logInfo0(*) 'from file ', TRIM(FileName_Tri)

    READ(IO%UNIT%FileIn,*)
    READ(IO%UNIT%FileIn,*)
    ! Read maximal degree of basis polynomials stored in the file.
    READ(IO%UNIT%FileIn,*) nMaxPoly
    IF(DISC%Galerkin%nPoly.GT.nMaxPoly) THEN
        logError(*) 'ERROR: Required polynomial for DG method is higher than the ones stored in file ', TRIM(FileName_Tri)
        STOP
    ENDIF

    MaxDegFr = (DISC%Galerkin%nPoly+1)*(DISC%Galerkin%nPoly+2)/2

    ALLOCATE(                                                                                        &
        DISC%Galerkin%cPoly_Tri(0:nMaxPoly, 0:nMaxPoly, 0:MaxDegFr-1, 0:nMaxPoly),                   &
        DISC%Galerkin%MassMatrix_Tri(MaxDegFr,MaxDegFr, 0:nMaxPoly),                                 &
        DISC%Galerkin%NonZeroCPoly_Tri(0:MaxDegFr-1,0:nMaxPoly),                                     &
        DISC%Galerkin%NonZeroCPolyIndex_Tri(3,1:nMaxPoly**3,0:MaxDegFr,0:nMaxPoly),                  &
        STAT = allocstat)
    IF(allocStat .NE. 0) THEN
        logError(*) 'ERROR: could not allocate all variables!'
        STOP
    END IF

    ! Read coefficients of basis functions and mass matrices up to degree nMaxPoly
    !DO iPoly = 0, nMaxPoly
    DO iPoly = 0, DISC%Galerkin%nPoly
        ! Read comment in front of the basis functions' coefficients
        logInfo0(*) 'Reading basis functions of order ', iPoly
        READ(IO%UNIT%FileIn,*)
        DegFr = (iPoly + 1)*(iPoly + 2)/2
        ! Read polynomial coefficients
        ! where the index of the degrees of freedom starts at zero
        DO iDegFr = 0, DegFr-1
            DO iEta = 0, iPoly
                DO iXi = 0, iPoly
                    READ(IO%UNIT%FileIn,*) DISC%Galerkin%cPoly_Tri(iXi,iEta,iDegFr,iPoly)
                ENDDO
            ENDDO
        ENDDO
        ! Read comment in front of the entries of the mass matrix
        READ(IO%UNIT%FileIn,*)
        logInfo0(*)  'Reading mass matrices   of order ', iPoly
        ! Read entries of the mass matrix
        DO k = 1, DegFr
            DO l = 1, DegFr
                READ(IO%UNIT%FileIn,*) DISC%Galerkin%MassMatrix_Tri(k,l,iPoly)
            ENDDO
        ENDDO
    ENDDO
    CLOSE(IO%UNIT%FileIn)

    DISC%Galerkin%NonZeroCPoly_Tri(:,:)          = 0
    DISC%Galerkin%NonZeroCPolyIndex_Tri(:,:,:,:) = -1
    !DO iPoly = 0, nMaxPoly
    DO iPoly = 0, DISC%Galerkin%nPoly
       DegFr = (iPoly + 1)*(iPoly + 2)/2
       DO iDegFr = 0, DegFr-1
            DO iEta = 0, iPoly
              DO iXi = 0, iPoly
                 IF(ABS(DISC%Galerkin%cPoly_Tri(iXi,iEta,iDegFr,iPoly)).GE.1e-6) THEN
                    DISC%Galerkin%NonZeroCPoly_Tri(iDegFr,iPoly) = DISC%Galerkin%NonZeroCPoly_Tri(iDegFr,iPoly) + 1
                    DISC%Galerkin%NonZeroCPolyIndex_Tri(1,DISC%Galerkin%NonZeroCPoly_Tri(iDegFr,iPoly),iDegFr,iPoly) = iXi
                    DISC%Galerkin%NonZeroCPolyIndex_Tri(2,DISC%Galerkin%NonZeroCPoly_Tri(iDegFr,iPoly),iDegFr,iPoly) = iEta
                 ENDIF
              ENDDO
            ENDDO
       ENDDO
    ENDDO
    ENDIF

        ! Compute and store surface gaussian integration points
!~         CALL TriangleQuadraturePoints(                         &
!~                  nIntGP     = DISC%Galerkin%nBndGP,            &
!~                  IntGaussP  = DISC%Galerkin%BndGaussP_Tet,     &
!~                  IntGaussW  = DISC%Galerkin%BndGaussW_Tet,     &
!~                  M          = DISC%Galerkin%nPoly+2,           &
!~                  IO         = IO,                              &
!~                  quiet      = .TRUE.                           )

  END SUBROUTINE Read2dGF

END MODULE dg_setup_mod
