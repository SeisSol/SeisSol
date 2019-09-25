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
       iniGalerkin3D_us_level1_new,        &
       iniGalerkin3D_us_level2_new,        &
       iniGalerkin3D_us_intern_new,        &
       icGalerkin3D_us_new,                &
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
    ! Interoperability with C needs continuous arrays in memory.
    deallocate(disc%galerkin%dgvar)

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
    use iso_c_binding, only: c_loc, c_null_char, c_bool
    use f_ftoc_bind_interoperability
    use calc_deltaT_mod
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
    real                            :: l_timeStepWidth
    real                            :: l_loads(3), l_scalings(3), l_cuts(2), l_timeScalings(2), l_gts
    integer                         :: iObject, iSide, iNeighbor, MPIIndex
    real, target                    :: materialVal(EQN%nBackgroundVar)
    logical(kind=c_bool)                        :: enableFreeSurfaceIntegration
    
    ! ------------------------------------------------------------------------!
    !
    CALL iniGalerkin3D_us_intern_new(EQN, DISC, MESH, BND, IC, SOURCE, OptionalFields, IO)
    !
    CALL BuildSpecialDGGeometry3D_new(OptionalFields%BackgroundValue,EQN,MESH,DISC,BND,MPI,IO)

    ! we have the material parameters and insphere diameters, let's get some time step widths out of this

    ! Set used initial conditions.
    call c_interoperability_setInitialConditionType(trim(IC%cICType) // c_null_char)

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
  do iElem = 1, MESH%nElem
    iSide = 0

    materialVal = OptionalFields%BackgroundValue(iElem,:)
    call c_interoperability_setMaterial( i_elem = iElem,                                         \
                                         i_side = iSide,                                         \
                                         i_materialVal = materialVal,                            \
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
      call c_interoperability_setupFSRMPointSources( momentTensor      = SOURCE%RP%MomentTensor,      &
                                                     velocityComponent = SOURCE%RP%VelocityComponent, &
                                                     numberOfSources   = SOURCE%RP%nSbfs(1),          &
                                                     centres           = SOURCE%RP%SpacePosition,     &
                                                     strikes           = SOURCE%RP%Strks,             &
                                                     dips              = SOURCE%RP%Dips,              &
                                                     rakes             = SOURCE%RP%Rake,              &
                                                     onsets            = SOURCE%RP%Tonset,            &
                                                     areas             = SOURCE%RP%Area,              &
                                                     timestep          = SOURCE%RP%t_samp,            &
                                                     numberOfSamples   = SOURCE%RP%nsteps,            &
                                                     timeHistories     = SOURCE%RP%TimeHist       )
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

    IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ALLOCATE( DISC%LocalIteration(MESH%nElem) )
        ALLOCATE( DISC%LocalTime(MESH%nElem)      )
        ALLOCATE( DISC%LocalDt(MESH%nElem)        )
        DISC%LocalIteration(:)  = 0.
        DISC%LocalTime(:)       = 0.
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
      DISC%DynRup%StateVar(:,:) = EQN%IniStateVar
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
      allocate(disc%DynRup%output_StateVar(DISC%Galerkin%nBndGP,MESH%Fault%nSide))

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
    use iso_c_binding, only: c_loc
    use f_ftoc_bind_interoperability
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
    REAL    :: iniGP_Plast(6)                                                 ! Initial stress loading for plastic calculations
    REAL    :: phi                                                            ! Value of the base function at GP      !
    REAL    :: Einv, v                                                        ! Inverse of Young's modulus, Poisson ratio v
    !
    REAL, POINTER :: IntGaussP(:,:)     =>NULL()
    REAL, POINTER :: IntGaussW(:)       =>NULL()
    REAL, POINTER :: IntGPBaseFunc(:,:) =>NULL()
    REAL, POINTER :: MassMatrix(:,:)    =>NULL()
    ! temporary degrees of freedom
    real    :: l_initialLoading( NUMBER_OF_BASIS_FUNCTIONS, 6 )
    real    :: l_plasticParameters(2)
    !-------------------------------------------------------------------------!
    !
    IF(.NOT.DISC%Galerkin%init) THEN
       logError(*) 'icGalerkin: SeisSol Interface not initialized!!'
       STOP
    ENDIF
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

    call c_interoperability_projectInitialField()
    !
    iPoly  = DISC%Galerkin%nPoly
    nIntGP = DISC%Galerkin%nIntGP
    nDegFr = DISC%Galerkin%nDegFr

    !$omp parallel do schedule(static) shared(eqn, disc, mesh, ic, source, io, iPoly, nIntGp, nDegFr) private(iElem, iIntGP, iDegFr, iVar, iVert, eType, locPoly, locDegFr, xi, eta, zeta, xGp, yGp, zGp, x, y, z, phi, intGaussP, intGaussW, intGPBaseFunc, massMatrix,l_initialLoading,l_plasticParameters, iniGP_plast)
    DO iElem = 1,MESH%nElem
        l_initialLoading=0
        l_plasticParameters=0
        
        IF(EQN%Plasticity.EQ.1 .AND. EQN%PlastMethod .EQ. 2) THEN
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
              IF(EQN%Plasticity.EQ.1 .AND. EQN%PlastMethod .EQ. 2) THEN !average approach for plasticity
              ! L2 projection of initial stress loading for the plastic calculations onto the DOFs
                iniGP_Plast(:) = EQN%IniStress(1:6,iElem)
                DO iDegFr = 1, nDegFr
                   phi = IntGPBaseFunc(iDegFr,iIntGP)
                   l_initialLoading(iDegFr,1:6) = l_initialLoading(iDegFr,1:6) + IntGaussW(iIntGP)*iniGP_plast(:)*phi
                ENDDO
             ENDIF
          ENDDO !iIntGP

          DO iDegFr = 1, nDegFr
            l_initialLoading(iDegFr, :) = l_initialLoading( iDegFr, : ) / massMatrix(iDegFr,iDegFr)
          ENDDO

          NULLIFY(intGaussP)
          NULLIFY(intGaussW)
          NULLIFY(IntGPBaseFunc)
          NULLIFY(MassMatrix)
        ENDIF

        IF(EQN%Plasticity.EQ.1 .AND. EQN%PlastMethod .EQ. 0) THEN !high-order points approach
        !elementwise assignement of the initial loading
           l_initialLoading(1,1:6) = EQN%IniStress(1:6,iElem)
        ENDIF

#ifdef USE_PLASTICITY
        ! initialize the element dependent plastic parameters
        l_plasticParameters(1) = EQN%PlastCo(iElem) !element-dependent plastic cohesion
        l_plasticParameters(2) = EQN%BulkFriction(iElem) !element-dependent bulk friction
        
        ! initialize loading in C
        call c_interoperability_setInitialLoading( i_meshId = iElem, \
                                                   i_initialLoading = pack( l_initialLoading, .true. ) )

        !initialize parameters in C
        call c_interoperability_setPlasticParameters( i_meshId            = iElem, \
                                                      i_plasticParameters = l_plasticParameters )
#endif
    ENDDO ! iElem

#ifdef USE_PLASTICITY
    call c_interoperability_setTv( tv = EQN%Tv )
#endif

#ifdef USE_PLASTICITY
    ! TODO: redundant (see iniGalerkin3D_us_level2_new) call to ensure correct intitial loading in copy layers.
    call c_interoperability_synchronizeCellLocalData();
#endif
    !
    logInfo0(*) 'DG initial condition projection done. '
    !
  END SUBROUTINE icGalerkin3D_us_new

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
    ALLOCATE( DISC%Galerkin%MaxWaveSpeed(MESH%nElem,MESH%nSideMax) )
    DO j=1,MESH%nSideMax
      DISC%Galerkin%WaveSpeed(:,j,1)=SQRT((MaterialVal(:,3)+2.*MaterialVal(:,2))/(MaterialVal(:,1)))
      DISC%Galerkin%WaveSpeed(:,j,2)=SQRT((MaterialVal(:,2))/(MaterialVal(:,1)))
      DISC%Galerkin%WaveSpeed(:,j,3)=SQRT((MaterialVal(:,2))/(MaterialVal(:,1)))
      DISC%Galerkin%MaxWaveSpeed(:,j)=SQRT((MaterialVal(:,3)+2.*MaterialVal(:,2))/(MaterialVal(:,1)))
    ENDDO
    !
    CONTINUE
    !
  END SUBROUTINE BuildSpecialDGGeometry3D_new

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
