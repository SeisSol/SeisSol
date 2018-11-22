!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
!!
!! @section LICENSE
!! Copyright (c) 2007-2016, SeisSol Group
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
!!
!! @section DESCRIPTION
!! Module containing Dynamic Rupture initial model setups
!! includes background stress and nucleation types
!!
!! Can be edited by users: Please add your models as a new subroutines

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE ini_model_DR_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE DGBasis_mod
  USE read_backgroundstress_mod
  !use StressReader
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  REAL, PARAMETER :: ZERO = 0.0D0
  !---------------------------------------------------------------------------!
  INTERFACE DR_setup
     MODULE PROCEDURE DR_setup
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: DR_setup
  PRIVATE :: DR_basic_ini
  private :: rotateStressToFaultCS
  !---------------------------------------------------------------------------!
  PRIVATE :: friction_RSF34
  PRIVATE :: friction_RSF7
  PRIVATE :: friction_RSF101
  PRIVATE :: friction_RSF103
  PRIVATE :: friction_LSW
  PRIVATE :: friction_LSW6
  !---------------------------------------------------------------------------!

  CONTAINS

  !> Interface to dynamic rupture initial models
  !<
  SUBROUTINE DR_setup(EQN,DISC,MESH,IO,BND)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations), target       :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tInputOutput)             :: IO
    TYPE (tBoundary)               :: BND
    !-------------------------------------------------------------------------!
    INTENT(IN)                      :: MESH, BND
    INTENT(INOUT)                   :: IO, EQN, DISC
    ! -------------------------------------------------------------------------

    ! Basic DR setup, valid for all models
    CALL DR_basic_ini(DISC,EQN,MESH,BND)
    !-------------------------------------------------------------------------!

    ! Initialize model dependent (space dependent) friction law parameters
    SELECT CASE(EQN%FL)
    CASE(2,16)
      ! Initialization of friction for linear slip weakening
      CALL friction_LSW(DISC,EQN,MESH,BND)
    CASE(3,4)
      ! Initialization of initial slip rate and friction for rate and state friction
      CALL friction_RSF34(DISC,EQN,MESH,BND)
    CASE(6)
      ! Initialization of friction and fault strength for bi-material linear slip weakening
      CALL friction_LSW6(DISC,EQN,MESH,BND)
    CASE(7)
      ! Initialization of initial slip rate and friction for fast velocity weakening friction
      CALL friction_RSF7(DISC,EQN,MESH,BND)
    CASE(101)
     ! Initialization of initial slip rate and friction for SCEC TPV103
     CALL friction_RSF101(DISC,EQN,MESH,BND)
    CASE(103)
     ! Initialization of initial slip rate and friction for SCEC TPV103
     CALL friction_RSF103(DISC,EQN,MESH,BND)
    END SELECT  ! Initialize model dependent rate-and-state friction law parameters type

  END SUBROUTINE DR_setup


  !> Initialization of basic dynamic rupture setup
  !<
  SUBROUTINE DR_basic_ini(DISC,EQN,MESH,BND)                 ! global variables
    use JacobiNormal_mod, only: RotationMatrix3D
    use f_ftoc_bind_interoperability
    use iso_c_binding, only: c_null_char, c_bool
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE (tBoundary)               :: BND
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    integer                             :: i
    real, allocatable, dimension(:,:)   :: nuc_xx,nuc_yy,nuc_zz,nuc_xy,nuc_yz,nuc_xz
    logical(kind=c_bool)                :: faultParameterizedByTraction
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: MESH, BND
    INTENT(INOUT) :: EQN,DISC
    !-------------------------------------------------------------------------!

    ! Allocation of DR fields
    ALLOCATE(  EQN%IniMu(DISC%Galerkin%nBndGP,MESH%Fault%nSide),            &
               EQN%IniBulk_xx(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniBulk_yy(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniBulk_zz(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniStateVar(DISC%Galerkin%nBndGP,MESH%Fault%nSide),      &
               EQN%IniShearXY(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniShearYZ(DISC%Galerkin%nBndGP,MESH%Fault%nSide),       &
               EQN%IniShearXZ(DISC%Galerkin%nBndGP,MESH%Fault%nSide)        )
    ALLOCATE(  DISC%DynRup%Strength(DISC%Galerkin%nBndGP,MESH%Fault%nSide)  )
    ALLOCATE(  DISC%DynRup%RF(DISC%Galerkin%nBndGP,MESH%Fault%nSide)        )
    ALLOCATE(  DISC%DynRup%DS(DISC%Galerkin%nBndGP,MESH%Fault%nSide)        )
    ALLOCATE(  DISC%DynRup%cohesion(DISC%Galerkin%nBndGP,MESH%Fault%nSide)  )

    ! Allocate and initialize magnitude output
    ALLOCATE(  DISC%DynRup%magnitude_out(MESH%Fault%nSide)                  )
    DISC%DynRup%magnitude_out(:) = .FALSE.

    IF (DISC%DynRup%magnitude_output_on.EQ.1) THEN
       ALLOCATE(  DISC%DynRup%averaged_Slip(MESH%Fault%nSide)        )
       !ini magnitude output
       DISC%DynRup%magnitude_out(:) = .TRUE.
       DISC%DynRup%averaged_Slip(:) = 0.0D0
    ENDIF

    ! ini rupture front output
    DISC%DynRup%RF = .FALSE.
    !ini dyn.stress ouput
    DISC%DynRup%DS = .FALSE.

    ! Initialize '+'side elements for RF und DS output
    IF ((DISC%DynRup%RFtime_on .EQ. 1) .AND. (DISC%DynRup%DS_output_on .EQ. 1) ) THEN
       ! Loop over every mesh element
           DO i = 1, MESH%Fault%nSide
              IF (MESH%FAULT%Face(i,1,1) .NE. 0) THEN
                 DISC%DynRup%RF(:,i) = .TRUE.
                 DISC%DynRup%DS(:,i) = .TRUE.
              ENDIF
           ENDDO
     ELSEIF ((DISC%DynRup%RFtime_on .EQ. 1) .AND. (DISC%DynRup%DS_output_on .EQ. 0 )) THEN
           DO i = 1, MESH%Fault%nSide
              IF (MESH%FAULT%Face(i,1,1) .NE. 0) THEN
                 DISC%DynRup%RF(:,i) = .TRUE.
              ENDIF
          ENDDO
    ENDIF

    faultParameterizedByTraction = c_interoperability_faultParameterizedByTraction(trim(DISC%DynRup%ModelFileName) // c_null_char)    
    
    if (faultParameterizedByTraction) then
      call c_interoperability_addFaultParameter("T_n" // c_null_char, EQN%IniBulk_xx)
      call c_interoperability_addFaultParameter("T_s" // c_null_char, EQN%IniShearXY)
      call c_interoperability_addFaultParameter("T_d" // c_null_char, EQN%IniShearXZ)
      EQN%IniBulk_yy(:,:) = 0.0d0
      EQN%IniBulk_zz(:,:) = 0.0d0
      EQN%IniShearYZ(:,:) = 0.0d0
    else
      call c_interoperability_addFaultParameter("s_xx" // c_null_char, EQN%IniBulk_xx)
      call c_interoperability_addFaultParameter("s_yy" // c_null_char, EQN%IniBulk_yy)
      call c_interoperability_addFaultParameter("s_zz" // c_null_char, EQN%IniBulk_zz)
      call c_interoperability_addFaultParameter("s_xy" // c_null_char, EQN%IniShearXY)
      call c_interoperability_addFaultParameter("s_yz" // c_null_char, EQN%IniShearYZ)
      call c_interoperability_addFaultParameter("s_xz" // c_null_char, EQN%IniShearXZ)
    endif

    !frictional parameter initialization
    SELECT CASE(EQN%FL)
    CASE(0)
       CONTINUE
    CASE(2,6,16)
       ALLOCATE(  DISC%DynRup%D_C(DISC%Galerkin%nBndGP,MESH%Fault%nSide)       )
       ALLOCATE(  DISC%DynRup%Mu_S(DISC%Galerkin%nBndGP,MESH%Fault%nSide)      )
       ALLOCATE(  DISC%DynRup%Mu_D(DISC%Galerkin%nBndGP,MESH%Fault%nSide)      )
       call c_interoperability_addFaultParameter("cohesion" // c_null_char, DISC%DynRup%cohesion)
       call c_interoperability_addFaultParameter("d_c" // c_null_char, DISC%DynRup%D_C)
       call c_interoperability_addFaultParameter("mu_s" // c_null_char, DISC%DynRup%Mu_S)
       call c_interoperability_addFaultParameter("mu_d" // c_null_char, DISC%DynRup%Mu_D)
       if (EQN%FL == 16) then
         ALLOCATE(  DISC%DynRup%forced_rupture_time(DISC%Galerkin%nBndGP,MESH%Fault%nSide))
         call c_interoperability_addFaultParameter("forced_rupture_time" // c_null_char, DISC%DynRup%forced_rupture_time)
       end if

    CASE(3,4,7,101,103)
      ALLOCATE(  DISC%DynRup%RS_a_array(DISC%Galerkin%nBndGP, MESH%Fault%nSide)        )
      call c_interoperability_addFaultParameter("rs_a" // c_null_char, DISC%DynRup%RS_a_array)
      if (EQN%FL == 103) then
        allocate( DISC%DynRup%RS_srW_array(DISC%Galerkin%nBndGP, MESH%Fault%nSide), &
                  DISC%DynRup%RS_sl0_array(DISC%Galerkin%nBndGP,MESH%Fault%nSide),  &
                  nuc_xx(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_yy(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_zz(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_xy(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_yz(DISC%Galerkin%nBndGP,MESH%Fault%nSide),                    &
                  nuc_xz(DISC%Galerkin%nBndGP,MESH%Fault%nSide)                     )
        call c_interoperability_addFaultParameter("rs_srW" // c_null_char, DISC%DynRup%RS_srW_array)
        call c_interoperability_addFaultParameter("RS_sl0" // c_null_char, DISC%DynRup%RS_sl0_array)
        if (faultParameterizedByTraction) then
          call c_interoperability_addFaultParameter("Tnuc_n" // c_null_char, nuc_xx)
          call c_interoperability_addFaultParameter("Tnuc_s" // c_null_char, nuc_xy)
          call c_interoperability_addFaultParameter("Tnuc_d" // c_null_char, nuc_xz)
          nuc_yy(:,:) = 0.0d0
          nuc_zz(:,:) = 0.0d0
          nuc_yz(:,:) = 0.0d0
        else
          call c_interoperability_addFaultParameter("nuc_xx" // c_null_char, nuc_xx)
          call c_interoperability_addFaultParameter("nuc_yy" // c_null_char, nuc_yy)
          call c_interoperability_addFaultParameter("nuc_zz" // c_null_char, nuc_zz)
          call c_interoperability_addFaultParameter("nuc_xy" // c_null_char, nuc_xy)
          call c_interoperability_addFaultParameter("nuc_yz" // c_null_char, nuc_yz)
          call c_interoperability_addFaultParameter("nuc_xz" // c_null_char, nuc_xz)
        endif
      end if
    END SELECT

    call c_interoperability_initializeFault(  trim(DISC%DynRup%ModelFileName) // c_null_char, &
                                              EQN%GPwise,                                     &
                                              MESH%ELEM%BndGP_Tri,                            &
                                              DISC%Galerkin%nBndGP                            )

    ! Rotate initial stresses to fault coordinate system
    allocate(EQN%InitialStressInFaultCS(DISC%Galerkin%nBndGP,6,MESH%Fault%nSide))
    call rotateStressToFaultCS(EQN,MESH,DISC%Galerkin%nBndGP,EQN%IniBulk_xx,EQN%IniBulk_yy,EQN%IniBulk_zz,EQN%IniShearXY,EQN%IniShearYZ,EQN%IniShearXZ,EQN%InitialStressInFaultCS,faultParameterizedByTraction)
    
    if (EQN%FL == 103) then
      allocate(EQN%NucleationStressInFaultCS(DISC%Galerkin%nBndGP,6,MESH%Fault%nSide))
      call rotateStressToFaultCS(EQN,MESH,DISC%Galerkin%nBndGP,nuc_xx,nuc_yy,nuc_zz,nuc_xy,nuc_yz,nuc_xz,EQN%NucleationStressInFaultCS,faultParameterizedByTraction)
      deallocate(nuc_xx,nuc_yy,nuc_zz,nuc_xy,nuc_yz,nuc_xz)
    end if
  END SUBROUTINE DR_basic_ini
  
  SUBROUTINE rotateStressToFaultCS(EQN,MESH,nBndGP,s_xx,s_yy,s_zz,s_xy,s_yz,s_xz,stressInFaultCS,faultParameterizedByTraction)
    use JacobiNormal_mod, only: RotationMatrix3D
    use iso_c_binding, only: c_bool
    use create_fault_rotationmatrix_mod, only: create_fault_rotationmatrix
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)                      :: EQN
    TYPE(tUnstructMesh)                   :: MESH
    integer                               :: nBndGP
    real, allocatable, dimension(:,:)     :: s_xx,s_yy,s_zz,s_xy,s_yz,s_xz
    real, allocatable, dimension(:,:,:)   :: stressInFaultCS
    logical(kind=c_bool)                  :: faultParameterizedByTraction
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    integer                             :: iFace, iBndGP
    real                                :: normal(3)
    real                                :: tangent1(3)
    real                                :: tangent2(3)
    real                                :: T(9,9)
    real                                :: iT(9,9)
    real                                :: rotmat(6,6)
    real                                :: iRotmat(6,6)
    real                                :: Stress(1:6,1:nBndGP)
    real                                :: StressinFaultCSTmp(6)
    !-------------------------------------------------------------------------!
    intent(in)                          :: EQN,MESH,nBndGP
    intent(inout)                       :: s_xx,s_yy,s_zz,s_xy,s_yz,s_xz
    intent(inout)                       :: stressInFaultCS
    
    do iFace = 1, MESH%Fault%nSide
      normal   = MESH%Fault%geoNormals( 1:3, iFace)
      tangent1 = MESH%Fault%geoTangent1(1:3, iFace)
      tangent2 = MESH%Fault%geoTangent2(1:3, iFace)
      CALL RotationMatrix3D(normal, tangent1, tangent2, T(:,:), iT(:,:), EQN)

      Stress(1,:) = s_xx(:,iFace)
      Stress(2,:) = s_yy(:,iFace)
      Stress(3,:) = s_zz(:,iFace)
      Stress(4,:) = s_xy(:,iFace)
      Stress(5,:) = s_yz(:,iFace)
      Stress(6,:) = s_xz(:,iFace)
      
      if (faultParameterizedByTraction) then
        call create_fault_rotationmatrix(rotmat,iFace,EQN,MESH,iRotmat)
        do iBndGP=1,nBndGP
          StressinFaultCSTmp = MATMUL(iRotmat, Stress(:,iBndGP))
          Stress(:,iBndGP) = StressinFaultCSTmp
        enddo
        s_xx(:,iFace) = Stress(1,:)
        s_yy(:,iFace) = Stress(2,:)
        s_zz(:,iFace) = Stress(3,:)
        s_xy(:,iFace) = Stress(4,:)
        s_yz(:,iFace) = Stress(5,:)
        s_xz(:,iFace) = Stress(6,:)
      endif

      do iBndGP=1,nBndGP
        StressinFaultCSTmp = MATMUL(iT(1:6,1:6), Stress(:,iBndGP))
        stressInFaultCS(iBndGP,:,iFace) = StressinFaultCSTmp
      enddo
    enddo
  END SUBROUTINE rotateStressToFaultCS


  !> Initialization of initial slip rate and friction for rate and state friction
  !<
  SUBROUTINE friction_RSF34(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate, X2
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

      ! element ID
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)

      ! get vertices of complete tet
      IF (MESH%Fault%Face(i,1,1) == 0) THEN
          ! iElem is in the neighbor domain
          ! The neighbor element belongs to a different MPI domain
          iNeighbor           = MESH%Fault%Face(i,1,2)          ! iNeighbor denotes "-" side
          iLocalNeighborSide  = MESH%Fault%Face(i,2,2)
          iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
          MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
          !
          xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
          yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
          zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
      ELSE
          !
          ! get vertices
          xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
          yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
          zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
      ENDIF
      !
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)
          !
          EQN%IniStateVar(i,iBndGP) = DISC%DynRup%NucRS_sv0
          !EQN%IniStateVar(i,iBndGP) = DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0*EXP((sstress/(nstress*DISC%DynRup%RS_b))-DISC%DynRup%RS_f0/DISC%DynRup%RS_b-DISC%DynRup%RS_a_array(iBndGP,i)/DISC%DynRup%RS_b*LOG(iniSlipRate/DISC%DynRup%RS_sr0))
          X2  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a)
          EQN%IniMu(iBndGP,i)=DISC%DynRup%RS_a * LOG(X2 + SQRT(X2**2 + 1.0))

      ENDDO ! iBndGP

  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF34      ! Initialization of initial slip rate and friction for rate and state friction

  !> Initialization of initial slip rate and friction for fast velocity weakening friction
  !<
  SUBROUTINE friction_RSF7(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

      ! element ID
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)

      ! get vertices of complete tet
      IF (MESH%Fault%Face(i,1,1) == 0) THEN
          ! iElem is in the neighbor domain
          ! The neighbor element belongs to a different MPI domain
          iNeighbor           = MESH%Fault%Face(i,1,2)          ! iNeighbor denotes "-" side
          iLocalNeighborSide  = MESH%Fault%Face(i,2,2)
          iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
          MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
          !
          xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
          yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
          zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
      ELSE
          !
          ! get vertices
          xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
          yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
          zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
      ENDIF
      !
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)
                   !
          EQN%IniStateVar(i,iBndGP)= (DISC%DynRup%RS_sl0*(DISC%DynRup%RS_a*EQN%IniBulk_yy(i,iBndGP)*iniSlipRate &
              + (iniSlipRate + DISC%DynRup%RS_sr0)*(DISC%DynRup%RS_f0*EQN%IniBulk_yy(i,iBndGP)-EQN%IniShearXY(i,iBndGP)))) &
              / (DISC%DynRup%RS_a*EQN%IniBulk_yy(i,iBndGP)*iniSlipRate-(iniSlipRate + DISC%DynRup%RS_sr0) &
              * (DISC%DynRup%RS_b*EQN%IniBulk_yy(i,iBndGP)-DISC%DynRup%RS_f0*EQN%IniBulk_yy(i,iBndGP)+EQN%IniShearXY(i,iBndGP)))
          EQN%IniMu(iBndGP,i)=DISC%DynRup%RS_f0+DISC%DynRup%RS_a*iniSlipRate/(iniSlipRate+DISC%DynRup%RS_sr0)-DISC%DynRup%RS_b*EQN%IniStateVar(i,iBndGP)/(EQN%IniStateVar(i,iBndGP)+DISC%DynRup%RS_sl0)
      ENDDO ! iBndGP

  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF7      ! Initialization of initial slip rate and friction for fast velocity weakening friction

  !> Initialization of initial slip rate and friction for SCEC TPV101
  !<
  SUBROUTINE friction_RSF101(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate, tmp
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

      ! element ID
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)

      ! get vertices of complete tet
      IF (MESH%Fault%Face(i,1,1) == 0) THEN
          ! iElem is in the neighbor domain
          ! The neighbor element belongs to a different MPI domain
          iNeighbor           = MESH%Fault%Face(i,1,2)          ! iNeighbor denotes "-" side
          iLocalNeighborSide  = MESH%Fault%Face(i,2,2)
          iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
          MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
          !
          xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
          yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
          zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
      ELSE
          !
          ! get vertices
          xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
          yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
          zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
      ENDIF
      !
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)
               !
          EQN%IniStateVar(i,iBndGP) = (DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0) * EXP((-EQN%IniShearXY(i,iBndGP)/EQN%IniBulk_yy(i,iBndGP)-DISC%DynRup%RS_f0-DISC%DynRup%RS_a_array(iBndGP,i)*LOG(iniSlipRate/DISC%DynRup%RS_sr0))/DISC%DynRup%RS_b)
          ! ASINH(X)=LOG(X+SQRT(X^2+1))
          tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a_array(iBndGP,i))
          EQN%IniMu(iBndGP,i)=DISC%DynRup%RS_a_array(iBndGP,i) * LOG(tmp + SQRT(tmp**2 + 1.0))

      ENDDO ! iBndGP

  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF101      ! Initialization of initial slip rate and friction for SCEC TPV101

  !> Initialization of initial slip rate and friction for SCEC TPV103
  !<
  SUBROUTINE friction_RSF103(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  USE JacobiNormal_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: iFace, iBndGP
  REAL                           :: iniSlipRate, tmp

  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)


  ! Loop over every mesh element
  DO iFace = 1, MESH%Fault%nSide
      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          tmp = ABS(SQRT(EQN%InitialStressInFaultCS(iBndGP,4,iFace)**2+EQN%InitialStressInFaultCS(iBndGP,6,iFace)**2)/(DISC%DynRup%RS_a_array(iBndGP,iFace)*EQN%InitialStressInFaultCS(iBndGP,1,iFace)))
          EQN%IniStateVar(iBndGP,iFace)=DISC%DynRup%RS_a_array(iBndGP,iFace)*LOG(2.0D0*DISC%DynRup%RS_sr0/iniSlipRate * (EXP(tmp)-EXP(-tmp))/2.0D0)
          ! ASINH(X)=LOG(X+SQRT(X^2+1))
          tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP(EQN%IniStateVar(iBndGP,iFace)/ DISC%DynRup%RS_a_array(iBndGP,iFace))
          EQN%IniMu(iBndGP,iFace)=DISC%DynRup%RS_a_array(iBndGP,iFace) * LOG(tmp + SQRT(tmp**2 + 1.0D0))
      ENDDO ! iBndGP
  ENDDO !    MESH%Fault%nSide

  END SUBROUTINE friction_RSF103      ! Initialization of initial slip rate and friction for SCEC TPV103

  !> Initialization of friction for linear slip weakening
  !<
  SUBROUTINE friction_LSW(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  EQN%IniMu(:,:) = DISC%DynRup%Mu_S(:,:)

  END SUBROUTINE friction_LSW

    SUBROUTINE friction_LSW6(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: iFace,iBndGP
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!
  
  EQN%IniMu(:,:) = DISC%DynRup%Mu_S(:,:)
  do iFace = 1, MESH%Fault%nSide
    do iBndGP = 1,DISC%Galerkin%nBndGP
      DISC%DynRup%Strength(iBndGP,iFace) = EQN%IniMu(iBndGP,iFace) * EQN%InitialStressInFaultCS(iBndGP,1,iFace)
    enddo
  enddo

  END SUBROUTINE friction_LSW6    ! Initialization of friction for bimaterial linear slip weakening



  END MODULE ini_model_DR_mod
