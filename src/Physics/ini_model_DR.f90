!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
!!
!! @section LICENSE
!! Copyright (c) 2014, SeisSol Group
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
  USE faultinput_mod
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
  !---------------------------------------------------------------------------!
  PRIVATE :: background_HOM
  PRIVATE :: background_TPV5
  PRIVATE :: background_STEP
  PRIVATE :: background_SMOOTH
  PRIVATE :: background_STEP2
  PRIVATE :: background_SMOOTH_GP
  PRIVATE :: background_TPV10
  PRIVATE :: background_TPV11
  PRIVATE :: background_DIP_RSF
  PRIVATE :: background_TPV1415
  PRIVATE :: background_TPV1617
  PRIVATE :: background_TOH1
  PRIVATE :: background_LAN1
  PRIVATE :: background_LAN2
  PRIVATE :: background_LAN3
  PRIVATE :: background_ALA
  PRIVATE :: background_NORTH
  PRIVATE :: background_TPV101
  PRIVATE :: background_TPV103
  !---------------------------------------------------------------------------!
  PRIVATE :: nucleation_STEP
  PRIVATE :: nucleation_SMOOTH_GP
  PRIVATE :: nucleation_ELLIPSE
  PRIVATE :: nucleation_TPV28_GP
  PRIVATE :: nucleation_TPV28

  PRIVATE :: background_TPV29
  PRIVATE :: background_SUMATRA
  PRIVATE :: background_TPV3132
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
                          
    ! Initialize background stress type
    SELECT CASE(DISC%DynRup%BackgroundType)
    CASE(0) 
       ! homogeneous case of background stress field
       CALL background_HOM(DISC,EQN,MESH)
    CASE(1)
       ! SCEC TPV5 test
       CALL background_TPV5(DISC,EQN,MESH)
    CASE(2)
       ! Depth dependence of stresses - step function
       CALL background_STEP(DISC,EQN,MESH)
    CASE(3)
       ! Depth dependence of stresses - smooth transition
       CALL background_SMOOTH(DISC,EQN,MESH)
    CASE(4)
       ! depth dependence of stresses and frictional LSW parameters - step function
       CALL background_STEP2(DISC,EQN,MESH)
    CASE(5)
       ! smooth depth dependence of stresses
       ! in contrast to type 3, in type 5 the stress is assigned to each Gaussian node       
       CALL background_SMOOTH_GP(DISC,EQN,MESH,BND)
    CASE(10,13)
       ! SCEC TPV10 test dipping fault subshear
       CALL background_TPV10(DISC,EQN,MESH,BND)
   !CASE(13)
       ! SCEC TPV12/13 with the same geometry but different depth-dependent stresses
       !CALL background_TPV10(DISC,EQN,MESH,BND)
    CASE(11)
       ! SCEC TPV11 test dipping fault supershear
       CALL background_TPV11(DISC,EQN,MESH,BND)
    CASE(12)
       ! Dipping fault with rate-and-state friction
       CALL background_DIP_RSF(DISC,EQN,MESH,BND)
    CASE(14,15)
       !  SCEC TPV14/15 test branching faults
       CALL background_TPV1415(DISC,EQN,MESH)
   CASE(16,17)
       !  SCEC TPV16/17 with heterogeneous initial stress field
       CALL background_TPV1617(EQN,MESH,IO,DISC,BND)
   CASE(29,30)
       !  SCEC TPV29 test case : rough fault (+plasticity = 30)
       CALL background_TPV29 (DISC,EQN,MESH,BND)
   CASE(31,32)
       !  SCEC TPV3132 test case : strike slip rupture in layered medium
       CALL background_TPV3132 (DISC,EQN,MESH,BND)
    CASE(50)
       ! Tohoku 1
       CALL background_TOH1(DISC,EQN,MESH)
    CASE(60)
       ! Landers full fault system
       CALL background_LAN1(DISC,EQN,MESH,BND)
    CASE(61)
       ! Landers segmented fault system 1
       CALL background_LAN2(DISC,EQN,MESH,BND)
    CASE(62)
       ! Landers segmented fault system 2
       CALL background_LAN3(DISC,EQN,MESH,BND)
    CASE(70)
       ! Alaska background stress model
       CALL background_ALA(DISC,EQN,MESH,BND)
    CASE(100)
       ! Northridge background stress model
       CALL background_NORTH(DISC,EQN,MESH,BND)   
    CASE(101)
       ! SCEC TPV101 test with rate-and-state friction (ageing law)
       CALL background_TPV101(DISC,EQN,MESH,BND)
    CASE(103)
       ! SCEC TPV103 test with velocity weakening friction (based on slip law)
       CALL background_TPV103(DISC,EQN,MESH,BND)
    CASE(120)
       CALL background_SUMATRA(DISC,EQN,MESH,BND)
    !
    ! Add your background stress model subroutine call here
    !   
    CASE DEFAULT
       logError(*) 'Chosen BackgroundType value ',DISC%DynRup%BackgroundType,' is not valid in present version of the code!'
       STOP
    END SELECT ! Initialize background stress type
    !-------------------------------------------------------------------------! 
    
    ! Initialize Nucleation type                  
    SELECT CASE(DISC%DynRup%Nucleation)
    CASE(0)
       ! 0=No nucleation zone
       CONTINUE
    CASE(1) 
       ! Nucleation by discontinuous jump on properties at [NucXMin,NucXMax] x [NucYMin,NucYMax]
       CALL nucleation_STEP(DISC,EQN,MESH)
    CASE(2)
       ! Nucleation by smooth jump on properties assigned to each Gaussian node
       CALL nucleation_SMOOTH_GP(DISC,EQN,MESH,BND)
    CASE(3)
       ! Nucleation by discontinuous elliptic nucleation zone
       CALL nucleation_ELLIPSE(DISC,EQN,MESH)
    CASE(28)
       ! Nucleation as in  SCEC TPV28 test assigned to each Gaussian node
       CALL nucleation_TPV28_GP(DISC,EQN,MESH,BND)
    CASE(29)
       ! Nucleation as in  SCEC TPV28 test assigned to each element
       CALL nucleation_TPV28(DISC,EQN,MESH)
    !
    ! Add your nucleation model subroutine call here
    ! 
    CASE DEFAULT  
       logError(*) 'Chosen DISC%DynRup%Nucleation type ',DISC%DynRup%Nucleation,' is not implemented in present version of the code!'
       STOP                    
    END SELECT  ! Initialize Nucleation type  
    !-------------------------------------------------------------------------! 
      
    ! Initialize model dependent (space dependent) friction law parameters
    SELECT CASE(EQN%FL)
    CASE(1,2,13,16,17)
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
    !-------------------------------------------------------------------------! 
    
    ! Read fault parameters from Par_file_faults
    if (DISC%DynRup%read_fault_file == 1) then
       call faultinput(disc,eqn,mesh,bnd,IO)
    end if
      
  END SUBROUTINE DR_setup


  !> Initialization of basic dynamic rupture setup
  !<
  SUBROUTINE DR_basic_ini(DISC,EQN,MESH,BND)                 ! global variables
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH   
    TYPE (tBoundary)               :: BND 
    !-------------------------------------------------------------------------!    
    ! Local variable declaration
    INTEGER			   :: i
    INTEGER                        :: iSide,iElem,iBndGP
    INTEGER                        :: iLocalNeighborSide,iNeighbor
    INTEGER                        :: MPIIndex, iObject
    REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
    REAL                           :: chi,tau
    REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
    REAL                           :: r, Vs, r_crit, hypox, hypoy, hypoz
    !-------------------------------------------------------------------------! 
    INTENT(IN)    :: MESH, BND
    INTENT(INOUT) :: EQN,DISC
    !-------------------------------------------------------------------------! 
  
    ! Allocation of DR fields
    ALLOCATE(  EQN%IniMu(MESH%Fault%nSide,DISC%Galerkin%nBndGP),            &
               EQN%IniBulk_xx(MESH%Fault%nSide,DISC%Galerkin%nBndGP),       &
               EQN%IniBulk_yy(MESH%Fault%nSide,DISC%Galerkin%nBndGP),       &
               EQN%IniBulk_zz(MESH%Fault%nSide,DISC%Galerkin%nBndGP),       &
               EQN%IniStateVar(MESH%Fault%nSide,DISC%Galerkin%nBndGP),      &
               EQN%IniShearXY(MESH%Fault%nSide,DISC%Galerkin%nBndGP),       &
               EQN%IniShearYZ(MESH%Fault%nSide,DISC%Galerkin%nBndGP),       &
               EQN%IniShearXZ(MESH%Fault%nSide,DISC%Galerkin%nBndGP)        )
    ALLOCATE(  DISC%DynRup%Strength(MESH%Fault%nSide,DISC%Galerkin%nBndGP)  )
    ALLOCATE(  DISC%DynRup%RF(MESH%Fault%nSide,DISC%Galerkin%nBndGP)        )
    ALLOCATE(  DISC%DynRup%cohesion(MESH%Fault%nSide,DISC%Galerkin%nBndGP)  )
      
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
    !T. Ulrich 8.2015 initial rupture time array (for Vr calculations)
    ALLOCATE(DISC%DynRup%rupture_time(MESH%Fault%nSide,DISC%Galerkin%nBndGP))
    DISC%DynRup%rupture_time(:,:)=0.
    
    !frictional parameter initialization
    SELECT CASE(EQN%FL)
    CASE(0)
       CONTINUE
    CASE(1,2,6,16,17)
       ! ini D_C and mu fields to constant on the entire fault (for LSW friction cases)
       ALLOCATE(  DISC%DynRup%D_C(MESH%Fault%nSide,DISC%Galerkin%nBndGP)       )
       ALLOCATE(  DISC%DynRup%Mu_S(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )
       ALLOCATE(  DISC%DynRup%Mu_D(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )
       DISC%DynRup%D_C(:,:)  = DISC%DynRup%D_C_ini
       DISC%DynRup%Mu_S(:,:) = DISC%DynRup%Mu_S_ini
       DISC%DynRup%Mu_D(:,:) = DISC%DynRup%Mu_D_ini
       EQN%IniMu(:,:)    =  DISC%DynRup%Mu_S_ini ! will be mapped to DISC%DynRup%Mu in dg_setup
       !
    CASE(13)!LSW with lower static coefficient of friction inside the nucleation zone needs additional initialisation for Mu_SNuc
       ! 
       ALLOCATE(  DISC%DynRup%D_C(MESH%Fault%nSide,DISC%Galerkin%nBndGP)       )
       ALLOCATE(  DISC%DynRup%Mu_S(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )
       ALLOCATE(  DISC%DynRup%Mu_D(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )
       DISC%DynRup%D_C(:,:)  = DISC%DynRup%D_C_ini       
       DISC%DynRup%Mu_D(:,:) = DISC%DynRup%Mu_D_ini
       
       !Mu_S is different inside a specified nucleation patch (patch is read in in nucleation case 13)
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
         
            xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
            yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
            zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
        ELSE 
           ! get vertices
            xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
            yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
            zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
        ENDIF
      
        DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
           ! Transformation of boundary GP's into XYZ coordinate system
            chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
            tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
            CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
            CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)

            IF ((xGP.LE. 1500.0D0) .AND. (xGP.GE.-1500.0D0)       &
                  .AND. (zGP.LE.-9093.26674D0) .AND. (zGP.GE.-11691.342951D0)) THEN
                DISC%DynRup%Mu_S(i,iBndGP) = DISC%DynRup%Mu_SNuc_ini
                EQN%IniMu(i,iBndGP) = DISC%DynRup%Mu_SNuc_ini ! will be mapped to DISC%DynRup%Mu in dg_setup
            ELSE
                DISC%DynRup%Mu_S(i,iBndGP) = DISC%DynRup%Mu_S_ini
                EQN%IniMu(i,iBndGP) = DISC%DynRup%Mu_S_ini ! will be mapped to DISC%DynRup%Mu in dg_setup
            ENDIF

         ENDDO ! iBndGP
      
     ENDDO !    MESH%Fault%nSide   

     CASE(30) !smooth forced rupture for benchmarks like TPV29/30 and TPV26/27

       ALLOCATE(  DISC%DynRup%D_C(MESH%Fault%nSide,DISC%Galerkin%nBndGP)       )
       ALLOCATE(  DISC%DynRup%Mu_S(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )
       ALLOCATE(  DISC%DynRup%Mu_D(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )
       ALLOCATE(DISC%DynRup%forced_rupture_time(MESH%Fault%nSide,DISC%Galerkin%nBndGP))

       DISC%DynRup%D_C(:,:)  = DISC%DynRup%D_C_ini
       DISC%DynRup%Mu_S(:,:) = DISC%DynRup%Mu_S_ini
       DISC%DynRup%Mu_D(:,:) = DISC%DynRup%Mu_D_ini
       EQN%IniMu(:,:)        =  DISC%DynRup%Mu_S_ini ! will be mapped to DISC%DynRup%Mu in dg_setup

       Vs = SQRT(EQN%mu/EQN%rho0)
       r_crit = DISC%DynRup%R_crit
       hypox = DISC%DynRup%XHypo
       hypoy = DISC%DynRup%YHypo
       hypoz = DISC%DynRup%ZHypo

       !calculate time of forced rupture for every BndGP, dependent of the distance to the hypocenter
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
         
            xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
            yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
            zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
        ELSE 
           ! get vertices
            xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
            yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
            zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
        ENDIF
      
        DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
           ! Transformation of boundary GP's into XYZ coordinate system
            chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
            tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
            CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
            CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)

           r = SQRT((xGP-hypox)**2+(yGP-hypoy)**2+(zGP-hypoz)**2)
           
           IF (r.LE.r_crit) THEN
              DISC%DynRup%forced_rupture_time(i,iBndGP) = r/(0.7d0*Vs)+(0.081d0*r_crit/(0.7d0*Vs))*(1d0/(1d0-(r/r_crit)*(r/r_crit))-1d0)
           ELSE
              DISC%DynRup%forced_rupture_time(i,iBndGP) = 1d9
           ENDIF

       ENDDO !iBndGP
      ENDDO !i

    CASE(3,4,7,12,101,103)
       ! ini initial slip rate fields to zero (for rate and state friction cases) 
       EQN%IniSlipRate1 = ZERO
       EQN%IniSlipRate2 = ZERO
       ! ini friction coefficient (will be copied on DISC%DynRup%Mu in dg_setup.f90 iniGalerkin3D_us_level2_new)
       EQN%IniMu(:,:)    =  DISC%DynRup%RS_f0 
    END SELECT
    
    ! ini of bimaterial case (simple planar case) 
    ! ALICE: The following line from original ini_model needs to be confirmed
    !DISC%DynRup%Strength(i,:) = EQN%IniMu(i,:)*EQN%IniBulk_yy(i,:)
    
  END SUBROUTINE DR_basic_ini
  
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  
  !> Homogeneous background stress field
  !<
  SUBROUTINE background_HOM(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! used for SCEC TPV3 test
  
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! constant background stress tensor and state variale
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
       
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_HOM ! Homogeneous background stress field
  
  !> SCEC TPV5 test case background stress field
  !<
  SUBROUTINE background_TPV5(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: LocX(3), LocY(3)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! SCEC TPV5 test
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/tpv5docs.html .
  ! center of nucleation patch is at 7.5km depth; size 3km x 3km
  ! fault reaches the surface at depth 0km
  ! units of mesh in meter

  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
      
      ! choose 2D fault plane for nucleation
      SELECT CASE(DISC%DynRup%NucDirX)
          CASE(1) !x direction
              LocX(:)=xp(1:3)
          CASE(2) !y direction
              LocX(:)=yp(1:3)
          CASE(3) !z direction  
              LocX(:)=zp(1:3)
      END SELECT
            
      SELECT CASE(DISC%DynRup%NucDirY)
          CASE(1) !x direction
              LocY(:)=xp(1:3)
          CASE(2) !y direction
              LocY(:)=yp(1:3)
          CASE(3) !z direction  
              LocY(:)=zp(1:3)
      END SELECT 
            
      ! right of the nucleation zone: square patch of lower initial shear stress centered at 7.5km depth
      IF(   MAXVAL(LocX(1:3)).LE.(9000.0D0) .AND. MINVAL(LocX(1:3)).GE.(6000.0D0)       &
          .AND. MAXVAL(LocY(1:3)).LE.(-6000.0D0) .AND. MINVAL(LocY(1:3)).GE.(-9000.0D0)) THEN
          EQN%IniShearXY(i,:)  = 62.0e6
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 0.0D0
          EQN%IniBulk_xx(i,:)  = DISC%DynRup%NucBulk_xx_0
          EQN%IniBulk_yy(i,:)  = DISC%DynRup%NucBulk_yy_0
          EQN%IniBulk_zz(i,:)  = DISC%DynRup%NucBulk_zz_0
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ENDIF
            
      ! left of the nucleation zone: square patch of higher initial shear stress centered at 7.5km depth
      IF(   MAXVAL(LocX(1:3)).LE.(-6000.0D0) .AND. MINVAL(LocX(1:3)).GE.(-9000.0D0)     &
          .AND. MAXVAL(LocY(1:3)).LE.(-6000.0D0) .AND. MINVAL(LocY(1:3)).GE.(-9000.0D0)) THEN
          EQN%IniShearXY(i,:)  = 78.0e6
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 0.0D0
          EQN%IniBulk_xx(i,:)  = DISC%DynRup%NucBulk_xx_0
          EQN%IniBulk_yy(i,:)  = DISC%DynRup%NucBulk_yy_0
          EQN%IniBulk_zz(i,:)  = DISC%DynRup%NucBulk_zz_0
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ENDIF
       
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_TPV5       ! SCEC TPV5 test  

  !> depth dependence of stresses - step function
  !<
  SUBROUTINE background_STEP(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! NOTE: y is depth, free surface is at y=+12500m, fault ends at y=-12500m
  ! add more/reduce nr of layers if necessary
      
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
      
      ! layer of stepwise different background stress
      IF ((sum(yp(:))/3.0D0) .LE. 12500.0D0 .AND. (sum(yp(:))/3.0D0) .GT. 12200.0D0 ) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 2.5e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 10.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ELSEIF ((sum(yp(:))/3.0D0) .LE. 12200.0D0 .AND. (sum(yp(:))/3.0D0) .GT. 10700.0D0) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 10.0e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 40.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ELSEIF ((sum(yp(:))/3.0D0) .LE. 10700.0D0 .AND. (sum(yp(:))/3.0D0) .GT. 9200.0D0) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 17.5e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 70.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ELSEIF ((sum(yp(:))/3.0D0) .LE. 9200.0D0) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 25.0e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 100.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0                              
      ENDIF
           
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_STEP       ! Depth dependent stress - step function

  !> depth dependence of stresses - smooth transition
  !<
  SUBROUTINE background_SMOOTH(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: average
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! NOTE: y is depth, free surface is at y=+12500m, fault ends at y=-12500m
  ! stress is align to a complete element!
      
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
      
      average = sum(yp(:))/3.0D0 - 12500.0D0
            
      IF (average .GT. -3000.0D0) THEN
          EQN%IniBulk_zz(i,:)  = 10.0e6 + 30.0e6*abs(average)/1000.0D0
          EQN%IniShearXZ(i,:)  =  2.5e6 + 7.5e6*abs(average)/1000.0D0
      ENDIF
      
     ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_SMOOTH       ! Depth dependent stress - smooth transition
  
  !> depth dependence of stresses and frictional LSW parameters - step function 
  !<
  SUBROUTINE background_STEP2(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! NOTE: y is depth, free surface is at y=+12500m, fault ends at y=-12500m
  ! add more/reduce nr of layers if necessary
      
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
     
      IF ((sum(yp(:))/3.0D0) .LE. 12500.0D0 .AND. (sum(yp(:))/3.0D0) .GT. 12200.0D0 ) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 2.5e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 10.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
           
          DISC%DynRup%D_C(i,:) = 6.0D0
          DISC%DynRup%Mu_S(i,:) = 3.0D0
          DISC%DynRup%Mu_D(i,:) = 3.0D0

      ELSEIF ((sum(yp(:))/3.0D0) .LE. 12200.0D0 .AND. (sum(yp(:))/3.0D0) .GT. 10700.0D0) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 10.0e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 40.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ELSEIF ((sum(yp(:))/3.0D0) .LE. 10700.0D0 .AND. (sum(yp(:))/3.0D0) .GT. 9200.0D0) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 17.5e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 70.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ELSEIF ((sum(yp(:))/3.0D0) .LE. 9200.0D0) THEN
          EQN%IniShearXY(i,:)  = 0.0D0
          EQN%IniShearYZ(i,:)  = 0.0D0
          EQN%IniShearXZ(i,:)  = 25.0e6
          EQN%IniBulk_xx(i,:)  = 0.0D0
          EQN%IniBulk_yy(i,:)  = 0.0D0
          EQN%IniBulk_zz(i,:)  = 100.0e6
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ENDIF 
     
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_STEP2       ! Depth dependent stress and LSW friction parameters - step function

  !> depth dependence of stresses - GP wise smooth transition
  !<
  SUBROUTINE background_SMOOTH_GP(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: average
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! NOTE: y is depth, free surface is at y=+12500m, fault ends at y=-12500m
  ! stress is align to a complete element!
      
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
      
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
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
          !
          !
          average = yGP - 12500.0D0 ! Note, average is not an average. I just recycle this variable.
          !
          IF (average .GT. -3000.0D0) THEN
              EQN%IniBulk_zz(i,iBndGP)  = 10.0e6 + 30.0e6*abs(average)/1000.0D0
              EQN%IniShearXZ(i,iBndGP)  = 2.5e6 + 7.5e6*abs(average)/1000.0D0
          ENDIF
          !
      ENDDO ! iBndGP
      
     ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_SMOOTH_GP       ! Depth dependent stress - GP wise smooth transition

  !> SCEC TPV10 test case
  !<
  SUBROUTINE background_TPV10(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: VertexSide(4,3)
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: average
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  LOGICAL                        :: nodewise=.FALSE.
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! SCEC TPV10  test with a dipping fault and static friction 0.760
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/download/TPV10_11_Description_v7.pdf
  ! fault is in the xz-plane, dip in negative y-direction
  ! depth is z negative
  ! center of nucleation patch: x=0.0 km, depth= -12.5km, size 3km x 3km
  ! fault reaches the surface at depth 0km
  ! units of mesh in meter

  ! SCEC TPV12/13 test with dipping fault with the same geometry as TPV10, different stress values
  ! Requires the correct mesh! For instructions see 
                
  ! switch for Gauss node wise stress assignment
  nodewise = .TRUE.
           
  ! cohesion 0.2MPa for SCEC TPV10/11 and TPV13/14
  DISC%DynRup%cohesion = -0.2e6
                      
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !
  
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
     
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
      
      ! Gauss node coordinate definition and stress assignment
      IF (nodewise) THEN
          
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
          DO iBndGP = 1,DISC%Galerkin%nBndGP
              
                ! Transformation of boundary GP's into XYZ coordinate system
                chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
                tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
                CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
                CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
                !
                !geometrically averaged depth (negative in z)
                average = zGP   ! Averaging not needed here
                !
                SELECT CASE(DISC%DynRup%BackgroundType)
                   CASE(10) !TPV10
                      ! depth dependent smooth initial stresses:
                      ! shear stress=0.55*normal stress and normal stress=7378Pa*down-dip distance
                      ! fault stresses rotated by 60 degree clockwise, shear stress positive
                      ! fault normal reference point at +y
                      !
                      EQN%IniBulk_yy(i,iBndGP)  = -2019.255513983125D0 * abs(average)/0.866025404D0
                      EQN%IniBulk_zz(i,iBndGP)  = -5358.744486016875D0 * abs(average)/0.866025404D0
                      EQN%IniShearYZ(i,iBndGP)  = 5223.717714560795D0 * abs(average)/0.866025404D0
                      ! Nucleation patch: set nucleation case to 0 in dyn file
                      ! shear stress=cohesion+(static fric + 0.0057)* initial normal stress
                      IF ((xGP.LE.1500.0D0) .AND. (xGP.GE.-1500.0D0)       &
                           .AND. (zGP.LE.-9093.26674D0) .AND. (zGP.GE.-11691.342951D0)) THEN
                         EQN%IniBulk_yy(i,iBndGP)  =  (-173205D0)+(-641.032721921598D0) * abs(average)/0.866025404D0
                         EQN%IniBulk_zz(i,iBndGP)  =  (173205D0)+(-6736.967278078402D0) * abs(average)/0.866025404D0
                         EQN%IniShearYZ(i,iBndGP)  =  (-100000D0)+(6019.435014560794D0) * abs(average)/0.866025404D0
                      ENDIF

                    CASE(13) !TPV12/13
                      ! depth dependent smooth initial stresses; fault normal reference point at +y
                      !
                      ! shear stress=0.549847*normal stress and normal stress=7390.01Pa*down-dip distance
                      ! fault stresses rotated by 30 degree clockwise, shear stress positive
                     IF (zGP.GE. -11951.15D0) THEN !down dip distance less than 13 800m
                         EQN%IniBulk_yy(i,iBndGP)  = -2023.521673446745D0 * abs(average)/0.866025404D0
                         EQN%IniBulk_zz(i,iBndGP)  = -5366.488326553255D0 * abs(average)/0.866025404D0
                         EQN%IniShearYZ(i,iBndGP)  = 5231.655611345521D0 * abs(average)/0.866025404D0
                       
                     
                      ! shear stress=0 and normal stress=14427.98Pa*down-dip distance
                      ! fault stresses rotated by 30 degree clockwise, shear stress positive
                     ELSE
                          EQN%IniBulk_yy(i,iBndGP)  =   -10820.985D0* abs(average)/0.866025404D0
                          EQN%IniBulk_zz(i,iBndGP)  =   -3606.995D0* abs(average)/0.866025404D0
                          EQN%IniShearYZ(i,iBndGP)  =   6247.49860264690D0* abs(average)/0.866025404D0
                     ENDIF
                ENDSELECT

            ENDDO ! iBndGP
      ! element wise stress assignment
      ELSE
          !
          !depth (negative in z)
          average = sum(zp(:))/3.0D0
 
          SELECT CASE(DISC%DynRup%BackgroundType)
              CASE(10) !TPV10
                ! depth dependent smooth initial stresses:
                ! shear stress=0.55*normal stress and normal stress=7378Pa*down-dip distance
                ! fault stresses rotated by 60 degree clockwise, shear stress positive
                ! fault normal reference point at +y
                !
                EQN%IniBulk_yy(i,:)  = -2019.255513983125D0 * abs(average)/0.866025404D0
                EQN%IniBulk_zz(i,:)  = -5358.744486016875D0 * abs(average)/0.866025404D0
                EQN%IniShearYZ(i,:)  = 5223.717714560795D0 * abs(average)/0.866025404D0
                ! Nucleation patch: set nucleation case to 0 in dyn file
                ! shear stress=cohesion+(static fric + 0.0057)* initial normal stress
                IF(   MAXVAL(xp(1:3)).LE.(1500.0D0) .AND. MINVAL(xp(1:3)).GE.(-1500.0D0)       &
                  .AND. MAXVAL(zp(1:3)).LE.(-9093.26674D0) .AND. MINVAL(zp(1:3)).GE.(-11691.342951D0)) THEN
                      EQN%IniBulk_yy(i,:)  =  (-173205D0)+(-641.032721921598D0) * abs(average)/0.866025404D0
                      EQN%IniBulk_zz(i,:)  =  (173205D0)+(-6736.967278078402D0) * abs(average)/0.866025404D0
                      EQN%IniShearYZ(i,:)  =  (-100000D0)+(6019.435014560794D0) * abs(average)/0.866025404D0
                ENDIF

              CASE(13) !TPV12/13
                ! depth dependent smooth initial stresses; fault normal reference point at +y
                ! shear stress=0.549847*normal stress and normal stress=7390.01Pa*down-dip distance
                ! fault stresses rotated by 30 degree clockwise, shear stress positive
                IF (zGP.GE. -11951.15D0) THEN !down dip distance less than 13 800m
                     EQN%IniBulk_yy(i,iBndGP)  = -2023.521673446745D0 * abs(average)/0.866025404D0
                     EQN%IniBulk_zz(i,iBndGP)  = -5366.488326553255D0 * abs(average)/0.866025404D0
                     EQN%IniShearYZ(i,iBndGP)  = 5231.655611345521D0 * abs(average)/0.866025404D0

                ! shear stress=0 and normal stress=14427.98Pa*down-dip distance
                ! fault stresses rotated by 30 degree clockwise, shear stress positive
               ELSE
                    EQN%IniBulk_yy(i,iBndGP)  =   -10820.985D0* abs(average)/0.866025404D0
                    EQN%IniBulk_zz(i,iBndGP)  =   -3606.995D0* abs(average)/0.866025404D0
                    EQN%IniShearYZ(i,iBndGP)  =   6247.49860264690D0* abs(average)/0.866025404D0
              ENDIF
          ENDSELECT
      ENDIF  
      
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_TPV10       ! SCEC TPV10  test

  !> SCEC TPV11 test case
  !<
  SUBROUTINE background_TPV11(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: VertexSide(4,3)
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: average
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  LOGICAL                        :: nodewise=.FALSE.
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! SCEC TPV11  test with a dipping fault and static friction 0.570
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/download/TPV10_11_Description_v7.pdf
  ! fault is in the xz-plane, dip in negative y-direction
  ! depth is z negative
  ! center of nucleation patch: x=0.0 km, depth= -12.5km, size 3km x 3km
  ! fault reaches the surface at depth 0km
  ! units of mesh in meter
              
  ! switch for Gauss node wise stress assignment
  nodewise = .TRUE.
           
  ! cohesion 0.2MPa for SCEC TPV10/11
  DISC%DynRup%cohesion = -0.2e6
                      
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !
  
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
     
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
      
      ! Gauss node coordinate definition and stress assignment
      IF (nodewise) THEN
          
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
          DO iBndGP = 1,DISC%Galerkin%nBndGP
              
              ! Transformation of boundary GP's into XYZ coordinate system
              chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
              tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
              CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
              CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
              !
              !geometrically averaged depth (negative in z)
              average = zGP   ! Averaging not needed here
              !
              ! depth dependent smooth initial stresses:
              ! shear stress=0.55*normal stress and normal stress=7378Pa*down-dip distance
              ! fault stresses rotated by 60 degree clockwise, shear stress positive
              ! fault normal reference point at +y
              !
              EQN%IniBulk_yy(i,iBndGP)  = -2019.255513983125D0 * abs(average)/0.866025404D0
              EQN%IniBulk_zz(i,iBndGP)  = -5358.744486016875D0 * abs(average)/0.866025404D0
              EQN%IniShearYZ(i,iBndGP)  = 5223.717714560795D0 * abs(average)/0.866025404D0
              ! Nucleation patch: set nucleation case to 0 in dyn file
              ! shear stress=cohesion+(static fric + 0.0057)* initial normal stress
              IF ((xGP.LE.1500.0D0) .AND. (xGP.GE.-1500.0D0)       &
                  .AND. (zGP.LE.-9093.26674D0) .AND. (zGP.GE.-11691.342951D0)) THEN
                  EQN%IniBulk_yy(i,iBndGP)  =  (-173205D0)+(-1855.044453454701D0) * abs(average)/0.866025404D0
                  EQN%IniBulk_zz(i,iBndGP)  =  (173205D0)+(-5522.955546545299D0) * abs(average)/0.866025404D0
                  EQN%IniShearYZ(i,iBndGP)  =  (-100000D0)+(5318.525014560794D0) * abs(average)/0.866025404D0
              ENDIF
          ENDDO ! iBndGP
      ! element wise stress assignment
      ELSE
          !
          !depth (negative in z)
          average = sum(zp(:))/3.0D0
          !
          ! depth dependent smooth initial stresses:
          ! shear stress=0.55*normal stress and normal stress=7378Pa*down-dip distance
          ! fault stresses rotated by 60 degree clockwise, shear stress positive
          ! fault normal reference point at +y
          !
          EQN%IniBulk_yy(i,:)  = -2019.255513983125D0 * abs(average)/0.866025404D0
          EQN%IniBulk_zz(i,:)  = -5358.744486016875D0 * abs(average)/0.866025404D0
          EQN%IniShearYZ(i,:)  = 5223.717714560795D0 * abs(average)/0.866025404D0
          ! Nucleation patch: set nucleation case to 0 in dyn file
          ! shear stress=cohesion+(static fric + 0.0057)* initial normal stress
          IF(   MAXVAL(xp(1:3)).LE.(1500.0D0) .AND. MINVAL(xp(1:3)).GE.(-1500.0D0)       &
              .AND. MAXVAL(zp(1:3)).LE.(-9093.26674D0) .AND. MINVAL(zp(1:3)).GE.(-11691.342951D0)) THEN
              EQN%IniBulk_yy(i,:)  =  (-173205D0)+(-1855.044453454701D0) * abs(average)/0.866025404D0
              EQN%IniBulk_zz(i,:)  =  (173205D0)+(-5522.955546545299D0) * abs(average)/0.866025404D0
              EQN%IniShearYZ(i,:)  =  (-100000D0)+(5318.525014560794D0) * abs(average)/0.866025404D0
          ENDIF
      ENDIF  
      
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE background_TPV11       ! SCEC TPV11  test

  !> Dipping fault with rate-and-state friction
  !<
  SUBROUTINE background_DIP_RSF(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: average,nstress,sstress,tmp
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Dipping fault and rate-and-state friction
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/download/TPV10_11_Description_v7.pdf
  ! fault is in the xz-plane, dip in negative y-direction
  ! depth is z negative
  ! center of nucleation patch: x=0.0 km, depth= -10.392km (12 km down-dip distance), size 3km x 3km
  ! Velocity strengthening layer at shallow depth (3km) where (a-b) abruptly changes from -0.004 to 0.004
  ! fault reaches the surface at depth 0km
  ! units of mesh in meter
  
  ALLOCATE(  DISC%DynRup%RS_a_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP)        )
  
  DISC%DynRup%cohesion = -1.0e6
   
  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)
             
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
      DISC%DynRup%RS_a_array(i,:) = DISC%DynRup%RS_a
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
      
         !geometrically averaged depth (negative in z)
          average = zGP  ! Averaging not needed here
          
          ! depth dependent smooth initial stresses:
          ! shear stress=0.55*normal stress and normal stress=7378Pa*down-dip distance
          ! fault stresses rotated by 60 degree clockwise, shear stress positive
          ! fault normal reference point at +y
          !
          EQN%IniBulk_yy(i,iBndGP)  = -2019.255513983125D0 * abs(average)/0.866025404D0
          EQN%IniBulk_zz(i,iBndGP)  = -5358.744486016875D0 * abs(average)/0.866025404D0
          EQN%IniShearYZ(i,iBndGP)  = 5223.717714560795D0 * abs(average)/0.866025404D0 

          ! effective normal stress
          nstress = 7378D0*abs(average)/0.866025404D0
                
          EQN%IniStateVar(i,iBndGP) = DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0*EXP((0.55D0*nstress/(nstress*DISC%DynRup%RS_b))-DISC%DynRup%RS_f0/DISC%DynRup%RS_b-DISC%DynRup%RS_a_array(i,iBndGP)/DISC%DynRup%RS_b*LOG(iniSlipRate/DISC%DynRup%RS_sr0))
          ! resulting changes in mu ini
          tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a_array(i,iBndGP))
          EQN%IniMu(i,iBndGP)=DISC%DynRup%RS_a_array(i,iBndGP) * LOG(tmp + SQRT(tmp**2 + 1.0))
          ! Nucleation patch: set nucleation case to 0 in dyn file
          ! shear stress=cohesion+(static fric + 0.7)* initial normal stress
          IF ((xGP.LE.1500.0D0) .AND. (xGP.GE.-1500.0D0)       &
              .AND. (zGP.LE.-9093.26674D0) .AND. (zGP.GE.-11691.342951D0)) THEN
              EQN%IniBulk_yy(i,iBndGP)  =  (-866025.4037844389D0)+(2772.89605785806D0) * abs(average)/0.866025404D0
              EQN%IniBulk_zz(i,iBndGP)  =  (866025.4037844389D0)+(-10150.89605785807D0) * abs(average)/0.866025404D0
              EQN%IniShearYZ(i,iBndGP)  =  (-499999.9999999999D0)+(7990.46771456079D0) * abs(average)/0.866025404D0

              ! effective normal stress
              nstress = 7378D0*abs(average)/0.866025404D0
              sstress = -1.0e6 + (0.6D0+0.8D0)*nstress
          
              ! resulting changes in SV_ini
              EQN%IniStateVar(i,iBndGP) = DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0*EXP((sstress/(nstress*DISC%DynRup%RS_b))-DISC%DynRup%RS_f0/DISC%DynRup%RS_b-DISC%DynRup%RS_a_array(i,iBndGP)/DISC%DynRup%RS_b*LOG(iniSlipRate/DISC%DynRup%RS_sr0))
                     
              !resulting changes in mu ini
              tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a_array(i,iBndGP))
              EQN%IniMu(i,iBndGP)=DISC%DynRup%RS_a_array(i,iBndGP) * LOG(tmp + SQRT(tmp**2 + 1.0))
          ENDIF

      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_DIP_RSF       ! Dipping fault with rate-and-state friction

  !> SCEC TPV14/15 test case branching faults
  !<
  SUBROUTINE background_TPV1415(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! SCEC TPV14/15 test with branching
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/download/TPV14_15_Description_v08.pdf .
  ! fault is in the xz-plane, branch is in '-y'
  ! center of nucleation patch is at x=0.0 (along strike) and  7.5km depth; size 3km x 3km
  ! fault reaches the surface at depth z=0km
  ! units of mesh in meter
  
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
  
      ! fault branch in direction '-y' has different ini stress 
      IF (sum(yp(1:3))/3.0D0 .LT. -0.1D0) THEN ! -0.1 to be safe
          ! normal stress as at main fault, but transformed (-120.0MPa)
          ! different sign convention than SCEC!
          ! TPV 14 right-lateral (sxx=0, sxy=+70MPa und syy=-120MPa)
          IF (DISC%DynRup%BackgroundType.EQ.14) THEN              
              EQN%IniBulk_xx(i,:)  =   30.6217782649107e6
              EQN%IniBulk_yy(i,:)  = -150.6217782649107e6
              EQN%IniShearXY(i,:)  =  -16.9615242270664e6
          ! TPV 15 left-lateral (sxx=0, sxy=-78 und syy=-120)
          ELSEIF (DISC%DynRup%BackgroundType.EQ.15) THEN  
              EQN%IniBulk_xx(i,:)  = -97.549981495186302e6
              EQN%IniBulk_yy(i,:)  = -22.450018504813706e6
              EQN%IniShearXY(i,:)  = -90.961524227066263e6              
          ENDIF
      ENDIF
       
  ENDDO !    MESH%Fault%nSide   
      
  END SUBROUTINE background_TPV1415       ! SCEC TPV14/15 test branching faults  

  !> SCEC TPV16/17 test case with heterogeneous initial stress field
  !<     
  SUBROUTINE background_TPV1617(EQN,MESH,IO,DISC,BND) 
  !-------------------------------------------------------------------------!
  USE read_backgroundstress_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tInputOutput)             :: IO
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND
  INTENT(INOUT) :: DISC,EQN, IO
  !-------------------------------------------------------------------------! 
  ! SCEC TPV16/17 test with heterogeneous stress background field
  ! Requires the correct mesh! 
  ! For instructions see http://scecdata.usc.edu/cvws/tpv16_17docs.html .
  ! fault is in the xz-plane
  ! fault reaches the surface at depth z = 0km (negative z is depth)
  ! origin at upper left corner of fault
  ! units of mesh in meter
  
  ! OPEN backgroundstress field
  CALL read_scec_stress(DISC,IO)          
  ALLOCATE(DISC%DynRup%forced_rupture_time(MESH%Fault%nSide,DISC%Galerkin%nBndGP))
 
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! constant background stress tensor and state variale
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
       
  ENDDO !    MESH%Fault%nSide   
 
  CALL extract_backgroundstress_2dplanefault(EQN,MESH,IO,DISC,BND)        
  
  END SUBROUTINE background_TPV1617 !  SCEC TPV16/17 with heterogeneous initial stress field     
   
  !> SCEC TPV29 test case : rough fault
  !> T. ULRICH 02.2015
  !> Landers1 used as a model
  !<
  SUBROUTINE background_TPV29 (DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: b11, b33, b13, Omega, g, Pf, zIncreasingCohesion
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! TPV29
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function (gravity)
  ! NOTE: z negative is depth, free surface is at z=0
  b11 = 1.025837D0
  b33 = 0.974162D0
  b13 =-0.158649D0
  g = 9.8D0    
  zIncreasingCohesion = -4000.
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      !EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
      
          ! for possible variation
          !DISC%DynRup%D_C(i,iBndGP)  = DISC%DynRup%D_C_ini
          !DISC%DynRup%Mu_S(i,iBndGP) = DISC%DynRup%Mu_S_ini
          !DISC%DynRup%Mu_D(i,iBndGP) = DISC%DynRup%Mu_D_ini
          !
          IF (zGP.GE.-17000.0D0) THEN
              Omega = 1D0
          ELSEIF (zGP.GE.-22000D0) THEN
              Omega = (zGP+22000D0)/5000D0
          ELSE
              Omega = 0D0
          ENDIF
          Pf = -1000D0 * g * zGP
          
          EQN%IniBulk_zz(i,iBndGP)  =  2670d0*g*zGP
          EQN%IniBulk_xx(i,iBndGP)  =  Omega*(b11*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniBulk_yy(i,iBndGP)  =  Omega*(b33*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniShearXY(i,iBndGP)  =  Omega*(b13*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniShearXZ(i,iBndGP)  =  0D0
          EQN%IniShearYZ(i,iBndGP)  =  0d0
          EQN%IniBulk_xx(i,iBndGP)  =  EQN%IniBulk_xx(i,iBndGP) + Pf
          EQN%IniBulk_yy(i,iBndGP)  =  EQN%IniBulk_yy(i,iBndGP) + Pf
          EQN%IniBulk_zz(i,iBndGP)  =  EQN%IniBulk_zz(i,iBndGP) + Pf
          !EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0

          ! manage cohesion
          IF (zGP.GE.zIncreasingCohesion) THEN
              ! higher cohesion near free surface
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6-0.0002d6*(zGP-zIncreasingCohesion)
          ELSE
              ! set cohesion
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6
          ENDIF
          
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_TPV29 

  !> SUMATRA test case
  !> T. ULRICH 06.2015
  !> tpv29 used as a model
  !<
  SUBROUTINE background_SUMATRA (DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: b11, b22, b12, b13, b23, Omega, g, Pf, zIncreasingCohesion
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! TPV29
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function (gravity)
  ! NOTE: z negative is depth, free surface is at z=0
  b11 = 1.1897
  b22 = 1.2179
  b12 = 0.0664
  b13 = 0.0764
  b23 = 0.0944

  g = 9.8D0    
  zIncreasingCohesion = -4000.
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      !EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
      
          ! for possible variation
          !DISC%DynRup%D_C(i,iBndGP)  = DISC%DynRup%D_C_ini
          !DISC%DynRup%Mu_S(i,iBndGP) = DISC%DynRup%Mu_S_ini
          !DISC%DynRup%Mu_D(i,iBndGP) = DISC%DynRup%Mu_D_ini
          !

          IF (zGP.GE.-45000.0D0) THEN
              Omega = 1D0
          ELSEIF (zGP.GE.-50000D0) THEN
              Omega = (zGP+50000D0)/5000D0
          ELSE
              Omega = 0D0
          ENDIF

          Pf = -1000D0 * g * zGP

          EQN%IniBulk_zz(i,iBndGP)  =  2670d0*g*zGP
          EQN%IniBulk_xx(i,iBndGP)  =  Omega*(b11*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniBulk_yy(i,iBndGP)  =  Omega*(b22*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniShearXY(i,iBndGP)  =  Omega*(b12*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniShearXZ(i,iBndGP)  =  Omega*(b13*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniShearYZ(i,iBndGP)  =  Omega*(b23*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniBulk_xx(i,iBndGP)  =  EQN%IniBulk_xx(i,iBndGP) + Pf
          EQN%IniBulk_yy(i,iBndGP)  =  EQN%IniBulk_yy(i,iBndGP) + Pf
          EQN%IniBulk_zz(i,iBndGP)  =  EQN%IniBulk_zz(i,iBndGP) + Pf

          ! manage cohesion
          IF (zGP.GE.zIncreasingCohesion) THEN
              ! higher cohesion near free surface
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6-0.0002d6*(zGP-zIncreasingCohesion)
          ELSE
              ! set cohesion
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6
          ENDIF
          
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_SUMATRA

  !> SCEC TPV3132 test case : strike slip rupture in layered medium
  !> T. ULRICH 02.2015
  !<
  SUBROUTINE background_TPV3132 (DISC,EQN,MESH,BND)
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
  LOGICAL                        :: FoundLayer
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  INTEGER                        :: iLayer
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: x, z, zb
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: mu, mu0
  REAL                           :: xHypo, zHypo, r
  REAL, parameter :: PI = 4 * atan (1.0d0)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! TPV29
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function (gravity)
  ! NOTE: z negative is depth, free surface is at z=0
  logError(*) 'initialization Gauss wise'
  mu0=32.03812032d9
  xHypo = 0d0
  zHypo = -7500d0  
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      !EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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
          !zb = SUM(zV(1:4))/4.0D0      

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
      
          z = zGP  
          x = xGP
          !get mu value at GP (more precise that getting directly mu element) 
           DO iLayer = 1, EQN%nLayers - 1
             IF( (z.LE.EQN%MODEL(iLayer,1)).AND.(z.GT.EQN%MODEL(iLayer+1,1)) ) THEN
               FoundLayer = .TRUE.
               mu = EQN%MODEL(iLayer,3) + (EQN%MODEL(iLayer+1,3) - EQN%MODEL(iLayer,3)) /  &
                                                                (EQN%MODEL(iLayer+1,1)   - EQN%MODEL(iLayer,1)  ) *  &
                                                                ( z - EQN%MODEL(iLayer,1) )
               EXIT
             ENDIF
           ENDDO   
           IF(.NOT.FoundLayer) THEN
             logError(*) 'Layered Medium Error: No layer found for depth', z
             STOP
           ENDIF      
          EQN%IniBulk_xx(i,iBndGP)  = -60d6*mu/mu0
          EQN%IniBulk_yy(i,iBndGP)  = -60d6*mu/mu0
          EQN%IniBulk_zz(i,iBndGP)  =  0d0
          EQN%IniShearXY(i,iBndGP)  =  30d6*mu/mu0
          EQN%IniShearXZ(i,iBndGP)  =  0D0
          EQN%IniShearYZ(i,iBndGP)  =  0d0

          ! distance to hypocenter (approx plane fault)
          r = sqrt( ((x-xHypo)*(x-xHypo))+((z-zHypo)*(z-zHypo)))

          IF (r.LE.1400D0) THEN
             EQN%IniShearXY(i,iBndGP)  =  EQN%IniShearXY(i,iBndGP)+4.95d6*mu/mu0
          ELSEIF (r.LE.2000D0) THEN
             EQN%IniShearXY(i,iBndGP)  =  EQN%IniShearXY(i,iBndGP)+2.475d6*(mu/mu0)*(1d0+dcos(PI*(r-1400d0)/600d0))
          ENDIF
          ! manage cohesion
          IF (z.GE.-2400.0D0) THEN
              ! higher cohesion near free surface
              DISC%DynRup%cohesion(i,iBndGP) = -0.000425d6*(2400d0+z)
          ELSE
              ! set cohesion
              DISC%DynRup%cohesion(i,iBndGP) = 0d0
          ENDIF
                
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_TPV3132       


  !> Tohoku1 backround stress model
  !<
  SUBROUTINE background_TOH1(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j,circle
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: dcpatch_radius, dcpatch_distance, dcpatch_center(2,4)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Requires the correct mesh!
  ! center of nucleation patch is at -52.67794752km/-409.5199856km depth; size 20km x 20km
  ! units of mesh are in kilometer but are scaled by 1000. -> all parameters here in m !
  ! topography and depth values: z positive with depth, z = 0 -> sea level?
  ! ini D_C is specified in the dyn file
  
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! constant background stress tensor and state variale 
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
      EQN%IniStateVar(i,:) =  EQN%RS_sv0 
  
      ! D_C depth dependence
      IF (MINVAL(zp(1:3)) .GT. 20000.0D0) THEN
          DISC%DynRup%D_C(:,:) = 2.0D0
      ENDIF
            
      ! 4 circles of lower D_C dependent on the lateral directions x and y
      dcpatch_radius = 10000.0D0
      dcpatch_center(1,1) = -205000.0D0
      dcpatch_center(2,1) = -480000.0D0
      dcpatch_center(1,2) = -205000.0D0
      dcpatch_center(2,2) = -472000.0D0
      dcpatch_center(1,3) = -147000.0D0
      dcpatch_center(2,3) = -409000.5D0
      dcpatch_center(1,4) = -137000.0D0
      dcpatch_center(2,4) = -326000.0D0
            
      DO circle = 1,4
          dcpatch_distance = SQRT( (MAXVAL(xp(1:3)) - dcpatch_center(1,circle))**2 + (MAXVAL(yp(1:3)) - dcpatch_center(2,circle))**2 )
          IF (dcpatch_distance .LT. dcpatch_radius) THEN
              DISC%DynRup%D_C(i,:) = 1.0D0
          ENDIF
      ENDDO
           
  ENDDO !    MESH%Fault%nSide   
  END SUBROUTINE background_TOH1       ! Tohoku 1

  
  !> Landers1 fault system backround stress model
  !<
  SUBROUTINE background_LAN1(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: z
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Landers
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function with 3 layers
  ! NOTE: z negative is depth, free surface is at z=+
             
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
      
          z = zGP  
          
          ! for possible variation
          !DISC%DynRup%D_C(i,iBndGP)  = DISC%DynRup%D_C_ini
          !DISC%DynRup%Mu_S(i,iBndGP) = DISC%DynRup%Mu_S_ini
          !DISC%DynRup%Mu_D(i,iBndGP) = DISC%DynRup%Mu_D_ini
          !
          ! N11 , grad2:
          !
          ! highest fault point is appr at z=1390m, offset 2000 is obtained by testing to shift the gradient up
          !
          EQN%IniBulk_xx(i,iBndGP)  =  -10.6723e6*(abs(z-2000.0D0))/1000.0D0
          EQN%IniBulk_yy(i,iBndGP)  =  -29.3277e6*(abs(z-2000.0D0))/1000.0D0
          EQN%IniBulk_zz(i,iBndGP)  =  -20.0000e6*(abs(z-2000.0D0))/1000.0D0
          EQN%IniShearXY(i,iBndGP)  =   -3.7687e6*(abs(z-2000.0D0))/1000.0D0
          EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
          EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
          EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0

          ! manage D_C
          IF (z.GT.-4000.0D0) THEN
              ! higher D_C to surpress supershear rupture at free surface
              DISC%DynRup%D_C(i,iBndGP) = DISC%DynRup%D_C_ini+0.6D0*(1.0D0+COS(4.0D0*ATAN(1.0D0) * abs(z)/4000.0D0))
          ELSEIF (z.LT.-12000.0D0) THEN
              ! higher D_C in depth mimic brittle ductile transition
              DISC%DynRup%D_C(i,iBndGP) = DISC%DynRup%D_C_ini+1.0D0*(1.0D0+COS(4.0D0*ATAN(1.0D0) * abs(z)/4000.0D0))
          ENDIF

          ! overwrite positive z area
          IF(z .GT. 0.0) THEN
              !EQN%IniBulk_xx(i,iBndGP)  =  0.0D0
              !EQN%IniBulk_yy(i,iBndGP)  =  0.0D0
              !EQN%IniBulk_zz(i,iBndGP)  =  0.0D0
              !EQN%IniShearXY(i,iBndGP)  =  0.0D0
              !EQN%IniShearYZ(i,iBndGP)  =  0.0D0
              !EQN%IniShearXZ(i,iBndGP)  =  0.0D0
              !EQN%IniStateVar(i,iBndGP) =  0.0D0

              DISC%DynRup%D_C(i,iBndGP) =  2.0D0
          ENDIF
          ! set cohesion
          DISC%DynRup%cohesion(i,iBndGP) = -2.0e6
      
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_LAN1       ! Landers 1

  !> Landers2 segmented fault system backround stress model
  !<  
  SUBROUTINE background_LAN2(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  USE COMMON_operators_mod, ONLY: OpenFile, XYinTriangle
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: x,y,z
  REAL                           :: x7(3),y7(3),x8(3),y8(3),x9(3),y9(3),x10(3),y10(3),x11(3),y11(3) ! stores coordinates of enclosing triangulars for fault segmentation
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: epsilon                                ! tolerance needed by function XYinTriangle
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Landers - with individual segments
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function with 3 layers
  ! NOTE: z negative is depth, free surface is at z=+

  ! Define segment triangulars:
  ! segment 7 (last point far right of the fault)
  x7 = (/-7900.0, -7680.0, -7000.0/)
  y7 = (/30010.0, 27070.0, 28500.0/)

  ! segment 8
  x8 = (/-7900.0, -2500.0, -2500.0/)
  y8 = (/30010.0, 22000.0, 30010.0/)

  ! segment 9
  x9 = (/-7900.0, -7900.0, -14500.0/)
  y9 = (/29000.0, 37720.0, 37720.0/)

  ! segment 10
  x10 = (/-13770.0, -19700.0, -13500.0/)
  y10 = (/37720.0, 43500.0, 43500.0/)

  ! segment 11
  x11 = (/-13770.0, -23300.0, -23300.0/)
  y11 = (/37720.0, 46700.0, 37720.0/)

  epsilon = 1.0e-3 ! to test if point lies in tri
             
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
         
          x = xGP
          y = yGP
          z = zGP  
          
          ! ini background
          ! N7E
          EQN%IniBulk_xx(i,iBndGP)  =  -8.6774e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniBulk_yy(i,iBndGP)  =  -24.8426e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniShearXY(i,iBndGP)  =  -2.0152e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
          EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
          EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0
          !
          ! ATTENTION: Keep order as is!
          !
          ! Segment #1
          IF (y .LE. 8600.0) THEN

              EQN%IniStateVar(i,iBndGP) =  1 ! check EQN%RS_sv0

          ! Segment #2
          ELSEIF ((x.LE. -1800.0).AND.(x.GE.-7500.0) &
              .AND.(y.LE.16500.0).AND.(y.GE.8600.0)) THEN

              EQN%IniStateVar(i,iBndGP)  =  2 ! check

          ! Segment #3
          ELSEIF ((x.LE. -1220.0).AND.(x.GE.-1800.0) &
              .AND.(y.LE.14500.0).AND.(y.GE.8600.0)) THEN

              EQN%IniStateVar(i,iBndGP)  = 3 ! check

          ! Segment #4
          ELSEIF ((x.LE. +1200.0).AND.(x.GE.-1220.0) &
              .AND.(y.LE.14500.0).AND.(y.GE.9800.0)) THEN

              EQN%IniStateVar(i,iBndGP)  = 4 ! check
           
          ! Segment #5 + #6
          ! first set large chunk to 5/6 and then overwrite it eventually with branch parameters 7,8,9
          ELSEIF ((x.LE.-1220.0).AND.(x.GE.-10800.0) &
              .AND.(y.LE.37720.0).AND.(y.GE.14500.0)) THEN
              
              EQN%IniBulk_xx(i,iBndGP)  =  -13.3719e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -20.148e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -7.6098e6*(abs(z-1390.0D0))/1000.0D0

    
              ! Segment #7
              IF (XYinTriangle(x,y,x7,y7,epsilon)) THEN
   
                  EQN%IniStateVar(i,iBndGP)  = 7 ! check

              ! Segment #8
              ELSEIF (XYinTriangle(x,y,x8,y8,epsilon)) THEN

                  EQN%IniStateVar(i,iBndGP)  = 8 ! check

              ! Segment #9
              ELSEIF (XYinTriangle(x,y,x9,y9,epsilon)) THEN

                  EQN%IniStateVar(i,iBndGP)  = 9 ! check
              ENDIF

          ! Segment #9
          ! missing 9 pieces
          ELSEIF (XYinTriangle(x,y,x9,y9,epsilon)) THEN

              EQN%IniStateVar(i,iBndGP)  = 9 ! check

          ! Segment #10
          ELSEIF (XYinTriangle(x,y,x10,y10,epsilon)) THEN

              EQN%IniStateVar(i,iBndGP)  = 10 ! check

          ! Segment #11
          ELSEIF (XYinTriangle(x,y,x11,y11,epsilon)) THEN

              EQN%IniStateVar(i,iBndGP)  = 11 ! check
          ! Segment #12
          ELSEIF (y.GE. 46700.0) THEN

              EQN%IniStateVar(i,iBndGP)  = 12 ! check

          ENDIF

                  ! manage D_C
          IF (z.GT.-4000.0D0) THEN
              ! higher D_C to surpress supershear rupture at free surface
              DISC%DynRup%D_C(i,iBndGP) = 0.8D0+0.6D0*(1.0D0+COS(4.0D0*ATAN(1.0D0) * abs(z)/4000.0D0))
          !ELSEIF (z.LT.-12000.0D0) THEN
                  !   ! higher D_C in depth mimic brittle ductile transition
                  !   DISC%DynRup%D_C(i,iBndGP) = 0.8D0+10e-3*(abs(z)-12000.0D0)
          ENDIF

          ! overwrite positive z area
          IF(z .GT. 0.0) THEN
              !EQN%IniBulk_xx(i,iBndGP)  =  0.0D0
              !EQN%IniBulk_yy(i,iBndGP)  =  0.0D0
              !EQN%IniBulk_zz(i,iBndGP)  =  0.0D0
              !EQN%IniShearXY(i,iBndGP)  =  0.0D0
              !EQN%IniShearYZ(i,iBndGP)  =  0.0D0
              !EQN%IniShearXZ(i,iBndGP)  =  0.0D0
              !EQN%IniStateVar(i,iBndGP) =  0.0D0

              DISC%DynRup%D_C(i,iBndGP) =  2.0D0
          ENDIF
                
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_LAN2       ! Landers 2 segmented fault system
  
  !> Landers3 segmented fault system backround stress model
  !<  
  SUBROUTINE background_LAN3(DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  USE COMMON_operators_mod, ONLY: OpenFile, XYinTriangle
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  TYPE (tBoundary)               :: BND
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: x,y,z
  REAL                           :: x7(3),y7(3),x8(3),y8(3),x9(3),y9(3),x10(3),y10(3),x11(3),y11(3) ! stores coordinates of enclosing triangulars for fault segmentation
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                         :: epsilon                                ! tolerance needed by function XYinTriangle
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Landers - with individual segments
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function with 3 layers
  ! NOTE: z negative is depth, free surface is at z=+

  ! difference to case 61: segment 5 is splitted into 2. splitting point is the curve, where the rupture usually stops.

  ! Define segment triangulars:
  ! segment 7 (last point far right of the fault)
  x7 = (/-7900.0, -7680.0, -7000.0/)
  y7 = (/30010.0, 27070.0, 28500.0/)

  ! segment 8
  x8 = (/-7900.0, -2500.0, -2500.0/)
  y8 = (/30010.0, 22000.0, 30010.0/)

  ! segment 9
  x9 = (/-7900.0, -7900.0, -14500.0/)
  y9 = (/29000.0, 37720.0, 37720.0/)

  ! segment 10
  x10 = (/-13770.0, -19700.0, -13500.0/)
  y10 = (/37720.0, 43500.0, 43500.0/)

  ! segment 11
  x11 = (/-13770.0, -23300.0, -23300.0/)
  y11 = (/37720.0, 46700.0, 37720.0/)

  epsilon = 1.0e-3 ! to test if point lies in tri
               
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
         
          x = xGP
          y = yGP
          z = zGP  
          
                    ! ini background
          ! N22
          EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
          EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
          EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
          EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0 

          ! ATTENTION: Keep order as is!
          !
          ! Segment #1
          IF (y .LE. 8600.0) THEN

              EQN%IniStateVar(i,iBndGP) =  1 ! check EQN%RS_sv0
              ! 22NE
              EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0

              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0 


          ! Segment #2
          ELSEIF ((x.LE. -1800.0).AND.(x.GE.-7500.0) &
              .AND.(y.LE.16500.0).AND.(y.GE.8600.0)) THEN

              EQN%IniStateVar(i,iBndGP)  =  2 ! check
              !N22
              EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0

              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0 


          ! Segment #3
          ELSEIF ((x.LE. -1220.0).AND.(x.GE.-1800.0) &
              .AND.(y.LE.14500.0).AND.(y.GE.8600.0)) THEN

              EQN%IniStateVar(i,iBndGP)  = 3 ! check
                       ! 11NE 
              EQN%IniBulk_xx(i,iBndGP)  =  -9.0366*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -24.4834*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -3.1205*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0 

          ! Segment #4
          ELSEIF ((x.LE. +1200.0).AND.(x.GE.-1220.0) &
              .AND.(y.LE.14500.0).AND.(y.GE.9800.0)) THEN

              EQN%IniStateVar(i,iBndGP)  = 4 ! check
                              ! 11NE 
              EQN%IniBulk_xx(i,iBndGP)  =  -9.0366*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -24.4834*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -3.1205*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0

                
          ! Segment #5 (small, until critical curve)
          ELSEIF ((x.LE.-1220.0).AND.(x.GE.-1800.0) &
              .AND.(y.LE.16600.0).AND.(y.GE.14500.0)) THEN

              EQN%IniStateVar(i,iBndGP)  = 5 ! check
                         ! 11NE 
              EQN%IniBulk_xx(i,iBndGP)  =  -9.0366*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -24.4834*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -3.1205*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0

          ! Segment #6
          ! first set large chunk to 5/6 and then overwrite it eventually with branch parameters 7,8,9
          ELSEIF ((x.LE.-1800.0).AND.(x.GE.-10800.0) &
              .AND.(y.LE.37720.0).AND.(y.GE.16600.0)) THEN

              EQN%IniStateVar(i,iBndGP)  = 6 ! check   
                                 ! 11NE 
              EQN%IniBulk_xx(i,iBndGP)  =  -9.0366*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -24.4834*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -3.1205*(abs(z)-1390.0D0)/1000.0D0
              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0


              ! Segment #7
              IF (XYinTriangle(x,y,x7,y7,epsilon)) THEN
   
                  EQN%IniStateVar(i,iBndGP)  = 7 ! check
                  ! N22
                  EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
                  EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
                  EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0 
              ! Segment #8
              ELSEIF (XYinTriangle(x,y,x8,y8,epsilon)) THEN

                  EQN%IniStateVar(i,iBndGP)  = 8 ! check             
                  ! N22
                  EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
                  EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0

              ! Segment #9
              ELSEIF (XYinTriangle(x,y,x9,y9,epsilon)) THEN

                  EQN%IniStateVar(i,iBndGP)  = 9 ! check
                  ! N22
                  EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
                  EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
                  EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
                  EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0
              ENDIF

          ! Segment #9
          ! missing 9 pieces
          ELSEIF (XYinTriangle(x,y,x9,y9,epsilon)) THEN

              EQN%IniStateVar(i,iBndGP)  = 9 ! check
              ! N22
              EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0


          ! Segment #10
          ELSEIF (XYinTriangle(x,y,x10,y10,epsilon)) THEN

              EQN%IniStateVar(i,iBndGP)  = 10 ! check
              ! N22
              EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0


          ! Segment #11
          ELSEIF (XYinTriangle(x,y,x11,y11,epsilon)) THEN

              EQN%IniStateVar(i,iBndGP)  = 11 ! check
              ! N22
              EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearYZ(i,iBndGP)  =  EQN%ShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0

    
   
          ! Segment #12
          ELSEIF (y.GE. 46700.0) THEN

              EQN%IniStateVar(i,iBndGP)  = 12 ! check
              ! N22
              EQN%IniBulk_xx(i,iBndGP)  =  -10.766e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_yy(i,iBndGP)  =  -22.754e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniBulk_zz(i,iBndGP)  =  -16.76e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXY(i,iBndGP)  =  -5.788e6*(abs(z-1390.0D0))/1000.0D0
              EQN%IniShearXZ(i,iBndGP)  =  EQN%ShearXZ_0
              EQN%IniStateVar(i,iBndGP) =  EQN%RS_sv0


          ENDIF

                  ! manage D_C
          IF (z.GT.-4000.0D0) THEN
              ! higher D_C to surpress supershear rupture at free surface
              DISC%DynRup%D_C(i,iBndGP) = 0.8D0+0.6D0*(1.0D0+COS(4.0D0*ATAN(1.0D0) * abs(z)/4000.0D0))
          !ELSEIF (z.LT.-12000.0D0) THEN
                  !   ! higher D_C in depth mimic brittle ductile transition
                  !   DISC%DynRup%D_C(i,iBndGP) = 0.8D0+10e-3*(abs(z)-12000.0D0)
          ENDIF

          ! overwrite positive z area
          IF(z .GT. 0.0) THEN
              !EQN%IniBulk_xx(i,iBndGP)  =  0.0D0
              !EQN%IniBulk_yy(i,iBndGP)  =  0.0D0
              !EQN%IniBulk_zz(i,iBndGP)  =  0.0D0
              !EQN%IniShearXY(i,iBndGP)  =  0.0D0
              !EQN%IniShearYZ(i,iBndGP)  =  0.0D0
              !EQN%IniShearXZ(i,iBndGP)  =  0.0D0
              !EQN%IniStateVar(i,iBndGP) =  0.0D0

              DISC%DynRup%D_C(i,iBndGP) = 2.0 ! 2.0D0
          ENDIF
          !
      ENDDO ! iBndGP
                          
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_LAN3       ! Landers 3 segmented fault system
  
  !> Alaska dipping fault backround stress model
  !<  
  SUBROUTINE background_ALA(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Alaska1
  ! Requires the correct mesh: alaska1.neu
  ! fault: dip at 10 deg, strike: parallel to x
  ! dipping fault points towards y, trench at y=z=0
  ! depth is z negative
  ! center of nucleation patch at hypocenter: 22.38km along strike, -89.26km along dip
  ! fault does reach the surface
  ! units of mesh in meter

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
         
          ! depth dependent smooth initial stresses:
          ! fault normal reference point at -y
          !
          EQN%IniBulk_yy(i,iBndGP)  = -2403.380382456552D0 * abs(zGP)/0.173648177D0
          EQN%IniBulk_zz(i,iBndGP)  = -2974.619617543449D0 * abs(zGP)/0.173648177D0
          EQN%IniShearYZ(i,iBndGP)  =  3269.892310776356D0 * abs(zGP)/0.173648177D0
          !
          ! Depth dependent stress values in nucleation patch -> set nucleation case to 0 in dyn file!
          ! shear stress=(static fric + 0.0057)* initial normal stress
          IF ((xGP.LE.(10000.0D0+22380.0D0)) .AND. (xGP.GE.(-10000.0D0+22380.0D0))       &
              .AND. (zGP.LE.(10000.0D0+87904.0D0)) .AND. (zGP.GE.(-10000.0D0+87904.0D0))) THEN
              EQN%IniBulk_yy(i,iBndGP)  =  (-641.032721921598D0) * abs(zGP)/0.173648177D0
              EQN%IniBulk_zz(i,iBndGP)  =  (-6736.967278078402D0) * abs(zGP)/0.173648177D0
              EQN%IniShearYZ(i,iBndGP)  =  (6019.435014560794D0) * abs(zGP)/0.173648177D0
          ENDIF
         
      ENDDO ! iBndGP
                          
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_ALA       ! Alaska dipping fault 

  !> Northridge dipping fault backround stress model
  !<  
  SUBROUTINE background_NORTH(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: average
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Northridge Thrust Fault Scenario based on SCEC TPV10
  ! Requires the correct mesh! Build from SCEC Community fault model
  ! fault: dip at 35.5 deg, strike: 117.82130071543423
  ! depth is z negative
  ! center of nucleation patch at hypocenter of 1994 EQ
  ! fault does not reach the surface
  ! units of mesh in meter

  ! higher cohesion than for SCEC TPV10/11
  DISC%DynRup%cohesion = -20e6
            
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
         
                    !geometrically averaged depth (negative in z)
          average = zGP   ! Averaging not needed here
          !
          ! depth dependent smooth initial stresses:
          ! fault normal reference point at +y
          !
          EQN%IniBulk_yy(i,iBndGP)  = -2403.380382456552D0 * abs(average)/0.380801179408141D0
          EQN%IniBulk_zz(i,iBndGP)  = -2974.619617543449D0 * abs(average)/0.380801179408141D0
          EQN%IniShearYZ(i,iBndGP)  =  3269.892310776356D0 * abs(average)/0.380801179408141D0
          ! Nucleation patch: set nucleation case to 0 in dyn file
          ! shear stress=cohesion+(static fric + 0.0057)* initial normal stress
          IF ((xGP.LE.(2500.0D0+357208.0D0)) .AND. (xGP.GE.(-2500.0D0+357208.0D0))       &
              .AND. (zGP.LE.-17250D0) .AND. (zGP.GE.-18250D0)) THEN
              EQN%IniBulk_yy(i,iBndGP)  =  (-17320500D0)+(-641.032721921598D0) * abs(average)/0.380801179408141D0
              EQN%IniBulk_zz(i,iBndGP)  =  (17320500D0)+(-6736.967278078402D0) * abs(average)/0.380801179408141D0
              EQN%IniShearYZ(i,iBndGP)  =  (-10000000D0)+(6019.435014560794D0) * abs(average)/0.380801179408141D0
          ENDIF
         
      ENDDO ! iBndGP
                          
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_NORTH       !  Northridge Thrust Fault Scenario based on SCEC TPV10

  !> SCEC TPV101 test with rate-and-state friction (ageing law) 
  !<  
  SUBROUTINE background_TPV101(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate, tmp, Boxx, Boxz
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! 101: SCEC TPV101 test with rate-and-state friction (ageing law) in full space
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/download/SCEC_validation_ageing_law.pdf
  ! vertical fault in the xz-plane,30 km long by 15 km deep, surrounded by a transition region 3 km thick
  ! depth is z negative, hypocenter centered along-strike and along-dip
  ! nucleation is time-dependent shear stress perturbation in a circular region with radius of 3 km
  ! fault reaches the surface at depth 0km
  ! units of mesh in meter
  ! Gauss node wise stress and friction properties assignment
  ! fault normal reference point at +y

  ALLOCATE(  DISC%DynRup%RS_a_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP)        )
 
  ! cohesion 0MPa
  DISC%DynRup%cohesion = 0.0
            
  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)
              
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
      DISC%DynRup%RS_a_array(i,:) = DISC%DynRup%RS_a
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
         
                              ! friction law changes in a
          ! smoothed Boxcar function in transition region (3km)
          IF ( ((xGP.GT.15000.0D0).AND.(xGP.LT.18000.0D0))      &
              .OR.((xGP.LT.-15000.0D0).AND.(xGP.GT.-18000.0D0)) &
              .OR.((zGP.LT.-15000.0D0).AND.(zGP.GT.-18000.0D0))) THEN
              ! TANH(X)=(1-EXP(-2X))/(1+EXP(-2X))
              tmp = 3000.0D0/(ABS(xGP)-15000.0D0-3000.0D0)+3000.0D0/(ABS(xGP)-15000.0D0)
              ! x
              Boxx = 0.5D0*(1.0D0+((1.0D0-EXP(-2.0D0*tmp))/(1.0D0+EXP(-2.0D0*tmp))))
              ! z
              tmp = 3000.0D0/(ABS(zGP)-3000.0D0)+3000.0D0/ABS(zGP)
              Boxz = 0.5D0*(1.0D0+((1.0D0-EXP(-2.0D0*tmp))/(1.0D0+EXP(-2.0D0*tmp))))
              ! smooth boxcar
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a+0.008D0*(1.0D0-(Boxx*Boxz))
          ELSEIF ((xGP.GE.18000.0D0) .OR. (xGP.LE.-18000.0D0)        &
              .OR. (zGP.LE.-18000.0D0)) THEN
              ! velocity strengthening in exterior region (3km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a+0.008D0
              !
          ELSE
              ! velocity-weakening in interior fault region (30 km * 15 km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a
              !
          ENDIF
          ! resulting changes in SV_ini
          EQN%IniStateVar(i,iBndGP) = (DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0) * EXP((-EQN%ShearXY_0/EQN%Bulk_yy_0-DISC%DynRup%RS_f0-DISC%DynRup%RS_a_array(i,iBndGP)*LOG(iniSlipRate/DISC%DynRup%RS_sr0))/DISC%DynRup%RS_b)
                ! Nucleation in Evaluate friction special case
                    
      ENDDO ! iBndGP
                          
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_TPV101       !  SCEC TPV101 test with rate-and-state friction (ageing law) 

  !> SCEC TPV103 test with velocity weakening friction (based on slip law) 
  !<  
  SUBROUTINE background_TPV103(DISC,EQN,MESH,BND)
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
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate, tmp, Boxx, Boxz
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! 103: SCEC TPV103 test with velocity weakening friction (based on slip law) in full space
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/download/SCEC_validation_slip_law.pdf
  ! vertical velocity-weakening fault in the xz-plane, 30 km long by 15 km deep,
  ! surrounded by a transition region 3 km thick, embedded in velocity-strengthening region
  ! depth is z negative, hypocenter centered along-strike and along-dip
  ! nucleation is time- and space dependent shear stress perturbation in a circular region with max. radius of 3 km
  ! fault reaches the surface at depth 0km
  ! units of mesh in meter
  ! Gauss node wise stress and friction properties assignment
  ! fault normal reference point at +y  

  ALLOCATE(  DISC%DynRup%RS_a_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP)        )
  ALLOCATE(  DISC%DynRup%RS_srW_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )

  ! cohesion 0MPa          
  DISC%DynRup%cohesion = 0.0
            
  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)
              
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! switch for rupture front output: RF
      IF (DISC%DynRup%RF_output_on == 1) THEN
          ! rupture front output just for + side elements!
          IF (MESH%FAULT%Face(i,1,1) .NE. 0) DISC%DynRup%RF(i,:) = .TRUE.
      ENDIF
      
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)  
      
      EQN%IniBulk_xx(i,:)  =  EQN%Bulk_xx_0
      EQN%IniBulk_yy(i,:)  =  EQN%Bulk_yy_0
      EQN%IniBulk_zz(i,:)  =  EQN%Bulk_zz_0
      EQN%IniShearXY(i,:)  =  EQN%ShearXY_0
      EQN%IniShearYZ(i,:)  =  EQN%ShearYZ_0
      EQN%IniShearXZ(i,:)  =  EQN%ShearXZ_0
            
      ! ini frictional parameters
      EQN%IniStateVar(i,:) =  EQN%RS_sv0
      DISC%DynRup%RS_a_array(i,:) = DISC%DynRup%RS_a
                
      ! Gauss node coordinate definition and stress assignment
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

      DO iBndGP = 1,DISC%Galerkin%nBndGP
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
          ! friction law changes in a
          ! smoothed Boxcar function in transition region (3km)
          IF ( ((xGP.GT.15000.0D0).AND.(xGP.LT.18000.0D0))      &
              .OR.((xGP.LT.-15000.0D0).AND.(xGP.GT.-18000.0D0)) &
              .OR.((zGP.LT.-15000.0D0).AND.(zGP.GT.-18000.0D0))) THEN
              ! TANH(X)=(1-EXP(-2X))/(1+EXP(-2X))
              tmp = 3000.0D0/(ABS(xGP)-15000.0D0-3000.0D0)+3000.0D0/(ABS(xGP)-15000.0D0)
              ! x
              Boxx = 0.5D0*(1.0D0+((1.0D0-EXP(-2.0D0*tmp))/(1.0D0+EXP(-2.0D0*tmp))))
              ! z
              tmp = 3000.0D0/(ABS(zGP)-3000.0D0)+3000.0D0/ABS(zGP)
              Boxz = 0.5D0*(1.0D0+((1.0D0-EXP(-2.0D0*tmp))/(1.0D0+EXP(-2.0D0*tmp))))
              ! smooth boxcar
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a+0.01D0*(1.0D0-(Boxx*Boxz))
              DISC%DynRup%RS_srW_array(i,iBndGP)=DISC%DynRup%RS_srW+0.9D0*(1.0D0-(Boxx*Boxz))
          ELSEIF ((xGP.GE.18000.0D0) .OR. (xGP.LE.-18000.0D0)        &
              .OR. (zGP.LE.-18000.0D0)) THEN
              ! velocity strengthening in exterior region (3km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a+0.01D0
              DISC%DynRup%RS_srW_array(i,iBndGP) = DISC%DynRup%RS_srW+0.9D0
              !
          ELSE
              ! velocity-weakening in interior fault region (30 km * 15 km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a
              DISC%DynRup%RS_srW_array(i,iBndGP) = DISC%DynRup%RS_srW
              !
          ENDIF
          ! resulting changes in SV_ini
          EQN%IniStateVar(i,iBndGP) = (DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0) * EXP((-EQN%ShearXY_0/EQN%Bulk_yy_0-DISC%DynRup%RS_f0-DISC%DynRup%RS_a_array(i,iBndGP)*LOG(iniSlipRate/DISC%DynRup%RS_sr0))/DISC%DynRup%RS_b)
      ! Nucleation in Evaluate friction special case  
      ENDDO ! iBndGP
                          
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_TPV103       !  SCEC TPV103 test with velocity weakening friction (based on slip law)     
    
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  
  !>  Nucleation by discontinuous jump on properties at [NucXMin,NucXMax] x [NucYMin,NucYMax]  
  !<
  SUBROUTINE nucleation_STEP(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: LocX(3), LocY(3)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Discontinuous jump on properties at [NucXMin,NucXMax] x [NucYMin,NucYMax]

  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
      
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
     
      ! choose 2D fault plane for nucleation
      SELECT CASE(DISC%DynRup%NucDirX)
          CASE(1) !x direction
              LocX(:)=xp(1:3)
          CASE(2) !y direction
              LocX(:)=yp(1:3)
          CASE(3) !z direction  
              LocX(:)=zp(1:3)
      END SELECT
            
      SELECT CASE(DISC%DynRup%NucDirY)
          CASE(1) !x direction
              LocY(:)=xp(1:3)
          CASE(2) !y direction
              LocY(:)=yp(1:3)
          CASE(3) !z direction  
              LocY(:)=zp(1:3)
      END SELECT 
      IF(   MAXVAL(LocX(1:3)).LE.(DISC%DynRup%NucXMax) .AND. MINVAL(LocX(1:3)).GE.(DISC%DynRup%NucXMin)     &
          .AND. MAXVAL(LocY(1:3)).LE.(DISC%DynRup%NucYMax) .AND. MINVAL(LocY(1:3)).GE.(DISC%DynRup%NucYMin)) THEN
          EQN%IniShearXY(i,:)  = DISC%DynRup%NucShearXY_0
          EQN%IniShearYZ(i,:)  = DISC%DynRup%NucShearYZ_0
          EQN%IniShearXZ(i,:)  = DISC%DynRup%NucShearXZ_0
          EQN%IniBulk_xx(i,:)  = DISC%DynRup%NucBulk_xx_0
          EQN%IniBulk_yy(i,:)  = DISC%DynRup%NucBulk_yy_0    
          EQN%IniBulk_zz(i,:)  = DISC%DynRup%NucBulk_zz_0
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ENDIF
            
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE nucleation_STEP       ! Nucleation by discontinuous jump on properties at [NucXMin,NucXMax] x [NucYMin,NucYMax]  
  
  !> Nucleation by smooth jump on properties assigned to each Gaussian node
  !<
  SUBROUTINE nucleation_SMOOTH_GP(DISC,EQN,MESH,BND)
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
  REAL                           :: r_i, r_a, r_b, r_s, x_hc, y_hc, d, sigma_n, sigma_i, sigma_y
  REAL                           :: LocX(3), LocY(3)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! - as in SPICE Code Validation
  ! - or in Galis et al. 2010, Fig. 3d. nomenclature is adapted from this paper.            
      
  ! compute/copy variables
  x_hc = DISC%DynRup%NucXmin ! x_hc = x coordinate of the hypocenter
  y_hc = DISC%DynRup%NucYmin ! y_hc = y coordinate of the hypocenter
  r_a  = DISC%DynRup%NucXmax ! r_a = major semi-axis
  r_b  = DISC%DynRup%NucYmax ! r_b = minor semi-axis
  r_s  = DISC%DynRup%r_s     ! width of the smooth transition
  sigma_n = 120.0e6 ! for the moment fix.
  sigma_i = 81.341e6 ! for the moment fix.
  sigma_y = DISC%DynRup%Mu_S_ini * sigma_n !  static traction sigma_y - assume constant Mu_S_ini along complete fault!
  r_i = r_a - r_s/EQN%PI * acos(2.0D0*(sigma_y - EQN%ShearXZ_0)/(sigma_i - EQN%ShearXZ_0)-1.0D0)
            !        
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
          CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
          !
          ! Select the 2 dimensions of the nucleation zone
          ! major axis
          SELECT CASE(DISC%DynRup%NucDirX)
              CASE(1) !x direction
                  LocX(1)=xGP
              CASE(2) !y direction
                  LocX(1)=yGP
              CASE(3) !z direction  
                  LocX(1)=zGP
          END SELECT
          !
          ! minor axis
          SELECT CASE(DISC%DynRup%NucDirY)
              CASE(1) !x direction
                  LocY(1)=xGP
              CASE(2) !y direction
                  LocY(1)=yGP
              CASE(3) !z direction  
                  LocY(1)=zGP
          END SELECT
          !
          d = SQRT((LocX(1)-x_hc)**2 + (r_a/r_b)**2*(LocY(1)-y_hc)**2)
          !
          IF (d.LE.r_i) THEN                   
              EQN%IniShearXY(i,iBndGP)  = DISC%DynRup%NucShearXY_0
              EQN%IniShearYZ(i,iBndGP)  = DISC%DynRup%NucShearYZ_0
              EQN%IniShearXZ(i,iBndGP)  = DISC%DynRup%NucShearXZ_0
              EQN%IniBulk_xx(i,iBndGP)  = DISC%DynRup%NucBulk_xx_0
              EQN%IniBulk_yy(i,iBndGP)  = DISC%DynRup%NucBulk_yy_0
              EQN%IniBulk_zz(i,iBndGP)  = DISC%DynRup%NucBulk_zz_0
              EQN%IniStateVar(i,iBndGP) = DISC%DynRup%NucRS_sv0
          !
          ELSEIF (d.LE.(r_i+r_s) .AND. d.GT.r_i) THEN
              EQN%IniShearXY(i,iBndGP)  = EQN%ShearXY_0 + 0.5D0*(DISC%DynRup%NucShearXY_0-EQN%ShearXY_0) * (1.0D0+COS(EQN%PI*(d-r_i)/r_s))
              EQN%IniShearYZ(i,iBndGP)  = EQN%ShearYZ_0 + 0.5D0*(DISC%DynRup%NucShearYZ_0-EQN%ShearYZ_0) * (1.0D0+COS(EQN%PI*(d-r_i)/r_s))
              EQN%IniShearXZ(i,iBndGP)  = EQN%ShearXZ_0 + 0.5D0*(DISC%DynRup%NucShearXZ_0-EQN%ShearXZ_0) * (1.0D0+COS(EQN%PI*(d-r_i)/r_s))
              EQN%IniBulk_xx(i,iBndGP)  = EQN%Bulk_xx_0 + 0.5D0*(DISC%DynRup%NucBulk_xx_0-EQN%Bulk_xx_0) * (1.0D0+COS(EQN%PI*(d-r_i)/r_s))
              EQN%IniBulk_yy(i,iBndGP)  = EQN%Bulk_yy_0 + 0.5D0*(DISC%DynRup%NucBulk_yy_0-EQN%Bulk_yy_0) * (1.0D0+COS(EQN%PI*(d-r_i)/r_s))
              EQN%IniBulk_zz(i,iBndGP)  = EQN%Bulk_zz_0 + 0.5D0*(DISC%DynRup%NucBulk_zz_0-EQN%Bulk_zz_0) * (1.0D0+COS(EQN%PI*(d-r_i)/r_s))
              EQN%IniStateVar(i,iBndGP) = EQN%RS_sv0    + 0.5D0*(DISC%DynRup%NucRS_sv0-EQN%RS_sv0)       * (1.0D0+COS(EQN%PI*(d-r_i)/r_s))
          ENDIF
      ENDDO ! iBndGP
      
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE nucleation_SMOOTH_GP       ! Nucleation by smooth jump on properties assigned to each Gaussian node 
  
  !> Nucleation by discontinuous elliptic nucleation zone
  !<
  SUBROUTINE nucleation_ELLIPSE(DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: LocX(3), LocY(3)
  REAL                           :: r_i, r_a, r_b, r_s, x_hc, y_hc, d, sigma_n, sigma_i, sigma_y
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! variables equal to nucleation by smooth propertie change (Case 2)            
  
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !
  
  ! compute/copy variables
  x_hc = DISC%DynRup%NucXmin ! x_hc = x coordinate of the hypocenter
  y_hc = DISC%DynRup%NucYmin ! y_hc = y coordinate of the hypocenter
  r_a  = DISC%DynRup%NucXmax ! r_a = major semi-axis
  r_b  = DISC%DynRup%NucYmax ! r_b = minor semi-axis
  r_s  = DISC%DynRup%r_s     ! width of the smooth transition
  r_i = r_a - r_s/EQN%PI * acos(2.0D0*(sigma_y - EQN%ShearXZ_0)/(sigma_i - EQN%ShearXZ_0)-1.0D0)
            !        
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
   
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
  
      ! get coordinates of fault face
      SELECT CASE(DISC%DynRup%NucDirX)
          CASE(1) !x direction
              LocX(:)=xp(1:3)
          CASE(2) !y direction
              LocX(:)=yp(1:3)
          CASE(3) !z direction  
              LocX(:)=zp(1:3)
      END SELECT
      !
      SELECT CASE(DISC%DynRup%NucDirY)
          CASE(1) !x direction
              LocY(:)=xp(1:3)
          CASE(2) !y direction 
              LocY(:)=yp(1:3)
          CASE(3) !z direction  
              LocY(:)=zp(1:3)
      END SELECT
            
      ! use barycenter to test if element is inside the ellipse
      d = SQRT((sum(LocX(1:3))/3.0D0-x_hc)**2 + (r_a/r_b)**2*(sum(LocY(1:3))/3.0D0-y_hc)**2)
  
      IF (d.LE.r_i) THEN                   
          EQN%IniShearXY(i,:)  = DISC%DynRup%NucShearXY_0
          EQN%IniShearYZ(i,:)  = DISC%DynRup%NucShearYZ_0
          EQN%IniShearXZ(i,:)  = DISC%DynRup%NucShearXZ_0
          EQN%IniBulk_xx(i,:)  = DISC%DynRup%NucBulk_xx_0
          EQN%IniBulk_yy(i,:)  = DISC%DynRup%NucBulk_yy_0
          EQN%IniBulk_zz(i,:)  = DISC%DynRup%NucBulk_zz_0
          EQN%IniStateVar(i,:) = DISC%DynRup%NucRS_sv0
      ENDIF
      
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE nucleation_ELLIPSE       ! Nucleation by discontinuous elliptic nucleation zone
  
  !> Nucleation as in  SCEC TPV28 test assigned to each Gaussian node
  !<
  SUBROUTINE nucleation_TPV28_GP(DISC,EQN,MESH,BND)
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
  REAL                           :: radius
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! Nucleation as in  SCEC TPV28 test with two hills (circular 2-stage-nucleation)
  ! Requires the correct mesh! For instructions see http://scecdata.usc.edu/cvws/download/TPV28_Description_v06.pdf
  ! vertical fault in the xz-plane, 30 km long by 15 km deep
  ! depth is z negative, hypocenter centered along-strike and along-dip
  ! nucleation is 2-stage stress perturbation in a circular regions with radii of 1.4 km and 2 km
  ! fault reaches the surface at depth 0km
  ! units of mesh in meter
  ! Gauss node wise nucleation properties assignment
  ! fault normal reference point at +y

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
             ! nucleation shear stress is added to the shear stress obtained by resolving 
          ! the stress tensor onto the fault surface
          !
          ! radial distance to hypocenter at x=0, z=-7500
          radius=SQRT(xGP**2+(zGP+7500.0D0)**2)
          ! 1. inner circular nucleation zone, constantly overstressed
          IF (radius.LE.1400.0D0) THEN
              EQN%IniShearXY(i,iBndGP)  =  11600000.0D0 + EQN%ShearXY_0
          ! 2. outer circular nucleation zone, smooth gradient
          ELSEIF ((radius.GT.1400.0D0) .AND. (radius.LE.2000.0D0)) THEN
              EQN%IniShearXY(i,iBndGP)  =  5800000.0D0 * (1.0D0+COS(EQN%PI*(radius-1400.0D0)/600.0D0)) + EQN%ShearXY_0 
          ENDIF
      ENDDO ! iBndGP
      
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE nucleation_TPV28_GP       ! Nucleation as in  SCEC TPV28 test assigned to each Gaussian node
  
  !> Nucleation as in  SCEC TPV28 test assigned to each element
  !<
  SUBROUTINE nucleation_TPV28 (DISC,EQN,MESH)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  !-------------------------------------------------------------------------!
  IMPLICIT NONE
  !-------------------------------------------------------------------------!
  TYPE(tDiscretization), target  :: DISC
  TYPE(tEquations)               :: EQN
  TYPE(tUnstructMesh)            :: MESH
  !-------------------------------------------------------------------------!
  ! Local variable declaration
  INTEGER                        :: i,j
  INTEGER                        :: iSide,iElem
  INTEGER                        :: iLocalNeighborSide
  INTEGER                        :: VertexSide(4,3)
  REAL                           :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
  REAL                           :: LocX(3), LocY(3)
  REAL                           :: r_i, r_a, r_b, r_s, chi, tau, x_hc, y_hc, d, radius
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! same as Case 28 but elementwise stress assignment
  ! radial distance to hypocenter at x=0, z=-7500
  ! based on barycenter 
  
  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !
  
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
       
      ! element ID    
      iElem = MESH%Fault%Face(i,1,1)
      iSide = MESH%Fault%Face(i,2,1)         
   
      ! get coordinates needed for special background types and nucleation zone
      IF (iElem .NE. 0) THEN
          !
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
          ENDDO           
      ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
          !
          iLocalNeighborSide = MESH%Fault%Face(i,2,2)
          DO j=1,3
              xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
              zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(i,1,2)))
          ENDDO            
      ENDIF
      radius=SQRT( ( SUM(xp(1:3))/3.0D0 )**2+ ( SUM(zp(1:3))/3.0D0+7500.0D0 )**2 )
      ! 1. inner circular nucleation zone, constantly overstressed
      IF (radius.LE.1400.0D0) THEN
          EQN%IniShearXY(i,:)  =  11600000.0D0 + EQN%ShearXY_0
      ! 2. outer circular nucleation zone, smooth gradient
      ELSEIF ((radius.GT.1400.0D0) .AND. (radius.LE.2000.0D0)) THEN
          EQN%IniShearXY(i,:)  =  5800000.0D0 * (1.0D0+COS(EQN%PI*(radius-1400.0D0)/600.0D0)) + EQN%ShearXY_0 
      ENDIF
                
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE nucleation_TPV28       ! Nucleation as in  SCEC TPV28 test assigned to each element
  
  
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
          !EQN%IniStateVar(i,iBndGP) = DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0*EXP((sstress/(nstress*DISC%DynRup%RS_b))-DISC%DynRup%RS_f0/DISC%DynRup%RS_b-DISC%DynRup%RS_a_array(i,iBndGP)/DISC%DynRup%RS_b*LOG(iniSlipRate/DISC%DynRup%RS_sr0))
          X2  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a)
          EQN%IniMu(i,iBndGP)=DISC%DynRup%RS_a * LOG(X2 + SQRT(X2**2 + 1.0))   
             
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
          EQN%IniMu(i,iBndGP)=DISC%DynRup%RS_f0+DISC%DynRup%RS_a*iniSlipRate/(iniSlipRate+DISC%DynRup%RS_sr0)-DISC%DynRup%RS_b*EQN%IniStateVar(i,iBndGP)/(EQN%IniStateVar(i,iBndGP)+DISC%DynRup%RS_sl0)             
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
          EQN%IniStateVar(i,iBndGP) = (DISC%DynRup%RS_sl0/DISC%DynRup%RS_sr0) * EXP((-EQN%IniShearXY(i,iBndGP)/EQN%IniBulk_yy(i,iBndGP)-DISC%DynRup%RS_f0-DISC%DynRup%RS_a_array(i,iBndGP)*LOG(iniSlipRate/DISC%DynRup%RS_sr0))/DISC%DynRup%RS_b)
          ! ASINH(X)=LOG(X+SQRT(X^2+1))
          tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP((DISC%DynRup%RS_f0 + DISC%DynRup%RS_b*LOG(DISC%DynRup%RS_sr0*EQN%IniStateVar(i,iBndGP)/DISC%DynRup%RS_sl0)) / DISC%DynRup%RS_a_array(i,iBndGP))
          EQN%IniMu(i,iBndGP)=DISC%DynRup%RS_a_array(i,iBndGP) * LOG(tmp + SQRT(tmp**2 + 1.0))

      ENDDO ! iBndGP
      
  ENDDO !    MESH%Fault%nSide   
              
  END SUBROUTINE friction_RSF101      ! Initialization of initial slip rate and friction for SCEC TPV101
  
  !> Initialization of initial slip rate and friction for SCEC TPV103 
  !<
  SUBROUTINE friction_RSF103(DISC,EQN,MESH,BND)
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
          ! SINH(X)=(EXP(X)-EXP(-X))/2
          tmp = ABS(EQN%IniShearXY(i,iBndGP)/(DISC%DynRup%RS_a_array(i,iBndGP)*EQN%IniBulk_yy(i,iBndGP)))
          EQN%IniStateVar(i,iBndGP)=DISC%DynRup%RS_a_array(i,iBndGP)*LOG(2.0D0*DISC%DynRup%RS_sr0/iniSlipRate * (EXP(tmp)-EXP(-tmp))/2.0D0)
          ! ASINH(X)=LOG(X+SQRT(X^2+1))
          tmp  = iniSlipRate*0.5/DISC%DynRup%RS_sr0 * EXP(EQN%IniStateVar(i,iBndGP)/ DISC%DynRup%RS_a_array(i,iBndGP))
          EQN%IniMu(i,iBndGP)=DISC%DynRup%RS_a_array(i,iBndGP) * LOG(tmp + SQRT(tmp**2 + 1.0D0))
          
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
  
          EQN%IniMu(i,iBndGP) = DISC%DynRup%Mu_S(i,iBndGP)
          
      ENDDO ! iBndGP
      
  ENDDO !    MESH%Fault%nSide   
              
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
  
          EQN%IniMu(i,iBndGP) = DISC%DynRup%Mu_S(i,iBndGP)
          !Attention: normal stress is here assumed in yy direction
          DISC%DynRup%Strength(i,iBndGP) = EQN%IniMu(i,iBndGP)*EQN%IniBulk_yy(i,iBndGP)
          
      ENDDO ! iBndGP
      
  ENDDO !    MESH%Fault%nSide
  
  END SUBROUTINE friction_LSW6    ! Initialization of friction for bimaterial linear slip weakening
  
  
  
  END MODULE ini_model_DR_mod
