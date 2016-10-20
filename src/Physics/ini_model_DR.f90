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
  PRIVATE :: background_TPV2627
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

  PRIVATE :: background_SUMATRA
  PRIVATE :: background_SUMATRA_RS
  PRIVATE :: background_SUMATRA_GEO

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
       ! SCEC TPV10 test dipping fault subshear, SCEC TPV12/TPV13 dipping fault
       CALL background_TPV10(DISC,EQN,MESH,BND)
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
   CASE(26)
       !  26 = SCEC TPV26/27 with heterogeneous/depth dependent initial stress field
       CALL background_TPV2627(EQN,MESH,DISC,BND)
   CASE(33)
       !  SCEC TPV3132 test case : strike slip rupture in layered medium
       CALL background_TPV33 (DISC,EQN,MESH,BND)
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
    CASE(1201)
       CALL background_SUMATRA_GEO(DISC,EQN,MESH,BND)
    CASE(1202)
       CALL background_SUMATRA_RS(DISC,EQN,MESH,BND)
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
    ALLOCATE(  DISC%DynRup%DS(MESH%Fault%nSide,DISC%Galerkin%nBndGP)        )
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
    !ini dyn.stress ouput
    DISC%DynRup%DS = .FALSE.

    ! Initialize '+'side elements for RF und DS output
    IF ((DISC%DynRup%RFtime_on .EQ. 1) .AND. (DISC%DynRup%DS_output_on .EQ. 1) ) THEN
       ! Loop over every mesh element
           DO i = 1, MESH%Fault%nSide
              IF (MESH%FAULT%Face(i,1,1) .NE. 0) THEN
                 DISC%DynRup%RF(i,:) = .TRUE.
                 DISC%DynRup%DS(i,:) = .TRUE.
              ENDIF
           ENDDO
     ELSEIF ((DISC%DynRup%RFtime_on .EQ. 1) .AND. (DISC%DynRup%DS_output_on .EQ. 0 )) THEN
           DO i = 1, MESH%Fault%nSide
              IF (MESH%FAULT%Face(i,1,1) .NE. 0) THEN
                 DISC%DynRup%RF(i,:) = .TRUE.
              ENDIF
          ENDDO
    ENDIF


    !T. Ulrich 8.2015 initial rupture time array (for Vr calculations)
    ALLOCATE(DISC%DynRup%rupture_time(MESH%Fault%nSide,DISC%Galerkin%nBndGP))
    DISC%DynRup%rupture_time(:,:)=0.

    !time at which the shear stress is equal the dynamic stress
    ALLOCATE(DISC%DynRup%dynStress_time(MESH%Fault%nSide,DISC%Galerkin%nBndGP))
    DISC%DynRup%dynStress_time(:,:)=0.

    
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
       
       Vs = DISC%DynRup%Vs_nucl
       IF (abs(Vs).LE.1d-6) THEN
          Vs = SQRT(EQN%mu/EQN%rho0)
       ENDIF

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

                    CASE(13) ! S.Wollherr 2014 TPV12/13
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

  !SCEC benchmark TPV26/TPV27
  SUBROUTINE background_TPV2627(EQN,MESH,DISC,BND)
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
  REAL                           :: b11, b33, b13, Pf
  REAL                           :: omega !depth dependent
  !-------------------------------------------------------------------------!
  INTENT(IN)    :: MESH, BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------!

  ! based on SCEC TPV26/27 test right-lateral strike-slip fault, z negative in depth
  ! 26 with elastic and 27 with viscoplastic material properties
  b11 = 0.926793
  b33 = 1.073206
  b13 = -0.169029

  VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
  VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
  VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
  VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !

  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

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


      IF (EQN%GPwise.EQ.1) THEN
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

              ! Transformation of boundary GP's into XYZ coordinate system
              chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
              tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
              CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
              CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
              !
              !depth, negative in depth
              !average = zGP   ! Averaging not needed here

              !

              Pf = 9800.0D0* abs(zGP) !fluid pressure, hydrostatic with water table at the surface
              IF (zGP.GE. -15000.0D0) THEN !depth less than 15000m
               omega = 1.0D0
              ELSEIF ((zGP.LT. -15000.0D0) .AND. (zGP .GE. -20000.0D0) ) THEN !depth between 15000 and 20000m
               omega = (20000.0D0-abs(zGP))/5000.0D0
              ELSE ! depth more than 20000m
               omega = 0.0D0
              ENDIF


             EQN%IniBulk_zz(i,iBndGP)  = -2670D0*9.8D0 * abs(zGP)
             EQN%IniBulk_xx(i,iBndGP)  = omega*(b11*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1-omega)*EQN%IniBulk_zz(i,iBndGP)
             EQN%IniBulk_yy(i,iBndGP)  = omega*(b33*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1-omega)*EQN%IniBulk_zz(i,iBndGP)
             EQN%IniShearXY(i,iBndGP)  = omega*(b13*(EQN%IniBulk_zz(i,iBndGP)+Pf))
             EQN%IniShearXZ(i,iBndGP)  = 0.0D0
             EQN%IniShearYZ(i,iBndGP)  = 0.0D0
             !add fluid pressure
             EQN%IniBulk_xx(i,iBndGP)  = EQN%IniBulk_xx(i,iBndGP)+Pf
             EQN%IniBulk_yy(i,iBndGP)  = EQN%IniBulk_yy(i,iBndGP)+Pf
             EQN%IniBulk_zz(i,iBndGP)  = EQN%IniBulk_zz(i,iBndGP)+Pf

             !depth dependent frictional cohesion, negative in seissol, in benchmark positive



            IF (zGP.GE.-5000.0D0) THEN
               DISC%DynRup%cohesion(i,iBndGP) = -0.4D6 - 0.00072D6*(5000D0-abs(zGP))
            ELSE
               DISC%DynRup%cohesion(i,iBndGP) = -0.4D6
            ENDIF

      ENDDO ! iBndGP

      ! element wise stress assignment
      ELSE
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
      !
      !depth (negative in z)
      zGP = sum(zp(:))/3.0D0 !average but named zGP as above


      Pf = 9800.0D0* abs(zGP) !fluid pressure, hydrostatic with water table at the surface

      IF (zGP.GE. -15000.0D0) THEN !depth less than 15000m
       omega = 1.0D0
      ELSEIF ((zGP.LT. -15000.0D0) .AND. (zGP .GE. -20000.0D0) ) THEN !depth between 15000 and 20000m
       omega = (20000.0D0-abs(zGP))/5000.0D0
      ELSE ! depth more than 20000m
       omega = 0.0D0
      ENDIF


     EQN%IniBulk_zz(i,:)  = -2670D0*9.8D0 * abs(zGP)
     EQN%IniBulk_xx(i,:)  = omega*(b11*(EQN%IniBulk_zz(i,:)+Pf)-Pf)+(1-omega)*EQN%IniBulk_zz(i,:)
     EQN%IniBulk_yy(i,:)  = omega*(b33*(EQN%IniBulk_zz(i,:)+Pf)-Pf)+(1-omega)*EQN%IniBulk_zz(i,:)
     EQN%IniShearXY(i,:)  = omega*(b13*(EQN%IniBulk_zz(i,:)+Pf))
     EQN%IniShearXZ(i,:)  = 0.0D0
     EQN%IniShearYZ(i,:)  = 0.0D0
     !add fluid pressure
     EQN%IniBulk_xx(i,:)  = EQN%IniBulk_xx(i,:)+Pf
     EQN%IniBulk_yy(i,:)  = EQN%IniBulk_yy(i,:)+Pf
     EQN%IniBulk_zz(i,:)  = EQN%IniBulk_zz(i,:)+Pf


     !depth dependent frictional cohesion, negative in seissol, in benchmark positive

    IF (zGP.GE.-5000.0D0) THEN
       DISC%DynRup%cohesion(i,:) = -0.4D6 - 0.00072D6*(5000D0-abs(zGP))
    ELSE
       DISC%DynRup%cohesion(i,:) = -0.4D6
    ENDIF



    ENDIF !node or elementwise

    ENDDO !    MESH%Fault%nSide
  END SUBROUTINE
   
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
  INTEGER                        :: k, nLayers, Laterally_homogenous_Stress
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: b11, b22, b12, b13, b23, Omega, g, Pf, zIncreasingCohesion
  REAL                           :: b11_N, b22_N, b12_N, b13_N, b23_N
  REAL                           :: b11_C, b22_C, b12_C, b13_C, b23_C
  REAL                           :: b11_S, b22_S, b12_S, b13_S, b23_S
  REAL                           :: yN1, yN2, yS1, yS2, alpha
  REAL                           :: sigzz, Rz, zLayers(20), rhoLayers(20)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! TPV29
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function (gravity)
  ! NOTE: z negative is depth, free surface is at z=0
  Laterally_homogenous_Stress = 1

  IF (Laterally_homogenous_Stress.EQ.1) THEN
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds
     b11 = 1.1854
     b22 = 1.3162
     b12 = 0.3076
     b13 = 0.1259
     b23 = 0.1555
  ELSE
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds,  South (strike = 25+90+180)
     b11_S = 1.0487
     b22_S = 1.4529
     b12_S = 0.2409
     b13_S = 0.0846
     b23_S = 0.1814
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds,  Center (strike = 40+90+180)
     b11_C = 1.1962
     b22_C = 1.3054
     b12_C = 0.3097
     b13_C = 0.1286
     b23_C = 0.1533
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds,  Center (strike = 75+90+180)
     b11_N = 1.5231
     b22_N = 0.9785
     b12_N = 0.1572
     b13_N = 0.1933
     b23_N = 0.0518

     ! 10.5/8.5/5/4
     yN2 = 1160695.0941260615
     yN1 = 939574.3060454715
     yS2 = 552664.2968779367
     yS1 = 442127.3902531094
  ENDIF

  g = 9.8D0    
  zIncreasingCohesion = -10000.
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide

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

          ! TO BE USED WITH 1d Layered medium
          !free surface assumed at z=-2000m
          !properties of continental crust
          nLayers = 6
          zLayers (1:6) = (/ 0d0,-2000d0, -6000d0, -12000d0, -23000d0,-600d6 /)
          rhoLayers (1:6) = (/ 1000d0, 2720d0, 2860d0, 3050d0, 3300d0, 3375d0 /)
          sigzz = 0d0

          DO k=2,nLayers
             IF (zGP.GT.zLayers(k)) THEN
                sigzz = sigzz + rhoLayers(k-1)*(zGP-zLayers(k-1))*g
                EXIT
             ELSE
                sigzz = sigzz + rhoLayers(k-1)*(zLayers(k)-zLayers(k-1))*g
             ENDIF
          ENDDO

          IF (zGP.LT.-25000D0) THEN
             Rz = (-zGp - 25000D0)/150e3
             !Rz = (-zGp - 25000D0)/50e3
             !Rz = 0d0
          ELSE
             Rz = 0.
          ENDIF
          !Rz = max(0D0,min(0.99999999999d0, Rz))
          !DISC%DynRup%Mu_D(i,iBndGP) = (1d0-Rz)*DISC%DynRup%Mu_D_ini + Rz*DISC%DynRup%Mu_S_ini
          !Rz = 0d0

          Omega = max(0D0,min(1d0, 1D0-Rz))

          IF (Laterally_homogenous_Stress.EQ.0) THEN
             ! The stress varies along y (along lon=cst)
             ! cst_N
             !**yN2
             ! lin from cst_N to cst_C
             !**yN1
             ! cst_C
             !**yS2
             ! lin from cst_C to cst_S
             !**yS1
             ! cst_S

             IF (yGP.GE.yN2) THEN
                b11 = b11_N
                b22 = b22_N
                b12 = b12_N
                b13 = b13_N
                b23 = b23_N
             ELSE IF ((yGP.GE.yN1).AND.(yGP.LT.yN2)) THEN
                alpha = (yGP-yN1)/(yN2-yN1)
                b11 = alpha * b11_N + (1d0-alpha)* b11_C
                b22 = alpha * b22_N + (1d0-alpha)* b22_C
                b12 = alpha * b12_N + (1d0-alpha)* b12_C
                b13 = alpha * b13_N + (1d0-alpha)* b13_C
                b23 = alpha * b23_N + (1d0-alpha)* b23_C
             ELSE IF ((yGP.GE.yS2).AND.(yGP.LT.yN1)) THEN
                b11 = b11_C
                b22 = b22_C
                b12 = b12_C
                b13 = b13_C
                b23 = b23_C
             ELSE IF ((yGP.GE.yS1).AND.(yGP.LT.yS2)) THEN
                alpha = (yGP-yS1)/(yS2-yS1)
                b11 = alpha * b11_C + (1d0-alpha)* b11_S
                b22 = alpha * b22_C + (1d0-alpha)* b22_S
                b12 = alpha * b12_C + (1d0-alpha)* b12_S
                b13 = alpha * b13_C + (1d0-alpha)* b13_S
                b23 = alpha * b23_C + (1d0-alpha)* b23_S
             ELSE
                b11 = b11_S
                b22 = b22_S
                b12 = b12_S
                b13 = b13_s
                b23 = b23_s
             ENDIF
          ENDIF

          Pf = -1000D0 * g * zGP
          EQN%IniBulk_zz(i,iBndGP)  =  sigzz
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
              !DISC%DynRup%cohesion(i,iBndGP) = -0.4d6-0.0002d6*(zGP-zIncreasingCohesion)
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6-1.0d6*(zGP-zIncreasingCohesion)/(-zIncreasingCohesion)
          ELSE
              ! set cohesion
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6
          ENDIF
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   

  END SUBROUTINE background_SUMATRA
                
!> SUMATRA test case with RS friction
  !> T. ULRICH 07.2016
  !<
  SUBROUTINE background_SUMATRA_RS (DISC,EQN,MESH,BND)
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
  INTEGER                        :: k, nLayers, Laterally_homogenous_Stress
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: b11, b22, b12, b13, b23, Omega, g, Pf, zIncreasingCohesion
  REAL                           :: b11_N, b22_N, b12_N, b13_N, b23_N
  REAL                           :: b11_C, b22_C, b12_C, b13_C, b23_C
  REAL                           :: b11_S, b22_S, b12_S, b13_S, b23_S
  REAL                           :: yN1, yN2, yS1, yS2, alpha
  REAL                           :: sigzz, Rz, zLayers(20), rhoLayers(20)
  REAL                           :: zBoStartTapering, zToStartTapering, zBoStopTapering, zToStopTapering
  REAL                           :: zTotaperingWidth, zBotaperingWidth, RS_a_inc, RS_srW_inc, Boxx, Boxz, tmp
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! TPV29
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function (gravity)
  ! NOTE: z negative is depth, free surface is at z=0

  ALLOCATE(  DISC%DynRup%RS_a_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP)        )
  ALLOCATE(  DISC%DynRup%RS_srW_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP)      )

  zToStartTapering = -10d3
  zBoStartTapering = -25d3
  zTotaperingWidth = 5d3
  zBotaperingWidth = 25d3
  RS_a_inc = 0.01d0
  RS_srW_inc = 0.9d0

  zBoStopTapering = zBoStartTapering - zBotaperingWidth
  zToStopTapering = zToStartTapering + zTotaperingWidth

  Laterally_homogenous_Stress = 1

  IF (Laterally_homogenous_Stress.EQ.1) THEN
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds
     b11 = 1.2644
     b22 = 1.4165
     b12 = 0.3577
     b13 = 0.0855
     b23 = 0.1055
  ELSE
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds,  South (strike = 25+90+180)
     b11_S = 1.0487
     b22_S = 1.4529
     b12_S = 0.2409
     b13_S = 0.0846
     b23_S = 0.1814
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds,  Center (strike = 40+90+180)
     b11_C = 1.1962
     b22_C = 1.3054
     b12_C = 0.3097
     b13_C = 0.1286
     b23_C = 0.1533
     !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds,  Center (strike = 75+90+180)
     b11_N = 1.5231
     b22_N = 0.9785
     b12_N = 0.1572
     b13_N = 0.1933
     b23_N = 0.0518

     ! 10.5/8.5/5/4
     yN2 = 1160695.0941260615
     yN1 = 939574.3060454715
     yS2 = 552664.2968779367
     yS1 = 442127.3902531094
  ENDIF

  g = 9.8D0    
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
      
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

          ! TO BE USED WITH 1d Layered medium
          !free surface assumed at z=-2000m
          !properties of continental crust
          nLayers = 6
          zLayers (1:6) = (/ 0d0,-2000d0, -6000d0, -12000d0, -23000d0,-600d6 /)
          rhoLayers (1:6) = (/ 1000d0, 2720d0, 2860d0, 3050d0, 3300d0, 3375d0 /)
          sigzz = 0d0

          DO k=2,nLayers
             IF (zGP.GT.zLayers(k)) THEN
                sigzz = sigzz + rhoLayers(k-1)*(zGP-zLayers(k-1))*g
                EXIT
             ELSE
                sigzz = sigzz + rhoLayers(k-1)*(zLayers(k)-zLayers(k-1))*g
             ENDIF
          ENDDO

          Omega = 1d0

          IF (Laterally_homogenous_Stress.EQ.0) THEN
             ! The stress varies along y (along lon=cst)
             ! cst_N
             !**yN2
             ! lin from cst_N to cst_C
             !**yN1
             ! cst_C
             !**yS2
             ! lin from cst_C to cst_S
             !**yS1
             ! cst_S

             IF (yGP.GE.yN2) THEN
                b11 = b11_N
                b22 = b22_N
                b12 = b12_N
                b13 = b13_N
                b23 = b23_N
             ELSE IF ((yGP.GE.yN1).AND.(yGP.LT.yN2)) THEN
                alpha = (yGP-yN1)/(yN2-yN1)
                b11 = alpha * b11_N + (1d0-alpha)* b11_C
                b22 = alpha * b22_N + (1d0-alpha)* b22_C
                b12 = alpha * b12_N + (1d0-alpha)* b12_C
                b13 = alpha * b13_N + (1d0-alpha)* b13_C
                b23 = alpha * b23_N + (1d0-alpha)* b23_C
             ELSE IF ((yGP.GE.yS2).AND.(yGP.LT.yN1)) THEN
                b11 = b11_C
                b22 = b22_C
                b12 = b12_C
                b13 = b13_C
                b23 = b23_C
             ELSE IF ((yGP.GE.yS1).AND.(yGP.LT.yS2)) THEN
                alpha = (yGP-yS1)/(yS2-yS1)
                b11 = alpha * b11_C + (1d0-alpha)* b11_S
                b22 = alpha * b22_C + (1d0-alpha)* b22_S
                b12 = alpha * b12_C + (1d0-alpha)* b12_S
                b13 = alpha * b13_C + (1d0-alpha)* b13_S
                b23 = alpha * b23_C + (1d0-alpha)* b23_S
             ELSE
                b11 = b11_S
                b22 = b22_S
                b12 = b12_S
                b13 = b13_s
                b23 = b23_s
             ENDIF
          ENDIF

          Pf = -1000D0 * g * zGP
          EQN%IniBulk_zz(i,iBndGP)  =  sigzz
          EQN%IniBulk_xx(i,iBndGP)  =  Omega*(b11*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniBulk_yy(i,iBndGP)  =  Omega*(b22*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniShearXY(i,iBndGP)  =  Omega*(b12*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniShearXZ(i,iBndGP)  =  Omega*(b13*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniShearYZ(i,iBndGP)  =  Omega*(b23*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniBulk_xx(i,iBndGP)  =  EQN%IniBulk_xx(i,iBndGP) + Pf
          EQN%IniBulk_yy(i,iBndGP)  =  EQN%IniBulk_yy(i,iBndGP) + Pf
          EQN%IniBulk_zz(i,iBndGP)  =  EQN%IniBulk_zz(i,iBndGP) + Pf
          

          IF ( ((zGP.GT.zToStartTapering).AND.(zGP.LT.zToStopTapering))      &
              .OR.((zGP.LT.zBoStartTapering).AND.(zGP.GT.zBoStopTapering))) THEN
              ! TANH(X)=(1-EXP(-2X))/(1+EXP(-2X))

              IF (zGP.LT.zBoStartTapering) THEN
                 ! z
                 tmp = -zBotaperingWidth/ABS(zGP-zBoStopTapering)+zBotaperingWidth/(ABS(zGP-zBoStartTapering))
                 tmp = EXP(-2.0D0*tmp)
                 if ((tmp-1).EQ.tmp) THEN
                    Boxz = 0d0
                 ELSE
                    Boxz = 0.5D0*(1.0D0+((1.0D0-tmp)/(1.0D0+tmp)))
                 ENDIF
              ELSEIF (zGP.GT.zToStartTapering) THEN
                 ! z
                 tmp = -zTotaperingWidth/ABS(zGP-zToStopTapering)+zTotaperingWidth/ABS(zGP-zToStartTapering)
                 tmp = EXP(-2.0D0*tmp)
                 if ((tmp-1).EQ.tmp) THEN
                    Boxz = 0d0
                 ELSE
                    Boxz = 0.5D0*(1.0D0+((1.0D0-tmp)/(1.0D0+tmp)))
                 ENDIF
              ELSE
                 Boxz =1d0
              ENDIF
              Boxx  =1d0
              ! smooth boxcar
              DISC%DynRup%RS_a_array(i,iBndGP) = DISC%DynRup%RS_a+RS_a_inc*(1.0D0-(Boxx*Boxz))
              DISC%DynRup%RS_srW_array(i,iBndGP) = DISC%DynRup%RS_srW+RS_srW_inc*(1.0D0-(Boxx*Boxz))
          ELSEIF ((zGP.GE.zToStopTapering) .OR. (zGP.LE.zBoStopTapering)) THEN
              ! velocity strengthening in exterior region (3km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a+RS_a_inc
              DISC%DynRup%RS_srW_array(i,iBndGP) = DISC%DynRup%RS_srW+RS_srW_inc
              !
          ELSE
              ! velocity-weakening in interior fault region (30 km * 15 km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a
              DISC%DynRup%RS_srW_array(i,iBndGP) = DISC%DynRup%RS_srW
              !
          ENDIF
          ! resulting changes in SV_ini done in friction_RSF103
          ! Nucleation in Evaluate friction special case

      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_SUMATRA_RS

  !> SUMATRA test case
  !> T. ULRICH 06.2015
  !> tpv29 used as a model
  !<
  SUBROUTINE background_SUMATRA_GEO (DISC,EQN,MESH,BND)
  !-------------------------------------------------------------------------!
  USE DGBasis_mod
  use JacobiNormal_mod, only: RotationMatrix3D
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
  INTEGER                        :: k, nLayers
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: b11, b22, b12, b13, b23, Omega, g, Pf, zIncreasingCohesion
  REAL                           :: sigzz, Rz, zLayers(20), rhoLayers(20) 
  REAL                           :: zLocal, ux(3),uy(3),uz(3),LocalStress(6),T(eqn%nVar,eqn%nVar), iT(eqn%nVar,eqn%nVar)
  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH, BND 
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  ! stress is assigned to each Gaussian node
  ! depth dependent stress function (gravity)

  !New parameters R=0.6, stress accounting for the 1d layered velocity, sii = sm - ds
  b11 = 1.1854
  b22 = 1.3162
  b12 = 0.3071
  b13 = 0.1259
  b23 = 0.1555

  g = 9.8D0    
  zIncreasingCohesion = -10000.
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
      
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
          ! R is taken at lon, lat  = (90,8)
          zLocal = sqrt(xGP**2+yGP**2+zGP**2)-6377726.19283
          ! TO BE USED WITH 1d Layered medium
          !free surface assumed at z=-2000m
          !properties of continental crust
          nLayers = 6
          zLayers (1:6) = (/ 0d0,-2000d0, -6000d0, -12000d0, -23000d0,-600d6 /)
          rhoLayers (1:6) = (/ 1000d0, 2720d0, 2860d0, 3050d0, 3300d0, 3375d0 /)
          sigzz = 0d0

          DO k=2,nLayers
             IF (zLocal.GT.zLayers(k)) THEN
                sigzz = sigzz + rhoLayers(k-1)*(zLocal-zLayers(k-1))*g
                EXIT
             ELSE
                sigzz = sigzz + rhoLayers(k-1)*(zLayers(k)-zLayers(k-1))*g
             ENDIF
          ENDDO

          IF (zLocal.LT.-25000D0) THEN
             Rz = (-zLocal - 25000D0)/150e3
          ELSE
             Rz = 0.
          ENDIF

          Omega = max(0D0,min(1d0, 1D0-Rz))
          
          Pf = -1000D0 * g * zLocal

          EQN%IniBulk_zz(i,iBndGP)  =  sigzz
          EQN%IniBulk_xx(i,iBndGP)  =  Omega*(b11*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniBulk_yy(i,iBndGP)  =  Omega*(b22*(EQN%IniBulk_zz(i,iBndGP)+Pf)-Pf)+(1d0-Omega)*EQN%IniBulk_zz(i,iBndGP)
          EQN%IniShearXY(i,iBndGP)  =  Omega*(b12*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniShearXZ(i,iBndGP)  =  Omega*(b13*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniShearYZ(i,iBndGP)  =  Omega*(b23*(EQN%IniBulk_zz(i,iBndGP)+Pf))
          EQN%IniBulk_xx(i,iBndGP)  =  EQN%IniBulk_xx(i,iBndGP) + Pf
          EQN%IniBulk_yy(i,iBndGP)  =  EQN%IniBulk_yy(i,iBndGP) + Pf
          EQN%IniBulk_zz(i,iBndGP)  =  EQN%IniBulk_zz(i,iBndGP) + Pf
          
          uz = (/xGP,yGP,zGP/)
          uz = uz/sqrt(uz(1)**2+uz(2)**2+uz(3)**2)

          ux(1) = - uz(2) 
          ux(2) =   uz(1)
          ux(3) =   0.
          ux = ux/sqrt(ux(1)**2+ux(2)**2+ux(3)**2)
          
          uy(1) = uz(2)*ux(3) - uz(3)*ux(2)
          uy(2) = uz(3)*ux(1) - uz(1)*ux(3)
          uy(3) = uz(1)*ux(2) - uz(2)*ux(1)

          localStress = (/EQN%IniBulk_xx(i,iBndGP), EQN%IniBulk_yy(i,iBndGP), EQN%IniBulk_zz(i,iBndGP),\
                          EQN%IniShearXY(i,iBndGP), EQN%IniShearYZ(i,iBndGP), EQN%IniShearXZ(i,iBndGP)/)

        ! compute & store rotation matrices:
        !   xyz to face-aligned coordinate system
        !   face-aligned coordinate system to xyz 
        call RotationMatrix3D( ux, uy, uz, T(:,:), iT(:,:),EQN )


          localStress=MATMUL(T(:6,:6),localStress)

          EQN%IniBulk_xx(i,iBndGP) = localStress(1)
          EQN%IniBulk_yy(i,iBndGP) = localStress(2)
          EQN%IniBulk_zz(i,iBndGP) = localStress(3)
          EQN%IniShearXY(i,iBndGP) = localStress(4)
          EQN%IniShearYZ(i,iBndGP) = localStress(5)
          EQN%IniShearXZ(i,iBndGP) = localStress(6)

          ! manage cohesion
          IF (zLocal.GE.zIncreasingCohesion) THEN
              ! higher cohesion near free surface
              !DISC%DynRup%cohesion(i,iBndGP) = -0.4d6-0.0002d6*(zGP-zIncreasingCohesion)
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6-1.0d6*(zLocal-zIncreasingCohesion)/(-zIncreasingCohesion)
          ELSE
              ! set cohesion
              DISC%DynRup%cohesion(i,iBndGP) = -0.4d6
          ENDIF
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_SUMATRA_GEO

  !> SCEC TPV33 test case : strike slip rupture in wave guide zone
  !> T. ULRICH 01.2016
  !<
  SUBROUTINE background_TPV33 (DISC,EQN,MESH,BND)
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
  REAL                           :: mu, Rx, Rz, Rt
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
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
      
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
          IF (xGP.LT.-9800D0) THEN
             Rx = (-xGp - 9800D0)/10e3
          ELSEIF (xGP.GT.1100D0) THEN
             Rx = (xGp - 1100D0)/10e3
          ELSE
             Rx = 0.
          ENDIF
          IF (zGP.LT.-8000D0) THEN
             Rz = (-zGp - 8000D0)/10e3
          ELSEIF (zGP.GT.-2300D0) THEN
             Rz = (zGp + 2300D0)/10e3
          ELSE
             Rz = 0.
          ENDIF
          Rt = min(1D0,sqrt(Rx**2+Rz**2))

          EQN%IniBulk_xx(i,iBndGP)  = -60d6
          EQN%IniBulk_yy(i,iBndGP)  = -60d6
          EQN%IniBulk_zz(i,iBndGP)  =  0d0
          EQN%IniShearXY(i,iBndGP)  =  30e6*(1d0-Rt)
          EQN%IniShearXZ(i,iBndGP)  =  0D0
          EQN%IniShearYZ(i,iBndGP)  =  0d0

          xHypo = -6D3
          zHypo = -6D3
          ! distance to hypocenter (approx plane fault)
          r = sqrt( ((x-xHypo)*(x-xHypo))+((z-zHypo)*(z-zHypo)))
          
          IF (r.LE.550D0) THEN
             EQN%IniShearXY(i,iBndGP)  =  EQN%IniShearXY(i,iBndGP)+3.150d6
          ELSEIF (r.LE.800D0) THEN
             EQN%IniShearXY(i,iBndGP)  =  EQN%IniShearXY(i,iBndGP)+1.575d6*(1d0+dcos(PI*(r-550d0)/250d0))
          ENDIF
                
      ENDDO ! iBndGP
                
  ENDDO !    MESH%Fault%nSide   
                
  END SUBROUTINE background_TPV33       

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
          ! original setup
          !EQN%IniBulk_xx(i,iBndGP)  =  -10.6723e6*(abs(z-2000.0D0))/1000.0D0
          !EQN%IniBulk_yy(i,iBndGP)  =  -29.3277e6*(abs(z-2000.0D0))/1000.0D0
          !EQN%IniBulk_zz(i,iBndGP)  =  -20.0000e6*(abs(z-2000.0D0))/1000.0D0
          !EQN%IniShearXY(i,iBndGP)  =   -3.7687e6*(abs(z-2000.0D0))/1000.0D0
          !can now be changed from the input
          EQN%IniBulk_xx(i,iBndGP)  =  EQN%Bulk_xx_0*(abs(z-2000.0D0))/1000.0D0
          EQN%IniBulk_yy(i,iBndGP)  =  EQN%Bulk_yy_0*(abs(z-2000.0D0))/1000.0D0
          EQN%IniBulk_zz(i,iBndGP)  =  EQN%Bulk_zz_0*(abs(z-2000.0D0))/1000.0D0
          EQN%IniShearXY(i,iBndGP)  =  EQN%ShearXY_0 *(abs(z-2000.0D0))/1000.0D0
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
          ! depth dependent, constant for cohesion_max = 0 (is 0 if not otherwise declared in the parameter file)
          IF (z.GT.-DISC%DynRup%cohesion_depth) THEN
              DISC%DynRup%cohesion(i,iBndGP) =  DISC%DynRup%cohesion_0 - DISC%DynRup%cohesion_max*(DISC%DynRup%cohesion_depth+zGP)/(DISC%DynRup%cohesion_depth+1500.0)
          ELSE
              DISC%DynRup%cohesion(i,iBndGP) = DISC%DynRup%cohesion_0
          ENDIF

      
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
  REAL                           :: xLeStartTapering, xRiStartTapering, xLeStopTapering, xRiStopTapering, xtaperingWidth
  REAL                           :: zStartTapering, zStopTapering, ztaperingWidth
  REAL                           :: RS_a_inc,RS_srW_inc
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

  xLeStartTapering = -15d3
  xRiStartTapering =  15d3
  zStartTapering = -15d3
  xtaperingWidth = 3d3
  ztaperingWidth = 3d3
  RS_a_inc = 0.01d0
  RS_srW_inc = 0.9d0

  xLeStopTapering = xLeStartTapering - xtaperingWidth
  xRiStopTapering = xRiStartTapering + xtaperingWidth
  zStopTapering = zStartTapering - ztaperingWidth



  ! cohesion 0MPa          
  DISC%DynRup%cohesion = 0.0
            
  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)
              
  ! Loop over every mesh element
  DO i = 1, MESH%Fault%nSide
      
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
          IF ( ((xGP.GT.xRiStartTapering).AND.(xGP.LT.xRiStopTapering))      &
              .OR.((xGP.LT.xLeStartTapering).AND.(xGP.GT.xLeStopTapering)) &
              .OR.((zGP.LT.zStartTapering).AND.(zGP.GT.zStopTapering))) THEN
              ! TANH(X)=(1-EXP(-2X))/(1+EXP(-2X))

              ! x left tapering
              IF (xGP.LT.xLeStartTapering) THEN
                 tmp = -xtaperingWidth/ABS(xGP-xLeStopTapering)+xtaperingWidth/ABS(xGP-xLeStartTapering)
                 tmp = EXP(-2.0D0*tmp)
                 if ((tmp-1).EQ.tmp) THEN
                    Boxx = 0d0
                 ELSE
                    ! x
                    Boxx = 0.5D0*(1.0D0+((1.0D0-tmp)/(1.0D0+tmp)))
                 ENDIF
              ELSE
              ! x right tapering
              IF (xGP.GT.xRiStartTapering) THEN
                 tmp = -xtaperingWidth/ABS(xGP-xRiStopTapering)+xtaperingWidth/ABS(xGP-xRiStartTapering)
                 tmp = EXP(-2.0D0*tmp)
                 if ((tmp-1).EQ.tmp) THEN
                    Boxx = 0d0
                 ELSE
                    ! x
                    Boxx = 0.5D0*(1.0D0+((1.0D0-tmp)/(1.0D0+tmp)))
                 ENDIF
              ELSE
                 Boxx =1d0
              ENDIF
              ENDIF

              IF (zGP.LT.zStartTapering) THEN
                 ! z
                 tmp = -ztaperingWidth/ABS(zGP-zStopTapering)+ztaperingWidth/ABS(zGP-zStartTapering)
                 tmp = EXP(-2.0D0*tmp)
                 if ((tmp-1).EQ.tmp) THEN
                    Boxz = 0d0
                 ELSE
                    Boxz = 0.5D0*(1.0D0+((1.0D0-tmp)/(1.0D0+tmp)))
                 ENDIF
              ELSE
                 Boxz =1d0
              ENDIF
              ! smooth boxcar
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a+RS_a_inc*(1.0D0-(Boxx*Boxz))
              DISC%DynRup%RS_srW_array(i,iBndGP)=DISC%DynRup%RS_srW+RS_srW_inc*(1.0D0-(Boxx*Boxz))
          ELSEIF ((xGP.GE.xRiStopTapering) .OR. (xGP.LE.xLeStopTapering)        &
              .OR. (zGP.LE.zStopTapering)) THEN
              ! velocity strengthening in exterior region (3km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a+RS_a_inc
              DISC%DynRup%RS_srW_array(i,iBndGP) = DISC%DynRup%RS_srW+RS_srW_inc
              !
          ELSE
              ! velocity-weakening in interior fault region (30 km * 15 km)
              DISC%DynRup%RS_a_array(i,iBndGP)=DISC%DynRup%RS_a
              DISC%DynRup%RS_srW_array(i,iBndGP) = DISC%DynRup%RS_srW
              !
          ENDIF
          ! resulting changes in SV_ini done in friction_RSF103
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
  INTEGER                        :: i
  INTEGER                        :: iSide,iElem,iBndGP
  INTEGER                        :: iLocalNeighborSide,iNeighbor
  INTEGER                        :: MPIIndex, iObject
  INTEGER                        :: allocstat            ! Allocation status return int.   !
  REAL                           :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
  REAL                           :: chi,tau
  REAL                           :: xi, eta, zeta, XGp, YGp, ZGp
  REAL                           :: iniSlipRate, tmp
  REAL                           :: Stress(6,1:DISC%Galerkin%nBndGP)
  REAL                           :: NormalVect_n(3)                                            ! Normal vector components         !
  REAL                           :: NormalVect_s(3)                                            ! Normal vector components         !
  REAL                           :: NormalVect_t(3)                                            ! Normal vector components         !
  REAL                           :: T(EQN%nVar,EQN%nVar)                                       ! Transformation matrix            !
  REAL                           :: iT(EQN%nVar,EQN%nVar)                                      ! inverse Transformation matrix    !

  !-------------------------------------------------------------------------! 
  INTENT(IN)    :: MESH,BND
  INTENT(INOUT) :: DISC,EQN
  !-------------------------------------------------------------------------! 
  
  EQN%IniSlipRate1 = DISC%DynRup%RS_iniSlipRate1
  EQN%IniSlipRate2 = DISC%DynRup%RS_iniSlipRate2
  iniSlipRate = SQRT(EQN%IniSlipRate1**2 + EQN%IniSlipRate2**2)

  IF(.NOT.ALLOCATED(DISC%DynRup%RS_a_array)) THEN
     logWarning(*) 'RS_a_array not allocated before call of friction_RSF103. Allocating now.'
     ALLOCATE(DISC%DynRup%RS_a_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP),STAT=allocstat)
     IF (allocStat .NE. 0) THEN
        logError(*) 'friction_RSF103: could not allocate all variables!'
        STOP
     END IF
     DISC%DynRup%RS_a_array(:,:) = DISC%DynRup%RS_a
  ENDIF
  IF(.NOT.ALLOCATED(DISC%DynRup%RS_srW_array)) THEN
     logWarning(*) 'RS_srW_array not allocated before call of friction_RSF103. Allocating now.'
     ALLOCATE(DISC%DynRup%RS_srW_array(MESH%Fault%nSide,DISC%Galerkin%nBndGP),STAT=allocstat)
     IF (allocStat .NE. 0) THEN
        logError(*) 'friction_RSF103: could not allocate all variables!'
        STOP
     END IF
     DISC%DynRup%RS_srW_array(:,:) = DISC%DynRup%RS_srW
  ENDIF


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

      !Background stress rotation to face's reference system
      !
      Stress(1,:)=EQN%IniBulk_xx(i,:)
      Stress(2,:)=EQN%IniBulk_yy(i,:)
      Stress(3,:)=EQN%IniBulk_zz(i,:)
      Stress(4,:)=EQN%IniShearXY(i,:)
      Stress(5,:)=EQN%IniShearYZ(i,:)
      Stress(6,:)=EQN%IniShearXZ(i,:)

      ! Local side's normal vector from "+" side = iElem
      NormalVect_n = MESH%Fault%geoNormals(1:3,i)
      NormalVect_s = MESH%Fault%geoTangent1(1:3,i)
      NormalVect_t = MESH%Fault%geoTangent2(1:3,i)
      !
      CALL RotationMatrix3D(NormalVect_n,NormalVect_s,NormalVect_t,T(:,:),iT(:,:),EQN)

      DO iBndGP=1, DISC%Galerkin%nBndGP
         Stress(:,iBndGP)=MATMUL(iT(1:6,1:6),Stress(:,iBndGP))
      ENDDO

      DO iBndGP = 1,DISC%Galerkin%nBndGP ! Loop over all Gauss integration points
          !
          ! Transformation of boundary GP's into XYZ coordinate system
          chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
          tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
          CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
          CALL TetraTrafoXiEtaZeta2XYZ(xGp,yGp,zGp,xi,eta,zeta,xV,yV,zV)
          ! SINH(X)=(EXP(X)-EXP(-X))/2
          !tmp = ABS(EQN%IniShearXY(i,iBndGP)/(DISC%DynRup%RS_a_array(i,iBndGP)*EQN%IniBulk_yy(i,iBndGP)))
          tmp = ABS(SQRT(Stress(4,iBndGP)**2+Stress(6,iBndGP)**2)/(DISC%DynRup%RS_a_array(i,iBndGP)*Stress(1,iBndGP)))
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
