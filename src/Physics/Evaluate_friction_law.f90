!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2013-2016, SeisSol Group
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
!! Module containing friction laws

#ifdef BG 
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE Eval_friction_law_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE DGBasis_mod
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  REAL, PARAMETER :: u_0  = 10e-14 ! slip rate is considered as being zero for instaneous healing
  REAL, PARAMETER :: ZERO = 0.0D0
  !---------------------------------------------------------------------------!
  INTERFACE Eval_friction_law
     MODULE PROCEDURE Eval_friction_law
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Eval_friction_law
  PRIVATE :: no_fault
  PRIVATE :: Linear_slip_weakening_bimaterial
  PRIVATE :: Linear_slip_weakening_TPV1617
  PRIVATE :: rate_and_state
  PRIVATE :: rate_and_state_vw
  PRIVATE :: rate_and_state_nuc101
  PRIVATE :: rate_and_state_nuc103
  !---------------------------------------------------------------------------!
  CONTAINS
  
  !> Interface to friction laws
  !<
  SUBROUTINE Eval_friction_law(    TractionGP_XY,TractionGP_XZ,        & ! OUT: updated Traction
                                   NorStressGP,XYStressGP,XZStressGP,  & ! IN: Godunov status
                                   iFace,iSide,iElem,time,timePoints,          & ! IN: element ID, time, inv Trafo
                                   rho,rho_neig,w_speed,w_speed_neig,  & ! IN: background values
                                   EQN,DISC,MESH,MPI,IO,BND)             ! global variables
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations), target       :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI
    TYPE(tInputOutput)             :: IO
    TYPE (tBoundary)               :: BND
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER     :: nBndGP,iTimeGP,nTimeGP
    INTEGER     :: iFace,iSide,iElem
    REAL        :: TractionGP_XY(:,:)
    REAL        :: TractionGP_XZ(:,:)
    REAL        :: NorStressGP(:,:)
    REAL        :: XYStressGP(:,:)
    REAL        :: XZStressGP(:,:)
    REAL        :: time
    real        :: timePoints(:)
    REAL        :: rho,rho_neig,w_speed(:),w_speed_neig(:)
    real        :: DeltaT(1:DISC%Galerkin%nTimeGP)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: MESH,MPI,IO,NorStressGP,XYStressGP,XZStressGP
    INTENT(IN)    :: iFace,iSide,iElem,rho,rho_neig,w_speed,w_speed_neig,time
    INTENT(INOUT) :: EQN,DISC,TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------! 

    ! load number of GP iterations
    nBndGP  = DISC%Galerkin%nBndGP
    nTimeGP = DISC%Galerkin%nTimeGP
    
    DeltaT(1)=timePoints(1)
    DO iTimeGP=2,nTimeGP
       DeltaT(iTimeGP)=timePoints(iTimeGP)-timePoints(iTimeGP-1)
    ENDDO
    DeltaT(nTimeGP) = DeltaT(nTimeGP) + DeltaT(1) ! to fill last segment of Gaussian integration
       
    ! Evaluate friction law GP-wise
    SELECT CASE(EQN%FL)
        CASE(0) ! No fault
        
           CALL no_fault(TractionGP_XY,TractionGP_XZ,XYStressGP,XZStressGP)
           
        CASE(2,16) ! Coulomb model for LSW

           CALL Linear_slip_weakening_TPV1617(                                     & !
                                TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
                                NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                                iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                                rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                                time,DeltaT,                               & ! IN: time
                                DISC,EQN,MESH,MPI,IO)                          
                                
        CASE(3,4) ! Rate-and-state friction
        
           CALL rate_and_state(                                            & !
                                TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
                                NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                                iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                                rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                                time,DeltaT,                               & ! IN: time
                                DISC,EQN,MESH,MPI,IO)
       
        CASE(6) ! Coulomb model for LSW and bimaterial
        
           CALL Linear_slip_weakening_bimaterial(                          & !
                                TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
                                NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                                iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                                rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                                time,DeltaT,                               & ! IN: time
                                DISC,EQN,MESH,MPI,IO)

        CASE(7) ! severe velocity weakening friction as in Ampuero&Ben-Zion2008

        CALL rate_and_state_vw(                                            & !
                                TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
                                NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                                iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                                rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                                time,DeltaT,                               & ! IN: time
                                DISC,EQN,MESH,MPI,IO)
                       
        CASE(101) ! Specific conditions for SCEC TPV101
                      ! as case 3 (rate-and-state friction) aging law
                      ! + time and space dependent nucleation
            
           logError(*) 'Currently disabled'
!~            CALL rate_and_state_nuc101(                                     & !
!~                                 TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
!~                                 NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
!~                                 iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
!~                                 rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
!~                                 time,DeltaT,iT,                            & ! IN: time, inv Trafo
!~                                 DISC,EQN,MESH,MPI,IO,BND)

        CASE(103) ! Specific conditions for SCEC TPV103
                      ! Fast velocity-weakening friction with slip law
                      ! + time and space dependent nucleation

            CALL rate_and_state_nuc103(                                     & !
                                 TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
                                 NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                                 iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                                 rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                                 time,DeltaT,                               & ! IN: time, inv Trafo
                                 DISC,EQN,MESH,MPI,IO,BND)


        CASE DEFAULT
          logError(*) 'ERROR in friction.f90: friction law case',EQN%FL,' not implemented!'
          STOP
    END SELECT    

  END SUBROUTINE Eval_friction_law

  !> case 0: no frictional sliding
  !<
  PURE SUBROUTINE no_fault(TractionGP_XY,TractionGP_XZ,XYStressGP,XZStressGP)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Local variable declaration 
    REAL          :: XYStressGP(:,:)
    REAL          :: XZStressGP(:,:)
    REAL          :: TractionGP_XY(:,:)
    REAL          :: TractionGP_XZ(:,:)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: XYStressGP,XZStressGP
    INTENT(INOUT) :: TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------! 

    ! LocTrac = ShTest
    TractionGP_XY(:,:) = XYStressGP(:,:)
    TractionGP_XZ(:,:) = XZStressGP(:,:)          
    
  END SUBROUTINE no_fault


!> Special friction case 6: linear slip weakening with Prakash-Clifton regularization
!<
  SUBROUTINE Linear_slip_weakening_bimaterial(                                & !
                                   TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
                                   NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                                   iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                                   rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                                   time,DeltaT,                               & ! IN: time
                                   DISC,EQN,MESH,MPI,IO)
    !-------------------------------------------------------------------------!
    USE prak_clif_mod 
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI
    TYPE(tInputOutput)             :: IO
    ! Local variable declaration 
    INTEGER     :: iBndGP,iTimeGP,nBndGP,nTimeGP
    INTEGER     :: iFace,iSide,iElem
    REAL        :: LocTracXY,LocTracXZ
    REAL        :: NorStressGP(nBndGP,nTimeGP)
    REAL        :: XYStressGP(nBndGP,nTimeGP)
    REAL        :: XZStressGP(nBndGP,nTimeGP)
    REAL        :: TractionGP_XY(nBndGP,nTimeGP)
    REAL        :: TractionGP_XZ(nBndGP,nTimeGP)
    REAL        :: Strength_exp
    REAL        :: sigma
    REAL        :: LocMu, LocD_C, LocSlip, LocSlip1, LocSlip2, LocP, P, LocSR, ShTest
    REAL        :: LocMu_S, LocMu_D
    REAL        :: LocSR1,LocSR2
    REAL        :: P_0,cohesion, Strength
    REAL        :: rho,rho_neig,w_speed(:),w_speed_neig(:)
    REAL        :: time_inc,time
    REAL        :: Deltat(1:nTimeGP)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: NorStressGP,XYStressGP,XZStressGP,iFace,iSide,iElem
    INTENT(IN)    :: rho,rho_neig,w_speed,w_speed_neig,time,nBndGP,nTimeGP,DeltaT
    INTENT(IN)    :: EQN,MESH,MPI,IO
    INTENT(INOUT) :: DISC,TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------! 

    DO iBndGP=1,nBndGP
     !
     LocMu     = DISC%DynRup%Mu(iBndGP,iFace)
     LocMu_S   = DISC%DynRup%Mu_S(iBndGP,iFace)
     LocMu_D   = DISC%DynRup%Mu_D(iBndGP,iFace)
     LocD_C    = DISC%DynRup%D_C(iBndGP,iFace)
     LocSlip   = DISC%DynRup%Slip(iBndGP,iFace)
     LocSlip1   = DISC%DynRup%Slip1(iBndGP,iFace)
     LocSlip2   = DISC%DynRup%Slip2(iBndGP,iFace)
     LocSR1    = DISC%DynRup%SlipRate1(iBndGP,iFace)
     LocSR2    = DISC%DynRup%SlipRate2(iBndGP,iFace)
     cohesion  = DISC%DynRup%cohesion(iBndGP,iFace)      ! cohesion is negative since negative normal stress is compression
     P_0       = EQN%InitialStressInFaultCS(iBndGP,1,iFace)
     Strength_exp = DISC%DynRup%Strength(iBndGP,iFace)
     !
     DO iTimeGP=1,nTimeGP
       LocP   = NorStressGP(iBndGP,iTimeGP)
       time_inc = DeltaT(iTimeGP)
       !
       ! modify strength according to prakash clifton
       LocSR = SQRT(LocSR1**2 + LocSR2**2)
       sigma = LocP+P_0
       CALL prakash_cliff_fric(Strength_exp,sigma,LocSR,DISC%DynRup%v_star,DISC%DynRup%L,LocMu,time_inc)
        
       ShTest = SQRT((EQN%InitialStressInFaultCS(iBndGP,4,iFace) + XYStressGP(iBndGP,iTimeGP))**2 + (EQN%InitialStressInFaultCS(iBndGP,6,iFace) + XZStressGP(iBndGP,iTimeGP))**2)

       !Coulomb's law (we use old mu value, as mu, S, SR and Traction are interdependent!)
       IF(ShTest.GT.Strength) THEN

         ! 1 evaluate friction
         LocTracXY = ((EQN%InitialStressInFaultCS(iBndGP,4,iFace) + XYStressGP(iBndGP,iTimeGP))/ShTest)*Strength
         LocTracXZ = ((EQN%InitialStressInFaultCS(iBndGP,6,iFace) + XZStressGP(iBndGP,iTimeGP))/ShTest)*Strength
           
         ! 2 update stress change
         LocTracXY = LocTracXY - EQN%InitialStressInFaultCS(iBndGP,4,iFace)
         LocTracXZ = LocTracXZ - EQN%InitialStressInFaultCS(iBndGP,6,iFace)
           
       ELSE
         LocTracXY = XYStressGP(iBndGP,iTimeGP)
         LocTracXZ = XZStressGP(iBndGP,iTimeGP)
       ENDIF
       !
       !Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
       LocSR1     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXY-XYStressGP(iBndGP,iTimeGP))
       LocSR2     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXZ-XZStressGP(iBndGP,iTimeGP))
       LocSR      = SQRT(LocSR1**2 + LocSR2**2)
       !
       ! Update slip
       LocSlip1 = LocSlip1 + LocSR1*time_inc
       LocSlip2 = LocSlip2 + LocSR2*time_inc
       LocSlip = LocSlip + LocSR*time_inc
       !
       IF(ABS(LocSlip).LT.LocD_C) THEN
         LocMu = LocMu_S - (LocMu_S-LocMu_D)/LocD_C*ABS(LocSlip)
       ELSE
         LocMu = LocMu_D
       ENDIF

       ! instantaneous healing
       IF (DISC%DynRup%inst_healing == 1) THEN
           IF (LocSR .LT. u_0) THEN
               LocMu = LocMu_S
               ! reset slip history for LSW
               LocSlip = 0.0D0
           ENDIF
       ENDIF           
       !
       !Save traction for flux computation
       TractionGP_XY(iBndGP,iTimeGP) = LocTracXY
       TractionGP_XZ(iBndGP,iTimeGP) = LocTracXZ           
       !
     ENDDO ! iTimeGP=1,DISC%Galerkin%nTimeGP
     !
     ! output rupture front 
     ! outside of iTimeGP loop in order to safe an 'if' in a loop
     ! this way, no subtimestep resolution possible
     IF (DISC%DynRup%RF(iBndGP,iFace) .AND. LocSR .GT. 0.001D0) THEN
        DISC%DynRup%rupture_time(iBndGP,iFace)=time
        DISC%DynRup%RF(iBndGP,iFace) = .FALSE.
     ENDIF

     !output time when shear stress is equal to the dynamic stress after rupture arrived
     !currently only for linear slip weakening
      IF ( (DISC%DynRup%rupture_time(iBndGP,iFace).GT.0.0) .AND. (DISC%DynRup%rupture_time(iBndGP,iFace) .LE. time)) THEN
          IF(DISC%DynRup%DS(iBndGP,iFace) .AND. ABS(LocSlip).GE.LocD_C) THEN
          DISC%DynRup%dynStress_time(iBndGP,iFace)=time
          DISC%DynRup%DS(iBndGP,iFace) = .FALSE.
          ENDIF
      ENDIF

     IF (LocSR.GT.DISC%DynRup%PeakSR(iBndGP,iFace)) THEN
        DISC%DynRup%PeakSR(iBndGP,iFace) = LocSR
     ENDIF
     !
     DISC%DynRup%Mu(iBndGP,iFace)        = LocMu
     DISC%DynRup%SlipRate1(iBndGP,iFace) = LocSR1
     DISC%DynRup%SlipRate2(iBndGP,iFace) = LocSR2
     DISC%DynRup%Slip(iBndGP,iFace)      = LocSlip
     DISC%DynRup%Slip1(iBndGP,iFace)     = LocSlip1
     DISC%DynRup%Slip2(iBndGP,iFace)     = LocSlip2
     DISC%DynRup%TracXY(iBndGP,iFace)    = LocTracXY
     DISC%DynRup%TracXZ(iBndGP,iFace)    = LocTracXZ
     DISC%DynRup%Strength(iBndGP,iFace)  = Strength_exp
     !
    ENDDO ! iBndGP=1,DISC%Galerkin%nBndGP

  END SUBROUTINE Linear_slip_weakening_bimaterial

  
  !> friction case 16,17
  !> Specific conditions for SCEC TPV16/17
  !> basically, introduction of a time dependent forced rupture
  !<
  SUBROUTINE Linear_slip_weakening_TPV1617(TractionGP_XY,TractionGP_XZ,       & ! OUT: traction
                                   NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                                   iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                                   rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                                   time,DeltaT,                               & ! IN: time
                                   DISC,EQN,MESH,MPI,IO)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI
    TYPE(tInputOutput)             :: IO    
    ! Local variable declaration
    INTEGER     :: iBndGP,iTimeGP,nBndGP,nTimeGP
    INTEGER     :: iFace,iSide,iElem
    REAL        :: time
    REAL        :: NorStressGP(nBndGP,nTimeGP)
    REAL        :: XYStressGP(nBndGP,nTimeGP)
    REAL        :: XZStressGP(nBndGP,nTimeGP)
    REAL        :: TractionGP_XY(nBndGP,nTimeGP)
    REAL        :: TractionGP_XZ(nBndGP,nTimeGP)
    real        :: LocTracXY(nBndGP)
    real        :: LocTracXZ(nBndGP)
    real        :: tmpSlip(nBndGP)
    real        :: P(nBndGP)
    real        :: Strength(nBndGP), ShTest(nBndGP)
    real        :: LocSR(nBndGP)
    real        :: srFactor
    REAL        :: rho,rho_neig,w_speed(:),w_speed_neig(:)
    REAL        :: time_inc
    REAL        :: Deltat(1:nTimeGP)
    REAL        :: t_0
    REAL        :: f1(nBndGP), f2(nBndGP)
    real        :: tn
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: NorStressGP,XYStressGP,XZStressGP,iFace,iSide,iElem
    INTENT(IN)    :: rho,rho_neig,w_speed,w_speed_neig,time,nBndGP,nTimeGP,DeltaT
    INTENT(IN)    :: EQN,MESH,MPI,IO
    INTENT(INOUT) :: DISC,TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------! 
    t_0 = DISC%DynRup%t_0
    tmpSlip = 0.0D0
    
    srFactor = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))
    tn = time
    
    do iTimeGP=1,nTimeGP
      time_inc = DeltaT(iTimeGP)
      tn=tn + time_inc
      
      P = EQN%InitialStressInFaultCS(:,1,iFace) + NorStressGP(:,iTimeGP)
      
      Strength = -DISC%DynRup%cohesion(:,iFace) - DISC%DynRup%Mu(:,iFace) * MIN(P,ZERO)      
      ShTest = SQRT((EQN%InitialStressInFaultCS(:,4,iFace) + XYStressGP(:,iTimeGP))**2 + (EQN%InitialStressInFaultCS(:,6,iFace) + XZStressGP(:,iTimeGP))**2)
      
      where (ShTest > Strength)
        ! 1 evaluate friction
        LocTracXY = (EQN%InitialStressInFaultCS(:,4,iFace) + XYStressGP(:,iTimeGP)) / ShTest(:) * Strength(:)
        LocTracXZ = (EQN%InitialStressInFaultCS(:,6,iFace) + XZStressGP(:,iTimeGP)) / ShTest(:) * Strength(:)
        
        ! 2 update stress change
        LocTracXY = LocTracXY(:) - EQN%InitialStressInFaultCS(:,4,iFace)
        LocTracXZ = LocTracXZ(:) - EQN%InitialStressInFaultCS(:,6,iFace)
      elsewhere
        LocTracXY = XYStressGP(:,iTimeGP)
        LocTracXZ = XZStressGP(:,iTimeGP)        
      end where
      
      ! Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
      DISC%DynRup%SlipRate1(:,iFace)     = srFactor*(LocTracXY(:)-XYStressGP(:,iTimeGP))
      DISC%DynRup%SlipRate2(:,iFace)     = srFactor*(LocTracXZ(:)-XZStressGP(:,iTimeGP))
      LocSR                              = SQRT(DISC%DynRup%SlipRate1(:,iFace)**2 + DISC%DynRup%SlipRate2(:,iFace)**2)
      
      ! Update slip
      DISC%DynRup%Slip1(:,iFace) = DISC%DynRup%Slip1(:,iFace) + DISC%DynRup%SlipRate1(:,iFace)*time_inc
      DISC%DynRup%Slip2(:,iFace) = DISC%DynRup%Slip2(:,iFace) + DISC%DynRup%SlipRate2(:,iFace)*time_inc
      DISC%DynRup%Slip(:,iFace)  = DISC%DynRup%Slip(:,iFace)  + LocSR(:)*time_inc      
      tmpSlip = tmpSlip(:) + LocSR(:)*time_inc
      
     ! Modif T. Ulrich-> generalisation of tpv16/17 to 30/31
     f1=dmin1(ABS(DISC%DynRup%Slip(:,iFace))/DISC%DynRup%D_C(:,iFace),1d0)    

     IF(EQN%FL.EQ.16) THEN 
        IF (t_0.eq.0) THEN
         where (tn >= DISC%DynRup%forced_rupture_time(:,iFace))
            f2=1.
         elsewhere
            f2=0.
         end where
        ELSE
           f2=dmax1(0d0,dmin1((time-DISC%DynRup%forced_rupture_time(:,iFace))/t_0,1d0))
        ENDIF
     ELSE !no forced time rupture
        f2=0.
     ENDIF

     DISC%DynRup%Mu(:,iFace) = DISC%DynRup%Mu_S(:,iFace) - (DISC%DynRup%Mu_S(:,iFace)-DISC%DynRup%Mu_D(:,iFace))*dmax1(f1,f2)

     ! instantaneous healing
     IF (DISC%DynRup%inst_healing == 1) THEN
        where (LocSR.LT. u_0)
           DISC%DynRup%Mu(:,iFace) = DISC%DynRup%Mu_S(:,iFace)
           DISC%DynRup%Slip(:,iFace)  = 0.0
        endwhere
     ENDIF
    
     TractionGP_XY(:,iTimeGP) = LocTracXY(:)
     TractionGP_XZ(:,iTimeGP) = LocTracXZ(:)      
    enddo

     ! output rupture front 
     ! outside of iTimeGP loop in order to safe an 'if' in a loop
     ! this way, no subtimestep resolution possible
    where (DISC%DynRup%RF(:,iFace) .AND. LocSR .GT. 0.001D0)
      DISC%DynRup%rupture_time(:,iFace)=time
      DISC%DynRup%RF(:,iFace) = .FALSE.
    end where

    !output time when shear stress is equal to the dynamic stress after rupture arrived
    !currently only for linear slip weakening
    where ( (DISC%DynRup%rupture_time(:,iFace) .GT. 0.0) .AND. &
            (DISC%DynRup%rupture_time(:,iFace) .LE. time) .AND. &
             DISC%DynRup%DS(:,iFace) .AND. &
            (ABS(DISC%DynRup%Slip(:,iFace)) .GE. DISC%DynRup%D_C(:,iFace)))
      DISC%DynRup%dynStress_time(:,iFace)=time
      DISC%DynRup%DS(:,iFace) = .FALSE.
    end where

    where (LocSR(:).GT.DISC%DynRup%PeakSR(:,iFace))
      DISC%DynRup%PeakSR(:,iFace) = LocSR
    end where

    DISC%DynRup%TracXY(:,iFace)    = LocTracXY
    DISC%DynRup%TracXZ(:,iFace)    = LocTracXZ

    !---compute and store slip to determine the magnitude of an earthquake ---
    !    to this end, here the slip is computed and averaged per element
    !    in calc_seissol.f90 this value will be multiplied by the element surface
    !    and an output happened once at the end of the simulation
    IF (DISC%DynRup%magnitude_out(iFace)) THEN
        DISC%DynRup%averaged_Slip(iFace) = DISC%DynRup%averaged_Slip(iFace) + sum(tmpSlip)/nBndGP
    ENDIF

  END SUBROUTINE Linear_slip_weakening_TPV1617


  !> friction case 3,4: rate and state friction
  !> aging (3) and slip law (4)
  !<
  SUBROUTINE rate_and_state(TractionGP_XY,TractionGP_XZ,               & ! OUT: traction
                            NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                            iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                            rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                            time,DeltaT,                               & ! IN: time
                            DISC,EQN,MESH,MPI,IO)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization)          :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI
    TYPE(tInputOutput)             :: IO    
    ! Local variable declaration
    INTEGER     :: i,j
    INTEGER     :: iBndGP,iTimeGP,nBndGP,nTimeGP
    INTEGER     :: iFace,iSide,iElem
    INTEGER     :: nSRupdates, nSVupdates, SignSR
    REAL        :: time  
    REAL        :: LocTracXY,LocTracXZ
    REAL        :: NorStressGP(nBndGP,nTimeGP)
    REAL        :: XYStressGP(nBndGP,nTimeGP)
    REAL        :: XZStressGP(nBndGP,nTimeGP)
    REAL        :: TractionGP_XY(nBndGP,nTimeGP)
    REAL        :: TractionGP_XZ(nBndGP,nTimeGP)
    REAL        :: LocMu, LocD_C, LocSlip, LocSlip1, LocSlip2, LocP, P, LocSR, ShTest
    REAL        :: LocMu_S, LocMu_D
    REAL        :: LocSR1,LocSR2
    REAL        :: P_0,Strength,cohesion
    REAL        :: rho,rho_neig,w_speed(:),w_speed_neig(:)
    REAL        :: time_inc
    REAL        :: Deltat(1:nTimeGP)
    REAL        :: SV0, tmp, tmp2, SRtest, NR, dNR
    REAL        :: LocSV
    REAL        :: RS_f0,RS_a,RS_b,RS_sl0,RS_sr0
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: NorStressGP,XYStressGP,XZStressGP,iFace,iSide,iElem
    INTENT(IN)    :: rho,rho_neig,w_speed,w_speed_neig,time,nBndGP,nTimeGP,DeltaT
    INTENT(IN)    :: EQN,MESH,MPI,IO
    INTENT(INOUT) :: DISC,TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------! 

    DO iBndGP=1,nBndGP
     !
     LocSlip   = DISC%DynRup%Slip(iBndGP,iFace)
     LocSlip1   = DISC%DynRup%Slip1(iBndGP,iFace)
     LocSlip2   = DISC%DynRup%Slip2(iBndGP,iFace)
     LocSR1    = DISC%DynRup%SlipRate1(iBndGP,iFace)
     LocSR2    = DISC%DynRup%SlipRate2(iBndGP,iFace)
     LocSV     = DISC%DynRup%StateVar(iBndGP,iFace)
     P_0       = EQN%InitialStressInFaultCS(iBndGP,1,iFace)
     cohesion  = DISC%DynRup%cohesion(iBndGP,iFace)      ! cohesion is negative since negative normal stress is compression
     !
     !logInfo(*) 'state variable evaluation', DISC%DynRup%StateVar(iBndGP,iFace), EQN%IniStateVar(iBndGP,iFace), 'iGP', iBndGP
     DO iTimeGP=1,nTimeGP
       LocP   = NorStressGP(iBndGP,iTimeGP)
       time_inc = DeltaT(iTimeGP)
       !
       RS_f0  = DISC%DynRup%RS_f0 
       RS_a   = DISC%DynRup%RS_a  
       RS_b   = DISC%DynRup%RS_b  
       RS_sl0 = DISC%DynRup%RS_sl0
       RS_sr0 = DISC%DynRup%RS_sr0
       !       
       !SignSR1   = SIGN(1.0,LocSR1)                    ! Gets the sign of the slip rate
       !SignSR2   = SIGN(1.0,LocSR2)                    ! Gets the sign of the slip rate
       !
       ! load traction and normal stress
       P      = LocP+P_0
       ShTest = SQRT((EQN%InitialStressInFaultCS(iBndGP,4,iFace) + XYStressGP(iBndGP,iTimeGP))**2 + (EQN%InitialStressInFaultCS(iBndGP,6,iFace) + XZStressGP(iBndGP,iTimeGP))**2)
       !
       ! We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
       ! ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )
       !           
       SV0=LocSV    ! Careful, the SV must always be corrected using SV0 and not LocSV!
       
       !
       ! The following process is adapted from that described by Kaneko et al. (2008)
       nSRupdates = 5
       nSVupdates = 2
       !
       LocSR      = SQRT(LocSR1**2 + LocSR2**2)
       tmp        = ABS(LocSR)
       !
       DO j=1,nSVupdates   !This loop corrects SV values
         !
         LocSR=ABS(LocSR)
         !
         IF(EQN%FL.EQ.3) THEN         ! aging law
           LocSV=SV0*EXP(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-EXP(-tmp*time_inc/RS_sl0))
         ELSEIF(EQN%FL.EQ.4) THEN     ! slip law
           LocSV=RS_sl0/tmp*(tmp*SV0/RS_sl0)**(EXP(-tmp*time_inc/RS_sl0))                 
         ENDIF
         !
         ! Newton-Raphson algorithm to determine the value of the slip rate.
         ! We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
         !  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
         ! In our case we equalize the values of the traction for two equations:
         !             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
         !             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
         !               where mu=a*asinh(SR/2/SR0*exp((F0+b*log(SR0*SV/L))/a (eq. 2a of Lapusta and Rice (2003))
         !
         SRtest=LocSR  ! We use as first guess the SR value of the previous time step
         !              
         DO i=1,nSRupdates  !This loop corrects SR values
           tmp          = 0.5/RS_sr0* EXP( (RS_f0+RS_b*LOG(RS_sr0*LocSV/RS_sl0) ) /RS_a)
           tmp2         = tmp*SRtest
           NR           = -(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) * &
                        (ABS(P)*RS_a*LOG(tmp2+SQRT(tmp2**2+1.0))-ShTest)-SRtest !not sure if ShTest should be + or -...        
           dNR          = -(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) * &
                        (ABS(P)*RS_a/SQRT(1+tmp2**2)*tmp)-1.0
           SRtest = ABS(SRtest-NR/dNR) ! no ABS needed around NR/dNR at least for aging law
         ENDDO
         tmp=0.5*(LocSR+ABS(SRtest))  ! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
         LocSR=ABS(SRtest)               
           
       ENDDO !  j=1,nSVupdates   !This loop corrects SV values
       
       !                    
       IF(EQN%FL.EQ.3) THEN         ! aging law
         LocSV=SV0*EXP(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-EXP(-tmp*time_inc/RS_sl0))
       ELSEIF(EQN%FL.EQ.4) THEN     ! slip law
         LocSV=RS_sl0/tmp*(tmp*SV0/RS_sl0)**(EXP(-tmp*time_inc/RS_sl0))                 
       ENDIF
       !               
       tmp  = 0.5 * (LocSR)/RS_sr0 * EXP((RS_f0 + RS_b*LOG(RS_sr0*LocSV/RS_sl0)) / RS_a)  
       !
       LocMu    = RS_a * LOG(tmp + SQRT(tmp**2 + 1.0))
       ! 2D:
       !LocTrac  = -(ABS(S_0)-LocMu*(LocP+P_0))*(S_0/ABS(S_0))
       !LocTrac  = ABS(LocTrac)*(-SignSR)  !!! line commented as it leads NOT to correct results
       ! update stress change
       LocTracXY = -((EQN%InitialStressInFaultCS(iBndGP,4,iFace) + XYStressGP(iBndGP,iTimeGP))/ShTest)*(LocMu*P+ABS(cohesion))
       LocTracXZ = -((EQN%InitialStressInFaultCS(iBndGP,6,iFace) + XZStressGP(iBndGP,iTimeGP))/ShTest)*(LocMu*P+ABS(cohesion))
       LocTracXY = LocTracXY - EQN%InitialStressInFaultCS(iBndGP,4,iFace)
       LocTracXZ = LocTracXZ - EQN%InitialStressInFaultCS(iBndGP,6,iFace)
       !
       ! Compute slip
       LocSlip   = LocSlip  + (LocSR)*time_inc ! ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening
       !
       !Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
       LocSR1     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXY-XYStressGP(iBndGP,iTimeGP))
       LocSR2     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXZ-XZStressGP(iBndGP,iTimeGP))

       LocSlip1   = LocSlip1  + (LocSR1)*time_inc 
       LocSlip2   = LocSlip2  + (LocSR2)*time_inc 

       !LocSR1     = SignSR1*ABS(LocSR1)
       !LocSR2     = SignSR2*ABS(LocSR2)           
       !
       !Save traction for flux computation
       TractionGP_XY(iBndGP,iTimeGP) = LocTracXY
       TractionGP_XZ(iBndGP,iTimeGP) = LocTracXZ           
       !
     ENDDO ! iTimeGP=1,DISC%Galerkin%nTimeGP
      
     !
     ! output rupture front 
     ! outside of iTimeGP loop in order to safe an 'if' in a loop
     ! this way, no subtimestep resolution possible
     IF (DISC%DynRup%RF(iBndGP,iFace) .AND. LocSR .GT. 0.001D0) THEN
        DISC%DynRup%rupture_time(iBndGP,iFace)=time
        DISC%DynRup%RF(iBndGP,iFace) = .FALSE.
     ENDIF
     IF (LocSR.GT.DISC%DynRup%PeakSR(iBndGP,iFace)) THEN
        DISC%DynRup%PeakSR(iBndGP,iFace) = LocSR
     ENDIF
     !
     DISC%DynRup%Mu(iBndGP,iFace)        = LocMu
     DISC%DynRup%SlipRate1(iBndGP,iFace) = LocSR1
     DISC%DynRup%SlipRate2(iBndGP,iFace) = LocSR2
     DISC%DynRup%Slip(iBndGP,iFace)      = LocSlip
     DISC%DynRup%Slip1(iBndGP,iFace)     = LocSlip1
     DISC%DynRup%Slip2(iBndGP,iFace)     = LocSlip2
     DISC%DynRup%StateVar(iBndGP,iFace)  = LocSV
     DISC%DynRup%TracXY(iBndGP,iFace)    = LocTracXY
     DISC%DynRup%TracXZ(iBndGP,iFace)    = LocTracXZ

     !
    ENDDO ! iBndGP=1,DISC%Galerkin%nBndGP
  END SUBROUTINE rate_and_state

  !> friction case 7: severe velocity weakening rate and state friction
  !< after Ampuero and Ben-Zion 2008
  !<
  SUBROUTINE rate_and_state_vw(TractionGP_XY,TractionGP_XZ,            & ! OUT: traction
                            NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                            iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                            rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                            time,DeltaT,                               & ! IN: time
                            DISC,EQN,MESH,MPI,IO)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization)          :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI
    TYPE(tInputOutput)             :: IO
    ! Local variable declaration
    INTEGER     :: i,j
    INTEGER     :: iBndGP,iTimeGP,nBndGP,nTimeGP
    INTEGER     :: iFace,iSide,iElem
    INTEGER     :: nSRupdates, nSVupdates, SignSR
    REAL        :: time
    REAL        :: LocTracXY,LocTracXZ
    REAL        :: NorStressGP(nBndGP,nTimeGP)
    REAL        :: XYStressGP(nBndGP,nTimeGP)
    REAL        :: XZStressGP(nBndGP,nTimeGP)
    REAL        :: TractionGP_XY(nBndGP,nTimeGP)
    REAL        :: TractionGP_XZ(nBndGP,nTimeGP)
    REAL        :: LocMu, LocD_C, LocSlip, LocSlip1, LocSlip2, LocP, P, LocSR, ShTest
    REAL        :: LocMu_S, LocMu_D
    REAL        :: LocSR1,LocSR2
    REAL        :: P_0,Strength,cohesion
    REAL        :: rho,rho_neig,w_speed(:),w_speed_neig(:)
    REAL        :: time_inc
    REAL        :: Deltat(1:nTimeGP)
    REAL        :: SV0, tmp, tmp2, SRtest, NR, dNR
    REAL        :: LocSV
    REAL        :: RS_f0,RS_a,RS_b,RS_sl0,RS_sr0, Tc, coeft
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: NorStressGP,XYStressGP,XZStressGP,iFace,iSide,iElem
    INTENT(IN)    :: rho,rho_neig,w_speed,w_speed_neig,time,nBndGP,nTimeGP,DeltaT
    INTENT(IN)    :: EQN,MESH,MPI,IO
    INTENT(INOUT) :: DISC,TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------!
    ! friction develops as                    mu = mu_s + a V/(V+Vc) - b SV/(SV + Dc)
    ! Note the typo in eq.1 of Ampuero&Ben-Zion, 2008
    ! state variable SV develops as     dSV / dt = (V-SV) / Tc
    ! parameters: static friction mu_s, char. velocity scale Vc, charact. timescale Tc,
    ! charact. length scale Dc, direct and evolution effect coeff. a,b
    ! Notice that Dc, a and b are recycled but not equivalent to cases 3 and 4
    ! steady-state friction value is:       mu = mu_s + (a - b) V/(V+Vc)
    ! dynamic friction value (if reached) mu_d = mu_s + (a - b)
    ! Tc tunes between slip-weakening and rate-weakening behavior
    !

    DO iBndGP=1,nBndGP
     !
     LocSlip   = DISC%DynRup%Slip(iBndGP,iFace)
     LocSlip1   = DISC%DynRup%Slip1(iBndGP,iFace)
     LocSlip2   = DISC%DynRup%Slip2(iBndGP,iFace)
     LocSR1    = DISC%DynRup%SlipRate1(iBndGP,iFace)
     LocSR2    = DISC%DynRup%SlipRate2(iBndGP,iFace)
     LocSV     = DISC%DynRup%StateVar(iBndGP,iFace)
     P_0       = EQN%InitialStressInFaultCS(iBndGP,1,iFace)
     !
     DO iTimeGP=1,nTimeGP
       LocP   = NorStressGP(iBndGP,iTimeGP)
       time_inc = DeltaT(iTimeGP)

       !
       RS_f0  = DISC%DynRup%RS_f0 ! equivalent to static friction coefficient
       RS_a   = DISC%DynRup%RS_a  ! direct effect
       RS_b   = DISC%DynRup%RS_b  ! evolution effect
       RS_sl0 = DISC%DynRup%RS_sl0 ! Dc, char. lengt scale
       RS_sr0 = DISC%DynRup%RS_sr0 ! Vc, char. velocity scale
       !
       ! load traction and normal stress
       P      = LocP+P_0
       ShTest = SQRT((EQN%InitialStressInFaultCS(iBndGP,4,iFace) + XYStressGP(iBndGP,iTimeGP))**2 + (EQN%InitialStressInFaultCS(iBndGP,6,iFace) + XZStressGP(iBndGP,iTimeGP))**2)
       !
       SV0=LocSV    ! Careful, the SV must always be corrected using SV0 and not LocSV!
       !
       ! The following process is adapted from that described by Kaneko et al. (2008)
       nSRupdates = 5
       nSVupdates = 2
       !
       LocSR      = SQRT(LocSR1**2 + LocSR2**2)
       !   charact. time scale Tc
       Tc = RS_sl0 / RS_sr0
       !   exponent
       coeft= EXP(-time_inc / Tc)
       !
       DO j=1,nSVupdates   !This loop corrects SV values
         !
         LocSR=ABS(LocSR)
         !   exact integration assuming constant V in this loop
         LocSV=Tc*LocSR*(1d0-coeft) + coeft*SV0
         !
         ! Newton-Raphson algorithm to determine the value of the slip rate.
         ! We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
         !  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
         ! In our case we equalize the values of the traction for two equations:
         !             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
         !             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
         !             where mu = mu_s + a V/(V+Vc) - b SV/(SV + Vc)
         !
         SRtest=LocSR  ! We use as first guess the SR value of the previous time step
         !
         DO i=1,nSRupdates  !This loop corrects SR values
           tmp          = RS_f0+RS_a*SRtest/(SRtest+RS_sr0)-RS_b*LocSV/(LocSV+RS_sl0)   !=mu
           NR           = -(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) * &
                        (ABS(P)*tmp-ShTest)-SRtest
           dNR          = -(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) * &
                        (ABS(P)*(RS_a/(SRtest+RS_sr0)-RS_a*SRtest/(SRtest+RS_sr0)**2)) -1.0
           SRtest = SRtest-NR/dNR
         ENDDO
         tmp=0.5*(LocSR+ABS(SRtest))  ! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
         LocSR=ABS(SRtest)

       ENDDO !  j=1,nSVupdates   !This loop corrects SV values
       !
       LocSV    = Tc*tmp*(1d0-coeft) + coeft*SV0
       !
       tmp  = 0.5 * (LocSR)/RS_sr0 * EXP((RS_f0 + RS_b*LOG(RS_sr0*LocSV/RS_sl0)) / RS_a)
       !
       LocMu    = RS_f0+RS_a*LocSR/(LocSR+RS_sr0)-RS_b*LocSV/(LocSV+RS_sl0)
       !
       ! update stress change
       LocTracXY = -((EQN%InitialStressInFaultCS(iBndGP,4,iFace) + XYStressGP(iBndGP,iTimeGP))/ShTest)*LocMu*P
       LocTracXZ = -((EQN%InitialStressInFaultCS(iBndGP,6,iFace) + XZStressGP(iBndGP,iTimeGP))/ShTest)*LocMu*P
       LocTracXY = LocTracXY - EQN%InitialStressInFaultCS(iBndGP,4,iFace)
       LocTracXZ = LocTracXZ - EQN%InitialStressInFaultCS(iBndGP,6,iFace)
       !
       ! Compute slip
       LocSlip   = LocSlip  + (LocSR)*time_inc ! ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening
       !
       !Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
       LocSR1     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXY-XYStressGP(iBndGP,iTimeGP))
       LocSR2     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXZ-XZStressGP(iBndGP,iTimeGP))

       LocSlip1   = LocSlip1  + (LocSR1)*time_inc 
       LocSlip2   = LocSlip2  + (LocSR2)*time_inc 
       !LocSR1     = SignSR1*ABS(LocSR1)
       !LocSR2     = SignSR2*ABS(LocSR2)
       !
       !Save traction for flux computation
       TractionGP_XY(iBndGP,iTimeGP) = LocTracXY
       TractionGP_XZ(iBndGP,iTimeGP) = LocTracXZ
       !
     ENDDO ! iTimeGP=1,DISC%Galerkin%nTimeGP
     !
     ! output rupture front
     ! outside of iTimeGP loop in order to safe an 'if' in a loop
     ! this way, no subtimestep resolution possible
     IF (DISC%DynRup%RF(iBndGP,iFace) .AND. LocSR .GT. 0.001D0) THEN
        DISC%DynRup%rupture_time(iBndGP,iFace)=time
        DISC%DynRup%RF(iBndGP,iFace) = .FALSE.
     ENDIF
     IF (LocSR.GT.DISC%DynRup%PeakSR(iBndGP,iFace)) THEN
        DISC%DynRup%PeakSR(iBndGP,iFace) = LocSR
     ENDIF
     !
     DISC%DynRup%Mu(iBndGP,iFace)        = LocMu
     DISC%DynRup%SlipRate1(iBndGP,iFace) = LocSR1
     DISC%DynRup%SlipRate2(iBndGP,iFace) = LocSR2
     DISC%DynRup%Slip(iBndGP,iFace)      = LocSlip
     DISC%DynRup%Slip1(iBndGP,iFace)     = LocSlip1
     DISC%DynRup%Slip2(iBndGP,iFace)     = LocSlip2
     DISC%DynRup%TracXY(iBndGP,iFace)    = LocTracXY
     DISC%DynRup%TracXZ(iBndGP,iFace)    = LocTracXZ
     DISC%DynRup%StateVar(iBndGP,iFace)  = LocSV
     !
    ENDDO ! iBndGP=1,DISC%Galerkin%nBndGP

  END SUBROUTINE rate_and_state_vw


  !> special friction case for SCEC TPV101: rate and state friction
  !> aging law
  !< with time and space dependent nucleation
  SUBROUTINE rate_and_state_nuc101(TractionGP_XY,TractionGP_XZ,        & ! OUT: traction
                            NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                            iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                            rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                            time,DeltaT,iT,                            & ! IN: time, inv Trafo
                            DISC,EQN,MESH,MPI,IO,BND)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization)          :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI
    TYPE(tInputOutput)             :: IO
    TYPE (tBoundary)               :: BND
    ! Local variable declaration
    INTEGER     :: i,j
    INTEGER     :: iBndGP,iTimeGP,nBndGP,nTimeGP
    INTEGER     :: iFace,iSide,iElem
    INTEGER     :: nSRupdates, nSVupdates, SignSR
    INTEGER     :: iNeighbor, iLocalNeighborSide
    INTEGER     :: MPIIndex, iObject
    REAL        :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
    REAL        :: time
    REAL        :: LocTracXY,LocTracXZ
    REAL        :: NorStressGP(nBndGP,nTimeGP)
    REAL        :: XYStressGP(nBndGP,nTimeGP)
    REAL        :: XZStressGP(nBndGP,nTimeGP)
    REAL        :: TractionGP_XY(nBndGP,nTimeGP)
    REAL        :: TractionGP_XZ(nBndGP,nTimeGP)
    REAL        :: iT(:,:)                                      ! inverse Transformation matrix    !
    REAL        :: Stress(6,1:nBndGP)
    REAL        :: LocMu, LocD_C, LocSlip, LocSlip1, LocSlip2, LocP, P, LocSR, ShTest
    REAL        :: LocMu_S, LocMu_D
    REAL        :: LocSR1,LocSR2
    REAL        :: P_0,Strength,cohesion
    REAL        :: rho,rho_neig,w_speed(:),w_speed_neig(:)
    REAL        :: time_inc
    REAL        :: Deltat(1:nTimeGP)
    REAL        :: SV0, tmp, tmp2, SRtest, NR, dNR
    REAL        :: LocSV
    REAL        :: RS_f0,RS_a,RS_b,RS_sl0,RS_sr0
    REAL        :: chi, tau, xi, eta, zeta, XGp, YGp, ZGp
    REAL        :: Rnuc, Tnuc, radius, Gnuc, Fnuc
    LOGICAL     :: nodewise=.FALSE.
    REAL        :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
    INTEGER     :: VertexSide(4,3)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: NorStressGP,XYStressGP,XZStressGP,iFace,iSide,iElem
    INTENT(IN)    :: rho,rho_neig,w_speed,w_speed_neig,time,nBndGP,nTimeGP,DeltaT
    INTENT(IN)    :: MESH,MPI,IO
    INTENT(INOUT) :: EQN,DISC,TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------!
    ! switch for Gauss node wise stress assignment
    nodewise = .TRUE.

    !Apply time dependent nucleation at global time step not sub time steps for simplicity
    !initialize time and space dependent nucleation
    Rnuc=3000.0D0
    Tnuc=1.0D0
    !
    IF (time.LE.Tnuc) THEN
       IF (nodewise) THEN        
    !
         ! Gauss node coordinate definition and stress assignment
         ! get vertices of complete tet
            IF (MESH%Fault%Face(iFace,1,1) == 0) THEN
                ! iElem is in the neighbor domain
                ! The neighbor element belongs to a different MPI domain
                iNeighbor           = MESH%Fault%Face(iFace,1,2)          ! iNeighbor denotes "-" side
                iLocalNeighborSide  = MESH%Fault%Face(iFace,2,2)
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
            DO iBndGP = 1,nBndGP
                !
                ! Transformation of boundary GP's into XYZ coordinate system
                chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
                tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
                CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
                CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
                !
                !radial distance to hypocenter
                radius=SQRT(xGP**2+(zGP+7500.0D0)**2)
                ! Inside nucleation patch add shear stress perturbation of 25 MPa along strike
                IF (radius.LT.Rnuc) THEN
                    Fnuc=EXP(radius**2/(radius**2-Rnuc**2))
                    IF (time.GT.0.0D0) THEN
                        Gnuc=EXP((time-Tnuc)**2/(time*(time-2.0D0*Tnuc)))
                    ELSE
                        Gnuc=0.0D0
                    ENDIF
                    EQN%IniShearXY(iFace,iBndGP)=EQN%ShearXY_0+25.0e6*Fnuc*Gnuc
                ENDIF ! Rnuc
                !
            ENDDO ! iBndGP
        !       
       ELSE
       ! get coordinates needed for nucleation zone
            IF (iElem .NE. 0) THEN
                !
                DO j=1,3
                    xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
                    yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
                    zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
                ENDDO           
            ELSEIF (iElem == 0) THEN ! in case "+" element is not present in the local domain
                !
                iLocalNeighborSide = MESH%Fault%Face(iFace,2,2)
                DO j=1,3
                    xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(iFace,1,2)))
                    yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(iFace,1,2)))
                    zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iLocalNeighborSide,j),MESH%Fault%Face(iFace,1,2)))
                ENDDO            
            ENDIF
            !max radial distance to hypocenter
            radius=SQRT( (MAXVAL(xp(1:3)))**2 + ((MAXVAL(zp(1:3))+7500.0D0))**2 ) 
                    ! Inside nucleation patch add shear stress perturbation of 25 MPa along strike
              IF (radius.LT.Rnuc) THEN
                    Fnuc=EXP(radius**2/(radius**2-Rnuc**2))
                    IF (time.GT.0.0D0) THEN
                        Gnuc=EXP((time-Tnuc)**2/(time*(time-2.0D0*Tnuc)))
                    ELSE
                        Gnuc=0.0D0
                    ENDIF
                    EQN%IniShearXY(iFace,:)=EQN%ShearXY_0+25.0e6*Fnuc*Gnuc
            ENDIF ! Rnuc
        ENDIF ! nodewise
    ENDIF ! Tnuc
    !
    !Background stress rotation to face's reference system
    !
    Stress(1,:)=EQN%IniBulk_xx(:,iFace)
    Stress(2,:)=EQN%IniBulk_yy(:,iFace)
    Stress(3,:)=EQN%IniBulk_zz(:,iFace)
    Stress(4,:)=EQN%IniShearXY(:,iFace)
    Stress(5,:)=EQN%IniShearYZ(:,iFace)
    Stress(6,:)=EQN%IniShearXZ(:,iFace)
    !
    DO iBndGP=1,nBndGP
       Stress(:,iBndGP)=MATMUL(iT(1:6,1:6),Stress(:,iBndGP))
    ENDDO
    !
    DO iBndGP=1,nBndGP
     !
     LocSlip   = DISC%DynRup%Slip(iBndGP,iFace)
     LocSlip1   = DISC%DynRup%Slip1(iBndGP,iFace)
     LocSlip2   = DISC%DynRup%Slip2(iBndGP,iFace)
     LocSR1    = DISC%DynRup%SlipRate1(iBndGP,iFace)
     LocSR2    = DISC%DynRup%SlipRate2(iBndGP,iFace)
     LocSV     = DISC%DynRup%StateVar(iBndGP,iFace)
     P_0       = Stress(1,iBndGP)
     !
     DO iTimeGP=1,nTimeGP
       LocP   = NorStressGP(iBndGP,iTimeGP)
       time_inc = DeltaT(iTimeGP)
       !
       RS_f0  = DISC%DynRup%RS_f0
       ! spatially variabel a	
       RS_a   = DISC%DynRup%RS_a_array(iFace,iBndGP)
       RS_b   = DISC%DynRup%RS_b
       RS_sl0 = DISC%DynRup%RS_sl0
       RS_sr0 = DISC%DynRup%RS_sr0
       !
       !SignSR1   = SIGN(1.0,LocSR1)                    ! Gets the sign of the slip rate
       !SignSR2   = SIGN(1.0,LocSR2)                    ! Gets the sign of the slip rate
       !
       ! load traction and normal stress
       P      = LocP+P_0
       ShTest = SQRT((Stress(4,iBndGP) + XYStressGP(iBndGP,iTimeGP))**2 + (Stress(6,iBndGP) + XZStressGP(iBndGP,iTimeGP))**2)
       !
       ! We use the regularized rate-and-state friction, after Rice & Ben-Zion (1996)
       ! ( Numerical note: ASINH(X)=LOG(X+SQRT(X^2+1)) )
       !
       SV0=LocSV    ! Careful, the SV must always be corrected using SV0 and not LocSV!
       !
       ! The following process is adapted from that described by Kaneko et al. (2008)
       nSRupdates = 5
       nSVupdates = 2
       !
       LocSR      = SQRT(LocSR1**2 + LocSR2**2)
       tmp        = ABS(LocSR)
       !
       DO j=1,nSVupdates   !This loop corrects SV values
         !
         LocSR=ABS(LocSR)
         !
         IF(EQN%FL.EQ.102) THEN         ! aging law
           LocSV=SV0*EXP(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-EXP(-tmp*time_inc/RS_sl0))
         ELSEIF(EQN%FL.EQ.104) THEN     ! slip law
           LocSV=RS_sl0/tmp*(tmp*SV0/RS_sl0)**(EXP(-tmp*time_inc/RS_sl0))
         ENDIF
         !
         ! Newton-Raphson algorithm to determine the value of the slip rate.
         ! We wish to find SR that fulfills g(SR)=f(SR), by building up the function NR=f-g , which has
         !  a derivative dNR = d(NR)/d(SR). We can then find SR by iterating SR_{i+1}=SR_i-( NR_i / dNR_i ).
         ! In our case we equalize the values of the traction for two equations:
         !             g = SR*mu/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
         !             f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
         !               where mu=a*asinh(SR/2/SR0*exp((F0+b*log(SR0*SV/L))/a (eq. 2a of Lapusta and Rice (2003))
         !
         SRtest=LocSR  ! We use as first guess the SR value of the previous time step
         !
         DO i=1,nSRupdates  !This loop corrects SR values
           tmp          = 0.5/RS_sr0* EXP( (RS_f0+RS_b*LOG(RS_sr0*LocSV/RS_sl0) ) /RS_a)
           tmp2         = tmp*SRtest
           NR           = -(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) * &
                        (ABS(P)*RS_a*LOG(tmp2+SQRT(tmp2**2+1.0))-ShTest)-SRtest !not sure if ShTest should be + or -...
           dNR          = -(1.0/w_speed(2)/rho+1.0/w_speed_neig(2)/rho_neig) * &
                        (ABS(P)*RS_a/SQRT(1+tmp2**2)*tmp) -1.0
           SRtest = SRtest-NR/dNR ! no ABS needed around NR/dNR at least for aging law
         ENDDO
         tmp=0.5*(LocSR+ABS(SRtest))  ! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
         LocSR=ABS(SRtest)

       ENDDO !  j=1,nSVupdates   !This loop corrects SV values
       !
       IF(EQN%FL.EQ.102) THEN         ! aging law
         LocSV=SV0*EXP(-tmp*time_inc/RS_sl0)+RS_sl0/tmp*(1.0-EXP(-tmp*time_inc/RS_sl0))
       ELSEIF(EQN%FL.EQ.104) THEN     ! slip law
         LocSV=RS_sl0/tmp*(tmp*SV0/RS_sl0)**(EXP(-tmp*time_inc/RS_sl0))
       ENDIF
       !
       tmp  = 0.5 * (LocSR)/RS_sr0 * EXP((RS_f0 + RS_b*LOG(RS_sr0*LocSV/RS_sl0)) / RS_a)
       !
       LocMu    = RS_a * LOG(tmp + SQRT(tmp**2 + 1.0))
       ! 2D:
       !LocTrac  = -(ABS(S_0)-LocMu*(LocP+P_0))*(S_0/ABS(S_0))
       !LocTrac  = ABS(LocTrac)*(-SignSR)  !!! line commented as it leads NOT to correct results
       ! update stress change
       LocTracXY = -((Stress(4,iBndGP) + XYStressGP(iBndGP,iTimeGP))/ShTest)*LocMu*P
       LocTracXZ = -((Stress(6,iBndGP) + XZStressGP(iBndGP,iTimeGP))/ShTest)*LocMu*P
       LocTracXY = LocTracXY - Stress(4,iBndGP)
       LocTracXZ = LocTracXZ - Stress(6,iBndGP)
       !
       ! Compute slip
       LocSlip   = LocSlip  + (LocSR)*time_inc ! ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening
       !
       !Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
       LocSR1     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXY-XYStressGP(iBndGP,iTimeGP))
       LocSR2     = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(LocTracXZ-XZStressGP(iBndGP,iTimeGP))

       LocSlip1   = LocSlip1  + (LocSR1)*time_inc 
       LocSlip2   = LocSlip2  + (LocSR2)*time_inc 
       !LocSR1     = SignSR1*ABS(LocSR1)
       !LocSR2     = SignSR2*ABS(LocSR2)
       !
       !Save traction for flux computation
       TractionGP_XY(iBndGP,iTimeGP) = LocTracXY
       TractionGP_XZ(iBndGP,iTimeGP) = LocTracXZ
       !
     ENDDO ! iTimeGP=1,DISC%Galerkin%nTimeGP
     !
     ! output rupture front
     ! outside of iTimeGP loop in order to safe an 'if' in a loop
     ! this way, no subtimestep resolution possible
     IF (DISC%DynRup%RF(iBndGP,iFace) .AND. LocSR .GT. 0.001D0) THEN
        DISC%DynRup%rupture_time(iBndGP,iFace)=time
        DISC%DynRup%RF(iBndGP,iFace) = .FALSE.
     ENDIF
     IF (LocSR.GT.DISC%DynRup%PeakSR(iBndGP,iFace)) THEN
        DISC%DynRup%PeakSR(iBndGP,iFace) = LocSR
     ENDIF
     !
     DISC%DynRup%Mu(iBndGP,iFace)        = LocMu
     DISC%DynRup%SlipRate1(iBndGP,iFace) = LocSR1
     DISC%DynRup%SlipRate2(iBndGP,iFace) = LocSR2
     DISC%DynRup%Slip(iBndGP,iFace)      = LocSlip
     DISC%DynRup%Slip1(iBndGP,iFace)     = LocSlip1
     DISC%DynRup%Slip2(iBndGP,iFace)     = LocSlip2
     DISC%DynRup%TracXY(iBndGP,iFace)    = LocTracXY
     DISC%DynRup%TracXZ(iBndGP,iFace)    = LocTracXZ
     DISC%DynRup%StateVar(iBndGP,iFace)  = LocSV
     !
    ENDDO ! iBndGP=1,DISC%Galerkin%nBndGP

  END SUBROUTINE rate_and_state_nuc101


  !> special friction case for SCEC TPV103: rate and state friction
  !> slip law with strong weakening
  !< with time and space dependent nucleation
  SUBROUTINE rate_and_state_nuc103(TractionGP_XY,TractionGP_XZ,        & ! OUT: traction
                            NorStressGP,XYStressGP,XZStressGP,         & ! IN: Godunov status
                            iFace,iSide,iElem,nBndGP,nTimeGP,          & ! IN: element ID and GP lengths
                            rho,rho_neig,w_speed,w_speed_neig,         & ! IN: background values
                            time,DeltaT,                               & ! IN: time
                            DISC,EQN,MESH,MPI,IO,BND)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization)          :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI
    TYPE(tInputOutput)             :: IO
    TYPE (tBoundary)               :: BND
    ! Local variable declaration
    INTEGER     :: i,j
    INTEGER     :: iBndGP,iTimeGP,nBndGP,nTimeGP
    INTEGER     :: iFace,iSide,iElem
    INTEGER     :: nSRupdates, nSVupdates, SignSR
    INTEGER     :: iNeighbor, iLocalNeighborSide
    INTEGER     :: MPIIndex, iObject
    REAL        :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
    REAL        :: time
    REAL        :: LocTracXY(nBndGP),LocTracXZ(nBndGP)
    REAL        :: NorStressGP(nBndGP,nTimeGP)
    REAL        :: XYStressGP(nBndGP,nTimeGP)
    REAL        :: XZStressGP(nBndGP,nTimeGP)
    REAL        :: TractionGP_XY(nBndGP,nTimeGP)
    REAL        :: TractionGP_XZ(nBndGP,nTimeGP)
    REAL        :: LocMu(nBndGP), LocD_C(nBndGP), LocSlip(nBndGP), LocSlip1(nBndGP), LocSlip2(nBndGP), LocP(nBndGP), P(nBndGP), LocSR(nBndGP), ShTest(nBndGP)
    REAL        :: LocMu_S, LocMu_D
    REAL        :: LocSR1(nBndGP),LocSR2(nBndGP)
    REAL        :: P_0(nBndGP),Strength(nBndGP),cohesion(nBndGP)
    REAL        :: rho,rho_neig,w_speed(:),w_speed_neig(:)
    REAL        :: time_inc
    REAL        :: Deltat(1:nTimeGP)
    REAL        :: SV0(nBndGP), tmp(nBndGP), tmp2(nBndGP), tmp3(nBndGP), SRtest(nBndGP), NR(nBndGP), dNR(nBndGP)
    REAL        :: LocSV(nBndGP)
    REAL        :: tmpSlip(nBndGP)
    REAL        :: RS_f0,RS_a(nBndGP),RS_b,RS_sl0(nBndGP),RS_sr0
    REAL        :: RS_fw,RS_srW(nBndGP),flv(nBndGP),fss(nBndGP),SVss(nBndGP)
    REAL        :: chi, tau, xi, eta, zeta, XGp, YGp, ZGp
    REAL        :: hypox, hypoy, hypoz
    REAL        :: Rnuc, Tnuc, radius, Gnuc, invZ, AlmostZero, aTolF
    REAL        :: prevtime,dt
    LOGICAL     :: has_converged
    LOGICAL     :: nodewise=.FALSE.
    REAL        :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
    INTEGER     :: VertexSide(4,3)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: NorStressGP,XYStressGP,XZStressGP,iFace,iSide,iElem
    INTENT(IN)    :: rho,rho_neig,w_speed,w_speed_neig,time,nBndGP,nTimeGP,DeltaT
    INTENT(IN)    :: MESH,MPI,IO
    INTENT(INOUT) :: EQN,DISC,TractionGP_XY,TractionGP_XZ
    !-------------------------------------------------------------------------!
    ! switch for Gauss node wise stress assignment
    nodewise = .TRUE.
    tmpSlip = 0.0D0

    !Apply time dependent nucleation at global time step not sub time steps for simplicity
    !initialize time and space dependent nucleation
    Tnuc = DISC%DynRup%t_0

    !TU 7.07.16: if the SR is too close to zero, we will have problems (NaN)
    !as a consequence, the SR is affected the AlmostZero value when too small
    AlmostZero = 1d-45
    !
    !PARAMETERS of THE optimisation loops
    !absolute tolerance on the function to be optimzed
    ! This value is quite arbitrary (a bit bigger as the expected numerical error) and may not be the most adapted
    !aTolF = 5e-15
    aTolF = 1e-8
    ! Number of iteration in the loops
    nSRupdates = 60
    nSVupdates = 2

    !dt = DISC%Galerkin%TimeGaussP(nTimeGP) + DeltaT(1)
    dt = sum(DeltaT(:))
    IF (time.LE.Tnuc) THEN
    IF (time.GT.0.0D0) THEN
        Gnuc=EXP((time-Tnuc)**2/(time*(time-2.0D0*Tnuc)))
        prevtime = time - dt
        IF (prevtime.GT.0.0D0) THEN
        Gnuc= Gnuc - EXP((prevtime-Tnuc)**2/(prevtime*(prevtime-2.0D0*Tnuc)))
        ENDIF
    ELSE
        Gnuc=0.0D0
    ENDIF

    !DISC%DynRup%NucBulk_** is already in fault coordinate system
    EQN%InitialStressInFaultCS(:,1,iFace)=EQN%InitialStressInFaultCS(:,1,iFace)+EQN%NucleationStressInFaultCS(:,1,iFace)*Gnuc
    EQN%InitialStressInFaultCS(:,2,iFace)=EQN%InitialStressInFaultCS(:,2,iFace)+EQN%NucleationStressInFaultCS(:,2,iFace)*Gnuc
    EQN%InitialStressInFaultCS(:,3,iFace)=EQN%InitialStressInFaultCS(:,3,iFace)+EQN%NucleationStressInFaultCS(:,3,iFace)*Gnuc
    EQN%InitialStressInFaultCS(:,4,iFace)=EQN%InitialStressInFaultCS(:,4,iFace)+EQN%NucleationStressInFaultCS(:,4,iFace)*Gnuc
    EQN%InitialStressInFaultCS(:,5,iFace)=EQN%InitialStressInFaultCS(:,5,iFace)+EQN%NucleationStressInFaultCS(:,5,iFace)*Gnuc
    EQN%InitialStressInFaultCS(:,6,iFace)=EQN%InitialStressInFaultCS(:,6,iFace)+EQN%NucleationStressInFaultCS(:,6,iFace)*Gnuc

    ENDIF ! Tnuc
    !
     !
     LocSlip   = DISC%DynRup%Slip(:,iFace)
     LocSlip1   = DISC%DynRup%Slip1(:,iFace)
     LocSlip2   = DISC%DynRup%Slip2(:,iFace)
     LocSR1    = DISC%DynRup%SlipRate1(:,iFace)
     LocSR2    = DISC%DynRup%SlipRate2(:,iFace)
     LocSV     = DISC%DynRup%StateVar(:,iFace)
     P_0       = EQN%InitialStressInFaultCS(:,1,iFace)
     !
     DO iTimeGP=1,nTimeGP
         !
         ! friction develops as                    mu = a * arcsinh[ V/(2*V0) * exp(SV/a) ]
         ! state variable SV develops as     dSV / dt = -(V - L) * (SV - SV_ss)
         !                                      SV_ss = a * ln[ 2*V0/V * sinh(mu_ss/a) ]
         !                                      mu_ss = mu_w + [mu_lv - mu_w] / [ 1 + (V/Vw)^8 ] ^ (1/8) ]
         !                                      mu_lv = mu_0 - (b-a) ln (V/V0)
         !
         LocP   = NorStressGP(:,iTimeGP)
         time_inc = DeltaT(iTimeGP)
         !
         RS_f0  = DISC%DynRup%RS_f0     ! mu_0, reference friction coefficient
         RS_sr0 = DISC%DynRup%RS_sr0    ! V0, reference velocity scale
         RS_fw  = DISC%DynRup%Mu_w      ! mu_w, weakening friction coefficient
         RS_srW = DISC%DynRup%RS_srW_array(:,iFace)    ! Vw, weakening sliding velocity, space dependent
         RS_a   = DISC%DynRup%RS_a_array(:,iFace) ! a, direct effect, space dependent
         RS_b   = DISC%DynRup%RS_b       ! b, evolution effect
         RS_sl0 = DISC%DynRup%RS_sl0_array(:,iFace)     ! L, char. length scale
         !
         ! load traction and normal stress
         P      = LocP+P_0
         ShTest = SQRT((EQN%InitialStressInFaultCS(:,4,iFace) + XYStressGP(:,iTimeGP))**2 + (EQN%InitialStressInFaultCS(:,6,iFace) + XZStressGP(:,iTimeGP))**2)
         !
         SV0=LocSV    ! Careful, the SV must always be corrected using SV0 and not LocSV!
         !
         ! The following process is adapted from that described by Kaneko et al. (2008)
         !
         LocSR      = SQRT(LocSR1**2 + LocSR2**2)
         LocSR = max(AlmostZero,LocSR)
         !
         tmp = LocSR
         invZ = (1.0d0/w_speed(2)/rho+1.0d0/w_speed_neig(2)/rho_neig)

         DO j=1,nSVupdates   !This loop corrects SV values
             !
             !1. update SV using Vold from the previous time step
             !   exact integration assuming constant V in this iteration
             !   low-velocity steady state friction coefficient
             flv = RS_f0 - (RS_b-RS_a)* LOG(tmp/RS_sr0)
             !   steady state friction coefficient
             fss = RS_fw + (flv - RS_fw)/(1.0D0+(tmp/RS_srW)**8d0)**(1.0D0/8.0D0)
             ! steady-state state variable with SINH(X)=(EXP(X)-EXP(-X))/2
             SVss = RS_a * LOG(2.0D0*RS_sr0/tmp * ( EXP(fss/RS_a)-EXP(-fss/RS_a))/2.0D0)
             !
             LocSV=SVss*(1.0d0-EXP(-tmp*time_inc/RS_sl0))+EXP(-tmp*time_inc/RS_sl0)*SV0

             !2. solve for Vnew , applying the Newton-Raphson algorithm as in Case 3 and 4
             !   but with different mu evolution
             ! SR fulfills g(SR)=f(SR), NR=f-g and dNR = d(NR)/d(SR)
             ! SR_{i+1}=SR_i-( NR_i / dNR_i )
             ! equalize:
             !         g = SR*MU/2/cs + T^G             (eq. 18 of de la Puente et al. (2009))
             !         f = (mu*P_0-|S_0|)*S_0/|S_0|     (Coulomb's model of friction)
             !  where mu = a * arcsinh[ V/(2*V0) * exp(SV/a) ]
             SRtest=LocSR  ! We use as first guess the SR value of the previous time step
             !
             tmp          = 0.5D0/RS_sr0* EXP(LocSV/RS_a)
             has_converged = .FALSE.

             DO i=1,nSRupdates  !This loop corrects SR values
                 ! for convenience
                 tmp2         = tmp*SRtest != X in ASINH(X) for mu calculation
                 NR           = -invZ * (ABS(P)*RS_a*LOG(tmp2+SQRT(tmp2**2+1.0))-ShTest)-SRtest
                 IF (maxval(abs(NR))<atolF) THEN
                    has_converged = .TRUE.
                    EXIT
                 ENDIF
                 dNR          = -invZ * (ABS(P)*RS_a/SQRT(1d0+tmp2**2)*tmp) -1.0
                 tmp3 = NR/dNR
                 SRtest = max(AlmostZero,SRtest - tmp3)
             ENDDO
             !
             ! 3. update theta, now using V=(Vnew+Vold)/2
             tmp=0.5d0*(LocSR+ABS(SRtest))  ! For the next SV update, use the mean slip rate between the initial guess and the one found (Kaneko 2008, step 6)
             !
             ! 4. solve again for Vnew
             LocSR=ABS(SRtest)
             !
         ENDDO !  j=1,nSVupdates   !This loop corrects SV values
         if (.NOT.has_converged) THEN
            !logError(*) 'nonConvergence RS Newton', time
            if (tmp(1).NE.tmp(1)) then
               logError(*) 'NaN detected', time
               STOP
            endif
         ENDIF
         !
         ! 5. get final theta, mu, traction and slip
         ! SV from mean slip rate in tmp
         flv = RS_f0 -(RS_b-RS_a)* LOG(tmp/RS_sr0)
         fss = RS_fw + (flv - RS_fw)/(1.0D0+(tmp/RS_srW)**8)**(1.0D0/8.0D0)
         SVss = RS_a * LOG(2.0D0*RS_sr0/tmp * ( EXP(fss/RS_a)-EXP(-fss/RS_a))/2.0D0)
         LocSV=Svss*(1.0d0-EXP(-tmp*time_inc/RS_sl0))+EXP(-tmp*time_inc/RS_sl0)*SV0
         !Mu from LocSR
         tmp = 0.5D0*(LocSR)/RS_sr0 * EXP(LocSV/RS_a)
         LocMu    = RS_a * LOG(tmp+SQRT(tmp**2+1.0D0))
         ! update stress change

         LocTracXY = -((EQN%InitialStressInFaultCS(:,4,iFace) + XYStressGP(:,iTimeGP))/ShTest)*LocMu*P
         LocTracXZ = -((EQN%InitialStressInFaultCS(:,6,iFace) + XZStressGP(:,iTimeGP))/ShTest)*LocMu*P
         LocTracXY = LocTracXY - EQN%InitialStressInFaultCS(:,4,iFace)
         LocTracXZ = LocTracXZ - EQN%InitialStressInFaultCS(:,6,iFace)
         !
         ! Compute slip
         LocSlip   = LocSlip  + (LocSR)*time_inc ! ABS of LocSR removed as it would be the accumulated slip that is usually not needed in the solver, see linear slip weakening
         !
         !Update slip rate (notice that LocSR(T=0)=-2c_s/mu*s_xy^{Godunov} is the slip rate caused by a free surface!)
         LocSR1     = -invZ*(LocTracXY-XYStressGP(:,iTimeGP))
         LocSR2     = -invZ*(LocTracXZ-XZStressGP(:,iTimeGP))

         !TU 07.07.16: correct LocSR1_2 to avoid numerical errors
         tmp = sqrt(LocSR1**2+LocSR2**2)
         where ( tmp.NE.0d0) 
            LocSR1 = LocSR*LocSR1/tmp
            LocSR2 = LocSR*LocSR2/tmp
         endwhere
         tmpSlip = tmpSlip(:) + tmp(:)*time_inc

         LocSlip1   = LocSlip1  + (LocSR1)*time_inc 
         LocSlip2   = LocSlip2  + (LocSR2)*time_inc 
         !LocSR1     = SignSR1*ABS(LocSR1)
         !LocSR2     = SignSR2*ABS(LocSR2)
         !
         !Save traction for flux computation
         TractionGP_XY(:,iTimeGP) = LocTracXY
         TractionGP_XZ(:,iTimeGP) = LocTracXZ
         !
     ENDDO ! iTimeGP=1,DISC%Galerkin%nTimeGP
     !
     ! output rupture front
     ! outside of iTimeGP loop in order to safe an 'if' in a loop
     ! this way, no subtimestep resolution possible
     where (DISC%DynRup%RF(:,iFace) .AND. LocSR .GT. 0.001D0)
        DISC%DynRup%rupture_time(:,iFace)=time
        DISC%DynRup%RF(:,iFace) = .FALSE.
     endwhere
     where (LocSR.GT.DISC%DynRup%PeakSR(:,iFace))
        DISC%DynRup%PeakSR(:,iFace) = LocSR
     endwhere
    !output time when shear stress is equal to the dynamic stress after rupture arrived
    !currently only for linear slip weakening
    where ( (DISC%DynRup%rupture_time(:,iFace) .GT. 0.0) .AND. &
            (DISC%DynRup%rupture_time(:,iFace) .LE. time) .AND. &
             DISC%DynRup%DS(:,iFace) .AND. &
             DISC%DynRup%Mu(:,iFace) .LE. (RS_fw+0.05*(RS_f0-RS_fw)))
      DISC%DynRup%dynStress_time(:,iFace)=time
      DISC%DynRup%DS(:,iFace) = .FALSE.
    end where
     !
     DISC%DynRup%Mu(:,iFace)        = LocMu
     DISC%DynRup%SlipRate1(:,iFace) = LocSR1
     DISC%DynRup%SlipRate2(:,iFace) = LocSR2
     DISC%DynRup%Slip(:,iFace)      = LocSlip
     DISC%DynRup%Slip1(:,iFace)     = LocSlip1
     DISC%DynRup%Slip2(:,iFace)     = LocSlip2
     DISC%DynRup%TracXY(:,iFace)    = LocTracXY
     DISC%DynRup%TracXZ(:,iFace)    = LocTracXZ
     DISC%DynRup%StateVar(:,iFace)  = LocSV

     IF (DISC%DynRup%magnitude_out(iFace)) THEN
        DISC%DynRup%averaged_Slip(iFace) = DISC%DynRup%averaged_Slip(iFace) + sum(tmpSlip)/nBndGP
     ENDIF
  !
 END SUBROUTINE rate_and_state_nuc103

 END MODULE
