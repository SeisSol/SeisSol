!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Atanas Atanasov (atanasoa AT in.tum.de, http://www5.in.tum.de/wiki/index.php/Atanas_Atanasov)
!! @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2013-2017, SeisSol Group
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
!! Routines handling fault output

#include <Initializer/preProcessorMacros.fpp>

  MODULE faultoutput_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  !
  INTERFACE faultoutput
     MODULE PROCEDURE faultoutput
  END INTERFACE
  !
  !---------------------------------------------------------------------------!
  PUBLIC   :: faultoutput
  public   :: calc_FaultOutput
  PRIVATE  :: write_FaultOutput_elementwise
  PRIVATE  :: write_FaultOutput_atPickpoint
  !---------------------------------------------------------------------------!

CONTAINS

!> Fault output interface to manage output calculation and type of format
!<
  SUBROUTINE faultoutput(EQN, DISC, MESH, IO, MPI, MaterialVal, BND, time, dt)

      !-------------------------------------------------------------------------!
      USE JacobiNormal_mod
      USE magnitude_output_mod
      !-------------------------------------------------------------------------!
      IMPLICIT NONE
      !-------------------------------------------------------------------------!
      ! Argument list declaration
      TYPE(tEquations)         :: EQN
      TYPE(tDiscretization)    :: DISC
      TYPE(tUnstructMesh)      :: MESH
      TYPE(tInputOutput)       :: IO                                            ! IO structure
      TYPE(tMPI)               :: MPI                                           ! MPI
      TYPE(tBoundary)          :: BND                                           ! BND    data structure
      REAL                     :: MaterialVal(MESH%nElem,EQN%nBackgroundVar)    ! Local Mean Values
      REAL                     :: dt, time                                      ! Timestep and time
	  ! local variable declaration
      LOGICAL                  :: isOnPickpoint
      LOGICAL                  :: isOnElementwise
      LOGICAL                  :: isOnFaultoutput
      INTEGER                  :: i, rank_int                                   ! Loop variables
      !-------------------------------------------------------------------------!
      isOnPickpoint = .FALSE.
      isOnElementwise = .FALSE.

      !-------------------------------------------------------------------------!
      !
      ! Calculate Fault Output of the previous timestep
      ! Here, because the complete updated dgvar value of the (MPI-)Neighbor is needed
      ! Note that this causes a dt timeshift in the DR output routines
      !
      !
      IF (DISC%DynRup%energy_rate_output_on.EQ.1) THEN
         IF ( MOD(DISC%iterationstep,DISC%DynRup%energy_rate_printtimeinterval).EQ.0 &
         .OR. (DISC%EndTime-time).LE.(dt*1.005d0) ) THEN
            CALL energy_rate_output(MaterialVal,time,DISC,MESH,MPI,IO)
         ENDIF
      ENDIF


      SELECT CASE(DISC%DynRup%OutputPointType)
       ! For historical reasons fault output DISC%DynRup%OutputPointType= 3 or 4 or 5
       ! Case 0 means no fault output
       !CASE(0)
       !  logInfo(*) 'Fault output case (0).'
       !  RETURN

       ! Fault Receiver Output
       CASE(3)
         ! Check if output is desired
         IF (.NOT. DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output ) THEN
             RETURN
         ENDIF
         !
         ! Check if output for this time step is desired
         ! (note that fault computation is always one behind due to MPI communication = -1!)
         IF ( MOD(DISC%iterationstep-1,DISC%DynRup%DynRup_out_atPickpoint%printtimeinterval).EQ.0 &
         .OR. (DISC%EndTime-time).LE.(dt*1.005d0) ) THEN
            IF (DISC%iterationstep.EQ.0) RETURN ! not the iteration 0
            CONTINUE
         ! print always first timestep
         ELSEIF(DISC%iterationstep .EQ. 1) THEN
            CONTINUE
         ELSE
            RETURN
         ENDIF
         CALL calc_FaultOutput(DISC%DynRup%DynRup_out_atPickpoint, DISC, EQN, MESH, MaterialVal, BND, time)
         CALL write_FaultOutput_atPickpoint(EQN, DISC, MESH, IO, MPI, MaterialVal, BND, time, dt)

       !Fault Element Output
       CASE(4)
         ! fault output iterations, including iteration 1 and last timestep
         IF (DISC%DynRup%DynRup_out_elementwise%printIntervalCriterion.EQ.1) THEN
            IF ( MOD(DISC%iterationstep-1,DISC%DynRup%DynRup_out_elementwise%printtimeinterval).EQ.0.OR. &
            ((min(DISC%EndTime,dt*DISC%MaxIteration)-time).LE.(dt*1.005d0))) THEN
               CONTINUE
!            ELSEIF(DISC%iterationstep .EQ. 0 ) THEN  ! uncomment to print setup at time 0 for intial stresses
!               CONTINUE
            ELSE  ! all timesteps without fault output
               RETURN
            ENDIF
         ELSE
            IF ((DISC%iterationstep.EQ.1).OR. ((DISC%iterationstep.GE.2) .AND. &
            ((MOD(time,DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec).LE.0.5*dt).OR. &
            ((ceiling(time/DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec)*DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec-time).LT.0.5*dt).OR. &
            ((min(DISC%EndTime,dt*DISC%MaxIteration)-time).LE.(dt*1.005d0))))) THEN
               CONTINUE
!            ELSEIF(DISC%iterationstep .EQ. 0 ) THEN  ! uncomment to print setup at time 0 for intial stresses
!               CONTINUE
            ELSE  ! all timesteps without fault output
               RETURN
            ENDIF
         ENDIF
         ! remember that fault output was written here
         isOnElementwise=.TRUE.
       ! combines option 3 and 4: output at individual stations and at complete fault
       CASE(5)
         ! check time for fault receiver output
         IF (.NOT. DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output ) THEN
             CONTINUE
         ELSEIF ( MOD(DISC%iterationstep-1,DISC%DynRup%DynRup_out_atPickpoint%printtimeinterval).EQ.0 &
         .OR. (min(DISC%EndTime,dt*DISC%MaxIteration)-time).LE.(dt*1.005d0) ) THEN
            isOnPickpoint = .TRUE.
         ENDIF
         ! check time for fault plane output
         IF (DISC%DynRup%DynRup_out_elementwise%printIntervalCriterion.EQ.1) THEN
            IF ( MOD(DISC%iterationstep-1,DISC%DynRup%DynRup_out_elementwise%printtimeinterval).EQ.0.OR. &
            ((min(DISC%EndTime,dt*DISC%MaxIteration)-time).LE.(dt*1.005d0))) THEN
               isOnElementwise=.TRUE.
            ENDIF
         ELSE
            IF ((DISC%iterationstep.EQ.1).OR. ((DISC%iterationstep.GE.2) .AND. &
            ((MOD(time,DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec).LE.0.5*dt).OR. &
            ((ceiling(time/DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec)*DISC%DynRup%DynRup_out_elementwise%printtimeinterval_sec-time).LT.0.5*dt).OR. &
            ((min(DISC%EndTime,dt*DISC%MaxIteration)-time).LE.(dt*1.005d0))))) THEN
               isOnElementwise=.TRUE.
            ENDIF
         ENDIF

         ! print always first timestep
         IF(DISC%iterationstep .EQ. 1) THEN
           isOnPickpoint = .TRUE.
           isOnElementwise=.TRUE.
         ENDIF
         !
         IF (isOnPickpoint) THEN
           CALL calc_FaultOutput(DISC%DynRup%DynRup_out_atPickpoint, DISC, EQN, MESH, MaterialVal, BND, time)
           CALL write_FaultOutput_atPickpoint(EQN, DISC, MESH, IO, MPI, MaterialVal, BND, time, dt)
         ENDIF
         !
       CASE DEFAULT
          ! no output
          CONTINUE
      END SELECT
      !
  END SUBROUTINE faultoutput
!
!> Fault output calculation at specific positions (receiver and elementwise)
!<
  SUBROUTINE calc_FaultOutput( DynRup_output, DISC, EQN, MESH, MaterialVal, BND, time )
    use  f_ftoc_bind_interoperability
    use iso_c_binding, only: c_loc

    !-------------------------------------------------------------------------!
    USE common_operators_mod
    USE JacobiNormal_mod
    USE DGBasis_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)              :: EQN                                      !< EQN global variable
    TYPE(tDiscretization), TARGET :: DISC                                     !< DISC global variable
    TYPE(tUnstructMesh)           :: MESH                                     !< MESH global variable
    TYPE(tInputOutput)            :: IO                                       !< IO structure
    TYPE(tMPI)                    :: MPI                                      !< MPI
    TYPE(tBoundary)               :: BND                                      !< BND    data structure
    REAL                          :: MaterialVal(MESH%nElem,EQN%nBackgroundVar) !< Local Mean Values
    REAL                          :: time                                     !< time
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    TYPE(tDynRup_output), target  :: DynRup_output                            !< Output data for Dynamic Rupture processes
    INTEGER :: iElem                                                          ! Element number                        !
    INTEGER :: iObject, MPIIndex, MPIIndex_DR
    INTEGER :: iDegFr                                                         ! Degree of freedom                     !
    INTEGER :: iVar                                                           ! Variable number                       !
    INTEGER :: iSide
    INTEGER :: iFace
    INTEGER :: iNeighbor, NeigBndGP
    INTEGER :: iLocalNeighborSide
    INTEGER :: NeighborPoly, LocPoly, MaxDegFr, iBndGP
    INTEGER :: NeighborDegFr, LocDegFr
    INTEGER :: OutVars, nOutPoints
    INTEGER :: i,j,m,k,iOutPoints                                             ! Loop variables                   !
    INTEGER :: SubElem, number_of_subtriangles                                ! elementwise fault refinement parameters
    REAL    :: NormalVect_n(3)                                                ! Normal vector components         !
    REAL    :: NormalVect_s(3)                                                ! Normal vector components         !
    REAL    :: NormalVect_t(3)                                                ! Normal vector components         !
    REAL    :: T(EQN%nVar,EQN%nVar)                                           ! Rotation matrix
    REAL    :: iT(EQN%nVar,EQN%nVar)                                          ! Rotation matrix
    REAL    :: V1(DISC%Galerkin%nDegFr,EQN%nVar)                              ! Reference state     (ref)
    REAL    :: V2(DISC%Galerkin%nDegFr,EQN%nVar)                              ! Reference state     (ref)
    REAL    :: Stress(6)                                                      ! The background stress tensor in vector form at a single BndGP
    REAL    :: SideVal(EQN%nVar), SideVal2(EQN%nVar)
    REAL    :: phi(2),rho, rho_neig, mu, mu_neig, lambda, lambda_neig
    REAL    :: LocXYStress, LocXZStress, TracXY, TracXZ, LocSRs, LocSRd
    REAL    :: w_speed(EQN%nNonZeroEV), TracEla, Trac, Strength, LocU, LocP, w_speed_neig(EQN%nNonZeroEV)
    REAL    :: S_0,P_0,S_XY,S_XZ
    REAL    :: MuVal, cohesion, LocSV
    REAL    :: LocYY, LocZZ, LocYZ                                            ! temporary stress values
    REAL    :: tmp_mat(1:6), LocMat(1:6), TracMat(1:6)                        ! temporary stress tensors
    REAL    :: rotmat(1:6,1:6)                                                ! Rotation matrix
    REAL    :: TmpMat(EQN%nBackgroundVar)                                     ! temporary material values
    REAL    :: NorDivisor,ShearDivisor,UVelDivisor
    REAL    :: strike_vector(1:3), dip_vector(1:3), crossprod(1:3) !for rotation of Slip from local to strike, dip coordinate
    REAL    :: cos1, sin1, scalarprod
    REAL, PARAMETER    :: ZERO = 0.0D0
    ! Parameters used for calculating Vr
    LOGICAL :: compute_Vr
    INTEGER :: nDegFr2d, jBndGP, i1, j1
    REAL    :: chi, tau, phiT, phi2T(2),Slowness, dt_dchi, dt_dtau, Vr
    REAL    :: dt_dx1, dt_dy1
    REAL    :: xV(4), yV(4), zV(4)
    REAL    :: xab(3), xac(3), grad2d(2,2), JacobiT2d(2,2)
    REAL, ALLOCATABLE  :: projected_RT(:)
    real, dimension( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES ) :: DOFiElem_ptr ! no: it's not a pointer..
    real, dimension( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES ) :: DOFiNeigh_ptr ! no pointer again
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: BND, DISC, EQN, MESH, MaterialVal, time
    INTENT(INOUT) :: DynRup_output
    !-------------------------------------------------------------------------!
    !
    ! register epik/scorep function
    EPIK_FUNC_REG("calc_faultoutput")
    SCOREP_USER_FUNC_DEFINE()
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN("calc_faultoutput")
    !
    nOutPoints = DynRup_output%nDR_pick                                        ! number of output receivers for this MPI domain
    ! in case of refined fault output elements, the rotation matrices are only stored for every mother tet
    ! in pickpoint case SubElem = 1, in refined elementwise case SubElem = number of subtets per element
    IF (DISC%DynRup%OutputPointType.EQ.4) THEN
        IF(DISC%DynRup%DynRup_out_elementwise%refinement_strategy.eq.1) THEN
            number_of_subtriangles=3
        ELSE IF (DISC%DynRup%DynRup_out_elementwise%refinement_strategy.eq.2) THEN
            number_of_subtriangles=4
        ELSE
            number_of_subtriangles=1
        ENDIF
        SubElem = number_of_subtriangles**DISC%DynRup%DynRup_out_elementwise%refinement
    ELSE
        SubElem = 1
    ENDIF
    DO iOutPoints = 1,nOutPoints                                               ! loop over number of output receivers for this domain
          !
          iFace               = DynRup_output%RecPoint(iOutPoints)%index       ! current receiver location
          iElem               = MESH%Fault%Face(iFace,1,1)
          iSide               = MESH%Fault%Face(iFace,2,1)
          !
          iNeighbor           = MESH%Fault%Face(iFace,1,2)
          iLocalNeighborSide  = MESH%Fault%Face(iFace,2,2)
          !
          w_speed(:)      = DISC%Galerkin%WaveSpeed(iElem,iSide,:)
          rho             = MaterialVal(iElem,1)
          !
          if( iElem == 0 ) then
            call c_interoperability_getNeighborDofsFromDerivatives( i_meshId = iNeighbor, \
                                                                    i_faceId = iLocalNeighborSide, \
                                                                    o_dofs   = dofiElem_ptr )
          else
            call c_interoperability_getDofsFromDerivatives( i_meshId = iElem, \
                                                            o_dofs   = DOFiElem_ptr)
          endif

          IF (iNeighbor == 0) THEN
            ! iNeighbor is in the neighbor domain
            ! The neighbor element belongs to a different MPI domain
            iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
            MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
            call c_interoperability_getNeighborDofsFromDerivatives( i_meshId = iElem, \
                                                                    i_faceId = iSide, \
                                                                    o_dofs   = DOFiNeigh_ptr )

            ! Bimaterial case only possible for elastic isotropic materials
            TmpMat(:)   = BND%ObjMPI(iObject)%NeighborBackground(:,MPIIndex)
            rho_neig    = TmpMat(1)
            mu_neig     = TmpMat(2)
            lambda_neig = TmpMat(3)
            w_speed_neig(1) = SQRT((lambda_neig+2.0D0*mu_neig)/rho_neig)  ! Will only work in elastic isotropic cases
            w_speed_neig(2) = SQRT(mu_neig/rho_neig)
            w_speed_neig(3) = w_speed_neig(2)
          ELSE
            ! normal case: iNeighbor present in local domain
            call c_interoperability_getDofsFromDerivatives( i_meshId = iNeighbor, \
                                                            o_dofs   = DOFiNeigh_ptr )
            w_speed_neig(:) = DISC%Galerkin%WaveSpeed(iNeighbor,iLocalNeighborSide,:)
            rho_neig        = MaterialVal(iNeighbor,1)
          ENDIF
          !
          V1(:,:)=0.
          V2(:,:)=0.
          !
          ! currently no p-adaptivity
          ! IF(DISC%Galerkin%pAdaptivity.GE.1) THEN
          !   LocPoly  = INT(DISC%Galerkin%LocPoly(iElem))
          ! ELSE
          !   LocPoly  = DISC%Galerkin%nPoly
          ! ENDIF
          LocPoly  = DISC%Galerkin%nPoly
          LocDegFr = (LocPoly+1)*(LocPoly+2)*(LocPoly+3)/6.0D0
          !
          ! Local side's normal and tangential vectors
          NormalVect_n = MESH%Fault%geoNormals(1:3,iFace)
          NormalVect_s = MESH%Fault%geoTangent1(1:3,iFace)
          NormalVect_t = MESH%Fault%geoTangent2(1:3,iFace)
          !
          ! Rotate DoF
          CALL RotationMatrix3D(NormalVect_n,NormalVect_s,NormalVect_t,T(:,:),iT(:,:),EQN)
          DO iDegFr=1,LocDegFr
            V1(iDegFr,:)=MATMUL(iT(:,:),dofiElem_ptr(iDegFr,1:EQN%nVar))
            V2(iDegFr,:)=MATMUL(iT(:,:),DOFiNeigh_ptr(iDegFr,1:EQN%nVar))
          ENDDO
          !
          ! load nearest boundary GP: iBndGP
          iBndGP = DynRup_output%OutInt(iOutPoints,1)

          !Background stress rotation to face's reference system
          Stress = EQN%InitialStressInFaultCS(iBndGP,:,iFace)
          !
          MuVal = DISC%DynRup%output_Mu(iBndGP,iFace)
          LocSV = DISC%DynRup%output_StateVar(iBndGP,iFace) ! load state variable of RS for output
          cohesion  = DISC%DynRup%cohesion(iBndGP,iFace)
          S_XY  = Stress(4)
          S_XZ  = Stress(6)
          P_0   = Stress(1)
          !
          ! Obtain values at output points
          SideVal  = 0.
          SideVal2 = 0.
          DO iDegFr = 1, LocDegFr
             ! Basis functions for the automatic evaluation of DOFs at fault output nodes
             phi(:) = DynRup_output%OutEval(iOutPoints,1,iDegFr,:)
             SideVal(:)  = SideVal(:)  + V1(iDegFr,:)*phi(1)
             SideVal2(:) = SideVal2(:) + V2(iDegFr,:)*phi(2)
          ENDDO
          !
          ! Compute divisors
          NorDivisor   = 1.0D0/(w_speed_neig(1)*rho_neig+w_speed(1)*rho)
          ShearDivisor = 1.0D0/(w_speed_neig(2)*rho_neig+w_speed(2)*rho)
          UVelDivisor  = 1.0D0/(w_speed(1)*rho)
          !
          ! Traction (elastic response, not at sliding fault!)
          LocXYStress     =  SideVal(4) + &
                            (((SideVal2(4)-SideVal(4)) + &
                            w_speed_neig(2)*rho_neig*(SideVal2(8)-SideVal(8)))* &
                            w_speed(2)*rho ) * ShearDivisor
          LocXZStress     =  SideVal(6) + &
                            (((SideVal2(6)-SideVal(6)) + &
                            w_speed_neig(2)*rho_neig*(SideVal2(9)-SideVal(9)))* &
                            w_speed(2)*rho ) * ShearDivisor
          ! Normal stress
          LocP = SideVal(1)+(((SideVal2(1)-SideVal(1))+w_speed_neig(1)*rho_neig*(SideVal2(7)-SideVal(7)))* &
                 w_speed(1)*rho) * NorDivisor
          ! Normal velocity
          LocU = SideVal(7)+ (LocP-SideVal(1)) * UVelDivisor

          ! compute missing sigma values for rotation
          LocYY = SideVal(2) + (LocP-SideVal(1)) * (1.0D0-2.0D0*w_speed(2)**2/w_speed(1)**2)
          LocZZ = SideVal(3) + (LocP-SideVal(1)) * (1.0D0-2.0D0*w_speed(2)**2/w_speed(1)**2)
          LocYZ = SideVal(5)

          ! absolute shear stress
          TracEla = SQRT((S_XY + LocXYStress)**2 + (S_XZ + LocXZStress)**2)

          ! compute fault strength
          SELECT CASE(EQN%FL)
          CASE DEFAULT
            ! linear slip weakening
            Strength = -MuVal*MIN(LocP+P_0,ZERO) - cohesion
          CASE(3,4)
             ! rate and state (once everything is tested and cohesion works for RS, this option could be merged to default)
             Strength = -MuVal*(LocP+P_0)
          CASE(6)
            ! exception for bimaterial with LSW case
            ! modify strength according to prakash clifton
            ! in the first step we just take the Strength values of the nearest GP
            Strength = DISC%DynRup%output_Strength(iBndGP,iFace)
          END SELECT
          !
          ! Coulomb's model for obtaining actual traction
          IF(TracEla.GT.ABS(Strength)) THEN
              TracXY = ((S_XY+LocXYStress)/TracEla)*Strength
              TracXZ = ((S_XZ+LocXZStress)/TracEla)*Strength
              !
              ! update stress change
              TracXY = TracXY - S_XY
              TracXZ = TracXZ - S_XZ
          ELSE
              Trac = TracEla
              TracXY = LocXYStress
              TracXZ = LocXZStress
          ENDIF
          !
          ! Rotate back to xyz-system
          ! create dummy matrix
          ! tmp_mat in xyz system:
          tmp_mat    = 0.0D0
          tmp_mat(1) = LocP
          tmp_mat(2) = LocYY
          tmp_mat(3) = LocZZ
          tmp_mat(4) = TracXY ! modified traction !
          tmp_mat(5) = LocYZ
          tmp_mat(6) = TracXZ ! modified traction !
          tmp_mat    = MATMUL(T(1:6,1:6),tmp_mat)

          ! rotate into fault system
          rotmat  = DynRup_output%rotmat((iOutPoints-1)/(SubElem)+1,1:6,1:6) ! in pickpoint case SubElem = 1, in refined elementwise case SubElem = number of subtets per element
          TracMat = MATMUL(rotmat,tmp_mat)

          ! Do it again but assume no slip occured
          tmp_mat    = 0.0D0
          tmp_mat(1) = LocP
          tmp_mat(2) = LocYY
          tmp_mat(3) = LocZZ
          tmp_mat(4) = LocXYStress ! traction of standard wave propagation without slip!
          tmp_mat(5) = LocYZ
          tmp_mat(6) = LocXZStress ! traction of standard wave propagation without slip!
          tmp_mat    = MATMUL(T(1:6,1:6),tmp_mat) ! tmp_mat in xyz system

          ! rotate into fault system
          LocMat = MATMUL(rotmat,tmp_mat)
          
          ! z must not be +/- (0,0,1) for the following to work (see also create_fault_rotationmatrix)
          strike_vector(1) = NormalVect_n(2)/sqrt(NormalVect_n(1)**2+NormalVect_n(2)**2)
          strike_vector(2) = -NormalVect_n(1)/sqrt(NormalVect_n(1)**2+NormalVect_n(2)**2)
          strike_vector(3) = 0.0D0
          dip_vector = NormalVect_n .x. strike_vector
          dip_vector = dip_vector / sqrt(dip_vector(1)**2+dip_vector(2)**2+dip_vector(3)**2)

          ! sliprate
          if (DISC%DynRup%SlipRateOutputType .eq. 1) then
            LocSRs = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(TracMat(4)-LocMat(4))
            LocSRd = -(1.0D0/(w_speed(2)*rho)+1.0D0/(w_speed_neig(2)*rho_neig))*(TracMat(6)-LocMat(6))
          else
            LocSRs = dot_product(SideVal2(8) * NormalVect_s + SideVal2(9) * NormalVect_t, strike_vector) - dot_product(SideVal(8) * NormalVect_s + SideVal(9) * NormalVect_t, strike_vector)
            LocSRd = dot_product(SideVal2(8) * NormalVect_s + SideVal2(9) * NormalVect_t, dip_vector   ) - dot_product(SideVal(8) * NormalVect_s + SideVal(9) * NormalVect_t, dip_vector   )
          end if


          ! TU 06.2015: average stress over the element in the considered direction
          ! (uncomment for rough geometry)
          !SideVal(1)  =  V1(1,1)*DynRup_output%OutEval(iOutPoints,1,1,1)
          !SideVal2(1) =  V2(1,1)*DynRup_output%OutEval(iOutPoints,1,1,2)
          !LocP = SideVal(1)+(((SideVal2(1)-SideVal(1))+w_speed_neig(1)*rho_neig*(SideVal2(7)-SideVal(7)))* &
          !       w_speed(1)*rho) * NorDivisor
          ! Store Values into Output vector OutVal

          OutVars = 0
          IF (DynRup_output%OutputMask(1).EQ.1) THEN
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = LocSRs !OutVars =1
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = LocSRd !OutVars =2
          ENDIF
          IF (DynRup_output%OutputMask(2).EQ.1) THEN
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = TracMat(4) !OutVars =3
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = TracMat(6) !OutVars =4
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = LocP !OutVars =5
          ENDIF
          IF (DynRup_output%OutputMask(3).EQ.1) THEN
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = LocU !OutVars =6
          ENDIF
          IF (DynRup_output%OutputMask(4).EQ.1) THEN
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = MuVal !OutVars =7
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = LocSV !OutVars =8
          ENDIF
          IF (DynRup_output%OutputMask(5).EQ.1) THEN
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = TracMat(4)+DISC%DynRup%DynRup_Constants(iOutPoints)%ts0 !OutVars =9
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = TracMat(6)+DISC%DynRup%DynRup_Constants(iOutPoints)%td0 !OutVars =10
              OutVars = OutVars + 1
              DynRup_output%OutVal(iOutPoints,1,OutVars) = LocP+DISC%DynRup%DynRup_Constants(iOutPoints)%p0 !OutVars =11
          ENDIF

          IF (DISC%DynRup%OutputPointType.EQ.4) THEN

              IF (DynRup_output%OutputMask(6).EQ.1) THEN
                  ! TU 07.15 rotate Slip from face reference coordinate to (strike,dip, normal) reference cordinate
                  cos1 = dot_product(strike_vector(:),NormalVect_s(:))
                  crossprod(:) = strike_vector(:) .x. NormalVect_s(:)

                  scalarprod = dot_product(crossprod(:),NormalVect_n(:))
                  !TU 2.11.15 :cos1**2 can be greater than 1 because of rounding errors -> min
                  IF (scalarprod.GT.0) THEN
                      sin1=sqrt(1-min(1d0,cos1**2))
                  ELSE
                      sin1=-sqrt(1-min(1d0,cos1**2))
                  ENDIF

                  OutVars = OutVars + 1
                  DynRup_output%OutVal(iOutPoints,1,OutVars)  = cos1 * DISC%DynRup%output_Slip1(iBndGP,iFace) - sin1* DISC%DynRup%output_Slip2(iBndGP,iFace)
                  OutVars = OutVars + 1
                  DynRup_output%OutVal(iOutPoints,1,OutVars)  = sin1 * DISC%DynRup%output_Slip1(iBndGP,iFace) + cos1 * DISC%DynRup%output_Slip2(iBndGP,iFace)
              ENDIF

              IF (DynRup_output%OutputMask(7).EQ.1) THEN
                ! Thomas Ulrich 3.08.2015
                ! Calculation of rupture velocity from rupture time array
                ! using the 2d basis functions and their derivatives

                ! To compute accurately dt/dxi, all rupture_time have to be defined on the element
                compute_VR = .TRUE.
                DO jBndGP = 1,DISC%Galerkin%nBndGP
                   IF (DISC%DynRup%output_rupture_time(jBndGP,iFace).EQ.0) THEN
                      compute_VR = .FALSE.
                      Vr = 0d0
                   ENDIF
                ENDDO

                IF (compute_VR.EQ..TRUE.) THEN

                ! The rupture velocity is based on the gradient of the rupture time
                ! The gradient is not expected to be accurate for the side gauss points
                !! We consider then the next internal gauss point if the closest gauss point is external
                i1 = int((iBndGP-1)/(DISC%Galerkin%nPoly + 2))+1
                j1 = modulo((iBndGP-1),(DISC%Galerkin%nPoly + 2))+1
                IF (i1.EQ.1) THEN
                   i1=i1+1
                ELSE IF (i1.EQ.(DISC%Galerkin%nPoly + 2)) THEN
                   i1=i1-1
                ENDIF
                IF (j1.EQ.1) THEN
                   j1=j1+1
                ELSE IF (j1.EQ.(DISC%Galerkin%nPoly + 2)) THEN
                   j1=j1-1
                ENDIF
                iBndGP = (i1-1)*(DISC%Galerkin%nPoly + 2)+j1

                !Uncomment for a selecting only the central GP
                !i1  = int((DISC%Galerkin%nPoly + 2 -1)/2)+1
                !iBndGP = (i1-1)*(DISC%Galerkin%nPoly + 2)+i1

                !project rupture_time onto the basis functions
                nDegFr2d = (DISC%Galerkin%nPoly + 1)*(DISC%Galerkin%nPoly + 2)/2
                !projection of the rupture time on the basis functions
                ALLOCATE(projected_RT(nDegFr2d))
                projected_RT(:)=0.

                DO jBndGP = 1,DISC%Galerkin%nBndGP
                  chi  = MESH%ELEM%BndGP_Tri(1,jBndGP)
                  tau  = MESH%ELEM%BndGP_Tri(2,jBndGP)
                  DO iDegFr = 1, nDegFr2d
                     call BaseFunc_Tri(phiT,iDegFr,chi,tau,DISC%Galerkin%nPoly,DISC)
                     projected_RT(iDegFr) = projected_RT(iDegFr) + &
                        DISC%Galerkin%BndGaussW_Tet(jBndGP)*DISC%DynRup%output_rupture_time(jBndGP,iFace)*phiT
                  ENDDO
                ENDDO
                DO iDegFr = 1,nDegFr2d
                   projected_RT(iDegFr) =  projected_RT(iDegFr)/DISC%Galerkin%MassMatrix_Tri(iDegFr,iDegFr,DISC%Galerkin%nPoly)
                ENDDO

                !calculation of the spatial derivatives of the rupture time
                chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
                tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
                dt_dchi=0d0
                dt_dtau=0d0
                DO iDegFr = 1, nDegFr2d
                   call BaseGrad_Tri(phi2T,iDegFr,chi,tau,DISC%Galerkin%nPoly,DISC)
                   dt_dchi = dt_dchi + projected_RT(iDegFr) * phi2T(1)
                   dt_dtau = dt_dtau + projected_RT(iDegFr) * phi2T(2)
                ENDDO

                !For calculating the 2d Jacobian matrix we will need the nodal coordinates
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

                SELECT CASE(iSide)
                 CASE(1)
                    xab = (/ xV(3)-xV(1), yV(3)-yV(1), zV(3)-zV(1) /)
                    xac = (/ xV(2)-xV(1), yV(2)-yV(1), zV(2)-zV(1) /)
                 CASE(2)
                    xab = (/ xV(2)-xV(1), yV(2)-yV(1), zV(2)-zV(1) /)
                    xac = (/ xV(4)-xV(1), yV(4)-yV(1), zV(4)-zV(1) /)
                 CASE(3)
                    xab = (/ xV(4)-xV(1), yV(4)-yV(1), zV(4)-zV(1) /)
                    xac = (/ xV(3)-xV(1), yV(3)-yV(1), zV(3)-zV(1) /)
                 CASE(4)
                    xab = (/ xV(3)-xV(2), yV(3)-yV(2), zV(3)-zV(2) /)
                    xac = (/ xV(4)-xV(2), yV(4)-yV(2), zV(4)-zV(2) /)
                ENDSELECT

                grad2d(1,1) = dot_product(NormalVect_s(:),xab)
                grad2d(1,2) = dot_product(NormalVect_t(:),xab)
                grad2d(2,1) = dot_product(NormalVect_s(:),xac)
                grad2d(2,2) = dot_product(NormalVect_t(:),xac)
                CALL MatrixInverse2x2(JacobiT2d,grad2d)

                dt_dx1 = JacobiT2d(1,1) * dt_dchi + JacobiT2d(1,2) * dt_dtau
                dt_dy1 = JacobiT2d(2,1) * dt_dchi + JacobiT2d(2,2) * dt_dtau
                Slowness = sqrt(dt_dx1**2 + dt_dy1**2)

                IF (Slowness.EQ.0d0) THEN
                  Vr=0d0
                ElSE
                  Vr=1d0/Slowness
                ENDIF
                DEALLOCATE(projected_RT)
                iBndGP = DynRup_output%OutInt(iOutPoints,1)
                ENDIF

                OutVars = OutVars + 1
                DynRup_output%OutVal(iOutPoints,1,OutVars)  = Vr

              ENDIF
              IF (DynRup_output%OutputMask(8).EQ.1) THEN
                  OutVars = OutVars + 1
                  DynRup_output%OutVal(iOutPoints,1,OutVars) = DISC%DynRup%output_Slip(iBndGP,iFace)
              ENDIF
              IF (DynRup_output%OutputMask(9).EQ.1) THEN
                  OutVars = OutVars + 1
                  DynRup_output%OutVal(iOutPoints,1,OutVars) = DISC%DynRup%output_PeakSR(iBndGP,iFace)
              ENDIF
              IF (DynRup_output%OutputMask(10).EQ.1) THEN
                  OutVars = OutVars + 1
                  DynRup_output%OutVal(iOutPoints,1,OutVars) = DISC%DynRup%output_rupture_time(iBndGP,iFace)
              ENDIF
              IF (DynRup_output%OutputMask(11).EQ.1) THEN
                  OutVars = OutVars + 1
                  DynRup_output%OutVal(iOutPoints,1,OutVars) = DISC%DynRup%output_dynStress_time(iBndGP,iFace)
              ENDIF
          ENDIF


          ! Store output
          IF (DISC%DynRup%OutputPointType.NE.4) THEN
          DynRup_output%CurrentPick(iOutPoints) = DynRup_output%CurrentPick(iOutPoints) +1
          !
          DynRup_output%TmpTime(DynRup_output%CurrentPick(iOutPoints)) = time
          !
          ELSE
          DynRup_output%CurrentPick(iOutPoints) = 1
          ENDIF
          DynRup_output%TmpState(iOutPoints,DynRup_output%CurrentPick(iOutPoints),:) = DynRup_output%OutVal(iOutPoints,1,:)
          !
    ENDDO ! iOutPoints = 1,nOutPoints
    !
    CONTINUE

    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()
    !
  END SUBROUTINE calc_FaultOutput

  !
  !> Subroutine initializing the fault output
  !<
  SUBROUTINE write_FaultOutput_elementwise(EQN, DISC, MESH, IO, MPI, MaterialVal, BND, time, dt)
      !-------------------------------------------------------------------------!
      USE JacobiNormal_mod
      use FaultWriter
      use, intrinsic :: iso_c_binding
      !-------------------------------------------------------------------------!
      IMPLICIT NONE
      !-------------------------------------------------------------------------!
      ! Argument list declaration
      TYPE(tEquations)         :: EQN
      TYPE(tDiscretization)    :: DISC
      TYPE(tUnstructMesh)      :: MESH
      TYPE(tInputOutput)       :: IO                                            ! IO structure
      TYPE(tMPI)               :: MPI                                           ! MPI
      TYPE(tBoundary)          :: BND                                           ! BND    data structure
      REAL                     :: MaterialVal(MESH%nElem,EQN%nBackgroundVar)    ! Local Mean Values
      REAL                     :: dt, time
      !-------------------------------------------------------------------------!
      call writeFault(time)
  END SUBROUTINE ! write_FaultOutput_elementwise
  !
  !> writing fault output at pickpoints to files
  !<
  SUBROUTINE write_FaultOutput_atPickpoint(EQN, DISC, MESH, IO, MPI, MaterialVal, BND, time, dt)
    !-------------------------------------------------------------------------!
    USE JacobiNormal_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)              :: EQN
    TYPE(tDiscretization)         :: DISC
    TYPE(tUnstructMesh)           :: MESH
    TYPE(tInputOutput)            :: IO                                         ! IO structure
    TYPE(tMPI)                    :: MPI                                        ! MPI
    TYPE(tBoundary)               :: BND                                        ! BND    data structure
    REAL                          :: MaterialVal(MESH%nElem,EQN%nBackgroundVar) ! Local Mean Values
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    REAL    :: dt, time                                                       ! Timestep and time
    INTEGER :: stat
    INTEGER :: nOutPoints
    INTEGER :: iOutPoints, k                                                  ! Loop variables                   !
    REAL, ALLOCATABLE       :: TmpStateIO(:)                                  ! local variable to avoid fortran runtime copy
    CHARACTER (LEN=5)       :: cmyrank
    CHARACTER (len=200)     :: ptsoutfile
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: EQN, MESH, MaterialVal, time
    INTENT(INOUT) :: DISC
    !-------------------------------------------------------------------------!
    !

    ALLOCATE ( TmpStateIO( UBOUND(DISC%DynRup%DynRup_out_atPickpoint%TmpState, 3) ) )
    nOutPoints = DISC%DynRup%DynRup_out_atPickpoint%nDR_pick                   ! number of output receivers for this MPI domain
    DO iOutPoints = 1,nOutPoints                                               ! loop over number of output receivers for this domain
          !
          IF(DISC%DynRup%DynRup_out_atPickpoint%CurrentPick(iOutPoints).GE.DISC%DynRup%DynRup_out_atPickpoint%MaxPickStore.OR.ABS(DISC%EndTime-time).LE.(dt*1.005d0)) THEN
#ifdef PARALLEL
            WRITE(cmyrank,'(I5.5)') MPI%myrank                                   ! myrank -> cmyrank
            WRITE(ptsoutfile, '(a,a15,i5.5,a1,a5,a4)') TRIM(IO%OutputFile),'-faultreceiver-',DISC%DynRup%DynRup_out_atPickpoint%RecPoint(iOutPoints)%globalreceiverindex,'-',TRIM(cmyrank),'.dat'!
#else
            WRITE(ptsoutfile, '(a,a15,i5.5,a4)') TRIM(IO%OutputFile),'-faultreceiver-',DISC%DynRup%DynRup_out_atPickpoint%RecPoint(iOutPoints)%globalreceiverindex,'.dat'
#endif
            !
            OPEN(UNIT     = DISC%DynRup%DynRup_out_atPickpoint%VFile(iOutPoints)                    , & !
                 FILE     = ptsoutfile                                       , & !
                 FORM     = 'FORMATTED'                                      , & !
                 STATUS   = 'OLD'                                            , & !
                 POSITION = 'APPEND'                                         , & !
                 RECL     = 80000                                            , & !
                 IOSTAT = stat                                                 ) !
            IF( stat.NE.0) THEN
               logError(*) 'cannot open ',ptsoutfile
               logError(*) 'Error status: ', stat
               STOP
            END IF
            !
            DO k=1,DISC%DynRup%DynRup_out_atPickpoint%CurrentPick(iOutPoints)
               TmpStateIO(:) = DISC%DynRup%DynRup_out_atPickpoint%TmpState(iOutPoints,k,:)
               WRITE(DISC%DynRup%DynRup_out_atPickpoint%VFile(iOutPoints),*) DISC%DynRup%DynRup_out_atPickpoint%TmpTime(k), TmpStateIO(:)
            ENDDO
            !
            DISC%DynRup%DynRup_out_atPickpoint%CurrentPick(iOutPoints) = 0
            !
            CLOSE( DISC%DynRup%DynRup_out_atPickpoint%VFile(iOutPoints) ) ! not the same numbers as in elementwise output
            !
          ENDIF
    ENDDO ! iOutPoints = 1,DISC%DynRup%nOutPoints
    DEALLOCATE ( TmpStateIO )
  END SUBROUTINE write_FaultOutput_atPickpoint
 !
 END MODULE faultoutput_mod
