!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2011, SeisSol Group
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
!! Dynamic Rupture Module: solves frictional sliding

#include <Initializer/preProcessorMacros.fpp>

MODULE Friction_mod
  !---------------------------------------------------------------------------!
  USE TypesDef

#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc
#endif
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE Friction
     MODULE PROCEDURE Friction
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Friction
  PRIVATE :: Get_Extrapolated_Boundary_Values
  PRIVATE :: space_time_integration
  PRIVATE :: get_godunov_state
  !---------------------------------------------------------------------------!
  CONTAINS


  !################################################################################################################
  !###                                             Unified Version                                              ###
  !################################################################################################################

  ! This subroutine computes the traction state imposed by a chosen friction law, given some slip, slip rate and
  ! normal pressure values at each space-time gaussian point
  ! The flux contribution is also handled by this routine for all variables at each space-time gaussian point using
  ! the Godunov state. Therefore, skip the flux contribution in galerkin3d_solver for a fault side (case (3)) and
  ! ensure that all flux calculations using the Godunov method
  !> main friction routine: calls all subroutines handling the computations
  !<
  SUBROUTINE Friction(EQN, DISC, MESH, MPI, IO, OptionalFields, BND, time, dt)
    !-------------------------------------------------------------------------!
    USE CauchyKovalewski_mod
    USE JacobiNormal_mod
    USE gauss_mod
    USE Eval_friction_law_mod
#ifdef GENERATEDKERNELS
    use f_ftoc_bind_interoperability
#endif

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    real                           :: dt
    TYPE(tEquations)               :: EQN
    TYPE(tDiscretization), target  :: DISC
    TYPE(tUnstructMesh)            :: MESH
    TYPE(tMPI)                     :: MPI                                     !
    TYPE(tUnstructOptionalFields)  :: OptionalFields                                ! OptionalFields data strucure
    TYPE(tBoundary)                :: BND                                           ! BND    data structure
    TYPE(tInputOutput)             :: IO
    REAL                           :: time
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: i                                                          ! Loop variable                    !
    INTEGER     :: iElem, LocElemType                               ! Element number                   !
    INTEGER     :: iNeighbor                                                  ! The element's neighbor           !
    INTEGER     :: iLocalNeighborSide                                         ! Local side in the neighbor       !
    INTEGER     :: iSide                                                      ! Local element side number        !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom       !
    INTEGER     :: iVar                                                       ! Counter for cons/prim variable   !
    INTEGER     :: iBndGP, iTimeGP, iFace
    INTEGER     :: LocPoly, LocDegFr
    INTEGER     :: ReactionTerm
    INTEGER     :: iPoly, iTimePoly
    INTEGER     :: LocnVar,LocDegFrMat
    INTEGER     :: MPIIndex, iObject, MPIIndex_DR
    REAL        :: NormalVect_n(3)                                            ! Normal vector components         !
    REAL        :: NormalVect_s(3)                                            ! Normal vector components         !
    REAL        :: NormalVect_t(3)                                            ! Normal vector components         !
    REAL        :: T(EQN%nVar,EQN%nVar)                                       ! Transformation matrix            !
    REAL        :: iT(EQN%nVar,EQN%nVar)                                      ! inverse Transformation matrix    !
    REAL        :: dudt(DISC%Galerkin%nDegFr,EQN%nVar+EQN%nAneFuncperMech)    ! Time derivative of degs. freedom !
    REAL        :: w_speed(EQN%nNonZeroEV),w_speed_neig(EQN%nNonZeroEV)
    REAL        :: godunov_state(DISC%Galerkin%nDegFr,EQN%nVar)               ! Auxilliary variable              !
    REAL        :: auxMatrix(EQN%nVar,EQN%nVar)                               ! Auxilliary matrix                !
    REAL        :: AniVec(3,3)
    REAL        :: TaylorDOF(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)  ! time - taylorseries for DOF
    REAL        :: Taylor1(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)
    REAL        :: Taylor2(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)
    REAL        :: BndVar1(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)
    REAL        :: BndVar2(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)    
    REAL        :: NorStressGP(1:DISC%Galerkin%nBndGP,1:DISC%Galerkin%nTimeGP)
    REAL        :: TractionGP_XY(1:DISC%Galerkin%nBndGP,1:DISC%Galerkin%nTimeGP)
    REAL        :: TractionGP_XZ(1:DISC%Galerkin%nBndGP,1:DISC%Galerkin%nTimeGP)
    REAL        :: UVelGP(1:DISC%Galerkin%nBndGP,1:DISC%Galerkin%nTimeGP)
    REAL        :: XYStressGP(1:DISC%Galerkin%nBndGP,1:DISC%Galerkin%nTimeGP)
    REAL        :: XZStressGP(1:DISC%Galerkin%nBndGP,1:DISC%Galerkin%nTimeGP)
    REAL        :: Trac_XY_DOF(DISC%Galerkin%nDegFr), Trac_XZ_DOF(DISC%Galerkin%nDegFr)
    REAL        :: Trac_XY_NgbDOF(DISC%Galerkin%nDegFr), Trac_XZ_NgbDOF(DISC%Galerkin%nDegFr)
    REAL        :: NorStressDOF(DISC%Galerkin%nDegFr),NorStress_NgbDOF(DISC%Galerkin%nDegFr)
    REAL        :: UVelDOF(DISC%Galerkin%nDegFr),UVel_NgbDOF(DISC%Galerkin%nDegFr)
    REAL        :: rho, rho_neig, mu, mu_neig, lambda, lambda_neig, rho_minus
    REAL        :: geoSurface
    REAL        :: TimeIntDof_iElem(DISC%Galerkin%nDegFr,EQN%nVar)            ! Time integrated dof
    REAL        :: TimeIntDof_iNeigh(DISC%Galerkin%nDegFr,EQN%nVar)           ! Time integrated dof
    REAL        :: FluxInt(DISC%Galerkin%nDegFr,DISC%Galerkin%nDegFr)         ! auxilary variable to store Flux Integration
    REAL        :: phi_array(DISC%Galerkin%nDegFr,DISC%Galerkin%nBndGP)
    !
#ifndef GENERATEDKERNELS
    REAL, POINTER                  :: DOFiElem_ptr(:,:) => NULL()             ! Actual dof
    REAL, POINTER                  :: DOFiNeigh_ptr(:,:) => NULL()            ! Actual dof
#endif

#ifndef GENERATEDKERNELS
    TYPE(tSparseTensor3), POINTER  :: AStar_Sp_ptr => NULL()                  ! Pointer to the star tensor
    TYPE(tSparseTensor3), POINTER  :: BStar_Sp_ptr => NULL()                  ! Pointer to the star tensor
    TYPE(tSparseTensor3), POINTER  :: CStar_Sp_ptr => NULL()                  ! Pointer to the star tensor
    TYPE(tSparseTensor3), POINTER  :: EStar_Sp_ptr => NULL()                  ! Pointer to the star tensor
    TYPE(tSparseTensor3), POINTER  :: AStar_Neighbor_Sp_ptr => NULL()         ! Pointer to the star tensor
    TYPE(tSparseTensor3), POINTER  :: BStar_Neighbor_Sp_ptr => NULL()         ! Pointer to the star tensor
    TYPE(tSparseTensor3), POINTER  :: CStar_Neighbor_Sp_ptr => NULL()         ! Pointer to the star tensor
    TYPE(tSparseTensor3), POINTER  :: EStar_Neighbor_Sp_ptr => NULL()         ! Pointer to the star tensor
    TYPE(tSparseTensor3b), POINTER :: Tens_xi_Sp_ptr   => NULL()              ! Pointer to tSparseTensor3b
    TYPE(tSparseTensor3b), POINTER :: Tens_eta_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_zeta_Sp_ptr => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_klm_Sp_ptr  => NULL()              !
#endif

    !
    REAL                           :: JacobiDet,JacobiDet_Neig, l_deltaTLower
    REAL                           :: TmpMat(EQN%nBackgroundVar)

    !-------------------------------------------------------------------------!
    INTENT(IN)    :: MESH, IO, BND, time, dt
    INTENT(INOUT) :: EQN, DISC

    ! register epik/scorep function Friction
    EPIK_FUNC_REG("Friction")
    SCOREP_USER_FUNC_DEFINE()
    !-------------------------------------------------------------------------!
    ! start the epik/scorep function Friction
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN("Friction")
    !
    ! set lower bound of the integration interval to zero
    l_deltaTLower = 0
    !
    ReactionTerm = 0
    !
    ! Determine Gaussian integration points in time
    call gausslegendre(0.,dt,DISC%Galerkin%TimeGaussP,DISC%Galerkin%TimeGaussW)
    !
    ! Precalculate the powers of the time GPs divided by the corresponding faculty of the Taylor series
    DO i = 0,DISC%Galerkin%nPoly
#ifndef NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
       DO iTimeGP = 1, DISC%Galerkin%nTimeGP
#else
       do iTimeGp = 1, NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
#endif
          DISC%Galerkin%dtPowerFactor(i,iTimeGP) = (DISC%Galerkin%TimeGaussP(iTimeGP)**i)/DISC%Galerkin%Faculty(i)
       ENDDO
    ENDDO
    !
    ! for p-adaptiviy + openMP the attribute of these variables must be changed in the openMP pragma and asigns moved into the iFace loop
    LocPoly     = DISC%Galerkin%nPoly
    LocDegFr    = DISC%Galerkin%nDegFr
    LocDegFrMat = DISC%Galerkin%nDegFrMat
    LocnVar     = EQN%nVar
    !
    !
#ifdef OMP
#ifndef GENERATEDKERNELS
     !$omp parallel private(iFace,iElem,iSide,iNeighbor,iLocalNeighborSide,LocElemType,geoSurface,NormalVect_n,NormalVect_s,NormalVect_t,T,iT,iObject,MPIIndex,MPIIndex_DR,DOFiElem_ptr,AStar_Sp_ptr,BStar_Sp_ptr,CStar_Sp_ptr,EStar_Sp_ptr,TmpMat,rho,mu,lambda,w_speed,AniVec,JacobiDet,DOFiNeigh_ptr,AStar_Neighbor_Sp_ptr,BStar_Neighbor_Sp_ptr,CStar_Neighbor_Sp_ptr,EStar_Neighbor_Sp_ptr,rho_neig,mu_neig,lambda_neig,w_speed_neig,rho_minus,FluxInt,phi_array,godunov_state,JacobiDet_neig,TimeIntDof_iElem,TimeIntDof_iNeigh,TaylorDOF,Taylor1,Taylor2,BndVar1,BndVar2,NorStressGP,XYStressGP,XZStressGP,UVelGP,iDegFr,iTimePoly,iTimeGP,iBndGP,TractionGP_XY,TractionGP_XZ,Trac_XY_DOF,Trac_XZ_DOF,NorStressDOF,UVelDOF,Trac_XY_NgbDOF,Trac_XZ_NgbDOF,NorStress_NgbDOF,UVel_NgbDOF,iVar,auxMatrix,dudt,iPoly)  private(Tens_xi_Sp_ptr,Tens_eta_Sp_ptr,Tens_zeta_Sp_ptr,Tens_klm_Sp_ptr)  shared(MESH,DISC,EQN,BND,IO,optionalFields,LocPoly,LocnVar,LocDegFr,LocDegFrMat,ReactionTerm,MPI) shared(dt,time) default(none)
     !$omp do schedule(static)
#else
    ! TODO, @breuera: what a mess..
    !$omp parallel private(iElem, iFace, iSide, iNeighbor, iLocalNeighborSide,locElemType, geoSurface, normalVect_n, normalVect_s, normalVect_t, t, iT, iObject, mpiIndex, mpiIndex_dr, tmpMat, rho, mu, lambda, w_speed, aniVec, jacobiDet, rho_neig, mu_neig, lambda_neig, w_speed_neig, rho_minus, fluxInt, phi_array, godunov_state, jacobiDet_neig, timeIntDof_iElem, timeIntDof_iNeigh, taylorDOF, taylor1, taylor2, bndVar1, bndVar2, norStressGP, xYStressGP, xZStressGP, uVelGP, iDegFr, iTimePoly, iTimeGP, iBndGP, tractionGP_XY, tractionGP_XZ, trac_XY_DOF, trac_XZ_DOF, norStressDOF, uVelDOF, trac_XY_NgbDOF, trac_XZ_NgbDOF, norStress_NgbDOF, uVel_NgbDOF, iVar, auxMatrix, dudt, iPoly ) shared( mesh, disc, eqn, bnd, io, optionalFields, locPoly, locnVar, locDegFr, locDegFrMat, reactionTerm, mpi, dt, time, l_deltaTLower) default( none ) 
    !$omp do schedule(static)
#endif
#endif
    DO iFace=1,mesh%fault%nSide
       !-----------------------------------------------------------------------------------------------------------
       ! STEP 1: Obtain time expansions of the elements at both sides of the fault, on the fault's reference system
       !-----------------------------------------------------------------------------------------------------------
       !
       !
       iElem               = MESH%Fault%Face(iFace,1,1)          ! Remark:
       iSide               = MESH%Fault%Face(iFace,2,1)          ! iElem denotes "+" side
       !
       iNeighbor           = MESH%Fault%Face(iFace,1,2)          ! iNeighbor denotes "-" side
       iLocalNeighborSide  = MESH%Fault%Face(iFace,2,2)
       !
       !
       IF (iElem == 0) THEN
         ! iElem is in the neighbor domain
         ! obtain parameters from "-" element

         LocElemType     = MESH%LocalElemType(iNeighbor)
         geoSurface      = 2.0D0*DISC%Galerkin%geoSurfaces(iLocalNeighborSide,iNeighbor)
         !
!        ! Local side's normal vector from "+" side = iElem
         NormalVect_n = MESH%Fault%geoNormals(1:3,iFace)
         NormalVect_s = MESH%Fault%geoTangent1(1:3,iFace)
         NormalVect_t = MESH%Fault%geoTangent2(1:3,iFace)
         !
         CALL RotationMatrix3D(NormalVect_n,NormalVect_s,NormalVect_t,T(:,:),iT(:,:),EQN)
         !
         ! IF (MESH%ELEM%MPIReference(iLocalNeighborSide,iNeighbor).NE.1) STOP 'MPI Reference error in friction' ! just for testing
         !
         ! The neighbor element belongs to a different MPI domain
         iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
         MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
         MPIIndex_DR = MESH%ELEM%MPINumber_DR(iLocalNeighborSide,iNeighbor)

         ! IF(MPIIndex.EQ.-1) THEN
         !   WRITE(IO%UNIT%errOut,*) 'Severe Error in friction. MPIIndex = -1 !', iNeighbor, iLocalNeighborSide
         !  STOP
         ! ENDIF
         ! IF(.NOT.BND%ObjMPI(iObject)%Init) THEN
         !    WRITE(IO%UNIT%errOut,*) ' SeisSol MPI Error in friction. Domain not initialized. '
         !    WRITE(IO%UNIT%errOut,*) iObject, MPIIndex
         !   STOP
         ! ENDIF

#ifndef GENERATEDKERNELS
         DOFiElem_ptr    => BND%ObjMPI(iObject)%MPI_DR_dgvar(:,:,MPIIndex_DR)

         ! Take star matrices from neighbor MPI element
         AStar_Sp_ptr    => BND%ObjMPI(iObject)%AStar_Sp(MPIIndex)                ! Star matrices in sparse version
         BStar_Sp_ptr    => BND%ObjMPI(iObject)%BStar_Sp(MPIIndex)                !
         CStar_Sp_ptr    => BND%ObjMPI(iObject)%CStar_Sp(MPIIndex)                      !
         EStar_Sp_ptr    => BND%ObjMPI(iObject)%EStar_Sp(MPIIndex)                !
#endif

         TmpMat(:)=BND%ObjMPI(iObject)%NeighborBackground(:,MPIIndex)
         rho    = TmpMat(1)
         mu     = TmpMat(2)
         lambda = TmpMat(3)
         w_speed(1) = SQRT((lambda+2*mu)/rho)  ! Will only work in elastic isotropic cases
         w_speed(2) = SQRT(mu/rho)
         w_speed(3) = w_speed(2)


       ELSE
         ! normal case: iElem present in local domain
         LocElemType     = MESH%LocalElemType(iElem)
         w_speed(:)      = DISC%Galerkin%WaveSpeed(iElem,iSide,:)
         rho             = OptionalFields%BackgroundValue(iElem,1)
         geoSurface      = 2.0D0*DISC%Galerkin%geoSurfaces(iSide,iElem)
         !
         ! Local side's normal and tangential vectors
         NormalVect_n = MESH%Fault%geoNormals(1:3,iFace)
         NormalVect_s = MESH%Fault%geoTangent1(1:3,iFace)
         NormalVect_t = MESH%Fault%geoTangent2(1:3,iFace)
         !
         CALL RotationMatrix3D(NormalVect_n,NormalVect_s,NormalVect_t,T(:,:),iT(:,:),EQN)
         !
         !
#ifndef GENERATEDKERNELS
         DOFiElem_ptr    => DISC%Galerkin%dgvar(1:LocDegFr,:,iElem,1)         ! Actual time degrees of freedom
         !

         ! Set Pointer for CK evaluation
         AStar_Sp_ptr    => DISC%Galerkin%AStar_Sp(iElem)                      ! Star matrices in sparse version
         BStar_Sp_ptr    => DISC%Galerkin%BStar_Sp(iElem)                      !
         CStar_Sp_ptr    => DISC%Galerkin%CStar_Sp(iElem)                      !
         EStar_Sp_ptr    => DISC%Galerkin%EStar_Sp(iElem)                      !
#endif

         JacobiDet       =  6.0d0*MESH%ELEM%Volume(iElem)

       ENDIF ! (iElem == 0)
       !
       !
       !
       IF (iNeighbor == 0) THEN
         ! iNeighbor is in the neighbor domain
         !
         ! IF (MESH%ELEM%MPIReference(iSide,iElem).NE.1) STOP 'MPI Reference error in friction' ! just for testing
         !
         ! The neighbor element belongs to a different MPI domain
         iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
         MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
         MPIIndex_DR = MESH%ELEM%MPINumber_DR(iSide,iElem)

#ifndef GENERATEDKERNELS
         ! IF(MPIIndex.EQ.-1) THEN
         !   WRITE(IO%UNIT%errOut,*) 'Severe Error in Galerkin3D in friction. MPIIndex = -1 for iNeighbor!', iElem, iSide
         !   STOP
         ! ENDIF
         ! IF(.NOT.BND%ObjMPI(iObject)%Init) THEN
         !    WRITE(IO%UNIT%errOut,*) ' SeisSol MPI Error. Domain not initialized. '
         !    WRITE(IO%UNIT%errOut,*) iObject, MPIIndex
         !    STOP
         ! ENDIF
         DOFiNeigh_ptr   => BND%ObjMPI(iObject)%MPI_DR_dgvar(:,:,MPIIndex_DR)

         ! Take star matrices from neighbor MPI element
         AStar_Neighbor_Sp_ptr    => BND%ObjMPI(iObject)%AStar_Sp(MPIIndex)                ! Star matrices in sparse version
         BStar_Neighbor_Sp_ptr    => BND%ObjMPI(iObject)%BStar_Sp(MPIIndex)                !
         CStar_Neighbor_Sp_ptr    => BND%ObjMPI(iObject)%CStar_Sp(MPIIndex)                !
         EStar_Neighbor_Sp_ptr    => BND%ObjMPI(iObject)%EStar_Sp(MPIIndex)                !
#endif

         TmpMat(:)=BND%ObjMPI(iObject)%NeighborBackground(:,MPIIndex)
         rho_neig    = TmpMat(1)
         mu_neig     = TmpMat(2)
         lambda_neig = TmpMat(3)
         w_speed_neig(1) = SQRT((lambda_neig+2*mu_neig)/rho_neig)  ! Will only work in elastic isotropic cases
         w_speed_neig(2) = SQRT(mu_neig/rho_neig)
         w_speed_neig(3) = w_speed_neig(2)

       ELSE
         ! normal case: iNeighbor present in local domain
         w_speed_neig(:) = DISC%Galerkin%WaveSpeed(iNeighbor,iLocalNeighborSide,:)
         rho_neig        = OptionalFields%BackgroundValue(iNeighbor,1)
         !

#ifndef GENERATEDKERNELS
         DOFiNeigh_ptr   => DISC%Galerkin%dgvar(1:LocDegFr,:,iNeighbor,1)

         !
         ! Set Pointer for CK evaluation
         AStar_Neighbor_Sp_ptr    => DISC%Galerkin%AStar_Sp(iNeighbor)         ! Star matrices in sparse version
         BStar_Neighbor_Sp_ptr    => DISC%Galerkin%BStar_Sp(iNeighbor)         !
         CStar_Neighbor_Sp_ptr    => DISC%Galerkin%CStar_Sp(iNeighbor)         !
         EStar_Neighbor_Sp_ptr    => DISC%Galerkin%EStar_Sp(iNeighbor)         !
#endif
         JacobiDet_neig           = 6.0d0*MESH%ELEM%Volume(iNeighbor)

       ENDIF ! (iNeighbor == 0)

#ifndef GENERATEDKERNELS

            Tens_xi_Sp_ptr   => DISC%Galerkin%ADGxi_Tet_Sp                    ! Stiff CK tensors for tetras
            Tens_eta_Sp_ptr  => DISC%Galerkin%ADGeta_Tet_Sp
            Tens_zeta_Sp_ptr => DISC%Galerkin%ADGzeta_Tet_Sp
            Tens_klm_sp_ptr  => DISC%Galerkin%ADGklm_Tet_Sp                   ! Tensor for space dependent reaction term

       ! Obtain time expansions from both sides    
        CALL CauchyKovalewski3D(                                            & !
                               TimeIntDOF    = TimeIntDof_iElem,            & ! Output
                               TimeDerDof    = TaylorDOF,                   & ! Output
                               Dof           = DOFiElem_ptr,                & ! Input
                               Dt            = Dt,                          & ! Input
                               A_Sp          = AStar_Sp_ptr,                & ! Input
                               B_Sp          = BStar_Sp_ptr,                & ! Input
                               C_Sp          = CStar_Sp_ptr,                & ! Input
                               E_Sp          = EStar_Sp_ptr,                & ! Input
                               ADGxi_Sp      = Tens_xi_Sp_ptr,              & ! Input
                               ADGeta_Sp     = Tens_eta_Sp_ptr,             & ! Input
                               ADGzeta_Sp    = Tens_zeta_Sp_ptr,            & ! Input
                               Mklm_Sp       = Tens_klm_sp_ptr,             & ! Input
                               ReactionTerm  = ReactionTerm,                & ! Input
                               LocDegFr      = LocDegFr,                    & ! Input
                               LocDegFrMat   = LocDegFrMat,                 & ! Input
                               LocPoly       = LocPoly,                     & ! Input
                               nVar          = LocnVar                      ) ! Input
       !
       Taylor1(1:LocDegFr,:,0:LocPoly)=TaylorDOF(1:LocDegFr,:,0:LocPoly)
        CALL CauchyKovalewski3D(                                            & !
                               TimeIntDOF    = TimeIntDof_iNeigh,           & ! Output
                               TimeDerDof    = TaylorDOF,                   & ! Output
                               Dof           = DOFiNeigh_ptr,               & ! Input
                               Dt            = Dt,                          & ! Input
                               A_Sp          = AStar_Neighbor_Sp_ptr,       & ! Input
                               B_Sp          = BStar_Neighbor_Sp_ptr,       & ! Input
                               C_Sp          = CStar_Neighbor_Sp_ptr,       & ! Input
                               E_Sp          = EStar_Neighbor_Sp_ptr,       & ! Input
                               ADGxi_Sp      = Tens_xi_Sp_ptr,              & ! Input
                               ADGeta_Sp     = Tens_eta_Sp_ptr,             & ! Input
                               ADGzeta_Sp    = Tens_zeta_Sp_ptr,            & ! Input
                               Mklm_Sp       = Tens_klm_sp_ptr,             & ! Input
                               ReactionTerm  = ReactionTerm,                & ! Input
                               LocDegFr      = LocDegFr,                    & ! Input
                               LocDegFrMat   = LocDegFrMat,                 & ! Input
                               LocPoly       = LocPoly,                     & ! Input
                               nVar          = LocnVar                      ) ! Input
       !
       Taylor2(1:LocDegFr,:,0:LocPoly)=TaylorDOF(1:LocDegFr,:,0:LocPoly)
       !
#else
       if( iElem .ne. 0 ) then
         call c_interoperability_getFaceDerInt( i_meshId                  = iElem             , \
                                                i_faceId                  = iSide             , \
                                                i_timeStepWidth           = dt                , \
                                                o_timeDerivativesCell     = taylor1(:,:,0)    , \
                                                o_timeDerivativesNeighbor = taylor2(:,:,0)    , \
                                                o_timeIntegratedCell      = timeIntDof_iElem  , \
                                                o_timeIntegratedNeighbor  = timeIntDof_iNeigh   \
                                              )
       else
         call c_interoperability_getFaceDerInt( i_meshId                  = iNeighbor          , \
                                                i_faceId                  = iLocalNeighborSide , \
                                                i_timeStepWidth           = dt                 , \
                                                o_timeDerivativesCell     = taylor2(:,:,0)     , \
                                                o_timeDerivativesNeighbor = taylor1(:,:,0)     , \
                                                o_timeIntegratedCell      = timeIntDof_iNeigh  , \
                                                o_timeIntegratedNeighbor  = timeIntDof_iElem     \
                                              )
       endif
#endif
       ! rotation of the quantities to face-aligned coordinate system
       do iPoly=0, LocPoly
           bndVar1( :, :, iPoly) = matmul( taylor1( :, :, iPoly), mesh%fault%forwardRotation( :, :, iFace) )
           bndVar2( :, :, iPoly) = matmul( taylor2( :, :, iPoly), mesh%fault%forwardRotation( :, :, iFace) )
       enddo

       !
       !
       !-----------------------------------------------------------------------------------------------------------
       ! STEP 2: Project space-time polynomials into gaussian points
       !-----------------------------------------------------------------------------------------------------------
       !
       !Discrete representation at space-time gaussian points of the variables at the boundary
       CALL Get_Extrapolated_Boundary_Values(             & !
                UVelGP,NorStressGP,XYStressGP,XZStressGP, & ! OUT: values at boundary GPs
                iFace,iSide,                              & ! IN: fault element ID
                LocDegFr,LocPoly,                         & ! IN: elements basis functions attributes
                BndVar1,BndVar2,                          & ! IN: rotated time expansions
                rho,rho_neig,w_speed,w_speed_neig,        & ! IN: material parameters
                MESH,DISC,EQN                             ) ! global variables

       !
       !-----------------------------------------------------------------------------------------------------------
       ! STEP 3: Impose traction (and corresponding velocity) given a certain friction law
       !-----------------------------------------------------------------------------------------------------------
       !
       CALL Eval_friction_law(                                & !
               TractionGP_XY,TractionGP_XZ,                   & ! OUT: updated Traction
               NorStressGP,XYStressGP,XZStressGP,             & ! IN: Godunov status
               iFace,iSide,iElem,time,iT,                     & ! IN: element ID, time, inv Trafo
               rho,rho_neig,w_speed(1:2),w_speed_neig(1:2),   & ! IN: background values
               EQN,DISC,MESH,MPI,IO,BND                       ) ! global variables

       !
       !-----------------------------------------------------------------------------------------------------------
       ! STEP 4: Compute exact fluxes for both elements and update element's state
       !         Including space and time integration of the traction, normal stress and fault perpendicular particle velocity
       !-----------------------------------------------------------------------------------------------------------
       !
       ! STEP 4a: For the element "+"
       IF (iElem .NE. 0) THEN

         ! space and time integration of the traction, normal stress and fault perpendicular particle velocity:
         phi_array(:,:) = MESH%ELEM%BndBF_GP_Tet(1:LocDegFr,1:DISC%Galerkin%nBndGP,iSide)
         CALL space_time_integration(                              &
              Trac_XY_DOF,Trac_XZ_DOF,NorStressDOF,UVelDOF,        & ! OUT: integrated variables
              TractionGP_XY,TractionGP_XZ,NorStressGP,UVelGP,      & ! IN: variables at GPs
              geoSurface,phi_array,                                & ! IN: geoSurface, phi: basis functions
              MESH,DISC                                            ) ! IN: global variables

         ! compute godunov state
         call get_godunov_state( godunov_state  = godunov_state, &
                                 timeIntDof     = timeIntDof_iElem, &
                                 trac_xy_dof    = trac_xy_dof, &
                                 trac_xz_dof    = trac_xz_dof, &
                                 norStressDof   = norStressDof, &
                                 uvelDof        = uvelDof, &
                                 fluxint        = disc%galerkin%fluxInt_tet( 1:locDegFr, 1:locDegFr, 1, 0, 1, iSide), &
                                 geoSurface     = geoSurface, &
                                 iT             = mesh%fault%forwardRotation( :, :, iFace), &
                                 w_speed        = w_speed(1:2), &
                                 rho            = rho, &
                                 disc           = disc )

         ! multiplication by face normal jacobian, inverse trafo-determinant and back rotation from face aligned to x-y-z coordinate system
         dudt(:,:) = matmul( godunov_state(:, :), mesh%fault%fluxSolver(:, :, 1, iFace) )

         ! multiply by inverse mass matrix
         do iDegFr = 1, LocDegFr
           dudt(iDegFr,:) = dudt(iDegFr, :) * (-disc%galerkin%iMassMatrix_tet(iDegFr, iDegFr, locPoly))
         enddo

         ! Write update buffer of fault elements, this is to avoid conflicts due to elements with more than one dynamic rupture face
         ! TODO: Once we are communicating time derivatives this step isn't required anymore.
         DISC%DynRup%DRupdates(1:LocDegFr,1:LocnVar,DISC%DynRup%DRupdatesPosition(iFace,1)) = dudt(1:LocDegFr,1:LocnVar)
       ENDIF

       ! STEP 4b: For the element "-"
       IF (iNeighbor .NE. 0) THEN

         ! space and time integration of the traction, normal stress and fault perpendicular particle velocity:
         phi_array(:,:) = DISC%DynRup%BndBF_GP_Tet(1:LocDegFr,1:DISC%Galerkin%nBndGP,iFace)
         CALL space_time_integration(                                        &
              Trac_XY_NgbDOF,Trac_XZ_NgbDOF,NorStress_NgbDOF,UVel_NgbDOF,    & ! OUT: integrated variables
              TractionGP_XY,TractionGP_XZ,NorStressGP,UVelGP,                & ! IN:  variables at GPs
              geoSurface,phi_array,                                          & ! IN:  geoSurface, phi: basis functions
              MESH,DISC                                                      ) ! IN: global variables

         ! compute godunov state
         rho_minus = -1.0D0*rho_neig ! ATTENTION: Introduce minus sign here to re-use get_godunov_state routine

         call get_godunov_state( godunov_state  = godunov_state, &
                                 timeIntDof     = timeIntDof_iNeigh, &
                                 trac_xy_dof    = trac_xy_ngbDof, &
                                 trac_xz_dof    = trac_xz_ngbDof, &
                                 norStressDof   = norStress_ngbDof, &
                                 uvelDof        = uvel_ngbDof, &
                                 fluxint        = disc%dynRup%fluxInt(1:locDegFr, 1:locDegFr, iFace), &
                                 geoSurface     = geoSurface, &
                                 iT             = mesh%fault%forwardRotation( :, :, iFace), &
                                 w_speed        = w_speed_neig(1:2), &
                                 rho            = rho_minus, &
                                 disc           = disc )

         ! multiplication by face normal jacobian, inverse trafo-determinant and back rotation from face aligned to x-y-z coordinate system
         dudt(:,:) = matmul( godunov_state(:, :), mesh%fault%fluxSolver(:, :, 2, iFace) )

         ! multiply by inverse mass matrix
         do iDegFr = 1, LocDegFr
           dudt(iDegFr,:) = dudt(iDegFr, :) * disc%galerkin%iMassMatrix_tet(iDegFr, iDegFr, locPoly)
         enddo
         
         ! Write update buffer of fault elements, this is to avoid conflicts due to elements with more than one dynamic rupture face
         ! TODO: Once we are communicating time derivatives this step isn't required anymore.
         DISC%DynRup%DRupdates(1:LocDegFr,1:LocnVar,DISC%DynRup%DRupdatesPosition(iFace,2)) = dudt(1:LocDegFr,1:LocnVar)
       ENDIF

#ifndef GENERATEDKERNELS
       ! Nullify pointers for next element
       NULLIFY(DOFiElem_ptr,DOFiNeigh_ptr)
       NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr, EStar_Sp_ptr)
       NULLIFY(AStar_Neighbor_Sp_ptr,BStar_Neighbor_Sp_ptr,CStar_Neighbor_Sp_ptr,EStar_Neighbor_Sp_ptr)
       NULLIFY(Tens_xi_Sp_ptr, Tens_eta_Sp_ptr, Tens_zeta_Sp_ptr, Tens_klm_sp_ptr)
#endif


    ENDDO ! iFace
#ifdef OMP
    !$omp end parallel   
#endif

#ifndef GENERATEDKERNELS
    ! Apply DR updates to dgvar
    ! TODO we might want to do this in parallel
    DO iFace=1,DISC%DynRup%nDRElems
      DISC%Galerkin%dgvar(1:DISC%Galerkin%nDegFr,1:EQN%nVar,DISC%DynRup%indicesOfDRElems(iFace),1) = DISC%Galerkin%dgvar(1:DISC%Galerkin%nDegFr,1:EQN%nVar,DISC%DynRup%indicesOfDRElems(iFace),1) + DISC%DynRup%DRupdates(1:DISC%Galerkin%nDegFr,1:EQN%nVar,iFace)
    ENDDO
#else
    do iFace=1,DISC%DynRup%nDRElems
      call c_interoperability_addToDofs( i_meshId           = DISC%DynRup%indicesOfDRElems(iFace), \
                                         i_update           = DISC%DynRup%DRupdates(1:DISC%Galerkin%nDegFr,1:EQN%nVar,iFace), \
                                         numberOfQuantities = EQN%nVar )
    enddo
#endif

    CONTINUE

  ! end the epik/scorep function Friction
  EPIK_FUNC_END()
  SCOREP_USER_FUNC_END()

  END SUBROUTINE Friction
  
  
  
  !> Obtain discrete representation at space-time gaussian points of the variables at the boundary
  !<
  PURE SUBROUTINE Get_Extrapolated_Boundary_Values(UVelGP,NorStressGP,XYStressGP,XZStressGP,iFace,iSide, &
                                                   LocDegFr,LocPoly,BndVar1,BndVar2,rho,rho_neig,        &
                                                   w_speed,w_speed_neig,MESH,DISC,EQN)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tUnstructMesh)            :: MESH                                    !
    TYPE(tDiscretization)          :: DISC                                    !
    TYPE(tEquations)               :: EQN                                     !
    !-------------------------------------------------------------------------!
    ! Local variable declaration 
    REAL        :: NorStressGP(DISC%Galerkin%nBndGP,DISC%Galerkin%nTimeGP)
    REAL        :: UVelGP(DISC%Galerkin%nBndGP,DISC%Galerkin%nTimeGP)
    REAL        :: XYStressGP(DISC%Galerkin%nBndGP,DISC%Galerkin%nTimeGP)
    REAL        :: XZStressGP(DISC%Galerkin%nBndGP,DISC%Galerkin%nTimeGP)
    REAL        :: BndVar1(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)
    REAL        :: BndVar2(DISC%Galerkin%nDegFr,EQN%nVarTotal,0:DISC%Galerkin%nPoly)
    REAL        :: LocNorStress, LocXYStress, LocXZStress, LocUVel
    REAL        :: rho,rho_neig,w_speed(2),w_speed_neig(2)
    REAL        :: phi1,phi2
    REAL        :: phi1_array(LocDegFr,DISC%Galerkin%nBndGP),phi2_array(LocDegFr,DISC%Galerkin%nBndGP)
    REAL        :: dtPowerFactor(0:LocPoly,1:DISC%Galerkin%nTimeGP)
    REAL        :: NorDivisor,ShearDivisor,UVelDivisor
    INTEGER     :: iFace,iSide
    INTEGER     :: iDegFr,iTimePoly,iTimeGP,iBndGP,nTimeGP,nBndGP
    INTEGER     :: LocDegFr,LocPoly
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: MESH,DISC,EQN,iFace,iSide,LocDegFr,LocPoly,BndVar1,BndVar2
    INTENT(IN)    :: rho,rho_neig,w_speed,w_speed_neig
    INTENT(INOUT) :: UVelGP,NorStressGP,XYStressGP,XZStressGP
    !-------------------------------------------------------------------------!
  
    NorStressGP(:,:) = 0.
    XYStressGP(:,:)  = 0.
    XZStressGP(:,:)  = 0.
    UVelGP(:,:)      = 0.
    
    ! auxiliary variables
#ifndef NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
    nTimeGP = DISC%Galerkin%nTimeGP
#else
    ntimeGP = NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
#endif
    nBndGP  = DISC%Galerkin%nBndGP
    phi1_array(:,:) = MESH%ELEM%BndBF_GP_Tet(1:LocDegFr,1:nBndGP,iSide)
    phi2_array(:,:) = DISC%DynRup%BndBF_GP_Tet(1:LocDegFr,1:nBndGP,iFace)
    dtPowerFactor(0:LocPoly,1:nTimeGP) = DISC%Galerkin%dtPowerFactor(0:LocPoly,1:nTimeGP)
    NorDivisor = 1.0D0/(w_speed_neig(1)*rho_neig+w_speed(1)*rho)
    ShearDivisor = 1.0D0/(w_speed_neig(2)*rho_neig+w_speed(2)*rho)
    UVelDivisor = 1.0D0/(w_speed(1)*rho) 
    !
    DO iDegFr=1,LocDegFr
      DO iTimePoly=0,LocPoly
        DO iTimeGP=1,nTimeGP
          DO iBndGP=1,nBndGP
            phi1 = phi1_array(iDegFr,iBndGP)
            phi2 = phi2_array(iDegFr,iBndGP)
            !
            LocNorStress    = phi1*BndVar1(iDegFr,1,iTimePoly) + &
                              ((phi2*BndVar2(iDegFr,1,iTimePoly)-phi1*BndVar1(iDegFr,1,iTimePoly) + &
                               w_speed_neig(1)*rho_neig*(phi2*BndVar2(iDegFr,7,iTimePoly)-phi1*BndVar1(iDegFr,7,iTimePoly))) * &
                               w_speed(1)*rho ) * NorDivisor
            LocXYStress     =  phi1*BndVar1(iDegFr,4,iTimePoly) + &
                              ((phi2*BndVar2(iDegFr,4,iTimePoly)-phi1*BndVar1(iDegFr,4,iTimePoly) + &
                               w_speed_neig(2)*rho_neig*(phi2*BndVar2(iDegFr,8,iTimePoly)-phi1*BndVar1(iDegFr,8,iTimePoly)))* &
                               w_speed(2)*rho ) * ShearDivisor
            LocXZStress     =  phi1*BndVar1(iDegFr,6,iTimePoly) + &
                              ((phi2*BndVar2(iDegFr,6,iTimePoly)-phi1*BndVar1(iDegFr,6,iTimePoly) + &
                               w_speed_neig(2)*rho_neig*(phi2*BndVar2(iDegFr,9,iTimePoly)-phi1*BndVar1(iDegFr,9,iTimePoly)))* &
                               w_speed(2)*rho ) * ShearDivisor

            LocUVel         = phi1*BndVar1(iDegFr,7,iTimePoly) + (LocNorStress-phi1*BndVar1(iDegFr,1,iTimePoly)) * UVelDivisor 

            UVelGP(iBndGP,iTimeGP)     = UVelGP(iBndGP,iTimeGP)     + LocUVel*dtPowerFactor(iTimePoly,iTimeGP)

            NorStressGP(iBndGP,iTimeGP) = NorStressGP(iBndGP,iTimeGP) + LocNorStress*dtPowerFactor(iTimePoly,iTimeGP)
            XYStressGP(iBndGP,iTimeGP)  = XYStressGP(iBndGP,iTimeGP)  + LocXYStress*dtPowerFactor(iTimePoly,iTimeGP)
            XZStressGP(iBndGP,iTimeGP)  = XZStressGP(iBndGP,iTimeGP)  + LocXZStress*dtPowerFactor(iTimePoly,iTimeGP)
            !
          ENDDO
        ENDDO
      ENDDO
    ENDDO
    !  
  END SUBROUTINE Get_Extrapolated_Boundary_Values

!> space time integration for traction, normal stress and fault perpendicular particle velocity
!<
  PURE SUBROUTINE space_time_integration(                              &
                  Trac_XY_DOF,Trac_XZ_DOF,NorStressDOF,UVelDOF,        & ! OUT: integrated variables
                  TractionGP_XY,TractionGP_XZ,NorStressGP,UVelGP,      & ! IN:  variables at GPs
                  geoSurface,phi_array,                                & ! IN:  geoSurface, phi: basis functions
                  MESH,DISC                                            ) ! IN:  global variables

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tUnstructMesh)            :: MESH                                    !
    TYPE(tDiscretization)          :: DISC                                    !
    !-------------------------------------------------------------------------!
    ! Local variable declaration 
    REAL        :: TractionGP_XY(:,:),TractionGP_XZ(:,:),NorStressGP(:,:),UVelGP(:,:)
    REAL        :: Trac_XY_DOF(:),Trac_XZ_DOF(:),NorStressDOF(:),UVelDOF(:)
    REAL        :: phi, phi_array(:,:)
    REAL        :: tmp_int
    REAL        :: geoSurface
    INTEGER     :: iBndGP,iDegFr,iTimeGP
    INTEGER     :: LocDegFr
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: MESH,DISC
    INTENT(IN)    :: TractionGP_XY,TractionGP_XZ,NorStressGP,UVelGP
    INTENT(IN)    :: geoSurface,phi_array
    INTENT(OUT)   :: Trac_XY_DOF,Trac_XZ_DOF,NorStressDOF,UVelDOF
    !-------------------------------------------------------------------------!
    !
    LocDegFr    = DISC%Galerkin%nDegFr
    Trac_XY_DOF  = 0.
    Trac_XZ_DOF  = 0.
    NorStressDOF = 0.
    UVelDOF      = 0.
    !
    DO iBndGP=1,DISC%Galerkin%nBndGP
      DO iDegFr=1,LocDegFr
        phi = phi_array(iDegFr,iBndGP)
#ifndef NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
        DO iTimeGP=1,DISC%Galerkin%nTimeGP
#else
        do iTimeGp=1,NUMBER_OF_TEMPORAL_INTEGRATION_POINTS
#endif
          tmp_int = phi*DISC%Galerkin%TimeGaussW(iTimeGP)*MESH%ELEM%BndGW_Tri(iBndGP)*geoSurface
          ! Traction:
          Trac_XY_DOF(iDegFr) = Trac_XY_DOF(iDegFr) + tmp_int * TractionGP_XY(iBndGP,iTimeGP)
          Trac_XZ_DOF(iDegFr) = Trac_XZ_DOF(iDegFr) + tmp_int * TractionGP_XZ(iBndGP,iTimeGP)
          ! NorStress
          NorStressDOF(iDegFr) = NorStressDOF(iDegFr) + tmp_int * NorStressGP(iBndGP,iTimeGP)
          ! UVel
          UVelDOF(iDegFr) = UVelDOF(iDegFr) + tmp_int * UVelGP(iBndGP,iTimeGP)
        ENDDO
      ENDDO
    ENDDO
    !
  END SUBROUTINE space_time_integration


!> get godunov state 
!<
  PURE SUBROUTINE get_godunov_state(                                        &
                  godunov_state,                                            & ! OUT: Godunov state
                  TimeIntDof,Trac_XY_DOF,Trac_XZ_DOF,NorStressDOF,UVelDOF,  & ! IN:  degrees of freedom
                  FluxInt,                                                  & ! IN:  Flux integration
                  geoSurface, iT,w_speed,rho,                               & ! IN:  geoSurface, Trafo, material properties
                  DISC                                                      ) ! IN:  global variables

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tDiscretization)          :: DISC                                    !< DISC global variable
    !-------------------------------------------------------------------------!
    ! Local variable declaration
    REAL        :: godunov_state(DISC%Galerkin%nDegFr,9)                      !< Godunov state
    REAL        :: TimeIntDof(DISC%Galerkin%nDegFr,9)                         !< Time integrated dof
    REAL        :: Trac_XY_DOF(:)                                             !< Time integrated dof
    REAL        :: Trac_XZ_DOF(:)                                             !< Time integrated dof
    REAL        :: NorStressDOF(:)                                            !< Time integrated dof
    REAL        :: UVelDOF(:)                                                 !< Time integrated dof
    REAL        :: FluxInt(:,:)                                               !< Flux integration matrix
    REAL        :: iT(:,:)                                                    !< inv Trafo
    REAL        :: geoSurface                                                 !< area of element surface
    REAL        :: rho                                                        !< densitiy rho
    REAL        :: w_speed(1:2)                                               !< wave speeds
    REAL        :: LocVarDOF(DISC%Galerkin%nDegFr,9)
    INTEGER     :: m,iDegFr,LocDegFr,iVar
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC,TimeIntDof,FluxInt
    INTENT(IN)    :: Trac_XY_DOF,Trac_XZ_DOF,NorStressDOF,UVelDOF 
    INTENT(IN)    :: geoSurface,rho,w_speed,iT
    INTENT(OUT)   :: godunov_state
    !-------------------------------------------------------------------------!
    !
    LocDegFr    = DISC%Galerkin%nDegFr
    !

    ! rotate time integrated DOFs from xyz- to fault-aligned coordinate system
    locVarDof = matmul( timeIntdof(:,:), iT )

    ! multiply flux matrix F^{-,i} with sub-matrix of time integrated dofs:
    !           _                 _
    !          | - - - * * * * * * |
    !          | - - - * * * * * * |
    !          | - - - * * * * * * |
    !          | - - - * * * * * * |
    !          | - - - * * * * * * |
    !          | - - - * * * * * * |
    !          | - - - * * * * * * |
    !          | - - - * * * * * * |
    !          |_- - - * * * * * *_|
    !                  | | | | | |
    !        stresses 1,2|1,3| | |
    !                   2,3  | | |
    !                        u v w particle
    !                              velocities
    ! 
    ! Remark: \f$\sigma^{2,3}\f$ and \f$u\f$ are replaced afterwards
    !
    godunov_state(:, 4:9) = matmul( fluxInt(:, :), locVarDof(:, 4:9) )
    !
    ! Multiplication with Surface area
    godunov_state(:,:) = godunov_state(:,:)*geoSurface
    !
    ! Create Godunov state at the interface
    godunov_state(:,8)  = godunov_state(:,8) + 1.0D0/(w_speed(2)*rho) * (Trac_XY_DOF(:)-godunov_state(:,4))
    godunov_state(:,9)  = godunov_state(:,9) + 1.0D0/(w_speed(2)*rho) * (Trac_XZ_DOF(:)-godunov_state(:,6))
    godunov_state(:,4)  = Trac_XY_DOF(:) !TractionDOF is already space-time integrated
    godunov_state(:,6)  = Trac_XZ_DOF(:) !TractionDOF is already space-time integrated
    godunov_state(:,1)  = NorStressDOF(:)
    godunov_state(:,2)  = 0.0D0 ! YYStressDOF(:)
    godunov_state(:,3)  = 0.0D0 ! ZZStressDOF(:)
    godunov_state(:,5)  = 0.0D0 ! YZStressDOF(:)
    godunov_state(:,7)  = UVelDOF(:)

  END SUBROUTINE get_godunov_state


END MODULE Friction_mod
