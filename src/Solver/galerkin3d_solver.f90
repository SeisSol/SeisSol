!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2007-2017, SeisSol Group
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
!! Sweeps over the mesh; in case of all non-classic code the used functionality is very limited.
!!

#include <Initializer/preProcessorMacros.fpp>

MODULE Galerkin3D_solver_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE COMMON_operators_mod
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  PUBLIC  :: ADERGalerkin3D_GTS
#ifndef GENERATEDKERNELS
  PUBLIC  :: ADERGalerkin3D_LTS
  PUBLIC  :: CKProcedureForEveryoneGTSUni
  PUBLIC  :: CKProcedureForEveryoneLTSUni
#endif
  !---------------------------------------------------------------------------!

CONTAINS

  !===========================================================================!
  !!                                                                         !!
  !!  ADERGalerkin3D evolves the ADER-DG solution for one timestep           !!
  !!                                                                         !!
  !===========================================================================!
  SUBROUTINE ADERGalerkin3D_GTS(time,dt,iteration,MaterialVal,OptionalFields,EQN,MESH,DISC,IC,SOURCE,BND,MPI,IO)
    !-------------------------------------------------------------------------!
    USE CauchyKovalewski_mod
    USE Galerkin_source_mod
    USE SP_MATMUL_mod
    USE gauss_mod
    USE Friction_mod
    USE faultoutput_mod
    USE Plasticity_mod
    use FaultWriter

#ifdef GENERATEDKERNELS
    use f_ftoc_bind_interoperability
    use iso_c_binding, only: c_loc
#endif
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)               :: EQN                                     !
    TYPE(tDiscretization), target  :: DISC                                    !
    TYPE(tUnstructMesh)            :: MESH                                    !
    TYPE(tUnstructOptionalFields)  :: OptionalFields                          !
    TYPE(tBoundary)                :: BND                                     !
    TYPE(tMPI)                     :: MPI                                     !
    TYPE(tInitialCondition)        :: IC                                      !
    TYPE(tSource)                  :: SOURCE                                  !
    TYPE(tInputOutput)             :: IO
    REAL                           :: time, dt                                ! Calculation time and timestep
    REAL                           :: MaterialVal(MESH%nElem,               & !
                                                  EQN%nBackgroundVar        ) ! Material Values
    INTEGER                        :: iteration                               ! Current iteration number
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !

    TYPE(tSparseTensor3), POINTER  :: AStar_Sp_ptr => NULL()                  ! Pointer to tSparseTensor3
    TYPE(tSparseTensor3), POINTER  :: BStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: CStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: EStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: FluxSolver_Sp_ptr => NULL()             !

    TYPE(tSparseTensor3b), POINTER :: Tens_xi_Sp_ptr   => NULL()              ! Pointer to tSparseTensor3b
    TYPE(tSparseTensor3b), POINTER :: Tens_eta_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_zeta_Sp_ptr => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_klm_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_k_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_k_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_k_Sp_ptr   => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_m_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_m_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_m_Sp_ptr   => NULL()              !

    TYPE(tSparseTensor3b), POINTER :: FluxInt_Sp_ptr    => NULL()             !

    REAL, POINTER                  :: MassMatrix_ptr(:,:) => NULL()           ! Pointer to the corresponding mass matrix
    REAL                           :: JacobiVol_MassMatrix(1:220)             ! 220 is dummy length good until O10!

    REAL, POINTER                  :: TimeIntDof_ptr(:,:)     => NULL()       ! Time integrated dof
    REAL, POINTER                  :: TimeDerDof_ptr(:,:,:)   => NULL()       ! Space-time polynomial pointer

    REAL                           :: dudt(DISC%Galerkin%nDegFr,            & !
                                           EQN%nVarTotal                    ) ! Auxilliary update dof
    REAL                           :: dudt_flux(DISC%Galerkin%nDegFr,       & !
                                                EQN%nVar                    ) ! Auxilliary update dof
    REAL                           :: dudt_src(DISC%Galerkin%nDegFr,        & !
                                               EQN%nVar                     ) ! Auxilliary update dof
    REAL                           :: dudt_plastic(DISC%Galerkin%nDegFr, 6  ) ! Auxilliary update dof for plastic calculations
    REAL                           :: auxvar(DISC%Galerkin%nDegFr,          & !
                                             EQN%nVar,                      & !
                                             DISC%Galerkin%nDegFrMat        ) ! Auxilliary variable
    REAL                           :: JacobiDetVolum                          ! Jacobian of the mapping
    REAL                           :: JacobiDetSurfa                          ! Jacobian of the mapping

    REAL                           :: CpuTime_ini                             ! Eval cpu time
    REAL                           :: CpuTime_end                             ! Eval cpu time
    REAL                           :: KineticEnergy_tmp                       ! kinetic energy
    INTEGER                        :: i
    INTEGER                        :: LocnVar                                 !
    INTEGER                        :: LocnPoly                                !
    INTEGER                        :: LocnDegFr                               !
    INTEGER                        :: LocnPolyMat                             !
    INTEGER                        :: LocnDegFrMat                            !
    INTEGER                        :: NeignVar                                !
    INTEGER                        :: NeignPoly                               !
    INTEGER                        :: NeignDegFr                              !
    INTEGER                        :: NeignPolyMat                            !
    INTEGER                        :: NeignDegFrMat                           !
    INTEGER                        :: iObject                                 !
    INTEGER                        :: MPIIndex                                !
    INTEGER                        :: iNeighbor                               !
    INTEGER                        :: iLocalNeighborSide                      !
    INTEGER                        :: iLocalNeighborVrtx                      !
    INTEGER                        :: iDegFr                                  ! Counter
    INTEGER                        :: iElem                                   ! Counter
    INTEGER                        :: rank_int                                ! Rank counter
#ifdef GENERATEDKERNELS
    INTEGER                        :: iElemST                                 ! Counter Source Terms
#endif
#ifndef GENERATEDKERNELS
    REAL, POINTER     :: DOFiElem_ptr(:,:)  => NULL()                         ! Actual dof
    REAL, POINTER     :: DOFiNeigh_ptr(:,:) => NULL()                         ! Actual dof
#else
    real, dimension( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES ) :: DOFiElem_ptr ! no: it's not a pointer..
    real, dimension( NUMBER_OF_BASIS_FUNCTIONS, NUMBER_OF_QUANTITIES ) :: DOFiNeigh_ptr ! no pointer again
#endif
    INTEGER                        :: iSide                                   ! Counter
    INTEGER                        :: iVar                                    ! Counter
    INTEGER                        :: nElem, nVarTotal
    REAL, POINTER :: IntGaussP(:,:)     =>NULL()
    REAL, POINTER :: IntGaussW(:)       =>NULL()
    REAL, POINTER :: IntGPBaseFunc(:,:) =>NULL()
    REAL, POINTER :: MassMatrix(:,:)    =>NULL()

    ! register ADERGalerkin3D_GTS function and epik/scorep region for volume integration, source terms and boundary integration
    EPIK_FUNC_REG("ADERGalerkin3D_GTS")
    SCOREP_USER_FUNC_DEFINE()
    ! register epik/scorep region for source terms
    EPIK_USER_REG(r_source, "sources")
    SCOREP_USER_REGION_DEFINE( r_source )
    ! register epik/scorep region dynamic rupture
    EPIK_USER_REG(r_dr, "dynamic_rupture")
    SCOREP_USER_REGION_DEFINE( r_dr )
    ! register epik/scorep region for dynamic rupture fault output
    EPIK_USER_REG(r_dr_output, "dynamic_rupture_fault_output")
    SCOREP_USER_REGION_DEFINE( r_dr_output )

    !-------------------------------------------------------------------------!

    ! start epik/scorep function ADERGalerkin3D_GTS
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN( "ADERGalerkin3D_GTS" )

    ! Determine Gaussian integration points in time
    call gausslegendre(0.,dt,DISC%Galerkin%TimeGaussP,DISC%Galerkin%TimeGaussW)
    !
    LocnPoly     = DISC%Galerkin%nPoly
    LocnDegFr    = (LocnPoly+1)*(LocnPoly+2)*(LocnPoly+3)/6
    LocnPolyMat  = DISC%Galerkin%nPolyMat
    LocnDegFrMat = (LocnPolyMat+1)*(LocnPolyMat+2)*(LocnPolyMat+3)/6
    LocnVar      = EQN%nVar
    !
    NULLIFY(TimeIntDof_ptr)
    NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr, EStar_Sp_ptr)
    NULLIFY(Tens_xi_Sp_ptr, Tens_eta_Sp_ptr, Tens_zeta_Sp_ptr, Tens_klm_sp_ptr)

    ! ==============================================================================
    ! DYNAMIC RUPTURE
    ! ==============================================================================
    EPIK_USER_START(r_dr)
    SCOREP_USER_REGION_BEGIN( r_dr, "dynamic_rupture", SCOREP_USER_REGION_TYPE_COMMON )
    ! Calculate Fault Output of the previous timestep
    ! Here, because the complete updated dgvar value of the (MPI-)Neighbor is needed
    ! Note that this causes a dt timeshift in the DR output routines
    EPIK_USER_START(r_dr_output)
    SCOREP_USER_REGION_BEGIN( r_dr_output, "fault_output", SCOREP_USER_REGION_TYPE_COMMON )
    CALL faultoutput(EQN, DISC, MESH, IO, MPI, MaterialVal, BND, time, dt)
    EPIK_USER_END(r_dr_output)
    SCOREP_USER_REGION_END( r_dr_output )

    IF (EQN%DR.EQ.1) THEN
       ! friction evaluation
       CALL Friction(EQN, DISC, MESH, MPI, IO, OptionalFields, BND, time, dt)
    ENDIF
    IF (DISC%DynRup%OutputPointType.EQ.4.OR.DISC%DynRup%OutputPointType.EQ.5) THEN
        ! close fault output
        IF ((min(DISC%EndTime,dt*DISC%MaxIteration)-time).LT.(dt*1.005d0)) THEN
            call closeFaultOutput()
        ENDIF
    ENDIF
    !
    EPIK_USER_END(r_dr)
    SCOREP_USER_REGION_END( r_dr )

    ! ==============================================================================
    ! Compute the volume integrals (stiffness terms) and face integrals (flux terms)
    ! with the time-integrated DOF
    ! ==============================================================================
    nElem = MESH%nElem

#if MEASURENODELEVELPERFORMANCE == 2
    ! initialize node level performance measurement
    call setNodeLevelDataStart( kernelType_volumeBoundary )
#endif

    ! start epik/scorep region for source terms and
    EPIK_USER_START(r_source)
    SCOREP_USER_REGION_BEGIN( r_source, "source_terms", SCOREP_USER_REGION_TYPE_COMMON )

#ifdef GENERATEDKERNELS

#else
    DO iElem = 1, nElem
        dudt(:,:) = 0.0d0

        TimeIntDof_ptr  => DISC%Galerkin%dgwork(:,:,iElem)                    ! Pointer to TimeIntDof(iElem)

        AStar_Sp_ptr    => DISC%Galerkin%AStar_Sp(iElem)                      ! Star matrices in sparse version
        BStar_Sp_ptr    => DISC%Galerkin%BStar_Sp(iElem)                      !
        CStar_Sp_ptr    => DISC%Galerkin%CStar_Sp(iElem)                      !
        EStar_Sp_ptr    => DISC%Galerkin%EStar_Sp(iElem)                      !

        Kxi_k_Sp_ptr   => DISC%Galerkin%Kxi_k_Tet_Sp                      ! Stiff tensors k for tetras
        Keta_k_Sp_ptr  => DISC%Galerkin%Keta_k_Tet_Sp
        Kzeta_k_Sp_ptr => DISC%Galerkin%Kzeta_k_Tet_Sp
        Kxi_m_Sp_ptr   => DISC%Galerkin%Kxi_m_Tet_Sp                      ! Stiff tensors m for tetras
        Keta_m_Sp_ptr  => DISC%Galerkin%Keta_m_Tet_Sp
        Kzeta_m_Sp_ptr => DISC%Galerkin%Kzeta_m_Tet_Sp

        MassMatrix_ptr   => DISC%Galerkin%MassMatrix_Tet(1:LocnDegFr,1:LocnDegFr,LocnPoly)
        JacobiDetVolum   =  6.0d0 * MESH%ELEM%Volume(iElem)

        ! Pre-compute combined JacobiDetVolum and mass matrix to avoid devision and large expression later
        JacobiVol_MassMatrix(:)  = 0.0D0
        DO iDegFr = 1, LocnDegFr
             JacobiVol_MassMatrix(iDegFr)  = 1.0D0 / ( JacobiDetVolum * MassMatrix_ptr(iDegFr,iDegFr) )
        ENDDO

        ! AStar
        CALL SPT_T3MUL_M(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),AStar_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat )
        CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kxi_k_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,AStar_Sp_ptr%nIndex3,AStar_Sp_ptr%Index3)
        CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kxi_m_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,AStar_Sp_ptr%nIndex3,AStar_Sp_ptr%Index3)
        ! BStar
        CALL SPT_T3MUL_M(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),BStar_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat )
        CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Keta_k_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,BStar_Sp_ptr%nIndex3,BStar_Sp_ptr%Index3)
        CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Keta_m_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,BStar_Sp_ptr%nIndex3,BStar_Sp_ptr%Index3)
        ! CStar
        CALL SPT_T3MUL_M(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),CStar_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat )
        CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kzeta_k_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,CStar_Sp_ptr%nIndex3,CStar_Sp_ptr%Index3)
        CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kzeta_m_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,CStar_Sp_ptr%nIndex3,CStar_Sp_ptr%Index3)
        !
        dudt = dudt * JacobiDetVolum
        !
        NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr)
        NULLIFY(Kxi_k_Sp_ptr,Keta_k_Sp_ptr,Kzeta_k_Sp_ptr)
        NULLIFY(Kxi_m_Sp_ptr,Keta_m_Sp_ptr,Kzeta_m_Sp_ptr)
        !================================================================================
        ! Compute the fluxes with the time-integrated DOF over each element boundary side
        !================================================================================
        DO iSide = 1,4
            !
            iNeighbor           = MESH%ELEM%SideNeighbor(iSide,iElem)
            iLocalNeighborSide  = MESH%ELEM%LocalNeighborSide(iSide,iElem)
            iLocalNeighborVrtx  = MESH%ELEM%LocalNeighborVrtx(iSide,iElem)
            !
            NeignPoly           = DISC%Galerkin%nPoly
            NeignDegFr          = (NeignPoly+1)*(NeignPoly+2)*(NeignPoly+3)/6
            NeignPolyMat        = DISC%Galerkin%nPolyMat
            NeignDegFrMat       = (NeignPolyMat+1)*(NeignPolyMat+2)*(NeignPolyMat+3)/6
            NeignVar            = EQN%nVar
            !
            ! Flux integral (element-side) dependent
            FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(0,1,iSide)
            JacobiDetSurfa = 2.0d0 * DISC%Galerkin%geoSurfaces(iSide,iElem)

            ! Flux solver
            FluxSolver_Sp_ptr => DISC%Galerkin%FLStar_Sp(iElem,iSide)
            !
            TimeIntDof_ptr    => DISC%Galerkin%DGwork(:,:,iElem)

            ! Flux contribution of the element itself
            auxvar(:,:,:)    = 0.0d0
            dudt_flux(:,:) = 0.0d0
            ! exception for the Dynamic Rupture case (3): all computations are done in friction
            IF (MESH%ELEM%Reference(iSide,iElem) .NE. 3) THEN
                CALL SPT_T3MUL_P(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),FluxSolver_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat)
                CALL SPL_T3MULB(dudt_flux,FluxInt_Sp_ptr,auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,FluxSolver_Sp_ptr%nIndex3,FluxSolver_Sp_ptr%Index3)
                dudt(1:LocnDegFr,1:EQN%nVar) = dudt(1:LocnDegFr,1:EQN%nVar) + dudt_flux * JacobiDetSurfa
            ENDIF

            NULLIFY(FluxInt_Sp_ptr)
            NULLIFY(FluxSolver_Sp_ptr)

            !====================================
            ! Flux contribution of the neighbor
            ! Boundary conditions, if necessary
            !====================================
            SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
            CASE(0,6) ! Standard case, not on a boundary
                !
                ! Flux integral (element-side) dependent
                FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(iLocalNeighborSide,iLocalNeighborVrtx,iSide)
                !
                ! Flux solver
                FluxSolver_Sp_ptr => DISC%Galerkin%FRStar_Sp(iElem,iSide)
                !
                !====================================
                ! If the neighbor element is an MPI
                ! 1.- construct an array with the MPI neighbor dof
                ! 2.- compute the time-integrated dof and assign them to NeigTimeIntDof
                !====================================
                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
                    iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
                    MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
                    ! debug information
                    !IF(MPIIndex.EQ.-1) THEN
                    !    WRITE(IO%UNIT%errOut,*) 'Severe Error in Galerkin3D. MPIIndex = -1 !', iElem, iSide
                    !    STOP
                    !ENDIF
                    !IF(.NOT.BND%ObjMPI(iObject)%Init) THEN
                    !    WRITE(IO%UNIT%errOut,*) ' SeisSol MPI Error. Domain not initialized. '
                    !    WRITE(IO%UNIT%errOut,*) iObject, MPIIndex
                    !    STOP
                    !ENDIF
                    !
                    ! Time-integrated DOF of the neighbor MPI element
                    TimeIntDof_ptr => BND%ObjMPI(iObject)%NeighborDOF(:,:,MPIIndex)
                    !
                    !
                ELSE
                    TimeIntDof_ptr  => DISC%Galerkin%DGwork(:,:,iNeighbor)
                    !
                    !
                ENDIF ! IF(MESH%ELEM%MPIReference(iSide,iElem).EQ.1)
                !
            CASE(1) ! Free surface
                !
                ! Flux integral (element-side) dependent
                FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(0,1,iSide)
                !
                ! Flux solver
                FluxSolver_Sp_ptr => DISC%Galerkin%FRStar_Sp(iElem,iSide)
                !
                TimeIntDof_ptr  => DISC%Galerkin%DGwork(:,:,iElem)

            CASE(3) ! Rupture Dynamics

               auxvar(:,:,:)    = 0.0d0
               dudt_flux(:,:)   = 0.0d0
               NULLIFY(FluxInt_Sp_ptr)
               NULLIFY(FluxSolver_Sp_ptr)
               ! all computations for the flux integration are done in friction for all variables
               CYCLE

            CASE(5) ! Absorbing boundary
                CYCLE
            END SELECT

            ! Flux contribution of the neighbor
            auxvar(:,:,:)     = 0.0d0
            dudt_flux(:,:)    = 0.0d0
            CALL SPT_T3MUL_P(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),FluxSolver_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat)
            CALL SPL_T3MULB(dudt_flux,FluxInt_Sp_ptr,auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,FluxSolver_Sp_ptr%nIndex3,FluxSolver_Sp_ptr%Index3)
            !
            dudt(1:LocnDegFr,1:EQN%nVar) = dudt(1:LocnDegFr,1:EQN%nVar) + dudt_flux * JacobiDetSurfa
            !
            NULLIFY(FluxInt_Sp_ptr)
            NULLIFY(FluxSolver_Sp_ptr)

        ENDDO ! iSide

        !==========================
        ! Source Terms
        !==========================
        ! GalerkinSource3d output in dudt_src the space and time integration of the corresponding source
        CALL GalerkinSource3D(dudt_src,iElem,time,dt,LocnDegFr,OptionalFields,EQN,DISC,MESH,SOURCE,IO)
        !
        dudt(1:LocnDegFr,1:EQN%nVar) = dudt(1:LocnDegFr,1:EQN%nVar) - dudt_src(1:LocnDegFr,1:EQN%nVar)
        !==========================
        ! Multiplication by inverse (diagonal) mass matrix
        !==========================
        nVarTotal = EQN%nVarTotal
        DO iDegFr = 1, LocnDegFr
            dudt(iDegFr,1:EQN%nVarTotal)  = - dudt(iDegFr,1:nVarTotal) * JacobiVol_MassMatrix(iDegFr)
        ENDDO
        !
        !==========================
        ! Update the values of the elastic variables for the next timestep
        !==========================
        DO iVar = 1, EQN%nVar
            DO iDegFr = 1, LocnDegFr
                DISC%Galerkin%dgvar(iDegFr,iVar,iElem,1) = DISC%Galerkin%dgvar(iDegFr,iVar,iElem,1) + dudt(iDegFr,iVar)
            ENDDO
        ENDDO
        !
        !Calculate the kinetic energy
        KineticEnergy_tmp = 0.0
        !average approach
        !KineticEnergy_tmp = EQN%Rho0*(DISC%Galerkin%dgvar(1,7, iElem,1)**2 + DISC%Galerkin%dgvar(1,8, iElem, 1)**2 + DISC%Galerkin%dgvar(1,9, iELem, 1)**2)
        !EQN%Energy(1,iElem) = 0.5*KineticEnergy_tmp*MESH%Elem%Volume(iElem)

        !using onyl the average of an element to calculate the energy
        !KineticEnergy_tmp = EQN%Rho0*(DISC%Galerkin%dgvar(1,7, iElem,1)**2 + DISC%Galerkin%dgvar(1,8, iElem, 1)**2 + DISC%Galerkin%dgvar(1,9, iELem, 1)**2)
        !using a dofs for calculating the energy
        DO iDegFr=1,LocnDegFr
           KineticEnergy_tmp = KineticEnergy_tmp + &
                               EQN%Rho0*MassMatrix_ptr(iDegFr,iDegFr)*(DISC%Galerkin%dgvar(iDegFr,7, iElem,1)**2 + &
                               DISC%Galerkin%dgvar(iDegFr,8, iElem, 1)**2 + DISC%Galerkin%dgvar(iDegFr,9, iELem, 1)**2)
        ENDDO

        EQN%Energy(1,iElem) = 0.5*KineticEnergy_tmp*6.0d0*MESH%Elem%Volume(iElem) !|J|=6*V  transformation from reference element

    ENDDO ! iElem
#endif



! ==============================================================================
! OFF-FAULT PLASTICITY
! ==============================================================================
        IF(EQN%Plasticity.EQ.1) THEN

           IF(EQN%PlastMethod .EQ. 0) THEN
              intGaussP     => DISC%Galerkin%intGaussP_Tet
              intGaussW     => DISC%Galerkin%intGaussW_Tet
           ENDIF
              DO iElem = 1, nElem !for every element
#ifndef GENERATEDKERNELS


                SELECT CASE(EQN%PlastMethod) !two different implementations
                  CASE(0) !high-order points implementation
                    CALL Plasticity_3D_high(DISC%Galerkin%dgvar(:,1:6,iElem,1), DISC%Galerkin%DOFStress(:,1:6,iElem), DISC%Galerkin%nDegFr, DISC%Galerkin%nDegFr, &
                                              EQN%BulkFriction, EQN%Tv, dt, EQN%mu, EQN%lambda, DISC%Galerkin%plasticParameters(1:3,iElem), &
                                              EQN%Energy(2:3,iElem), DISC%Galerkin%pstrain(1:7,iElem), intGaussP, intGaussW, &
                                              DISC, EQN%nVar, DISC%Galerkin%nIntGP)

                  CASE(2) !average approach with first dof
                    CALL Plasticity_3D_avg(DISC,DISC%Galerkin%dgvar(:,1:6,iElem,1), DISC%Galerkin%DOFStress(:,1:6,iElem), DISC%Galerkin%nDegFr, &
                                            DISC%Galerkin%nDegFr, EQN%BulkFriction, EQN%Tv, dt, EQN%mu, EQN%lambda,DISC%Galerkin%plasticParameters(1:3,iElem), &
                                            EQN%Energy(2:3,iElem), DISC%Galerkin%pstrain(1:7,iElem) )
                END SELECT


#endif
!for the GK version the plasticity call is moved to Interoperability.cpp
              ENDDO

        ENDIF
! ==============================================================================
    ! end epik/scorep region for source terms
    EPIK_USER_END(r_source)
    SCOREP_USER_REGION_END( r_source )

    NULLIFY(TimeIntDof_ptr)
    NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr, EStar_Sp_ptr)
    NULLIFY(Kxi_k_Sp_ptr, Keta_k_Sp_ptr, Kzeta_k_Sp_ptr)
    NULLIFY(Kxi_m_Sp_ptr, Keta_m_Sp_ptr, Kzeta_m_Sp_ptr)

    ! end epik function ADERGalerkin3D_GTS
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()

  END SUBROUTINE ADERGalerkin3D_GTS

#ifndef GENERATEDKERNELS
  !===========================================================================!
  !!                                                                         !!
  !!  ADERGalerkin3D LTS evolves the ADER-DG solution for one timestep       !!
  !!                                                                         !!
  !===========================================================================!

  SUBROUTINE ADERGalerkin3D_LTS(time,dt,iteration,MaterialVal,OptionalFields,EQN,MESH,DISC,IC,SOURCE,BND,MPI,IO)
    !-------------------------------------------------------------------------!
    USE CauchyKovalewski_mod
    USE Galerkin_source_mod
    USE SP_MATMUL_mod
    USE gauss_mod
#ifdef HDF
    USE receiver_hdf_mod
#else
    USE receiver_mod
#endif
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)               :: EQN                                     !
    TYPE(tDiscretization), target  :: DISC                                    !
    TYPE(tUnstructMesh)            :: MESH                                    !
    TYPE(tUnstructOptionalFields)  :: OptionalFields                          !
    TYPE(tBoundary)                :: BND                                     !
    TYPE(tMPI)                     :: MPI                                     !
    TYPE(tInitialCondition)        :: IC                                      !
    TYPE(tSource)                  :: SOURCE                                  !
    TYPE(tInputOutput)             :: IO                                      !
    REAL                           :: time, dt                                ! Calculation time and timestep
    REAL                           :: MaterialVal(MESH%nElem,               & !
                                                  EQN%nBackgroundVar        ) ! Material Values
    INTEGER                        :: iteration                               ! Current iteration number
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !

    TYPE(tSparseTensor3), POINTER  :: AStar_Sp_ptr => NULL()                  ! Pointer to tSparseTensor3
    TYPE(tSparseTensor3), POINTER  :: BStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: CStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: EStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: FluxSolver_Sp_ptr => NULL()             !

    TYPE(tSparseTensor3b), POINTER :: Tens_xi_Sp_ptr   => NULL()              ! Pointer to tSparseTensor3b
    TYPE(tSparseTensor3b), POINTER :: Tens_eta_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_zeta_Sp_ptr => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_klm_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_k_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_k_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_k_Sp_ptr   => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_m_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_m_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_m_Sp_ptr   => NULL()              !

    TYPE(tSparseTensor3b), POINTER :: FluxInt_Sp_ptr    => NULL()             !

    REAL, POINTER                  :: MassMatrix_ptr(:,:) => NULL()           ! Pointer to the corresponding mass matrix

    REAL, POINTER                  :: Dof_ptr(:,:)            => NULL()       ! Dof of element iElem at time tau=0
    REAL, POINTER                  :: TimeIntDof_ptr(:,:)     => NULL()       ! Time integrate dof
    REAL, POINTER                  :: TimeDerDof_ptr(:,:,:)   => NULL()       ! Time derivative dof
    !
    REAL, TARGET                   :: TimeIntDofElem(DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal          ) !
    REAL, TARGET                   :: TimeIntDofSide(DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal,         & !
                                                     MESH%nSideMax          ) !
    REAL, TARGET                   :: NeigTimeIntDofElem(                   & !
                                                     DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal          ) !
    REAL, TARGET                   :: NeigTimeIntDofSide(                   & !
                                                     DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal,         & !
                                                     MESH%nSideMax          ) !
    REAL, TARGET                   :: NeigTimeDerDof(                       & !
                                                     DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal,         & !
                                                     0:DISC%Galerkin%nPoly  ) !

    REAL                           :: dudt(DISC%Galerkin%nDegFr,            & !
                                           EQN%nVarTotal                    ) ! Auxilliary update dof
    REAL                           :: dudt_flux(DISC%Galerkin%nDegFr,       & !
                                                EQN%nVarTotal               ) ! Auxilliary update dof
    REAL                           :: dudt_src(DISC%Galerkin%nDegFr,        & !
                                               EQN%nVarTotal                ) ! Auxilliary update dof
    REAL                           :: auxvar(DISC%Galerkin%nDegFr,          & !
                                             EQN%nVarTotal,                 & !
                                             DISC%Galerkin%nDegFrMat        ) ! Auxilliary variable
    REAL                           :: JacobiDetVolum                          ! Jacobian of the mapping
    REAL                           :: JacobiDetSurfa                          ! Jacobian of the mapping

    REAL, SAVE                     :: PercentComplete
    REAL, SAVE                     :: NextPercentComplete
    REAL                           :: tmin(0:MESH%nSideMax)
    REAL                           :: tmax(0:MESH%nSideMax)
    REAL                           :: NeigTime(0:MESH%nSideMax)

    REAL                           :: CpuTime_ini                             ! Eval cpu time
    REAL                           :: CpuTime_end                             ! Eval cpu time

    INTEGER                        :: LocnVar                                 !
    INTEGER                        :: LocnPoly                                !
    INTEGER                        :: LocnDegFr                               !
    INTEGER                        :: LocnPolyMat                             !
    INTEGER                        :: LocnDegFrMat                            !
    INTEGER                        :: NeignVar                                !
    INTEGER                        :: NeignPoly                               !
    INTEGER                        :: NeignDegFr                              !
    INTEGER                        :: NeignPolyMat                            !
    INTEGER                        :: NeignDegFrMat                           !
    INTEGER                        :: ReactionTerm                            !

    INTEGER                        :: ncycles
    INTEGER                        :: SumElementUpdates
    INTEGER                        :: ElementsUpdated(MESH%nElem)
    INTEGER                        :: TotalUpdates(2)

    LOGICAL                        :: IterationCompleted
    LOGICAL                        :: UpdateElement(MESH%nElem)

    INTEGER                        :: iObject                                 !
    INTEGER                        :: MPIIndex                                !
    INTEGER                        :: iNeighbor                               !
    INTEGER                        :: iLocalNeighborSide                      !
    INTEGER                        :: iLocalNeighborVrtx                      !
    INTEGER                        :: counter, iPoly
    INTEGER                        :: iDegFr                                  ! Counter
    INTEGER                        :: iElem                                   ! Counter
    INTEGER                        :: iSide                                   ! Counter
    INTEGER                        :: iVar                                    ! Counter

    !-------------------------------------------------------------------------!

    IF(.NOT.DISC%Galerkin%init) THEN
        logError(*) 'ERROR stepQuadFreeGalerkin: SeisSol Interface not initialized!!'
        STOP
    ENDIF
    !
    dt   = MAXVAL(DISC%LocalDt(:))
    time = MAXVAL(DISC%LocalTime(:))
    !
    nCycles             = 0
    SumElementUpdates   = 0
    ElementsUpdated(:)  = 0
    IterationCompleted  = .FALSE.
    !
    ! At first timestep, compute approximate speedup
    ! and do the Cauchy-Kovalewski procedure for all elements
    IF(DISC%iterationstep.EQ.0) THEN
        NextPercentComplete = 0.
        TotalUpdates(1)     = 0
        TotalUpdates(2)     = 0
        tmin = MINVAL(DISC%LocalDt(:))
        DO iElem = 1, MESH%nElem
            TotalUpdates(1) = TotalUpdates(1) + CEILING( (DISC%EndTime-DISC%LocalTime(iElem))/DISC%LocalDt(iElem) )
            TotalUpdates(2) = TotalUpdates(2) + CEILING( (DISC%EndTime-DISC%LocalTime(iElem))/tmin(0))
        ENDDO
        DISC%TotalUpdates = TotalUpdates(1)
        PRINT *, ' Total Element updates with local  timestepping : ', TotalUpdates(1)
        PRINT *, ' Total Element updates with global timestepping : ', TotalUpdates(2)
        PRINT *, ' Estimated speedup : ', REAl(TotalUpdates(2))/REAL(TotalUpdates(1))
        !
    ENDIF ! DISC%iterationstep.EQ.0

    NULLIFY(TimeDerDof_ptr)
    NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr, EStar_Sp_ptr)
    NULLIFY(Kxi_k_Sp_ptr, Keta_k_Sp_ptr, Kzeta_k_Sp_ptr)
    NULLIFY(Kxi_m_Sp_ptr, Keta_m_Sp_ptr, Kzeta_m_Sp_ptr)

    DO WHILE( .NOT. IterationCompleted )
#ifdef HDF
        CALL receiver_hdf(EQN,MESH,DISC,MPI,IO)
#else
        CALL receiver(EQN,MESH,DISC,MPI,IO)
#endif
        ! Loop over all elements: Set update flag based on local timestepping
        DO iElem = 1, MESH%nElem
            !
            ! Control element update for local timestepping
            !
            ! Match printtime if necessary
            IF(DISC%LocalTime(iElem)+DISC%LocalDt(iElem).GE.DISC%PrintTime) THEN
                DISC%LocalDt(iElem) = DISC%printtime - DISC%LocalTime(iElem)
            ENDIF
            ! Match endtime if necessary
            IF(DISC%LocalTime(iElem)+DISC%LocalDt(iElem).GE.DISC%EndTime) THEN
                DISC%LocalDt(iElem) = DISC%EndTime - DISC%LocalTime(iElem)
            ENDIF
            !
            ! Find the minimum of t+dt over all neighboring elements
            tmin(0) = 1e20
            DO iSide = 1,4
                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
                    iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
                    MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
                    IF(MPIIndex.EQ.-1) THEN
                        PRINT *, 'Severe Error in Galerkin3D. MPIIndex = -1 !', iElem, iSide
                        STOP
                    ENDIF
                    NeigTime(0) = BND%ObjMPI(iObject)%NeighborTime(MPIIndex) + BND%ObjMPI(iObject)%NeighborDt(MPIIndex)
                    tmin(0) = MIN( tmin(0),NeigTime(0))
                ELSE
                    SELECT CASE (MESH%ELEM%Reference(iSide,iElem))
                    CASE(0,6)
                        iNeighbor   = MESH%ELEM%SideNeighbor(iSide,iElem)
                        NeigTime(0) = DISC%LocalTime(iNeighbor)+DISC%LocalDt(iNeighbor)
                        tmin(0) = MIN( tmin(0),NeigTime(0))
                    ENDSELECT
                ENDIF
            ENDDO ! iSide
            !
            ! Regular update criterion: Update unless some neighbor has to update first
            ! ========================
            !
            IF( (DISC%LocalTime(iElem)+DISC%LocalDt(iElem)).LE.tmin(0) ) THEN
                UpdateElement(iElem) = .TRUE.
            ELSE
                UpdateElement(iElem) = .FALSE.
            ENDIF
            !
            ! Exception: If printtime or endtime have been reached, then do not update the element
            !
            IF(DISC%LocalTime(iElem).GE.DISC%PrintTime) THEN
                UpdateElement(iElem)   = .FALSE.
                ElementsUpdated(iElem) = 1
            ENDIF
            IF(DISC%IterationCriterion.EQ.2) THEN
                ! If outer synchronization time is reached, do not update any more.
                IF(DISC%LocalTime(iElem).GE.time+dt) THEN
                    UpdateElement(iElem)   = .FALSE.
                    ElementsUpdated(iElem) = 1
                ENDIF
            ENDIF
            ! If endtime is reached, do not update any more.
            IF(DISC%LocalTime(iElem).GE.DISC%EndTime) THEN
                UpdateElement(iElem)   = .FALSE.
                ElementsUpdated(iElem) = 1
            ENDIF
            !
        ENDDO ! iElem
        !
        ! Loop over all elements: Update elements as necessary
        !
        DO iElem = 1, MESH%nElem
            !
            IF(.NOT.UpdateElement(iElem)) THEN
                CYCLE
            ENDIF
            !
            IF(DISC%Galerkin%pAdaptivity.GE.1) THEN
                LocnPoly  = DISC%Galerkin%LocPoly(iElem)
            ELSE
                LocnPoly  = DISC%Galerkin%nPoly
            ENDIF
            LocnDegFr    = (LocnPoly+1)*(LocnPoly+2)*(LocnPoly+3)/6
            LocnPolyMat  = DISC%Galerkin%nPolyMat
            LocnDegFrMat = (LocnPolyMat+1)*(LocnPolyMat+2)*(LocnPolyMat+3)/6
            LocnVar      = EQN%nVar
            !
            TimeDerDof_ptr  => DISC%Galerkin%DGTaylor(1:LocnDegFr,1:LocnVar,0:LocnPoly,iElem) ! Pointer to TimeDerDof(iElem)
            !
            tmin(0) = DISC%LocalTime(iElem)
            tmax(0) = DISC%LocalTime(iElem)+DISC%LocalDt(iElem)
            !
            CALL TimeIntTaylor( TimeIntDof = TimeIntDofElem,                            &
                                TimeDerDof = TimeDerDof_ptr,                            &
                                t0         = tmin(0),                                   &
                                t1         = tmin(0),                                   &
                                t2         = tmax(0),                                   &
                                LocnPoly   = LocnPoly,                                  &
                                LocnVar    = LocnVar,                                   &
                                LocnDegFr  = LocnDegFr                                  )
            !
            dudt = 0.0d0
            !
            AStar_Sp_ptr    => DISC%Galerkin%AStar_Sp(iElem)                      ! Star matrices in sparse version
            BStar_Sp_ptr    => DISC%Galerkin%BStar_Sp(iElem)                      !
            CStar_Sp_ptr    => DISC%Galerkin%CStar_Sp(iElem)                      !
            EStar_Sp_ptr    => DISC%Galerkin%EStar_Sp(iElem)                      !
            !
            Kxi_k_Sp_ptr   => DISC%Galerkin%Kxi_k_Tet_Sp                      ! Stiff tensors k for tetras
            Keta_k_Sp_ptr  => DISC%Galerkin%Keta_k_Tet_Sp
            Kzeta_k_Sp_ptr => DISC%Galerkin%Kzeta_k_Tet_Sp
            Kxi_m_Sp_ptr   => DISC%Galerkin%Kxi_m_Tet_Sp                      ! Stiff tensors m for tetras
            Keta_m_Sp_ptr  => DISC%Galerkin%Keta_m_Tet_Sp
            Kzeta_m_Sp_ptr => DISC%Galerkin%Kzeta_m_Tet_Sp
            !
            MassMatrix_ptr   => DISC%Galerkin%MassMatrix_Tet(1:LocnDegFr,1:LocnDegFr,LocnPoly)
            JacobiDetVolum   =  6.0d0 * MESH%ELEM%Volume(iElem)
            !
            ! AStar
            CALL SPT_T3MUL_M(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),AStar_Sp_ptr,TimeIntDofElem,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat )
            CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kxi_k_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,AStar_Sp_ptr%nIndex3,AStar_Sp_ptr%Index3)
            CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kxi_m_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,AStar_Sp_ptr%nIndex3,AStar_Sp_ptr%Index3)
            ! BStar
            CALL SPT_T3MUL_M(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),BStar_Sp_ptr,TimeIntDofElem,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat )
            CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Keta_k_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,BStar_Sp_ptr%nIndex3,BStar_Sp_ptr%Index3)
            CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Keta_m_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,BStar_Sp_ptr%nIndex3,BStar_Sp_ptr%Index3)
            ! CStar
            CALL SPT_T3MUL_M(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),CStar_Sp_ptr,TimeIntDofElem,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat )
            CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kzeta_k_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,CStar_Sp_ptr%nIndex3,CStar_Sp_ptr%Index3)
            CALL SPL_T3MULB(dudt(1:LocnDegFr,1:LocnVar),Kzeta_m_Sp_ptr, auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,CStar_Sp_ptr%nIndex3,CStar_Sp_ptr%Index3)
            !
            dudt = dudt * JacobiDetVolum
            !
            NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr)
            NULLIFY(Kxi_k_Sp_ptr,Keta_k_Sp_ptr,Kzeta_k_Sp_ptr)
            NULLIFY(Kxi_m_Sp_ptr,Keta_m_Sp_ptr,Kzeta_m_Sp_ptr)
            !
            !================================================================================
            ! Compute the fluxes with the time-integrated DOF over each element boundary side
            !================================================================================
            DO iSide = 1,4
                !
                iNeighbor           = MESH%ELEM%SideNeighbor(iSide,iElem)
                iLocalNeighborSide  = MESH%ELEM%LocalNeighborSide(iSide,iElem)
                iLocalNeighborVrtx  = MESH%ELEM%LocalNeighborVrtx(iSide,iElem)
                !
                NeignPoly           = DISC%Galerkin%nPoly
                NeignDegFr          = (NeignPoly+1)*(NeignPoly+2)*(NeignPoly+3)/6
                NeignPolyMat        = DISC%Galerkin%nPolyMat
                NeignDegFrMat       = (NeignPolyMat+1)*(NeignPolyMat+2)*(NeignPolyMat+3)/6
                NeignVar            = EQN%nVar
                !
                TimeDerDof_ptr  => DISC%Galerkin%DGTaylor(1:LocnDegFr,1:LocnVar,0:LocnPoly,iElem) ! Pointer to TimeDerDof(iElem)
                !
                !====================================
                ! If the neighbor element is an MPI
                ! 1.- construct an array with the MPI neighbor dof
                ! 2.- compute the time-integrated dof and assign them to NeigTimeIntDof
                !====================================
                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1)THEN
                    SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
                    CASE(0,6) ! Internal or periodic boundary
                        ! The neighbor element belongs to a different MPI domain
                        iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
                        MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
                        IF(MPIIndex.EQ.-1) THEN
                            logError(*) 'Severe Error in Galerkin3D. MPIIndex = -1 !', iElem, iSide
                            STOP
                        ENDIF
                        IF(.NOT.BND%ObjMPI(iObject)%Init) THEN
                            logError(*) 'SeisSol MPI Error. Domain not initialized. '
                            logError(*) iObject, MPIIndex
                            STOP
                        ENDIF
                        !
                        tmin(iSide) = MAX( DISC%LocalTime(iElem),                                                                 &
                                           BND%ObjMPI(iObject)%NeighborTime(MPIIndex)                                             )
                        tmax(iSide) = MIN( DISC%LocalTime(iElem) + DISC%LocalDt(iElem),                                           &
                                           BND%ObjMPI(iObject)%NeighborTime(MPIIndex) + BND%ObjMPI(iObject)%NeighborDt(MPIIndex)  )
                        IF(ABS(tmax(iSide)-tmin(iSide)).LE.(1e-10*dt)) THEN
                            CYCLE
                        ENDIF
                        !
                        NeigTime(iSide) = BND%ObjMPI(iObject)%NeighborTime(MPIIndex)
                        !
                        ! Get Time-Taylor series from MPI neighbor element
                        !
                        DO iVar = 1, EQN%nVarTotal
                         counter = 0
                         DO iDegFr = 1, DISC%Galerkin%nDegFr
                          DO iPoly = 0, DISC%Galerkin%nPoly
                            counter = counter + 1
                            NeigTimeDerDof(iDegFr,iVar,iPoly) = BND%ObjMPI(iObject)%NeighborDOF(counter,iVar,MPIIndex)
                          ENDDO
                         ENDDO
                        ENDDO
                        !
                    END SELECT ! MESH%ELEM%Reference(iSide,iElem)
                ELSE
                    SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
                    CASE(0,6)
                        tmin(iSide) = MAX( DISC%LocalTime(iElem), DISC%LocalTime(iNeighbor)    )
                        tmax(iSide) = MIN( DISC%LocalTime(iElem) + DISC%LocalDt(iElem),        &
                                           DISC%LocalTime(iNeighbor) + DISC%LocalDt(iNeighbor) )
                        IF(ABS(tmax(iSide)-tmin(iSide)).LE.(1e-10*dt)) THEN
                            CYCLE
                        ENDIF
                        NeigTimeDerDof(1:NeignDegFr,1:NeignVar,0:NeignPoly) = DISC%Galerkin%DGTaylor(1:NeignDegFr,1:NeignVar,0:LocnPoly,iNeighbor)
                        NeigTime(iSide) = DISC%LocalTime(iNeighbor)
                    CASE DEFAULT
                        tmin(iSide) = DISC%LocalTime(iElem)
                        tmax(iSide) = DISC%LocalTime(iElem) + DISC%LocalDt(iElem)
                    END SELECT
                ENDIF ! MESH%ELEM%MPIReference(iSide,iElem).EQ.1

                CALL TimeIntTaylor( TimeIntDof = TimeIntDofSide(1:LocnDegFr,                  &
                                                                1:LocnVar,                    &
                                                                iSide),                       &
                                    TimeDerDof = TimeDerDof_ptr,                              &
                                    t0         = DISC%LocalTime(iElem),                       &
                                    t1         = tmin(iSide),                                 &
                                    t2         = tmax(iSide),                                 &
                                    LocnPoly   = LocnPoly,                                    &
                                    LocnVar    = LocnVar,                                     &
                                    LocnDegFr  = LocnDegFr                                    )

                SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
                CASE(0,6)
                    CALL TimeIntTaylor( TimeIntDof = NeigTimeIntDofSide(1:NeignDegFr,         &
                                                                        1:NeignVar,           &
                                                                        iSide),               &
                                        TimeDerDof = NeigTimeDerDof,                          &
                                        t0         = NeigTime(iSide),                         &
                                        t1         = tmin(iSide),                             &
                                        t2         = tmax(iSide),                             &
                                        LocnPoly   = NeignPoly,                               &
                                        LocnVar    = NeignVar,                                &
                                        LocnDegFr  = NeignDegFr                               )
                ENDSELECT

                ! Flux integral (element-side) dependent
                FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(0,1,iSide)
                JacobiDetSurfa = 2.0d0 * DISC%Galerkin%geoSurfaces(iSide,iElem)
                ! Flux solver
                FluxSolver_Sp_ptr => DISC%Galerkin%FLStar_Sp(iElem,iSide)
                !
                TimeIntDof_ptr    => TimeIntDofSide(1:LocnDegFr,1:LocnVar,iSide)
                !
                ! Flux contribution of the element itself
                auxvar    = 0.0d0
                dudt_flux = 0.0d0
                CALL SPT_T3MUL_P(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),FluxSolver_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat)
                CALL SPL_T3MULB(dudt_flux,FluxInt_Sp_ptr,auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,FluxSolver_Sp_ptr%nIndex3,FluxSolver_Sp_ptr%Index3)
                dudt = dudt + dudt_flux * JacobiDetSurfa

                NULLIFY(FluxInt_Sp_ptr)
                NULLIFY(FluxSolver_Sp_ptr)

                !====================================
                ! Flux contribution of the neighbor
                ! Boundary conditions, if necessary
                !====================================
                SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
                CASE(0,6) ! Standard case, not on a boundary
                    !
                    ! Flux integral (element-side) dependent
                    FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(iLocalNeighborSide,iLocalNeighborVrtx,iSide)
                    !
                    ! Flux solver
                    FluxSolver_Sp_ptr => DISC%Galerkin%FRStar_Sp(iElem,iSide)
                    !
                    TimeIntDof_ptr    => NeigTimeIntDofSide(1:NeignDegFr,1:NeignVar,iSide)
                    !
                CASE(1) ! Free surface
                    !
                    ! Flux integral (element-side) dependent
                    FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(0,1,iSide)
                    !
                    ! Flux solver
                    FluxSolver_Sp_ptr => DISC%Galerkin%FRStar_Sp(iElem,iSide)
                    !
                    TimeIntDof_ptr    => TimeIntDofSide(1:LocnDegFr,1:LocnVar,iSide)
                CASE(5) ! Absorbing boundary
                    CYCLE
                END SELECT

                ! Flux contribution of the element itself
                auxvar    = 0.0d0
                dudt_flux = 0.0d0
                CALL SPT_T3MUL_P(auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),FluxSolver_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat)
                CALL SPL_T3MULB(dudt_flux,FluxInt_Sp_ptr,auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,FluxSolver_Sp_ptr%nIndex3,FluxSolver_Sp_ptr%Index3)
                !
                dudt = dudt + dudt_flux * JacobiDetSurfa
                !
                NULLIFY(FluxInt_Sp_ptr)
                NULLIFY(FluxSolver_Sp_ptr)

            ENDDO ! iSide
            !
            !==========================
            ! Source Terms
            !==========================
            !
            ! Determine Gaussian integration points in time
            CALL gausslegendre(0.,DISC%LocalDt(iElem),DISC%Galerkin%TimeGaussP,DISC%Galerkin%TimeGaussW)
            !
            ! GalerkinSource3d output in dudt_src the space and time integration of the corresponding source
            CALL GalerkinSource3D(dudt_src,iElem,DISC%LocalTime(iElem),DISC%LocalDt(iElem),LocnDegFr,OptionalFields,EQN,DISC,MESH,SOURCE,IO)
            !
            dudt = dudt - dudt_src
            !
            ! Add the update contribution calculated by other elements
            ! ========================================================
            !
            dudt = dudt + DISC%Galerkin%dgwork(1:LocnDegFr,1:LocnVar,iElem)
            !
            !==========================
            ! Multiplication by inverse (diagonal) mass matrix
            !==========================
            DO iDegFr = 1, LocnDegFr
                dudt(iDegFr,1:EQN%nVarTotal)  = - dudt(iDegFr,1:EQN%nVarTotal) / ( JacobiDetVolum * MassMatrix_ptr(iDegFr,iDegFr) )
            ENDDO
            !
            !==========================
            ! Update the values of the elastic variables for the next timestep
            !==========================
            DO iVar = 1, EQN%nVarTotal
                DO iDegFr = 1, LocnDegFr
                    DISC%Galerkin%dgvar(iDegFr,iVar,iElem,1) = DISC%Galerkin%dgvar(iDegFr,iVar,iElem,1) + dudt(iDegFr,iVar)
                ENDDO
            ENDDO
            !
            ! Reset contributions after having updated the element
            !
            DISC%Galerkin%dgwork(:,:,iElem)  = 0.
            DISC%LocalTime(iElem)            = DISC%LocalTime(iElem) + DISC%LocalDt(iElem)
            DISC%LocalIteration(iElem)       = DISC%LocalIteration(iElem) + 1
            ElementsUpdated(iElem)           = 1
            SumElementUpdates                = SumElementUpdates + 1
            !
            ! Do the Cauchy-Kovalewski procedure to provide new time derivatives for the Taylor series
            !
            Dof_ptr         => DISC%Galerkin%dgvar(:,:,iElem,1)                               ! Pointer to Dof(iElem)
            TimeDerDof_ptr  => DISC%Galerkin%DGTaylor(1:LocnDegFr,1:LocnVar,0:LocnPoly,iElem) ! Pointer to TimeDerDof(iElem)
            !
            AStar_Sp_ptr    => DISC%Galerkin%AStar_Sp(iElem)                      ! Star matrices in sparse version
            BStar_Sp_ptr    => DISC%Galerkin%BStar_Sp(iElem)                      !
            CStar_Sp_ptr    => DISC%Galerkin%CStar_Sp(iElem)                      !
            EStar_Sp_ptr    => DISC%Galerkin%EStar_Sp(iElem)                      !
            !
            Tens_xi_Sp_ptr   => DISC%Galerkin%ADGxi_Tet_Sp                    ! Stiff CK tensors for tetras
            Tens_eta_Sp_ptr  => DISC%Galerkin%ADGeta_Tet_Sp
            Tens_zeta_Sp_ptr => DISC%Galerkin%ADGzeta_Tet_Sp
            Tens_klm_sp_ptr  => DISC%Galerkin%ADGklm_Tet_Sp                   ! Tensor for space dependent reaction term

            CALL CauchyKovalewski3D(                                            & !
                                   TimeDerDof    = TimeDerDof_ptr,              & ! Output
                                   Dof           = Dof_ptr,                     & ! Input
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
                                   LocDegFr      = LocnDegFr,                   & ! Input
                                   LocDegFrMat   = LocnDegFrMat,                & ! Input
                                   LocPoly       = LocnPoly,                    & ! Input
                                   nVar          = LocnVar                      ) ! Input
            !
            NULLIFY(TimeDerDof_ptr)
            NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr, EStar_Sp_ptr)
            NULLIFY(Kxi_k_Sp_ptr, Keta_k_Sp_ptr, Kzeta_k_Sp_ptr)
            NULLIFY(Kxi_m_Sp_ptr, Keta_m_Sp_ptr, Kzeta_m_Sp_ptr)
            !
            !==========================
            ! Now send the corresponding flux contributions to the neighbors
            ! Compute the fluxes with the time-integrated DOF over each element boundary side
            !==========================
            DO iSide = 1,4

                IF(ABS(tmax(iSide)-tmin(iSide)).LE.(1e-10*dt)) THEN
                    CYCLE
                ENDIF

                dudt = 0.0d0

                iLocalNeighborSide  = MESH%ELEM%LocalNeighborSide(iSide,iElem)
                iLocalNeighborVrtx  = MESH%ELEM%LocalNeighborVrtx(iSide,iElem)

                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
                    iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
                    MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
                    ! Flux solver
                    FluxSolver_Sp_ptr => BND%ObjMPI(iObject)%FLStar_Sp(MPIIndex,iLocalNeighborSide)
                ELSE
                    SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
                    CASE(0,6)
                        iNeighbor = MESH%ELEM%SideNeighbor(iSide,iElem)
                        ! Flux solver
                        FluxSolver_Sp_ptr => DISC%Galerkin%FLStar_Sp(iNeighbor,iLocalNeighborSide)
                    CASE DEFAULT
                        CYCLE
                    END SELECT
                    !
                ENDIF

                ! Flux integral (element-side) dependent
                FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(0,1,iLocalNeighborSide)
                JacobiDetSurfa = 2.0d0 * DISC%Galerkin%geoSurfaces(iSide,iElem)
                !
                ! Time integrated dof of neighbor
                TimeIntDof_ptr    => NeigTimeIntDofSide(1:NeignDegFr,1:NeignVar,iSide)
                !
                ! Flux contribution of the element itself
                auxvar    = 0.0d0
                dudt_flux = 0.0d0
                CALL SPT_T3MUL_P(auxvar(1:NeignDegFr,1:NeignVar,1:NeignDegFrMat),FluxSolver_Sp_ptr,TimeIntDof_ptr,NeignVar,NeignVar,NeignDegFr,NeignDegFrMat)
                CALL SPL_T3MULB(dudt_flux,FluxInt_Sp_ptr,auxvar(1:NeignDegFr,1:NeignVar,1:NeignDegFrMat),NeignVar,NeignDegFr,NeignDegFr,NeignDegFrMat,FluxSolver_Sp_ptr%nIndex3,FluxSolver_Sp_ptr%Index3)
                dudt = dudt + dudt_flux * JacobiDetSurfa

                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
                    ! Flux solver
                    FluxSolver_Sp_ptr => BND%ObjMPI(iObject)%FRStar_Sp(MPIIndex,iLocalNeighborSide)
                ELSE
                    ! Flux solver
                    FluxSolver_Sp_ptr => DISC%Galerkin%FRStar_Sp(iNeighbor,iLocalNeighborSide)
                    !
                ENDIF

                ! Flux integral (element-side) dependent
                FluxInt_Sp_ptr => DISC%Galerkin%FluxInt_Tet_Sp(iSide,iLocalNeighborVrtx,iLocalNeighborSide)
                JacobiDetSurfa = 2.0d0 * DISC%Galerkin%geoSurfaces(iSide,iElem)
                !
                ! Time integrated dof of neighbor
                TimeIntDof_ptr    => TimeIntDofSide(1:LocnDegFr,1:LocnVar,iSide)
                !
                ! Flux contribution of the element itself
                auxvar    = 0.0d0
                dudt_flux = 0.0d0
                CALL SPT_T3MUL_P(auxvar(1:LocnDegFr,1:NeignVar,1:LocnDegFrMat),FluxSolver_Sp_ptr,TimeIntDof_ptr,LocnVar,LocnVar,LocnDegFr,LocnDegFrMat)
                CALL SPL_T3MULB(dudt_flux,FluxInt_Sp_ptr,auxvar(1:LocnDegFr,1:LocnVar,1:LocnDegFrMat),LocnVar,LocnDegFr,LocnDegFr,LocnDegFrMat,FluxSolver_Sp_ptr%nIndex3,FluxSolver_Sp_ptr%Index3)
                dudt = dudt + dudt_flux * JacobiDetSurfa

                !
                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
                    iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
                    MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
                    IF(MPIIndex.EQ.-1) THEN
                        PRINT *, 'Severe Error in Galerkin3D. MPIIndex = -1 !', iElem, iSide
                        STOP
                    ENDIF
                    SELECT CASE( BND%ObjMPI(iObject)%NeighborUpdate(MPIIndex) )
                    CASE(0)
                        BND%ObjMPI(iObject)%NeighborDuDt(:,:,MPIIndex) = BND%ObjMPI(iObject)%NeighborDuDt(:,:,MPIIndex) + dudt(:,:)
                    CASE(1)
                        BND%ObjMPI(iObject)%NeighborDuDt(:,:,MPIIndex) = BND%ObjMPI(iObject)%NeighborDuDt(:,:,MPIIndex) + dudt(:,:)
                    CASE DEFAULT
                        PRINT *, 'Severe Error in Galerkin3D. MPINeighborUpdate must be 0 or 1 but is ', BND%ObjMPI(iObject)%NeighborUpdate(MPIIndex)
                        STOP
                    END SELECT
                ELSE
                    DISC%Galerkin%dgwork(1:NeignDegFr,1:NeignVar,iNeighbor) = DISC%Galerkin%dgwork(1:NeignDegFr,1:NeignVar,iNeighbor) + dudt(1:NeignDegFr,1:NeignVar)
                ENDIF
                !
                CONTINUE
                !
            ENDDO ! iSide
            !
            CONTINUE
            !
        ENDDO ! iElem
        !
        nCycles = nCycles + 1
        !
        SELECT CASE(DISC%IterationCriterion)
        CASE(1)
            IterationCompleted = .TRUE.
        CASE(2)
            IF( MINVAL(ElementsUpdated(:)).EQ.1 ) THEN
                IterationCompleted = .TRUE.
            ENDIF
        CASE(3)
            IterationCompleted = .TRUE.
        END SELECT
        !
    ENDDO ! WHILE( .NOT. IterationCompleted )

    IF(MPI%myrank.EQ.0) THEN
        tmin(0) = MINVAL(DISC%LocalDt(:))
        tmax(0) = MAXVAL(DISC%LocalDt(:))
        IF(DISC%TotalUpdates.GT.20) THEN
            PercentComplete = INT(100*SUM( DISC%LocalIteration(:) )/DISC%TotalUpdates)
            IF( PercentComplete.GE.NextPercentComplete ) THEN
                PRINT *, '| ', iteration+1, ' iterations done. nCycles = ', nCycles
                PRINT *, '| Min / max dt          : ', tmin(0),tmax(0)
                IF(tmin(0).GT.0.) THEN
                    PRINT *, '| Max / min dt ratio    : ', tmax(0)/tmin(0)
                ENDIF
                PRINT *, '| Min / max t           : ', MINVAL(DISC%LocalTime(:)),MAXVAL(DISC%LocalTime(:))
                PRINT *, '| Min / max t elements  : ', MINLOC(DISC%LocalTime(:)),MAXLOC(DISC%LocalTime(:))
                PRINT *, '| Element updates       : ', SumElementUpdates, ' out of ', MESH%nElem
                PRINT *, '| Total element updates : ', SUM( DISC%LocalIteration(:) ), ' ( ', PercentComplete, ' % ) '
                PRINT *, '| ----- '
                NextPercentComplete = NextPercentComplete + 5
            ENDIF
        ENDIF
    ENDIF

  END SUBROUTINE ADERGalerkin3D_LTS


  SUBROUTINE CKProcedureForEveryoneGTSUni(time,dt,iteration,MaterialVal,OptionalFields,EQN,MESH,DISC,IC,SOURCE,BND,MPI,IO)
    !-------------------------------------------------------------------------!
    USE CauchyKovalewski_mod
    USE Galerkin_source_mod
    USE SP_MATMUL_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)               :: EQN                                     !
    TYPE(tDiscretization), target  :: DISC                                    !
    TYPE(tUnstructMesh)            :: MESH                                    !
    TYPE(tUnstructOptionalFields)  :: OptionalFields                          !
    TYPE(tBoundary)                :: BND                                     !
    TYPE(tMPI)                     :: MPI                                     !
    TYPE(tInitialCondition)        :: IC                                      !
    TYPE(tSource)                  :: SOURCE                                  !
    TYPE(tInputOutput)             :: IO                                      !
    REAL                           :: time, dt                                ! Calculation time and timestep
    REAL                           :: MaterialVal(MESH%nElem,               & !
                                                  EQN%nBackgroundVar        ) ! Material Values
    INTEGER                        :: iteration                               ! Current iteration number
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !

    TYPE(tSparseTensor3), POINTER  :: AStar_Sp_ptr => NULL()                  ! Pointer to tSparseTensor3
    TYPE(tSparseTensor3), POINTER  :: BStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: CStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: EStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: FluxSolver_Sp_ptr => NULL()             !

    TYPE(tSparseTensor3b), POINTER :: Tens_xi_Sp_ptr   => NULL()              ! Pointer to tSparseTensor3b
    TYPE(tSparseTensor3b), POINTER :: Tens_eta_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_zeta_Sp_ptr => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_klm_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_k_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_k_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_k_Sp_ptr   => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_m_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_m_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_m_Sp_ptr   => NULL()              !

    TYPE(tSparseTensor3b), POINTER :: FluxInt_Sp_ptr    => NULL()             !

    REAL, POINTER                  :: MassMatrix_ptr(:,:) => NULL()           ! Pointer to the corresponding mass matrix

    REAL, POINTER                  :: Dof_ptr(:,:)            => NULL()       ! Dof of element iElem at time tau=0
    REAL, POINTER                  :: TimeIntDof_ptr(:,:)     => NULL()       ! Time integrated dof
    REAL, POINTER                  :: TimeDerDof_ptr(:,:,:)   => NULL()       ! Space-time polynomial pointer

    REAL, TARGET                   :: NeigTimeIntDof(DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal          ) !

    REAL                           :: dudt(DISC%Galerkin%nDegFr,            & !
                                           EQN%nVarTotal                    ) ! Auxilliary update dof
    REAL                           :: dudt_flux(DISC%Galerkin%nDegFr,       & !
                                                EQN%nVarTotal               ) ! Auxilliary update dof
    REAL                           :: dudt_src(DISC%Galerkin%nDegFr,        & !
                                               EQN%nVarTotal                ) ! Auxilliary update dof
    REAL                           :: auxvar(DISC%Galerkin%nDegFr,          & !
                                             EQN%nVarTotal,                 & !
                                             DISC%Galerkin%nDegFrMat        ) ! Auxilliary variable
    REAL                           :: JacobiDetVolum                          ! Jacobian of the mapping
    REAL                           :: JacobiDetSurfa                          ! Jacobian of the mapping

    REAL                           :: CpuTime_ini                             ! Eval cpu time
    REAL                           :: CpuTime_end                             ! Eval cpu time

    INTEGER                        :: LocnVar                                 !
    INTEGER                        :: LocnPoly                                !
    INTEGER                        :: LocnDegFr                               !
    INTEGER                        :: LocnPolyMat                             !
    INTEGER                        :: LocnDegFrMat                            !
    INTEGER                        :: NeignVar                                !
    INTEGER                        :: NeignPoly                               !
    INTEGER                        :: NeignDegFr                              !
    INTEGER                        :: NeignPolyMat                            !
    INTEGER                        :: NeignDegFrMat                           !
    INTEGER                        :: ReactionTerm                            !
    INTEGER                        :: iObject                                 !
    INTEGER                        :: MPIIndex                                !
    INTEGER                        :: iNeighbor                               !
    INTEGER                        :: iLocalNeighborSide                      !
    INTEGER                        :: iLocalNeighborVrtx                      !
    INTEGER                        :: iDegFr                                  ! Counter
    INTEGER                        :: iElem                                   ! Counter
    INTEGER                        :: iSide                                   ! Counter
    INTEGER                        :: iVar                                    ! Counter
    INTEGER                        :: nElem

    !-------------------------------------------------------------------------!

    !
    LocnPoly     = DISC%Galerkin%nPoly
    LocnDegFr    = (LocnPoly+1)*(LocnPoly+2)*(LocnPoly+3)/6
    LocnPolyMat  = DISC%Galerkin%nPolyMat
    LocnDegFrMat = (LocnPolyMat+1)*(LocnPolyMat+2)*(LocnPolyMat+3)/6
    LocnVar      = EQN%nVarTotal
    !
    ! =================================================================
    ! Perform the ADER time-integration of the degrees of freedom (DOF)
    ! =================================================================

    nElem = MESH%nElem

    DO iElem = 1,nElem
#ifdef GENERATEDKERNELS
        ! we are out of here: no more time integration in galerkin3d
#else

        ReactionTerm = 0

        Dof_ptr         => DISC%Galerkin%dgvar(:,:,iElem,1)                   ! Pointer to Dof(iElem)
        TimeIntDof_ptr  => DISC%Galerkin%dgwork(:,:,iElem)                    ! Pointer to TimeIntDof(iElem)

        AStar_Sp_ptr    => DISC%Galerkin%AStar_Sp(iElem)                      ! Star matrices in sparse version
        BStar_Sp_ptr    => DISC%Galerkin%BStar_Sp(iElem)                      !
        CStar_Sp_ptr    => DISC%Galerkin%CStar_Sp(iElem)                      !
        EStar_Sp_ptr    => DISC%Galerkin%EStar_Sp(iElem)                      !

        Tens_xi_Sp_ptr   => DISC%Galerkin%ADGxi_Tet_Sp                    ! Stiff CK tensors for tetras
        Tens_eta_Sp_ptr  => DISC%Galerkin%ADGeta_Tet_Sp
        Tens_zeta_Sp_ptr => DISC%Galerkin%ADGzeta_Tet_Sp
        Tens_klm_sp_ptr  => DISC%Galerkin%ADGklm_Tet_Sp                   ! Tensor for space dependent reaction term

        CALL CauchyKovalewski3D(                                            & !
                               TimeIntDOF    = TimeIntDOF_ptr,              & ! Output
                               Dof           = Dof_ptr,                     & ! Input
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
                               LocDegFr      = LocnDegFr,                   & ! Input
                               LocDegFrMat   = LocnDegFrMat,                & ! Input
                               LocPoly       = LocnPoly,                    & ! Input
                               nVar          = LocnVar                      ) ! Input
#endif
    ENDDO ! iElem

  END SUBROUTINE CKProcedureForEveryoneGTSUni


  SUBROUTINE CKProcedureForEveryoneLTSUni(time,dt,iteration,MaterialVal,OptionalFields,EQN,MESH,DISC,IC,SOURCE,BND,MPI,IO)
    !-------------------------------------------------------------------------!
    USE CauchyKovalewski_mod
    USE Galerkin_source_mod
    USE SP_MATMUL_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)               :: EQN                                     !
    TYPE(tDiscretization), target  :: DISC                                    !
    TYPE(tUnstructMesh)            :: MESH                                    !
    TYPE(tUnstructOptionalFields)  :: OptionalFields                          !
    TYPE(tBoundary)                :: BND                                     !
    TYPE(tMPI)                     :: MPI                                     !
    TYPE(tInitialCondition)        :: IC                                      !
    TYPE(tSource)                  :: SOURCE                                  !
    TYPE(tInputOutput)             :: IO                                      !
    REAL                           :: time, dt                                ! Calculation time and timestep
    REAL                           :: MaterialVal(MESH%nElem,               & !
                                                  EQN%nBackgroundVar        ) ! Material Values
    INTEGER                        :: iteration                               ! Current iteration number
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !

    TYPE(tSparseTensor3), POINTER  :: AStar_Sp_ptr => NULL()                  ! Pointer to tSparseTensor3
    TYPE(tSparseTensor3), POINTER  :: BStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: CStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: EStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: FluxSolver_Sp_ptr => NULL()             !

    TYPE(tSparseTensor3b), POINTER :: Tens_xi_Sp_ptr   => NULL()              ! Pointer to tSparseTensor3b
    TYPE(tSparseTensor3b), POINTER :: Tens_eta_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_zeta_Sp_ptr => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_klm_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_k_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_k_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_k_Sp_ptr   => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kxi_m_Sp_ptr     => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Keta_m_Sp_ptr    => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Kzeta_m_Sp_ptr   => NULL()              !

    TYPE(tSparseTensor3b), POINTER :: FluxInt_Sp_ptr    => NULL()             !

    REAL, POINTER                  :: MassMatrix_ptr(:,:) => NULL()           ! Pointer to the corresponding mass matrix

    REAL, POINTER                  :: Dof_ptr(:,:)            => NULL()       ! Dof of element iElem at time tau=0
    REAL, POINTER                  :: TimeIntDof_ptr(:,:)     => NULL()       ! Time integrate dof
    REAL, POINTER                  :: TimeDerDof_ptr(:,:,:)   => NULL()       ! Time derivative dof
    !
    REAL, TARGET                   :: TimeIntDofElem(DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal          ) !
    REAL, TARGET                   :: TimeIntDofSide(DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal,         & !
                                                     MESH%nSideMax          ) !
    REAL, TARGET                   :: NeigTimeIntDofElem(                   & !
                                                     DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal          ) !
    REAL, TARGET                   :: NeigTimeIntDofSide(                   & !
                                                     DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal,         & !
                                                     MESH%nSideMax          ) !
    REAL, TARGET                   :: NeigTimeDerDof(                       & !
                                                     DISC%Galerkin%nDegFr,  & !
                                                     EQN%nVarTotal,         & !
                                                     0:DISC%Galerkin%nPoly  ) !

    REAL                           :: dudt(DISC%Galerkin%nDegFr,            & !
                                           EQN%nVarTotal                    ) ! Auxilliary update dof
    REAL                           :: dudt_flux(DISC%Galerkin%nDegFr,       & !
                                                EQN%nVarTotal               ) ! Auxilliary update dof
    REAL                           :: dudt_src(DISC%Galerkin%nDegFr,        & !
                                               EQN%nVarTotal                ) ! Auxilliary update dof
    REAL                           :: auxvar(DISC%Galerkin%nDegFr,          & !
                                             EQN%nVarTotal,                 & !
                                             DISC%Galerkin%nDegFrMat        ) ! Auxilliary variable
    REAL                           :: JacobiDetVolum                          ! Jacobian of the mapping
    REAL                           :: JacobiDetSurfa                          ! Jacobian of the mapping

    REAL, SAVE                     :: PercentComplete
    REAL, SAVE                     :: NextPercentComplete
    REAL                           :: tmin(0:MESH%nSideMax)
    REAL                           :: tmax(0:MESH%nSideMax)
    REAL                           :: NeigTime(0:MESH%nSideMax)

    REAL                           :: CpuTime_ini                             ! Eval cpu time
    REAL                           :: CpuTime_end                             ! Eval cpu time

    INTEGER                        :: LocnVar                                 !
    INTEGER                        :: LocnPoly                                !
    INTEGER                        :: LocnDegFr                               !
    INTEGER                        :: LocnPolyMat                             !
    INTEGER                        :: LocnDegFrMat                            !
    INTEGER                        :: NeignVar                                !
    INTEGER                        :: NeignPoly                               !
    INTEGER                        :: NeignDegFr                              !
    INTEGER                        :: NeignPolyMat                            !
    INTEGER                        :: NeignDegFrMat                           !
    INTEGER                        :: ReactionTerm                            !

    INTEGER                        :: ncycles
    INTEGER                        :: SumElementUpdates
    INTEGER                        :: ElementsUpdated(MESH%nElem)
    INTEGER                        :: TotalUpdates(2)

    LOGICAL                        :: IterationCompleted
    LOGICAL                        :: UpdateElement(MESH%nElem)

    INTEGER                        :: iObject                                 !
    INTEGER                        :: MPIIndex                                !
    INTEGER                        :: iNeighbor                               !
    INTEGER                        :: iLocalNeighborSide                      !
    INTEGER                        :: iLocalNeighborVrtx                      !
    INTEGER                        :: iDegFr                                  ! Counter
    INTEGER                        :: iElem                                   ! Counter
    INTEGER                        :: iSide                                   ! Counter
    INTEGER                        :: iVar                                    ! Counter
    !-------------------------------------------------------------------------!

    IF(.NOT.DISC%Galerkin%init) THEN
        logError(*) 'ERROR stepQuadFreeGalerkin: SeisSol Interface not initialized!!'
        STOP
    ENDIF
    !
    dt   = MAXVAL(DISC%LocalDt(:))
    time = MAXVAL(DISC%LocalTime(:))
    !
    nCycles             = 0
    SumElementUpdates   = 0
    ElementsUpdated(:)  = 0
    IterationCompleted  = .FALSE.
    !
    ! At first timestep, compute approximate speedup
    ! and do the Cauchy-Kovalewski procedure for all elements
    IF(DISC%iterationstep.EQ.0) THEN
        !
        DO iElem = 1, MESH%nElem
        !
            IF(DISC%Galerkin%pAdaptivity.GE.1) THEN
                LocnPoly  = DISC%Galerkin%LocPoly(iElem)
            ELSE
                LocnPoly  = DISC%Galerkin%nPoly
            ENDIF
            LocnDegFr    = (LocnPoly+1)*(LocnPoly+2)*(LocnPoly+3)/6
            LocnPolyMat  = DISC%Galerkin%nPolyMat
            LocnDegFrMat = (LocnPolyMat+1)*(LocnPolyMat+2)*(LocnPolyMat+3)/6
            LocnVar      = EQN%nVar
            ReactionTerm = 0
            !
            DISC%Galerkin%DGTaylor(:,:,:,iElem) = 0.
            !
            !
            Dof_ptr         => DISC%Galerkin%dgvar(:,:,iElem,1)                               ! Pointer to Dof(iElem)
            TimeDerDof_ptr  => DISC%Galerkin%DGTaylor(1:LocnDegFr,1:LocnVar,0:LocnPoly,iElem) ! Pointer to TimeDerDof(iElem)

            AStar_Sp_ptr    => DISC%Galerkin%AStar_Sp(iElem)                      ! Star matrices in sparse version
            BStar_Sp_ptr    => DISC%Galerkin%BStar_Sp(iElem)                      !
            CStar_Sp_ptr    => DISC%Galerkin%CStar_Sp(iElem)                      !
            EStar_Sp_ptr    => DISC%Galerkin%EStar_Sp(iElem)                      !

            Tens_xi_Sp_ptr   => DISC%Galerkin%ADGxi_Tet_Sp                    ! Stiff CK tensors for tetras
            Tens_eta_Sp_ptr  => DISC%Galerkin%ADGeta_Tet_Sp
            Tens_zeta_Sp_ptr => DISC%Galerkin%ADGzeta_Tet_Sp
            Tens_klm_sp_ptr  => DISC%Galerkin%ADGklm_Tet_Sp                   ! Tensor for space dependent reaction term

            CALL CauchyKovalewski3D(                                            & !
                                   TimeDerDof    = TimeDerDof_ptr,              & ! Output
                                   Dof           = Dof_ptr,                     & ! Input
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
                                   LocDegFr      = LocnDegFr,                   & ! Input
                                   LocDegFrMat   = LocnDegFrMat,                & ! Input
                                   LocPoly       = LocnPoly,                    & ! Input
                                   nVar          = LocnVar                      ) ! Input
        !
        ENDDO
        !
    ENDIF ! DISC%iterationstep.EQ.0
    !
 END SUBROUTINE CKProcedureForEveryoneLTSUni
! GENERATEDKERNELS
#endif

END MODULE Galerkin3D_solver_mod
