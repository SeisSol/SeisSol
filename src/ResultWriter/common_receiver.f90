!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Stefan Wenk (wenk AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/wenk)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE common_receiver_mod
  !--------------------------------------------------------------------------
  USE TypesDef
#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc
    use f_ftoc_bind_interoperability
#endif
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------

  INTERFACE common_ini_receiver
    MODULE PROCEDURE common_ini_receiver
  END INTERFACE

  INTERFACE common_receiver_ck
    MODULE PROCEDURE common_receiver_ck
  END INTERFACE

  INTERFACE common_receiver_interp
    MODULE PROCEDURE common_receiver_interp
  END INTERFACE

  INTERFACE common_receiver_close
    MODULE PROCEDURE common_receiver_close
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC  :: common_ini_receiver
  PUBLIC  :: common_receiver_ck
  PUBLIC  :: common_receiver_interp
  PUBLIC  :: common_receiver_close
  !----------------------------------------------------------------------------

CONTAINS
  !
  !< provide receiver locations for regular and hdf5 output
  !
  SUBROUTINE common_ini_receiver(EQN,MESH,DISC,SOURCE,IO,MPI)
    !--------------------------------------------------------------------------
    USE DGBasis_mod
#ifdef GENERATEDKERNELS
    use f_ftoc_bind_interoperability
    use iso_c_binding, only: c_ptr
#endif
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------
    TYPE (tEquations)        :: EQN
    TYPE (tUnstructMesh)     :: MESH
    TYPE (tDiscretization)   :: DISC
    TYPE (tSource)           :: SOURCE
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    !--------------------------------------------------------------------------
    ! local Variables
    !--------------------------------------------------------------------------
    INTEGER                  :: i, iElem, iPick, iCPU
    INTEGER,ALLOCATABLE      :: MPI_receiver_Element(:,:)
    INTEGER,ALLOCATABLE      :: MPI_receiver_Index(:)
    REAL                     :: io_x, io_y, io_z, t
    REAL                     :: xmin, xmax, ymin, ymax, zmin, zmax
    REAL                     :: xV(MESH%nVertexMax), yV(MESH%nVertexMax), zV(MESH%nVertexMax)
    REAL                     :: xi,eta,zeta
    !--------------------------------------------------------------------------
    INTENT(IN)               :: MESH, DISC, SOURCE
    INTENT(INOUT)            :: IO
    !--------------------------------------------------------------------------
    !
    ! Maximum extend of computational domain
    !
    xmin = MINVAL( MESH%VRTX%xyNode(1,:) )
    xmax = MAXVAL( MESH%VRTX%xyNode(1,:) )
    ymin = MINVAL( MESH%VRTX%xyNode(2,:) )
    ymax = MAXVAL( MESH%VRTX%xyNode(2,:) )
    zmin = MINVAL( MESH%VRTX%xyNode(3,:) )
    zmax = MAXVAL( MESH%VRTX%xyNode(3,:) )
    !
    ! Assume that the point is not inside the domain
    !
    IO%UnstructRecPoint(:)%inside =.false.
    IO%UnstructRecPoint(:)%index  = -1
    !
    ! Loop over all record and PGM points
    !
    IO%nlocalRecordPoint = 0
    DO i=1, IO%ntotalRecordPoint
      !
      io_x = IO%UnstructRecpoint(i)%X
      io_y = IO%UnstructRecpoint(i)%Y
      io_z = IO%UnstructRecpoint(i)%Z
      !
      IF ( (io_x.GE.xmin).and.(io_x.LE.xmax).and.(io_y.GE.ymin).and.(io_y.LE.ymax) .and.(io_z.GE.zmin).and.(io_z.LE.zmax)) THEN
        !
#ifdef OMP
  !$omp parallel private(iElem)
  !$omp do schedule(static)
#endif
        DO iElem = 1, MESH%nElem
          IF(XYZInElement(io_x, io_y, io_z, iElem,1e-5,MESH,DISC)) THEN
#ifdef OMP
  !$omp critical
#endif
            IO%UnstructRecpoint(i)%inside=.true.  ! Point is localized in an element
            IO%UnstructRecPoint(i)%index = iElem  ! Element number
            SELECT CASE(MESH%LocalElemType(iElem))
            CASE(4)
              xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
              yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
              zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
            CASE(6)
              xV(:) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(:,iElem))
              yV(:) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(:,iElem))
              zV(:) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(:,iElem))
            END SELECT
            CALL TrafoXYZ2XiEtaZeta(xi,eta,zeta,io_x,io_y,io_z,xV,yV,zV,MESH%LocalVrtxType(iElem))
            IO%UnstructRecpoint(i)%xi   = xi
            IO%UnstructRecpoint(i)%eta  = eta
            IO%UnstructRecpoint(i)%zeta = zeta
            IO%nlocalRecordPoint = IO%nlocalRecordPoint + 1
#ifdef OMP
  !$omp end critical
#else
            EXIT
#endif
          ENDIF
        ENDDO
#ifdef OMP
    !$omp end parallel   
#endif
      ENDIF
    ENDDO
    !
    !  Clean double receivers that are possibly appearing at MPI boundaries
    !
#ifdef PARALLEL
    IF(MPI%nCPU.GT.1) THEN
      IF(IO%MPIPickCleaningDone.EQ.0) THEN
        ! log info
        logInfo(*) 'Cleaning possible double receiver locations for MPI... '
        !
        ALLOCATE( MPI_receiver_Element(IO%ntotalRecordPoint,0:MPI%nCPU-1) )
        ALLOCATE( MPI_receiver_Index( IO%ntotalRecordPoint ) )

        MPI_receiver_Index(:) = IO%UnstructRecPoint(:)%index

        CALL MPI_ALLGATHER(MPI_receiver_Index,     IO%ntotalRecordPoint,MPI_INTEGER, &
                           MPI_receiver_Element,   IO%ntotalRecordPoint,MPI_INTEGER, &
                           MPI%commWorld, MPI%iErr                                          )

        DO iPick = 1, IO%ntotalRecordPoint
          ! stop execution if one receiver is outside the domain
          IF (SUM(MPI_receiver_Element(iPick,:)).EQ.-MPI%nCPU) THEN
            logError(*) 'Receiver ',iPick,' is not located in the computational domain!'
            STOP
          ENDIF
          !
          IF(IO%UnstructRecpoint(iPick)%inside) THEN
            DO iCPU = 0, MPI%myrank-1
            ! check IO%UnstructRecPoint(iPick)%index .ge. 0 in order to avoid double cleaning
              IF(MPI_receiver_Element(iPick,iCPU).NE.-1 .and. IO%UnstructRecPoint(iPick)%index .ge. 0) THEN
                logInfo(*) ' '
                logInfo(*) 'receiver number ', iPick, ' cleaned in element ', IO%UnstructRecPoint(iPick)%index
                logInfo(*) 'It was already found in element ', MPI_receiver_Element(iPick,iCPU), ' in CPU ', iCPU
                logInfo(*) ' '
                IO%UnstructRecpoint(iPick)%inside = .false.
                IO%UnstructRecPoint(iPick)%index  = -1
                IO%nlocalRecordPoint = IO%nlocalRecordPoint-1
              ENDIF
            ENDDO
          ENDIF
        ENDDO
        !
        DEALLOCATE( MPI_receiver_Element )
        DEALLOCATE( MPI_receiver_Index )
        logInfo(*) 'MPI receiver cleaning done.  '
        IO%MPIPickCleaningDone = 1
      ENDIF
    ENDIF
#else
    DO i = 1, IO%ntotalRecordPoint
      IF(.NOT.IO%UnstructRecpoint(i)%inside) THEN
        logError(*) 'Receiver ',i,' is not located in the computational domain!'
        STOP
      ENDIF
    ENDDO
#endif

#ifdef GENERATEDKERNELS
  do i=1, io%ntotalRecordPoint
    if( io%unstructRecPoint(i)%index .ge. 0 ) then
      call c_interoperability_addReceiver( i, io%unstructRecPoint(i)%index )
    endif
  enddo

  if (IO%pickDtType .eq. 2) then
    logError(*), "no support for time step width dependent receivers."
  else
    call c_interoperability_setReceiverSampling( io%pickdt )
  endif
#endif

  END SUBROUTINE common_ini_receiver
  !
  !< cauchy kovalewski 3d for receiver output
  !
  SUBROUTINE common_receiver_ck(EQN,MESH,DISC,IO,j,TaylorDof,dt,time,localpicktime,dt_op,time_op)
    !--------------------------------------------------------------------------
    USE DGBasis_mod
    USE CauchyKovalewski_mod
    !--------------------------------------------------------------------------
#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_ptr
#endif
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)              :: EQN
    TYPE (tUnstructMesh)           :: MESH
    TYPE (tDiscretization)         :: DISC
    TYPE (tInputOutput)            :: IO
    REAL                           :: dt, time, localpicktime
    REAL                           :: TaylorDOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal,0:DISC%Galerkin%nPolyRec)
    REAL,OPTIONAL                  :: dt_op, time_op
    INTEGER                        :: j
    !--------------------------------------------------------------------------
    REAL                           :: DOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal)
    !--------------------------------------------------------------------------
    TYPE(tSparseTensor3), POINTER  :: AStar_Sp_ptr => NULL()                  ! Pointer to tSparseTensor3
    TYPE(tSparseTensor3), POINTER  :: BStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: CStar_Sp_ptr => NULL()                  !
    TYPE(tSparseTensor3), POINTER  :: EStar_Sp_ptr => NULL()                  !
    !--------------------------------------------------------------------------
    TYPE(tSparseTensor3b), POINTER :: Tens_xi_Sp_ptr   => NULL()              ! Pointer to tSparseTensor3b
    TYPE(tSparseTensor3b), POINTER :: Tens_eta_Sp_ptr  => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_zeta_Sp_ptr => NULL()              !
    TYPE(tSparseTensor3b), POINTER :: Tens_klm_Sp_ptr  => NULL()              !
    !--------------------------------------------------------------------------
    INTEGER                        :: LocElemType                             !
    INTEGER                        :: LocnVar                                 !
    INTEGER                        :: LocnPoly                                !
    INTEGER                        :: LocnDegFr                               !
    INTEGER                        :: LocnPolyMat                             !
    INTEGER                        :: LocnDegFrMat                            !
    INTEGER                        :: ReactionTerm                            !
    INTEGER                        :: iElem                                   !
    !--------------------------------------------------------------------------
    INTENT(IN)                     :: EQN,MESH,DISC,j,time_op,dt_op
#ifdef GENERATEDKERNELS
    intent(out)                    :: taylorDOF
#else
    INTENT(OUT)                    :: TaylorDOF,time,dt,localpicktime
#endif
    !--------------------------------------------------------------------------
    !
    iElem = IO%UnstructRecPoint(j)%index
    !
#ifndef GENERATEDKERNELS
    SELECT CASE(DISC%Galerkin%DGMethod)
    CASE(3)
      localpicktime = IO%LocalPickTime(j)
      time          = DISC%LocalTime(iElem)
      dt            = DISC%LocalDt(iElem)
    CASE DEFAULT
      localpicktime = IO%picktime
      time          = time_op
      dt            = dt_op
    END SELECT
    !
!?  ! usually the receiver are in the refined mesh area and all have equally small timesteps.
!?  ! are there any cases where receivers are located in the coarse part of the mesh, where the timesteps are larger than the
!?  ! sampling rate required in the refined part of the mesh?
    ! exception handling
    IF(IO%pickdt.LT.dt) THEN
      logError(*) ' '
      logError(*) '|   Receiver sampling ', IO%pickdt,'s smaller than timestep ', dt,'s'
      logError(*) '|   Please increase receiver sampling rate'
      logError(*) ' '
      stop
    ENDIF
    !
    IF(EQN%LocAnelastic(iElem).EQ.1.OR.EQN%LocPoroelastic(iElem).EQ.1.OR.EQN%LocPoroelastic(iElem).EQ.2)THEN
      IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ReactionTerm = 0
      ELSE
        ReactionTerm = 1
      ENDIF
    ELSE
      ReactionTerm = 0
    ENDIF
    !
    ! GTS, LTS ADER-DG
    DOF = DISC%Galerkin%dgvar(:,:,iElem,1)
    !
    LocnPoly     = DISC%Galerkin%nPoly
    LocnDegFr    = (LocnPoly+1)*(LocnPoly+2)*(LocnPoly+3)/6
    LocnPolyMat  = DISC%Galerkin%nPolyMat
    LocnDegFrMat = (LocnPolyMat+1)*(LocnPolyMat+2)*(LocnPolyMat+3)/6
    LocnVar      = EQN%nVar
    LocElemType  = MESH%LocalElemType(iElem)
    !
    AStar_Sp_ptr    => DISC%Galerkin%AStar_Sp(iElem)                      ! Star matrices in sparse version
    BStar_Sp_ptr    => DISC%Galerkin%BStar_Sp(iElem)                      !
    CStar_Sp_ptr    => DISC%Galerkin%CStar_Sp(iElem)                      !
    EStar_Sp_ptr    => DISC%Galerkin%EStar_Sp(iElem)                      !

    SELECT CASE(LocElemType)
    CASE(4) ! Tetra
        Tens_xi_Sp_ptr   => DISC%Galerkin%ADGxi_Tet_Sp                    ! Stiff CK tensors for tetras
        Tens_eta_Sp_ptr  => DISC%Galerkin%ADGeta_Tet_Sp
        Tens_zeta_Sp_ptr => DISC%Galerkin%ADGzeta_Tet_Sp
        Tens_klm_sp_ptr  => DISC%Galerkin%ADGklm_Tet_Sp                   ! Tensor for space dependent reaction term
    CASE(6) ! Hexa
        Tens_xi_Sp_ptr   => DISC%Galerkin%ADGxi_Hex_Sp                    ! Stiff CK tensors for hexas
        Tens_eta_Sp_ptr  => DISC%Galerkin%ADGeta_Hex_Sp
        Tens_zeta_Sp_ptr => DISC%Galerkin%ADGzeta_Hex_Sp
        Tens_klm_sp_ptr  => DISC%Galerkin%ADGzeta_Hex_Sp !Temporary assignment We need to create the correct tensor. Used for space dependent reaction term
    CASE DEFAULT
        logError('(" ERROR in ADERGalerkin3D: MESH%LocalElemType(",I5,")=",I2," unknown!!!")') iElem, MESH%LocalElemType(iElem)
        STOP
    END SELECT ! LocElemType

    IF(DISC%Galerkin%CKMethod.EQ.1.AND.EQN%LocPoroelastic(iElem).NE.0) THEN
      ! NOT IMPLEMENTED
      !CALL FastSpaceTimeDG_TimeExpansion(                                                            &
      !                         TaylorDOF     = TaylorDOF,                                            &
      !                         IC_DOF        = DOF,                                                  &
      !                         SpInvSystem   = DISC%Galerkin%InvSystemMatrix(iElem),                 &
      !                         EQN           = EQN,                                                  &
      !                         DISC          = DISC                                                  )
      !
      logError(*) 'Error in SUBROUTINE receiver'
      logError(*) 'FastSpaceTimeDG_TimeExpansion not implemented in unified version'
      STOP
    ELSE
      CALL CauchyKovalewski3D(                                            & !
                             TimeDerDof    = TaylorDof,                   & ! Output
                             Dof           = Dof,                         & ! Input
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
      NULLIFY(AStar_Sp_ptr, BStar_Sp_ptr, CStar_Sp_ptr, EStar_Sp_ptr)
      NULLIFY(Tens_xi_Sp_ptr, Tens_eta_Sp_ptr, Tens_zeta_Sp_ptr, Tens_klm_sp_ptr)
    ENDIF
#else
     call c_interoperability_getTimeDerivatives( i_meshId         = iElem, \
                                                 o_timeDerivatives = taylorDof(:,:,0) )
#endif
    !
  END SUBROUTINE common_receiver_ck
  !
  !< get field values and interpolate to receiver sampling rate
  !
  SUBROUTINE common_receiver_interp(EQN,MESH,DISC,IO,j,TaylorDof,time,localpicktime,state,state_rot)
    !
    USE DGBasis_mod
    USE ConvertXiEta2XY_mod
    USE COMMON_operators_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)         :: EQN
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tDiscretization)    :: DISC
    TYPE (tInputOutput)       :: IO
    REAL                      :: time, localpicktime
    REAL                      :: state(EQN%nVarTotal), state_rot(9)
    REAL                      :: TaylorDOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal,0:DISC%Galerkin%nPolyRec)
    INTEGER                   :: j
    !--------------------------------------------------------------------------
    REAL, POINTER             :: cPoly3D(:,:,:,:,:)         => NULL()         ! Pointer to basis functions coeff
    INTEGER, POINTER          :: NonZeroCPoly(:,:)          => NULL()         !
    INTEGER, POINTER          :: NonZeroCPolyIndex(:,:,:,:) => NULL()         !
    !--------------------------------------------------------------------------
    REAL                      :: x(MESH%nVertexMax), y(MESH%nVertexMax), z(MESH%nVertexMax)
    REAL                      :: xP, yP, zP, xi, eta, zeta
    REAL                      :: DDX(3), DDY(3), DDZ(3)
    REAL                      :: InterpDOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal)
    REAL                      :: tau
    REAL                      :: dtk(0:DISC%Galerkin%nPoly), dtk_tmp(0:DISC%Galerkin%nPoly)
    REAL                      :: XYZcPoly(EQN%nVar,0:DISC%Galerkin%nPolyRec,0:DISC%Galerkin%nPolyRec,0:DISC%Galerkin%nPolyRec)
    REAL                      :: JacobiT(3,3)
    REAL                      :: phigrad(3), grad(3,3), tmp(3,3)
    INTEGER                   :: i, k, h, n, l, iElem, iDegFr
    INTEGER                   :: LocElemType
    real                      :: cPoly3DTmp(0:DISC%Galerkin%nPolyRec, 0:DISC%Galerkin%nPolyRec, 0:DISC%Galerkin%nPolyRec, 0:DISC%Galerkin%nDegFrRec-1)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,MESH,DISC,j,TaylorDof,time,localpicktime
    INTENT(OUT)               :: state, state_rot
    !--------------------------------------------------------------------------
    !
    ! Relative time
    tau = localpicktime - time
    !
    iElem = IO%UnstructRecPoint(j)%index
    LocElemType = MESH%LocalElemType(iElem)
    !
    ! Interpolate in time using Taylor series (or Legendre basis, for ST-DG)
    InterpDOF(:,:) = 0.
    !
    IF(EQN%LocPoroelastic(iElem).NE.2) THEN
      DO i = 0, DISC%Galerkin%nPoly
          InterpDOF(:,:) = InterpDOF(:,:) + TaylorDOF(:,:,i)*tau**i/DISC%Galerkin%Faculty(i)
      ENDDO
    ELSE
      dtk(:) = 0.d0
      DO i=0,DISC%Galerkin%nPoly
        DO k=0,DISC%Galerkin%nPoly
          dtk(i) = dtk(i) + (tau**k)*DISC%Galerkin%cTimePoly(k,i,DISC%Galerkin%nPoly)
        ENDDO
      ENDDO
      !
      DO i = 0, DISC%Galerkin%nPoly
           InterpDOF(:,:) = InterpDOF(:,:) + TaylorDOF(:,:,i)*dtk(i)
      ENDDO
    ENDIF
    !
    SELECT CASE(LocElemType)
    CASE(4) ! Tetra
      cPoly3D           => DISC%Galerkin%cPoly3D_Tet
      NonZeroCPoly      => DISC%Galerkin%NonZeroCPoly_Tet
      NonZeroCPolyIndex => DISC%Galerkin%NonZeroCPolyIndex_Tet
    CASE(6) ! Hexa
      cPoly3D           => DISC%Galerkin%cPoly3D_Hex
      NonZeroCPoly      => DISC%Galerkin%NonZeroCPoly_Hex
      NonZeroCPolyIndex => DISC%Galerkin%NonZeroCPolyIndex_Hex
    END SELECT
    !
    CALL GetStateXiEtaZeta(                                     &
        state             = state,                              &
        xi                = IO%UnstructRecPoint(j)%xi,          &
        eta               = IO%UnstructRecPoint(j)%eta,         &
        zeta              = IO%UnstructRecPoint(j)%zeta,        &
        nDegFr            = DISC%Galerkin%nDegFrRec,            &
        nVar              = EQN%nVar,                           &
        nPoly             = DISC%Galerkin%nPoly,                &
        DOF               = InterpDOF,                          &
        CPoly3D           = CPoly3D,                            &
        NonZeroCPoly      = NonZeroCPoly,                       &
        NonZeroCPolyIndex = NonZeroCPolyIndex                    )
    !
    IF(IO%Rotation.NE.0) THEN
      SELECT CASE(LocElemType)
        CASE(4)
          ! Tetrahedrons
          DO i=1,4
            x(i) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(i,iElem))
            y(i) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(i,iElem))
            z(i) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(i,iElem))
          ENDDO
          !
          !Getting the polynomial of the variables in the X,Y,Z coordinate system
          cPoly3DTmp = cPoly3D(0:DISC%Galerkin%nPolyRec, 0:DISC%Galerkin%nPolyRec, 0:DISC%Galerkin%nPolyRec, 0:DISC%Galerkin%nDegFrRec-1, DISC%Galerkin%nPolyRec)
          CALL ConvertXiEtaZeta2XYZ( &
                   XYZcPoly   = XYZcPoly,                                                &
                   u_hat      = InterpDOF,                                               &
                   cpoly      = cPoly3DTmp,                                              &
                   x          = x,                                                       &
                   y          = y,                                                       &
                   z          = z,                                                       &
                   nVar       = EQN%nVar,                                                &
                   nDegFr     = DISC%Galerkin%nDegFrRec,                                 &
                   nOrdPoly   = DISC%Galerkin%nPolyRec                                   )
          !
          XP = IO%UnstructRecPoint(j)%X
          YP = IO%UnstructRecPoint(j)%Y
          ZP = IO%UnstructRecPoint(j)%Z
          !
          DDX=0.
          DDY=0.
          DDZ=0.
          !Computation of space derivatives of the variables at the given receiver
          DO h=0,DISC%Galerkin%nPolyRec
            DO n=0,DISC%Galerkin%nPolyRec
              DO l=0,DISC%Galerkin%nPolyRec
                IF(h.gt.0) DDX(1:3)=DDX(1:3) + h * XP**(h-1) * YP**(n)   * ZP**(l)   * XYZcPoly(7:9,h,n,l)
                IF(n.gt.0) DDY(1:3)=DDY(1:3) + n * XP**(h)   * YP**(n-1) * ZP**(l)   * XYZcPoly(7:9,h,n,l)
                IF(l.gt.0) DDZ(1:3)=DDZ(1:3) + l * XP**(h)   * YP**(n)   * ZP**(l-1) * XYZcPoly(7:9,h,n,l)
              ENDDO
            ENDDO
          ENDDO
          !
          IF(IO%Rotation.EQ.1)THEN
             !Rotation computation (x 0.5)
             State_rot(1)=0.5*(DDY(3)-DDZ(2))
             State_rot(2)=0.5*(DDZ(1)-DDX(3))
             State_rot(3)=0.5*(DDX(2)-DDY(1))
          ELSEIF(IO%Rotation.EQ.2)THEN
             !Seismic Moment Tensor Contribution (see R.Graves & D.Wald (20019, BSSA, Vol.106, No.B5, p.8745-8766)
             State_rot(1)= DDX(1)
             State_rot(2)= DDY(1)
             State_rot(3)= DDZ(1)
             State_rot(4)= DDX(2)
             State_rot(5)= DDY(2)
             State_rot(6)= DDZ(2)
             State_rot(7)= DDX(3)
             State_rot(8)= DDY(3)
             State_rot(9)= DDZ(3)
          ELSEIF(IO%Rotation.EQ.3)THEN
             !Curl and Divergence computation
             State_rot(1)=DDY(3)-DDZ(2)
             State_rot(2)=DDZ(1)-DDX(3)
             State_rot(3)=DDX(2)-DDY(1)
             State_rot(4)=DDX(1)+DDY(2)+DDZ(3)
          ENDIF
          !
        CASE(6)
          ! Hexahedrons
          DO i=1,8
             x(i) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(i,iElem))
             y(i) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(i,iElem))
             z(i) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(i,iElem))
          ENDDO
          xi   = IO%UnstructRecPoint(j)%xi
          eta  = IO%UnstructRecPoint(j)%eta
          zeta = IO%UnstructRecPoint(j)%zeta
          ! Get the inverse of the Jacobian of the Element-Mapping
          CALL HexaTrafoXiEtaZetaGrad(grad,xi,eta,zeta,x,y,z)
          CALL MatrixInverse3x3(JacobiT,grad)
          ! Get the derivatives of the three velocity components in the xi,eta,zeta-reference system
          grad(:,:) = 0.
          DO iDegFr = 1, DISC%Galerkin%nDegFrRec
               CALL BaseGrad3D(phigrad,iDegFr,xi,eta,zeta,DISC%Galerkin%nPolyRec,CPoly3D,NonZeroCPoly,NonZeroCPolyIndex)
               grad(1,:) = grad(1,:) + phigrad(1)*InterpDOF(iDegFr,7:9)
               grad(2,:) = grad(2,:) + phigrad(2)*InterpDOF(iDegFr,7:9)
               grad(3,:) = grad(3,:) + phigrad(3)*InterpDOF(iDegFr,7:9)
          ENDDO
          ! Get the derivatives of the three velocity components in the x,y,z-physical system
          DDX(:) = 0.
          DDY(:) = 0.
          DDZ(:) = 0.
          !
          tmp = MATMUL(TRANSPOSE(JacobiT),grad)
          !
          DDX(:) = tmp(1,:)
          DDY(:) = tmp(2,:)
          DDZ(:) = tmp(3,:)
          !
          IF(IO%Rotation.EQ.1)THEN
             !Curl computation (x 0.5)
             State_rot(1)=0.5*(DDY(3)-DDZ(2))
             State_rot(2)=0.5*(DDZ(1)-DDX(3))
             State_rot(3)=0.5*(DDX(2)-DDY(1))
          ELSEIF(IO%Rotation.EQ.2)THEN
             !Seismic Moment Tensor Contribution (see R.Graves & D.Wald (20019, BSSA, Vol.106, No.B5, p.8745-8766)
             State_rot(1)= DDX(1)
             State_rot(2)= DDY(1)
             State_rot(3)= DDZ(1)
             State_rot(4)= DDX(2)
             State_rot(5)= DDY(2)
             State_rot(6)= DDZ(2)
             State_rot(7)= DDX(3)
             State_rot(8)= DDY(3)
             State_rot(9)= DDZ(3)
          ELSEIF(IO%Rotation.EQ.3)THEN
             !Curl and Divergence computation
             State_rot(1)=DDY(3)-DDZ(2)
             State_rot(2)=DDZ(1)-DDX(3)
             State_rot(3)=DDX(2)-DDY(1)
             State_rot(4)=DDX(1)+DDY(2)+DDZ(3)
          ENDIF
          !
        END SELECT ! LocElemType
        !
      WHERE (abs(State_rot(:)).LE.1e-40)
           State_rot(:) = 0.0
      END WHERE
      !
    ENDIF ! IO%Rotation.NE.0
    !
    WHERE (abs(State(:)).LE.1e-40)
         State(:) = 0.0
    END WHERE
    !
    NULLIFY(CPoly3D)
    NULLIFY(NonZeroCPoly)
    NULLIFY(NonZeroCPolyIndex)
    !
  END SUBROUTINE common_receiver_interp
  !
  !< close receiver
  !< deallocate derived datatypes
  SUBROUTINE common_receiver_close(DISC,IO,MPI)
    !--------------------------------------------------------------------------
#ifdef HDF
    USE hdf_output_utils_mod
#endif
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tDiscretization)   :: DISC
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    ! local Variables
    INTEGER                  :: i
    !--------------------------------------------------------------------------
    INTENT(IN)               :: DISC, MPI
    INTENT(INOUT)            :: IO
    !--------------------------------------------------------------------------
    !
    IF(ASSOCIATED(IO%UNIT%VFile)) THEN
      DEALLOCATE( IO%UNIT%VFile)
    ENDIF
    IF(ASSOCIATED(IO%LocalPickTime)) THEN
      DEALLOCATE( IO%LocalPickTime )
    ENDIF
#ifdef HDF
    IF(IO%nlocalRecordPoint.NE.0) THEN
      CALL close_hdf5_file(IO%hd_rec%file_id)
      ! mpi destructor
      CALL MPI_Group_free(IO%hd_rec%mpi_grp, MPI%iErr)
      CALL MPI_Comm_free(IO%hd_rec%mpi_comm, MPI%iErr)
      DEALLOCATE(IO%hd_rec)
    ENDIF
#endif
    !
  END SUBROUTINE common_receiver_close

END MODULE common_receiver_mod
