!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2008-2016, SeisSol Group
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

#include <Initializer/preProcessorMacros.fpp>

MODULE MPIExchangeValues_mod
  !----------------------------------------------------------------------------
  USE typesDef
  USE COMMON_operators_mod
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
#ifndef GENERATEDKERNELS
  !----------------------------------------------------------------------------
  INTERFACE MPIExchangeInvSystems
     MODULE PROCEDURE MPIExchangeInvSystems
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC :: MPIExchangeValues_GTS_init
  PUBLIC :: MPIExchangeValues_GTS
  PUBLIC :: MPIExchangeValues_LTS
  PUBLIC :: MPIExchangeJacobians_new
  PUBLIC :: MPIExchangeInvSystems
#endif
  PUBLIC :: MPIExchangeBackground
  !----------------------------------------------------------------------------

  !> Data structures for optimized GTS
  TYPE(tRealMessage), ALLOCATABLE :: send_message_gts(:)
  TYPE(tRealMessage), ALLOCATABLE :: recv_message_gts(:)
  INTEGER, ALLOCATABLE            :: send_request_gts(:)
  INTEGER, ALLOCATABLE            :: recv_request_gts(:)
  INTEGER, ALLOCATABLE            :: send_status_list_gts(:,:)
  INTEGER, ALLOCATABLE            :: recv_status_list_gts(:,:)
  integer, ALLOCATABLE :: bndDomainElementList(:,:)
  integer :: nBndDomainElementList

CONTAINS

#ifndef GENERATEDKERNELS
  !> Optimized version of MPIExchangeValues for GTS only
  !! Initialization, only done once
  !! Attention: Anelasticity not supported yet!!
  SUBROUTINE MPIExchangeValues_GTS_init(DISC,EQN,BND,MPI)
    IMPLICIT NONE

#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif

    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tDiscretization)          :: DISC
    TYPE (tEquations)               :: EQN
    TYPE (tBoundary)                :: BND
    TYPE (tMPI)                     :: MPI
    ! Local variable declaration
    INTEGER                         :: iDomain, MsgLength, LoopDegFr
    INTEGER                         :: i, ierr, iCPU, counter
    INTEGER, PARAMETER              :: MIN_MSG_LENGTH = 64

    IF(MPI%nCPU.EQ.1) THEN
        RETURN
    ENDIF

    LoopDegFr = DISC%Galerkin%nDegFrRec

#ifdef PARALLEL

    ALLOCATE(send_message_gts(BND%NoMPIDomains),        &
             recv_message_gts(BND%NoMPIDomains),        &
             send_request_gts(BND%NoMPIDomains),        &
             recv_request_gts(BND%NoMPIDomains),        &
             send_status_list_gts(MPI_STATUS_SIZE,BND%NoMPIDomains),    &
             recv_status_list_gts(MPI_STATUS_SIZE,BND%noMPIDomains))

    nBndDomainElementList = 0
    DO iDomain = 1, BND%NoMPIDomains
      ! Compute total length of message for each MPI domain interface
      ! Length = No. of MPI boundary Elements * ( No. of DOF * No. of Variables + No. of background values )
      ! In the anelastic case, also the anelastic functions have to be exchanged!
      ! Length = No. of MPI boundary Elements * ( No. of DOF * (No. of Variables + No. of Mechanisms * No. of anel. functions per Mechanism) + No. of background values )
      IF (EQN%DR.EQ.1) THEN
        ! exchange additionally the dgvar values of the fault elements needed for the Godunov state
        MsgLength = BND%ObjMPI(iDomain)%nElem*(LoopDegFr*EQN%nVarTotal) &
                  + BND%ObjMPI(iDomain)%nFault_MPI*LoopDegFr*EQN%nVarTotal
      ELSE
        MsgLength = BND%ObjMPI(iDomain)%nElem*(LoopDegFr*EQN%nVarTotal)
      ENDIF
      ! Allocate space for message
      MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
      ALLOCATE(send_message_gts(iDomain)%Content(MsgLength))
      ALLOCATE(recv_message_gts(iDomain)%Content(MsgLength))

      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU

      ! Initialize persistant MPI requests
      call MPI_SEND_INIT(send_message_gts(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                     send_request_gts(iDomain), iErr)
      call MPI_RECV_INIT(recv_message_gts(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                     recv_request_gts(iDomain), iErr)

      ! Count total number of boundary elements
      nBndDomainElementList = nBndDomainElementList + BND%ObjMPI(iDomain)%nElem
    ENDDO

    ALLOCATE(bndDomainElementList(2, nBndDomainElementList))

    ! Initialize bndDomainElementList
    counter = 0
    DO iDomain = 1, BND%NoMPIDomains
      DO i = 1, BND%ObjMPI(iDomain)%nElem
        counter = counter + 1
        bndDomainElementList(1,counter) = iDomain
        bndDomainElementList(2,counter) = i
      ENDDO
    ENDDO

  ! TODO: Deallocation of communication data structures is currently missing
#endif

  END SUBROUTINE MPIExchangeValues_GTS_init

  !> Optimized version of MPIExchangeValues for GTS only
  !! Attention: Anelasticity not supported yet!!
  SUBROUTINE MPIExchangeValues_GTS(DISC,EQN,BND)
    IMPLICIT NONE

#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif

    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tDiscretization)          :: DISC
    TYPE (tEquations)               :: EQN
    TYPE (tBoundary)                :: BND
    ! Local variable declaration
    INTEGER                         :: i, iDomain, counter, counter2, iElem, LoopDegFr
    INTEGER                         :: iDegFr, iVar, iErr

    !--------------------------------------------------------------------------

#ifdef PARALLEL

    LoopDegFr = DISC%Galerkin%nDegFrRec
   
#ifdef OMP
!$omp parallel
!$omp do schedule(static) private(counter, iDomain, i, iElem, counter2, iDegFr, iVar)
#endif
    DO counter = 1,nBndDomainElementList
      iDomain = bndDomainElementList(1, counter)
      i = bndDomainElementList(2, counter)

      iElem = BND%ObjMPI(iDomain)%DomainElements(i)

        ! Send DG degrees of freedom to neighbor domain
        counter2 = (i-1) * LoopDegFr*EQN%nVarTotal
        ! TODO: Use Fotran syntax
        DO iVar = 1, EQN%nVarTotal
          DO iDegFr = 1, LoopDegFr
            counter2 = counter2 + 1
            !Set up list of all DOF of the ELASTIC variables
            send_message_gts(iDomain)%Content(counter2) = DISC%Galerkin%dgwork(iDegFr,iVar,iElem)
          ENDDO
        ENDDO
!        ! Send elastic parameters (rho, mu, lambda) and anelastic coefficients (omega_n, YP_n, YS_n) to neighbor domain
!        DO iVar = 1, EQN%nBackgroundVar
!            counter = counter + 1
!            IF(counter.GT.MsgLength) THEN
!              logError(*) 'Sending message: Counter > MsgLength', counter, MsgLength
!              STOP
!            ENDIF
!            send_message(iDomain)%Content(counter) = OptionalFields%BackgroundValue(iElem,iVar)
!        ENDDO
    ENDDO
#ifdef OMP
!$omp end parallel
#endif

      !
      ! For Dynamic Rupture simulation add dgvar values now at the end of send_message
      ! -> don't modify variable counter!
      !
    IF (EQN%DR.EQ.1) THEN
#ifdef OMP
!$omp parallel
!$omp do schedule(static) private(counter, iDomain, i, iElem, iDegFr, iVar)
#endif
      DO iDomain = 1, BND%NoMPIDomains
         counter = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal) ! Correct starting point (after dgwork)
         DO i = 1, BND%ObjMPI(iDomain)%nFault_MPI
            iElem        = BND%ObjMPI(iDomain)%Domain_Fault_Elem(i)
            !
            DO iVar = 1, EQN%nVarTotal
              DO iDegFr = 1, LoopDegFr
                counter = counter + 1
                !Set up list of all DOF of the ELASTIC variables
                send_message_gts(iDomain)%Content(counter) = DISC%Galerkin%dgvar(iDegFr,iVar,iElem,1)
              ENDDO
            ENDDO
            !
         ENDDO ! nFault_MPI
      ENDDO
#ifdef OMP
!$omp end parallel
#endif
    ENDIF ! (EQN%DR.EQ.1)

    ! Start all communication
    call MPI_STARTALL(BND%NoMPIDomains,send_request_gts, ierr)
    call MPI_STARTALL(BND%NoMPIDomains,recv_request_gts, ierr)

    ! Wait until all communication has finished
    CALL MPI_WAITALL(BND%NoMPIDomains,send_request_gts,send_status_list_gts,ierr)
    CALL MPI_WAITALL(BND%NoMPIDomains,recv_request_gts,recv_status_list_gts,ierr)

#ifdef OMP
!$omp parallel
!$omp do schedule(static) private(counter, iDomain, i, iElem, counter2, iDegFr, iVar)
#endif
    ! Decode 1D information array obtained from neighbor CPUs
    DO counter = 1,nBndDomainElementList
      iDomain = bndDomainElementList(1, counter)
      i = bndDomainElementList(2, counter)

      counter2 = (i-1) * LoopDegFr*EQN%nVarTotal
      DO iVar = 1, EQN%nVarTotal
        DO iDegFr = 1, LoopDegFr
          counter2 = counter2 + 1
          BND%ObjMPI(iDomain)%NeighborDOF(iDegFr,iVar,i) = recv_message_gts(iDomain)%content(counter2)
        ENDDO
      ENDDO
        ! Decode elastic parameters (rho, mu, lambda) and anelastic coefficients (omega_n, YP_n, YS_n)
!        DO iVar = 1, EQN%nBackgroundVar
!            counter = counter + 1
!            IF(counter.GT.MsgLength) THEN
!              !logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
!              STOP
!            ENDIF
!            BND%ObjMPI(iDomain)%NeighborBackground(iVar,i) = recv_message(iDomain)%Content(counter)
!        ENDDO
        ! Get additional information for the local time stepping scheme
    ENDDO
#ifdef OMP
!$omp end parallel
#endif

    ! Decode 1D fault information (dgvar)
    IF (EQN%DR.EQ.1) THEN
#ifdef OMP
!$omp parallel
!$omp do schedule(static) private(counter, iDomain, i, iElem, iDegFr, iVar)
#endif
      DO iDomain = 1, BND%NoMPIDomains
         counter = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal) ! Correct starting point (after dgwork)
         DO i = 1, BND%ObjMPI(iDomain)%nFault_MPI
            !
            DO iVar = 1, EQN%nVarTotal
              DO iDegFr = 1, LoopDegFr
                counter = counter + 1

                BND%ObjMPI(iDomain)%MPI_DR_dgvar(iDegFr,iVar,i) = recv_message_gts(iDomain)%Content(counter)! DISC%Galerkin%dgvar(iDegFr,iVar,iElem)
              ENDDO
            ENDDO
            !
         ENDDO ! nFault_MPI
      ENDDO
#ifdef OMP
!$omp end parallel
#endif
    ENDIF ! (EQN%DR.EQ.1)

#endif
  END SUBROUTINE MPIExchangeValues_GTS

  !> Local timesteping version based on the old code
  SUBROUTINE MPIExchangeValues_LTS(DISC,EQN,BND,MESH,IO,OptionalFields,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#ifdef VTRACE
!    INCLUDE 'VT.inc'
#endif
#endif
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tDiscretization)          :: DISC
    TYPE (tEquations)               :: EQN
    TYPE (tBoundary)                :: BND
    TYPE (tUnstructMesh)            :: MESH
    TYPE (tInputOutput)             :: IO
    TYPE (tUnstructOptionalFields)  :: OptionalFields
    TYPE (tMPI)                     :: MPI
    ! Local variable declaration
    INTEGER                         :: iDomain, iCPU, iElem, iSide, iDegFr
    INTEGER                         :: iNeighbor, iObject, MPIIndex
    INTEGER                         :: iVar, iPoly, iMech, i,j,k
    INTEGER                         :: counter, cnt, MsgLength
    INTEGER                         :: UpdateElement
    INTEGER                         :: iErr
    TYPE(tRealMessage), ALLOCATABLE :: send_message(:)
    TYPE(tRealMessage), ALLOCATABLE :: recv_message(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: send_imessage(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: recv_imessage(:)
    INTEGER, ALLOCATABLE            :: send_request(:)
    INTEGER, ALLOCATABLE            :: recv_request(:)
    INTEGER, ALLOCATABLE            :: send_status_list(:,:)
    INTEGER, ALLOCATABLE            :: recv_status_list(:,:)
    REAL, ALLOCATABLE               :: SendDOF(:,:) 
    REAL                            :: DuDt(DISC%Galerkin%nDegFr,EQN%nVarTotal)
    REAL                            :: tmin,timeMPI,dtmax,dtmaxMPI
    INTEGER, PARAMETER              :: MIN_MSG_LENGTH = 64
    INTEGER                         :: LocElemType, LoopDegFr 
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,MESH,IO,MPI
    INTENT(INOUT)             :: BND,DISC
    !--------------------------------------------------------------------------

    IF(MPI%nCPU.EQ.1) THEN
        RETURN
    ENDIF

    IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        LoopDegFr = DISC%Galerkin%nDegFrST 
    ELSE
        LoopDegFr = DISC%Galerkin%nDegFrRec 
    ENDIF

#ifdef PARALLEL      
!#ifdef VTRACE
!      CALL VTBEGIN(MPI%Trace%hFunctions(3), MPI%iErr)
!#endif

    ALLOCATE(send_message(BND%NoMPIDomains),        &
             recv_message(BND%NoMPIDomains),        &
             send_imessage(BND%NoMPIDomains),       &
             recv_imessage(BND%NoMPIDomains),       &             
             send_request(BND%NoMPIDomains),        &
             recv_request(BND%NoMPIDomains),        &
             send_status_list(MPI_STATUS_SIZE,BND%NoMPIDomains),    &
             recv_status_list(MPI_STATUS_SIZE,BND%noMPIDomains),    &
             SendDOF(LoopDegFr,EQN%nVarTotal)       )

    DO iDomain = 1, BND%NoMPIDomains
      ! Compute total length of message for each MPI domain interface
      ! Length = No. of MPI boundary Elements * ( No. of DOF * No. of Variables + No. of background values )
      ! In the anelastic case, also the anelastic functions have to be exchanged!
      ! Length = No. of MPI boundary Elements * ( No. of DOF * (No. of Variables + No. of Mechanisms * No. of anel. functions per Mechanism) + No. of background values )
      IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrST*EQN%nVarTotal+EQN%nBackgroundVar)
        ! For method with local timestep, also exchange the 
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*(EQN%nVar+EQN%nAneFuncperMech)+2)
      ELSEIF (EQN%DR.EQ.1) THEN
        ! exchange additionally the dgvar values of the fault elements needed for the Godunov state
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal+EQN%nBackgroundVar) &
                  + BND%ObjMPI(iDomain)%nFault_MPI*DISC%Galerkin%nDegFrRec*EQN%nVarTotal
      ELSE
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal+EQN%nBackgroundVar) 
      ENDIF
      ! Allocate space for message
      MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
      ALLOCATE(send_message(iDomain)%Content(MsgLength))
      ALLOCATE(recv_message(iDomain)%Content(MsgLength))

      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
      counter = 0
      ! Put all information into a 1D array
      DO i = 1, BND%ObjMPI(iDomain)%nElem
        iElem = BND%ObjMPI(iDomain)%DomainElements(i)
        IF(DISC%Galerkin%DGMethod.EQ.3) THEN
            DO iVar = 1, EQN%nVarTotal 
             cnt = 0 
             DO iDegFr = 1, DISC%Galerkin%nDegFrRec 
              DO iPoly = 0, DISC%Galerkin%nPolyRec 
                cnt = cnt + 1 
                SendDOF(cnt,iVar) = DISC%Galerkin%dgtaylor(iDegFr,iVar,iPoly,iElem) 
              ENDDO
             ENDDO
            ENDDO
        ELSE
            SendDOF(:,:) = DISC%Galerkin%dgwork(:,:,iElem)
        ENDIF 
        ! Send DG degrees of freedom to neighbor domain
        DO iDegFr = 1, LoopDegFr 
          DO iVar = 1, EQN%nVarTotal
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Sending message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            !Set up list of all DOF of the ELASTIC variables
            send_message(iDomain)%Content(counter) = SendDOF(iDegFr,iVar)  ! DISC%Galerkin%dgwork(iDegFr,iVar,iElem) 
          ENDDO
        ENDDO
        ! Send elastic parameters (rho, mu, lambda) and anelastic coefficients (omega_n, YP_n, YS_n) to neighbor domain
        DO iVar = 1, EQN%nBackgroundVar 
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Sending message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            send_message(iDomain)%Content(counter) = OptionalFields%BackgroundValue(iElem,iVar)
        ENDDO
        ! Send additional information for the local time stepping scheme
        IF(DISC%Galerkin%DGMethod.EQ.3) THEn
            DO iDegFr = 1, DISC%Galerkin%nDegFr 
              DO iVar = 1, EQN%nVar+EQN%nAneFuncperMech
                counter = counter + 1
                IF(counter.GT.MsgLength) THEN
                  logError(*) 'Sending message: Counter > MsgLength', counter, MsgLength
                  STOP
                ENDIF
                ! For the local timestepping scheme, communicate all updates as well as the local time and timestep
                send_message(iDomain)%Content(counter) = BND%ObjMPI(iDomain)%NeighborDuDt(iDegFr,iVar,i) 
              ENDDO
            ENDDO
            BND%ObjMPI(iDomain)%NeighborDuDt(:,:,i) = 0.
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Sending message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            send_message(iDomain)%Content(counter) = DISC%LocalTime(iElem)
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Sending message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            send_message(iDomain)%Content(counter) = DISC%LocalDt(iElem)
        ENDIF
      ENDDO
      !
      ! For Dynamic Rupture simulation add dgvar values now at the end of send_message
      ! -> don't modify variable counter!
      !
      IF (EQN%DR.EQ.1) THEN
         DO i = 1, BND%ObjMPI(iDomain)%nFault_MPI
            iElem        = BND%ObjMPI(iDomain)%Domain_Fault_Elem(i)
            SendDOF(:,:) = DISC%Galerkin%dgvar(1:LoopDegFr,:,iElem,1)
            !
            ! copy SendDOF to a 1D vector
            DO iDegFr = 1, LoopDegFr 
              DO iVar = 1, EQN%nVarTotal
                counter = counter + 1
                IF(counter.GT.MsgLength) THEN
                  logError(*) 'Sending message: Counter > MsgLength', counter, MsgLength
                  STOP
                ENDIF
                !Set up list of all DOF of the ELASTIC variables
                send_message(iDomain)%Content(counter) = SendDOF(iDegFr,iVar)  ! DISC%Galerkin%dgvar(iDegFr,iVar,iElem) 
              ENDDO
            ENDDO
            !
         ENDDO ! nFault_MPI
      
      ENDIF ! (EQN%DR.EQ.1)
      !
      ! Post a send request to neighbor CPU
      CALL MPI_ISEND(send_message(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                     send_request(iDomain), iErr)
    ENDDO

    
    DO iDomain = 1, BND%NoMPIDomains
      ! Compute total length of message for each MPI domain interface
      ! Length = No. of MPI boundary Elements * ( No. of DOF * No. of Variables + No. of background values )
      ! In the anelastic case, also the anelastic functions have to be exchanged!
      ! Length = No. of MPI boundary Elements * ( No. of DOF * (No. of Variables + No. of Mechanisms * No. of anel. functions per Mechanism) + No. of background values )
      IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ! For method with local timestep, also exchange the 
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrST*EQN%nVarTotal+EQN%nBackgroundVar)
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*(EQN%nVar+EQN%nAneFuncperMech)+2)
      ELSEIF (EQN%DR.EQ.1) THEN
        ! exchange additionally the dgvar values of the fault elements needed for the Godunov state
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal+EQN%nBackgroundVar) &
                  + BND%ObjMPI(iDomain)%nFault_MPI*DISC%Galerkin%nDegFrRec*EQN%nVarTotal          
      ELSE 
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal+EQN%nBackgroundVar)
      ENDIF
      MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
      ! Post a receive request from neighbor CPU
      CALL MPI_IRECV(recv_message(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                     recv_request(iDomain), iErr)
    ENDDO

    ! Wait until all communication has finished
    CALL MPI_WAITALL(BND%NoMPIDomains,send_request,send_status_list,ierr)
    CALL MPI_WAITALL(BND%NoMPIDomains,recv_request,recv_status_list,ierr)

    ! Decode 1D information array obtained from neighbor CPUs
    DO iDomain = 1, BND%NoMPIDomains
      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
      IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        ! For method with local timestep, also exchange the 
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrST*EQN%nVarTotal+EQN%nBackgroundVar)
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*(EQN%nVar+EQN%nAneFuncperMech)+2)
      ELSEIF (EQN%DR.EQ.1) THEN
        ! exchange additionally the dgvar values of the fault elements needed for the Godunov state
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal+EQN%nBackgroundVar) &
                  + BND%ObjMPI(iDomain)%nFault_MPI*DISC%Galerkin%nDegFrRec*EQN%nVarTotal        
      ELSE
        MsgLength = BND%ObjMPI(iDomain)%nElem*(DISC%Galerkin%nDegFrRec*EQN%nVarTotal+EQN%nBackgroundVar) 
      ENDIF
      counter = 0
      ! Get all information from the 1D array
      DO i = 1, BND%ObjMPI(iDomain)%nElem
        iElem = BND%ObjMPI(iDomain)%DomainElements(i)
        DO iDegFr = 1, LoopDegFr  
         DO iVar = 1, EQN%nVarTotal
            counter = counter + 1 
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            BND%ObjMPI(iDomain)%NeighborDOF(iDegFr,iVar,i) = recv_message(iDomain)%content(counter)
         ENDDO
        ENDDO
        ! Decode elastic parameters (rho, mu, lambda) and anelastic coefficients (omega_n, YP_n, YS_n)  
        DO iVar = 1, EQN%nBackgroundVar 
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              !logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            BND%ObjMPI(iDomain)%NeighborBackground(iVar,i) = recv_message(iDomain)%Content(counter) 
        ENDDO
        ! Get additional information for the local time stepping scheme
        IF(DISC%Galerkin%DGMethod.EQ.3) THEN
            DO iDegFr = 1, DISC%Galerkin%nDegFr 
              DO iVar = 1, EQN%nVar+EQN%nAneFuncperMech
                counter = counter + 1
                IF(counter.GT.MsgLength) THEN
                  logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
                  STOP
                ENDIF
                ! For the local timestepping scheme, communicate all updates as well as the local time and timestep
                DuDt(iDegFr,iVar) = recv_message(iDomain)%content(counter)
              ENDDO
            ENDDO
            DISC%Galerkin%dgwork(:,:,iElem) = DISC%Galerkin%dgwork(:,:,iElem) + DuDt(:,:)
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            BND%ObjMPI(iDomain)%NeighborTime(i) = recv_message(iDomain)%content(counter)
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
              STOP
            ENDIF
            BND%ObjMPI(iDomain)%NeighborDt(i) = recv_message(iDomain)%content(counter)
        ENDIF
      ENDDO
      !
      ! Decode 1D fault information (dgvar)
      IF (EQN%DR.EQ.1) THEN
         DO i = 1, BND%ObjMPI(iDomain)%nFault_MPI
            !
            DO iDegFr = 1, LoopDegFr 
              DO iVar = 1, EQN%nVarTotal
                counter = counter + 1
                IF(counter.GT.MsgLength) THEN
                  logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
                  STOP
                ENDIF

                BND%ObjMPI(iDomain)%MPI_DR_dgvar(iDegFr,iVar,i) = recv_message(iDomain)%Content(counter)! DISC%Galerkin%dgvar(iDegFr,iVar,iElem) 
              ENDDO
            ENDDO
            !
         ENDDO ! nFault_MPI  
      ENDIF ! (EQN%DR.EQ.1)      
      !
    ENDDO
    
    ! Deallocate all communication variables
    DO iDomain = 1, BND%NoMPIDomains
      DEALLOCATE(send_message(iDomain)%Content)
      DEALLOCATE(recv_message(iDomain)%Content)
    ENDDO
    !
    ! Can MPI neighbors do an update? This is important for communication of the dudt variable
    !
    IF(DISC%Galerkin%DGMethod.EQ.3) THEN

        DO iDomain = 1, BND%NoMPIDomains
          ! Compute total length of message for each MPI domain interface
          MsgLength = BND%ObjMPI(iDomain)%nElem
          MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
          ! Allocate space for message
          ALLOCATE(send_imessage(iDomain)%Content(MsgLength))
          ALLOCATE(recv_imessage(iDomain)%Content(MsgLength))
          ! Number of neighbor CPU
          iCPU = BND%ObjMPI(iDomain)%CPU
          counter = 0
          ! Put all information into a 1D array
          DO i = 1, BND%ObjMPI(iDomain)%nElem
            iElem = BND%ObjMPI(iDomain)%DomainElements(i)
            ! Send 1 if element can do an update within this local timestep or send 0 if it cannot.
            tmin = 1e20
            LocElemType = MESH%LocalElemType(iElem)
            DO iSide = 1, LocElemType
                IF (MESH%ELEM%MPIReference(iSide,iElem).EQ.1) THEN
                    iObject  = MESH%ELEM%BoundaryToObject(iSide,iElem)
                    MPIIndex = MESH%ELEM%MPINumber(iSide,iElem)
                    IF(MPIIndex.EQ.-1) THEN
                        PRINT *, 'Severe Error in Galerkin3D. MPIIndex = -1 !', iElem, iSide
                        STOP
                    ENDIF
                    timeMPI = BND%ObjMPI(iObject)%NeighborTime(MPIIndex) + BND%ObjMPI(iObject)%NeighborDt(MPIIndex)
                    tmin    = MIN( tmin,timeMPI)
                ELSE
                    SELECT CASE (MESH%ELEM%Reference(iSide,iElem))
                    CASE(0,6)
                        iNeighbor = MESH%ELEM%SideNeighbor(iSide,iElem)
                        timeMPI   = DISC%LocalTime(iNeighbor)+DISC%LocalDt(iNeighbor)
                        tmin      = MIN( tmin,timeMPI)
                    ENDSELECT
                ENDIF
            ENDDO ! iSide
            !
            IF( (DISC%LocalTime(iElem)+DISC%LocalDt(iElem)).LE.(tmin)) THEN
                UpdateElement = 1
            ELSE
                UpdateElement = 0
            ENDIF
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              PRINT *, 'Fatal error: Counter > MsgLength', counter, MsgLength, ' in CPU : ', MPI%myrank
              STOP
            ENDIF
            send_imessage(iDomain)%Content(counter) = UpdateElement
          ENDDO
          ! Post a send request to neighbor CPU
          CALL MPI_ISEND(send_imessage(iDomain)%Content, MsgLength, MPI_INTEGER, iCPU, 1, MPI%commWorld, &
                         send_request(iDomain), iErr)
        ENDDO

        DO iDomain = 1, BND%NoMPIDomains
          ! Compute total length of message for each MPI domain interface
          MsgLength = BND%ObjMPI(iDomain)%nElem
          MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
          ! Number of neighbor CPU
          iCPU = BND%ObjMPI(iDomain)%CPU
          ! Post a receive request from neighbor CPU
          CALL MPI_IRECV(recv_imessage(iDomain)%Content, MsgLength, MPI_INTEGER, iCPU, 1, MPI%commWorld, &
                         recv_request(iDomain), iErr)
        ENDDO

        ! Wait until all communication has finished
        CALL MPI_WAITALL(BND%NoMPIDomains,send_request,send_status_list,ierr)
        CALL MPI_WAITALL(BND%NoMPIDomains,recv_request,recv_status_list,ierr)

        ! Decode 1D information array obtained from neighbor CPUs
        DO iDomain = 1, BND%NoMPIDomains
          ! Number of neighbor CPU
          iCPU = BND%ObjMPI(iDomain)%CPU
          MsgLength = BND%ObjMPI(iDomain)%nElem
          counter = 0
          ! Get all information from the 1D array
          DO i = 1, BND%ObjMPI(iDomain)%nElem
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              logError(*) 'Fatal error: Counter > MsgLength', counter, MsgLength, ' in CPU : ', MPI%myrank
              STOP
            ENDIF
            BND%ObjMPI(iDomain)%NeighborUpdate(i) = recv_imessage(iDomain)%Content(counter) 
          ENDDO
        ENDDO
    
        ! Deallocate all communication variables
        DO iDomain = 1, BND%NoMPIDomains
          DEALLOCATE(send_imessage(iDomain)%Content)
          DEALLOCATE(recv_imessage(iDomain)%Content)
        ENDDO

    ENDIF
    !
    DEALLOCATE(send_message,       &
               recv_message,       &
               send_imessage,      &
               recv_imessage,      &
               send_request,       &
               recv_request,       &
               send_status_list,   &
               recv_status_list,   &
               SendDOF             )

#ifdef VTRACE
!      CALL VTEND(MPI%Trace%hFunctions(3), MPI%iErr)
#endif
#endif
  END SUBROUTINE MPIExchangeValues_LTS

!> Communicates Star matrices only once before the time loop starts
!<
  SUBROUTINE MPIExchangeJacobians_new(DISC,EQN,BND,MESH,IO,OptionalFields,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#ifdef VTRACE
!    INCLUDE 'VT.inc'
#endif
#endif
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tDiscretization)          :: DISC
    TYPE (tEquations)               :: EQN
    TYPE (tBoundary)                :: BND
    TYPE (tUnstructMesh)            :: MESH
    TYPE (tInputOutput)             :: IO
    TYPE (tUnstructOptionalFields)  :: OptionalFields
    TYPE (tMPI)                     :: MPI
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,MESH,IO,MPI
    INTENT(INOUT)             :: BND,DISC
    !--------------------------------------------------------------------------
    ! Local variable declaration
    INTEGER                         :: iDomain, iCPU, iElem, iSide, iDegFr, iLocElem
    INTEGER                         :: iNeighbor, iObject, MPIIndex
    INTEGER                         :: iVar, iMech, i,j,k
    INTEGER                         :: counter, MsgLength
    INTEGER                         :: UpdateElement
    INTEGER                         :: iErr
    INTEGER                         :: p,q
    TYPE(tRealMessage), ALLOCATABLE :: send_message(:)
    TYPE(tRealMessage), ALLOCATABLE :: recv_message(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: send_imessage(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: recv_imessage(:)
    INTEGER, ALLOCATABLE            :: send_request(:)
    INTEGER, ALLOCATABLE            :: recv_request(:)
    INTEGER, ALLOCATABLE            :: send_status_list(:,:)
    INTEGER, ALLOCATABLE            :: recv_status_list(:,:)
    REAL                            :: DuDt(DISC%Galerkin%nDegFr,EQN%nVarTotal)
    REAL                            :: tmin,timeMPI,dtmax,dtmaxMPI
    INTEGER, PARAMETER              :: MIN_MSG_LENGTH = 64
    REAL, POINTER                   :: AStar(:,:,:)
    REAL, POINTER                   :: BStar(:,:,:)
    REAL, POINTER                   :: CStar(:,:,:)
    REAL, POINTER                   :: EStar(:,:,:)
    REAL, POINTER                   :: FLStar(:,:,:)
    REAL, POINTER                   :: FRStar(:,:,:)
    !
    INTEGER                         :: LocElemType
    !
    !--------------------------------------------------------------------------
    !
    IF(MPI%nCPU.EQ.1) THEN
        RETURN
    ENDIF

    ALLOCATE( AStar(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) )
    ALLOCATE( BStar(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) )
    ALLOCATE( CStar(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) )
    ALLOCATE( EStar(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) )
    ALLOCATE( FLStar(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) )
    ALLOCATE( FRStar(EQN%nVarTotal,EQN%nVarTotal,DISC%Galerkin%nDegFrMat) )

#ifdef PARALLEL      
#ifdef VTRACE
!      CALL VTBEGIN(MPI%Trace%hFunctions(3), MPI%iErr)
#endif

    ALLOCATE(send_message(BND%NoMPIDomains),        &
             recv_message(BND%NoMPIDomains),        &
             send_imessage(BND%NoMPIDomains),       &
             recv_imessage(BND%NoMPIDomains),       &             
             send_request(BND%NoMPIDomains),        &
             recv_request(BND%NoMPIDomains),        &
             send_status_list(MPI_STATUS_SIZE,BND%NoMPIDomains),    &
             recv_status_list(MPI_STATUS_SIZE,BND%noMPIDomains)     )

    DO iDomain = 1, BND%NoMPIDomains
        ! Compute total length of message for each MPI domain interface
        ! The factor 4 is for matrices A, B, C and E
        MsgLength = BND%ObjMPI(iDomain)%nElem * 4 * EQN%nVarTotal * EQN%nVarTotal * DISC%Galerkin%nDegFrMat
        ! Add Element type
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem
        ! Add the flux matrices (solution of Riemann problem)
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem * 2 * EQN%nVarTotal * EQN%nVarTotal * DISC%Galerkin%nDegFrMat *&
		MESH%nSideMax
        ! Allocate space for message
        MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
        ALLOCATE(send_message(iDomain)%Content(MsgLength))
        ALLOCATE(recv_message(iDomain)%Content(MsgLength))
        ! Number of neighbor CPU
        iCPU = BND%ObjMPI(iDomain)%CPU
        counter = 0
        ! Put all information into a 1D array
        DO i = 1, BND%ObjMPI(iDomain)%nElem
            iElem = BND%ObjMPI(iDomain)%DomainElements(i)
            CALL UnpackSparseTensor3(AStar, DISC%Galerkin%AStar_Sp(iElem), EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat)
            CALL UnpackSparseTensor3(BStar, DISC%Galerkin%BStar_Sp(iElem), EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat)
            CALL UnpackSparseTensor3(CStar, DISC%Galerkin%CStar_Sp(iElem), EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat)
            CALL UnpackSparseTensor3(EStar, DISC%Galerkin%EStar_Sp(iElem), EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat)
            ! Send DG degrees of freedom to neighbor domain
            DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                DO q = 1, EQN%nVarTotal
                    DO p = 1, EQN%nVarTotal
                        counter = counter + 1
                        IF(counter.GT.MsgLength) THEN
                          logError(*) 'Send message: Counter > MsgLength', counter, MsgLength
                          STOP
                        ENDIF
                        ! Setup list of all Jacobian matrices
                        send_message(iDomain)%Content(counter) = AStar(p,q,iDegFr)
                        counter = counter + 1
                        send_message(iDomain)%Content(counter) = BStar(p,q,iDegFr)
                        counter = counter + 1
                        send_message(iDomain)%Content(counter) = CStar(p,q,iDegFr)
                        counter = counter + 1
                        send_message(iDomain)%Content(counter) = EStar(p,q,iDegFr)
                    ENDDO
                ENDDO
            ENDDO
            LocElemType = MESH%LocalElemType(iElem)
            ! Send the element type
            counter = counter + 1
            send_message(iDomain)%Content(counter) = LocElemType
            DO iSide=1,LocElemType
                CALL UnpackSparseTensor3(FLStar, DISC%Galerkin%FLStar_Sp(iElem,iSide), &
			EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat)
                CALL UnpackSparseTensor3(FRStar, DISC%Galerkin%FRStar_Sp(iElem,iSide), &
			EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat)
                DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                    DO q = 1, EQN%nVarTotal
                        DO p = 1, EQN%nVarTotal
                            counter = counter + 1
                            IF(counter.GT.MsgLength) THEN
                              logError(*) 'Send message: Counter > MsgLength', counter, MsgLength
                              STOP
                            ENDIF
                            ! Setup list of all Jacobian matrices
                            send_message(iDomain)%Content(counter) = FLStar(p,q,iDegFr)
                            counter = counter + 1
                            send_message(iDomain)%Content(counter) = FRStar(p,q,iDegFr)
                        ENDDO ! p
                    ENDDO ! q
                ENDDO ! iDegFr
            ENDDO ! iSide
        ENDDO
        ! Post a send request to neighbor CPU
        CALL MPI_ISEND(send_message(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                     send_request(iDomain), iErr)
        !
    ENDDO ! iDomain
    
    DO iDomain = 1, BND%NoMPIDomains
        ! Compute total length of message for each MPI domain interface
        ! The factor 4 is for matrices A, B, C and E
        MsgLength = BND%ObjMPI(iDomain)%nElem * 4 * EQN%nVarTotal * EQN%nVarTotal * DISC%Galerkin%nDegFrMat
        ! Add Element type
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem
        ! Add the flux matrices (solution of Riemann problem)
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem * 2 * EQN%nVarTotal * EQN%nVarTotal * DISC%Galerkin%nDegFrMat * &
		MESH%nSideMax
        !
        MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
        ! Number of neighbor CPU
        iCPU = BND%ObjMPI(iDomain)%CPU
        ! Post a receive request from neighbor CPU
        CALL MPI_IRECV(recv_message(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                       recv_request(iDomain), iErr)
    ENDDO ! iDomain

    ! Wait until all communication has finished
    CALL MPI_WAITALL(BND%NoMPIDomains,send_request,send_status_list,ierr)
    CALL MPI_WAITALL(BND%NoMPIDomains,recv_request,recv_status_list,ierr)

    ! Decode 1D information array obtained from neighbor CPUs
    DO iDomain = 1, BND%NoMPIDomains
      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
        ! The factor 4 is for matrices A, B, C and E
        MsgLength = BND%ObjMPI(iDomain)%nElem * 4 * EQN%nVarTotal * EQN%nVarTotal * DISC%Galerkin%nDegFrMat
        ! Add Element type
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem
        ! Add the flux matrices (solution of Riemann problem)
        MsgLength = MsgLength + BND%ObjMPI(iDomain)%nElem * 2 * EQN%nVarTotal * EQN%nVarTotal * DISC%Galerkin%nDegFrMat * &
		MESH%nSideMax
        !
        counter = 0
        ! Get all information from the 1D array
        DO i = 1, BND%ObjMPI(iDomain)%nElem
            iElem = BND%ObjMPI(iDomain)%DomainElements(i)
            ! Send DG degrees of freedom to neighbor domain
            DO iDegFr = 1, DISC%Galerkin%nDegFrMat
              DO q = 1, EQN%nVarTotal
                DO p = 1, EQN%nVarTotal
                  counter = counter + 1
                  IF(counter.GT.MsgLength) THEN
                    logError(*) 'Recv message: Counter > MsgLength', counter, MsgLength
                    STOP
                  ENDIF
                  ! Decode list of all Jacobian matrices
                  AStar(p,q,iDegFr) = recv_message(iDomain)%content(counter)
                  counter = counter + 1
                  BStar(p,q,iDegFr) = recv_message(iDomain)%content(counter)
                  counter = counter + 1
                  CStar(p,q,iDegFr) = recv_message(iDomain)%content(counter)
                  counter = counter + 1
                  EStar(p,q,iDegFr) = recv_message(iDomain)%content(counter)
                ENDDO
              ENDDO
            ENDDO
            CALL IniSparseTensor3(BND%ObjMPI(iDomain)%AStar_Sp(i), AStar, EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat) 
            CALL IniSparseTensor3(BND%ObjMPI(iDomain)%BStar_Sp(i), BStar, EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat) 
            CALL IniSparseTensor3(BND%ObjMPI(iDomain)%CStar_Sp(i), CStar, EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat) 
            CALL IniSparseTensor3(BND%ObjMPI(iDomain)%EStar_Sp(i), EStar, EQN%nVarTotal, EQN%nVarTotal, DISC%Galerkin%nDegFrMat) 
            !
            counter = counter + 1
            LocElemType = recv_message(iDomain)%content(counter)
            DO iSide=1,LocElemType
                DO iDegFr = 1, DISC%Galerkin%nDegFrMat
                    DO q = 1, EQN%nVarTotal
                        DO p = 1, EQN%nVarTotal
                            counter = counter + 1
                            IF(counter.GT.MsgLength) THEN
                                !logError(*), 'Recv message: Counter > MsgLength', counter, MsgLength
                                STOP
                            ENDIF
                            ! Decode list of all flux matrices
                            FLStar(p,q,iDegFr) = recv_message(iDomain)%content(counter)
                            counter = counter + 1
                            FRStar(p,q,iDegFr) = recv_message(iDomain)%content(counter)
                        ENDDO ! p
                    ENDDO ! q
                ENDDO ! iDegFr
                CALL IniSparseTensor3(BND%ObjMPI(iDomain)%FLStar_Sp(i,iSide), FLStar, EQN%nVarTotal, &
			EQN%nVarTotal, DISC%Galerkin%nDegFrMat) 
                CALL IniSparseTensor3(BND%ObjMPI(iDomain)%FRStar_Sp(i,iSide), FRStar, EQN%nVarTotal, &
			EQN%nVarTotal, DISC%Galerkin%nDegFrMat) 
            ENDDO ! iSide
        ENDDO ! i
        BND%ObjMPI(iDomain)%Init = .TRUE.
    ENDDO ! iDomain
    
    ! Deallocate all communication variables
    DO iDomain = 1, BND%NoMPIDomains
        DEALLOCATE(send_message(iDomain)%Content)
        DEALLOCATE(recv_message(iDomain)%Content)
    ENDDO ! iDomain
    !
    DEALLOCATE(send_message,       &
               recv_message,       &
               send_imessage,      &
               recv_imessage,      &
               send_request,       &
               recv_request,       &
               send_status_list,   &
               recv_status_list    )

#ifdef VTRACE
!      CALL VTEND(MPI%Trace%hFunctions(3), MPI%iErr)
#endif
#endif

    DEALLOCATE( AStar )
    DEALLOCATE( BStar )
    DEALLOCATE( CStar )
    DEALLOCATE( EStar )
    DEALLOCATE( FLStar )
    DEALLOCATE( FRStar )

    CONTINUE

  END SUBROUTINE MPIExchangeJacobians_new


#ifdef PARALLEL
  SUBROUTINE MPIExchangeInvSystems(DISC,EQN,BND,MESH,IO,OptionalFields,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------

    INCLUDE 'mpif.h'
#ifdef VTRACE
!    INCLUDE 'VT.inc'
#endif

    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tDiscretization)          :: DISC
    TYPE (tEquations)               :: EQN
    TYPE (tBoundary)                :: BND
    TYPE (tUnstructMesh)            :: MESH
    TYPE (tInputOutput)             :: IO
    TYPE (tUnstructOptionalFields)  :: OptionalFields
    TYPE (tMPI)                     :: MPI
    ! Local variable declaration
    INTEGER                         :: iDomain, iCPU, iElem
    INTEGER                         :: iNeighbor, iObject, MPIIndex
    INTEGER                         :: iVar, iMech, i, j, NZ
    INTEGER                         :: count_in, count_re, MsgLengthInt, MsgLengthReal
    INTEGER                         :: iErr
    TYPE(tRealMessage), ALLOCATABLE :: send_message(:)
    TYPE(tRealMessage), ALLOCATABLE :: recv_message(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: send_imessage(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: recv_imessage(:)
    INTEGER, ALLOCATABLE            :: send_request(:)
    INTEGER, ALLOCATABLE            :: recv_request(:)
    INTEGER, ALLOCATABLE            :: send_status_list(:,:)
    INTEGER, ALLOCATABLE            :: recv_status_list(:,:)
    INTEGER, PARAMETER              :: MIN_MSG_LENGTH = 64
    INTEGER                         :: nNonZero(DISC%Galerkin%nDegFrST*EQN%nVarTotal)
    REAL                            :: NonZero(DISC%Galerkin%InvSyst_MaxnNonZeros)
    INTEGER                         :: NonZeroIndex1(DISC%Galerkin%InvSyst_MaxnNonZeros)
    INTEGER                         :: NonZeroIndex2(DISC%Galerkin%InvSyst_MaxnNonZeros)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,MESH,IO,MPI
    INTENT(INOUT)             :: BND,DISC
    !--------------------------------------------------------------------------

    ! The InvSystemMatrices will be passed to the neighbours. Two messages will be sent: one real and one integer.
    ! If a non-porous element is neighbour, then no InvSystemMatrix is built, although a full message is sent.
    IF(MPI%nCPU.EQ.1) THEN
        RETURN
    ENDIF
    IF(DISC%Galerkin%CKMethod.NE.1) THEN
        RETURN
    ENDIF

    ALLOCATE(send_message(BND%NoMPIDomains),         &
             recv_message(BND%NoMPIDomains),         &
             send_imessage(BND%NoMPIDomains),        &
             recv_imessage(BND%NoMPIDomains),        &
             send_request(BND%NoMPIDomains),         &
             recv_request(BND%NoMPIDomains),         &
             send_status_list(MPI_STATUS_SIZE,BND%NoMPIDomains),    &
             recv_status_list(MPI_STATUS_SIZE,BND%noMPIDomains)     )

    logInfo(*)  'Exchanging Inverse System values between processors... '

    DO iDomain = 1, BND%NoMPIDomains
      ! Compute total length of messages for each MPI domain interface
      MsgLengthInt   = BND%ObjMPI(iDomain)%nElem *(DISC%Galerkin%InvSyst_MaxnNonZeros * 2 + DISC%Galerkin%nDegFrST*EQN%nVarTotal)
      MsgLengthReal  = BND%ObjMPI(iDomain)%nElem * DISC%Galerkin%InvSyst_MaxnNonZeros
      ALLOCATE(send_message(iDomain)%Content(MsgLengthReal))
      ALLOCATE(recv_message(iDomain)%Content(MsgLengthReal))
      ALLOCATE(send_imessage(iDomain)%Content(MsgLengthInt))
      ALLOCATE(recv_imessage(iDomain)%Content(MsgLengthInt))

      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
  
      count_in = 0
      count_re = 0
      ! Put all information into a 1D array
      DO i = 1, BND%ObjMPI(iDomain)%nElem
        iElem = BND%ObjMPI(iDomain)%DomainElements(i)
        NonZero(:)       = 0.
        NonZeroIndex1(:) = 0
        NonZeroIndex2(:) = 0
        IF(EQN%LocPoroelastic(iElem).NE.0) THEN
          nNonZero(:)         = DISC%Galerkin%InvSystemMatrix(iElem)%nNonZero(:)
          nZ                  = nNonZero(DISC%Galerkin%nDegFrST*EQN%nVarTotal)
          NonZero(1:nZ)       = DISC%Galerkin%InvSystemMatrix(iElem)%NonZero(1:nZ)
          NonZeroIndex1(1:nZ) = DISC%Galerkin%InvSystemMatrix(iElem)%NonZeroIndex1(1:nZ)
          NonZeroIndex2(1:nZ) = DISC%Galerkin%InvSystemMatrix(iElem)%NonZeroIndex2(1:nZ)
        ELSE
          nNonZero(:)      = 0
          nZ               = 1
          NonZero(:)       = 0.
          NonZeroIndex1(:) = 0
          NonZeroIndex2(:) = 0
        ENDIF
        !
        ! Pack inverse system values in 1D arrays to be sent to neighbor domain
        DO j=1,DISC%Galerkin%InvSyst_MaxnNonZeros
          count_re = count_re+1
          send_message(iDomain)%Content(count_re)   = NonZero(j)
        ENDDO
        DO j=1,DISC%Galerkin%nDegFrST*EQN%nVarTotal
          count_in=count_in+1
          send_imessage(iDomain)%Content(count_in)  = nNonZero(j)
        ENDDO
        DO j=1,DISC%Galerkin%InvSyst_MaxnNonZeros
          count_in=count_in+1
          send_imessage(iDomain)%Content(count_in)  = NonZeroIndex1(j)
        ENDDO
        DO j=1,DISC%Galerkin%InvSyst_MaxnNonZeros
          count_in=count_in+1
          send_imessage(iDomain)%Content(count_in)  = NonZeroIndex2(j)
        ENDDO
      ENDDO
      ! Post a send request to neighbor CPU (INTEGER VALUES ONLY: nNonZeros + NonZeroIndex1 + NonZeroIndex2)
      CALL MPI_ISEND(send_imessage(iDomain)%Content, MsgLengthInt,  MPI_INTEGER,   iCPU, 1, &
		MPI%commWorld,send_request(iDomain), iErr)
      !
    ENDDO

    DO iDomain = 1, BND%NoMPIDomains
      iCPU = BND%ObjMPI(iDomain)%CPU
      MsgLengthInt   = BND%ObjMPI(iDomain)%nElem *(DISC%Galerkin%InvSyst_MaxnNonZeros * 2 + DISC%Galerkin%nDegFrST*EQN%nVarTotal)
      CALL MPI_IRECV(recv_imessage(iDomain)%Content, MsgLengthInt, MPI_INTEGER, iCPU, 1, MPI%commWorld, &
		recv_request(iDomain), iErr)
    ENDDO

    ! Wait until all communication has finished
    CALL MPI_WAITALL(BND%NoMPIDomains,send_request,send_status_list,ierr)
    CALL MPI_WAITALL(BND%NoMPIDomains,recv_request,recv_status_list,ierr)

    DO iDomain = 1, BND%NoMPIDomains
      iCPU = BND%ObjMPI(iDomain)%CPU
      MsgLengthReal  = BND%ObjMPI(iDomain)%nElem * DISC%Galerkin%InvSyst_MaxnNonZeros
      CALL MPI_ISEND(send_message(iDomain)%Content,  MsgLengthReal,  MPI%MPI_AUTO_REAL, iCPU, 1, &
		MPI%commWorld,send_request(iDomain), iErr)
    ENDDO

    DO iDomain = 1, BND%NoMPIDomains
      iCPU = BND%ObjMPI(iDomain)%CPU
      MsgLengthReal  = BND%ObjMPI(iDomain)%nElem * DISC%Galerkin%InvSyst_MaxnNonZeros
      CALL MPI_IRECV(recv_message(iDomain)%Content, MsgLengthReal, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
		recv_request(iDomain), iErr)
    ENDDO

    ! Wait until all communication has finished
    CALL MPI_WAITALL(BND%NoMPIDomains,send_request,send_status_list,ierr)
    CALL MPI_WAITALL(BND%NoMPIDomains,recv_request,recv_status_list,ierr)

    ! Decode 1D information array obtained from neighbor CPUs
    DO iDomain = 1, BND%NoMPIDomains
      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
      count_in = 0
      count_re = 0
      ! Get all information from the 1D array
      DO i = 1, BND%ObjMPI(iDomain)%nElem
        iElem = BND%ObjMPI(iDomain)%DomainElements(i)
        ! Read in all information for one element
        DO j=1,DISC%Galerkin%nDegFrST*EQN%nVarTotal
          count_in = count_in+1
          nNonZero(j) = recv_imessage(iDomain)%Content(count_in)
        ENDDO
        DO j=1,DISC%Galerkin%InvSyst_MaxnNonZeros
          count_in = count_in+1
          NonZeroIndex1(j) = recv_imessage(iDomain)%Content(count_in)
        ENDDO
        DO j=1,DISC%Galerkin%InvSyst_MaxnNonZeros
          count_in = count_in+1
          NonZeroIndex2(j) = recv_imessage(iDomain)%Content(count_in)
        ENDDO
        DO j=1,DISC%Galerkin%InvSyst_MaxnNonZeros
          count_re = count_re+1
          NonZero(j) = recv_message(iDomain)%Content(count_re)
        ENDDO

        !Number of non-zeros for particular element
        NZ=nNonZero(DISC%Galerkin%nDegFrST*EQN%nVarTotal)

        IF(NZ.EQ.0) THEN !Is a non-porous element, no InvSystemMatrix has to be built!

          CONTINUE

        ELSE !Is a porous element, initialize and save InvSystemMatrix

          ALLOCATE( BND%ObjMPI(iDomain)%InvSystemMatrix(i)%nNonZero(DISC%Galerkin%nDegFrST*EQN%nVarTotal) )       
          ALLOCATE( BND%ObjMPI(iDomain)%InvSystemMatrix(i)%NonZero(      NZ)  )
          ALLOCATE( BND%ObjMPI(iDomain)%InvSystemMatrix(i)%NonZeroIndex1(NZ)  )
          ALLOCATE( BND%ObjMPI(iDomain)%InvSystemMatrix(i)%NonZeroIndex2(NZ)  )

          BND%ObjMPI(iDomain)%InvSystemMatrix(i)%nNonZero(:)       = nNonZero(:)
          BND%ObjMPI(iDomain)%InvSystemMatrix(i)%NonZero(:)        = NonZero(1:NZ)
          BND%ObjMPI(iDomain)%InvSystemMatrix(i)%NonZeroIndex1(:)  = NonZeroIndex1(1:NZ) 
          BND%ObjMPI(iDomain)%InvSystemMatrix(i)%NonZeroIndex2(:)  = NonZeroIndex2(1:NZ)

        ENDIF

      ENDDO
    ENDDO

    logInfo(*)  'All Inverse Systems successfully exchanged and decoded! '

    ! Deallocate all communication variables
    DO iDomain = 1, BND%NoMPIDomains
      DEALLOCATE(send_message(iDomain)%Content)
      DEALLOCATE(recv_message(iDomain)%Content)
    ENDDO
    !
    DEALLOCATE(send_message,       &
               recv_message,       &
               send_imessage,      &
               recv_imessage,      &
               send_request,       &
               recv_request,       &
               send_status_list,   &
               recv_status_list    )


  END SUBROUTINE MPIExchangeInvSystems

#else
  SUBROUTINE MPIExchangeInvSystems
  END SUBROUTINE MPIExchangeInvSystems
#endif

! GENERATEDKERNELS
#endif

  ! Exchange background values
  SUBROUTINE MPIExchangeBackground(DISC,EQN,BND,MESH,IO,OptionalFields,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#ifdef VTRACE
!    INCLUDE 'VT.inc'
#endif
#endif
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tDiscretization)          :: DISC
    TYPE (tEquations)               :: EQN
    TYPE (tBoundary)                :: BND
    TYPE (tUnstructMesh)            :: MESH
    TYPE (tInputOutput)             :: IO
    TYPE (tUnstructOptionalFields)  :: OptionalFields
    TYPE (tMPI)                     :: MPI
    ! Local variable declaration
    INTEGER                         :: iDomain, iCPU, iElem, iSide, iDegFr
    INTEGER                         :: iNeighbor, iObject, MPIIndex
    INTEGER                         :: iVar, iMech, i,j,k
    INTEGER                         :: counter, MsgLength
    INTEGER                         :: UpdateElement
    INTEGER                         :: iErr
    TYPE(tRealMessage), ALLOCATABLE :: send_message(:)
    TYPE(tRealMessage), ALLOCATABLE :: recv_message(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: send_imessage(:)
    TYPE(tIntegerMessage), ALLOCATABLE :: recv_imessage(:)
    INTEGER, ALLOCATABLE            :: send_request(:)
    INTEGER, ALLOCATABLE            :: recv_request(:)
    INTEGER, ALLOCATABLE            :: send_status_list(:,:)
    INTEGER, ALLOCATABLE            :: recv_status_list(:,:)
    REAL                            :: DuDt(DISC%Galerkin%nDegFr,EQN%nVarTotal)
    REAL                            :: tmin,timeMPI,dtmax,dtmaxMPI
    INTEGER, PARAMETER              :: MIN_MSG_LENGTH = 64
    !--------------------------------------------------------------------------
    INTENT(IN)                :: EQN,MESH,IO,MPI
    INTENT(INOUT)             :: BND,DISC
    !--------------------------------------------------------------------------
    IF(MPI%nCPU.EQ.1) THEN
        RETURN
    ENDIF

#ifdef PARALLEL      
#ifdef VTRACE
!      CALL VTBEGIN(MPI%Trace%hFunctions(3), MPI%iErr)
#endif

    ALLOCATE(send_message(BND%NoMPIDomains),        &
             recv_message(BND%NoMPIDomains),        &
             send_imessage(BND%NoMPIDomains),       &
             recv_imessage(BND%NoMPIDomains),       &             
             send_request(BND%NoMPIDomains),        &
             recv_request(BND%NoMPIDomains),        &
             send_status_list(MPI_STATUS_SIZE,BND%NoMPIDomains),    &
             recv_status_list(MPI_STATUS_SIZE,BND%noMPIDomains)     )

    DO iDomain = 1, BND%NoMPIDomains
      ! Compute total length of message for each MPI domain interface
      ! Length = No. of MPI boundary Elements * ( No. of background values )
      MsgLength = BND%ObjMPI(iDomain)%nElem*(EQN%nBackgroundVar)
      ! Allocate space for message
      MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
      ALLOCATE(send_message(iDomain)%Content(MsgLength))
      ALLOCATE(recv_message(iDomain)%Content(MsgLength))

      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
      counter = 0
      ! Put all information into a 1D array
      DO i = 1, BND%ObjMPI(iDomain)%nElem
        iElem = BND%ObjMPI(iDomain)%DomainElements(i)
        ! Send elastic parameters (rho, mu, lambda) and anelastic coefficients (omega_n, YP_n, YS_n) to neighbor domain
        DO iVar = 1, EQN%nBackgroundVar 
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              PRINT *, 'Fatal error: Counter > MsgLength', counter, MsgLength, ' in CPU : ', MPI%myrank
              STOP
            ENDIF
            send_message(iDomain)%Content(counter) = OptionalFields%BackgroundValue(iElem,iVar)
        ENDDO
      ENDDO
      ! Post a send request to neighbor CPU
      CALL MPI_ISEND(send_message(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                     send_request(iDomain), iErr)
    ENDDO

    
    DO iDomain = 1, BND%NoMPIDomains
      ! Compute total length of message for each MPI domain interface
      ! Length = No. of MPI boundary Elements * ( No. of background values )
      MsgLength = BND%ObjMPI(iDomain)%nElem*(EQN%nBackgroundVar)
      MsgLength = MAX(MIN_MSG_LENGTH,MsgLength)
      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
      ! Post a receive request from neighbor CPU
      CALL MPI_IRECV(recv_message(iDomain)%Content, MsgLength, MPI%MPI_AUTO_REAL, iCPU, 1, MPI%commWorld, &
                     recv_request(iDomain), iErr)
    ENDDO

    ! Wait until all communication has finished
    CALL MPI_WAITALL(BND%NoMPIDomains,send_request,send_status_list,ierr)
    CALL MPI_WAITALL(BND%NoMPIDomains,recv_request,recv_status_list,ierr)

    ! Decode 1D information array obtained from neighbor CPUs
    DO iDomain = 1, BND%NoMPIDomains
      ! Number of neighbor CPU
      iCPU = BND%ObjMPI(iDomain)%CPU
      MsgLength = BND%ObjMPI(iDomain)%nElem*(EQN%nBackgroundVar)
      counter = 0
      ! Get all information from the 1D array
      DO i = 1, BND%ObjMPI(iDomain)%nElem
        iElem = BND%ObjMPI(iDomain)%DomainElements(i)
        ! Decode elastic parameters (rho, mu, lambda) and anelastic coefficients (omega_n, YP_n, YS_n)  
        DO iVar = 1, EQN%nBackgroundVar 
            counter = counter + 1
            IF(counter.GT.MsgLength) THEN
              PRINT *, 'Fatal error: Counter > MsgLength', counter, MsgLength, ' in CPU : ', MPI%myrank
              STOP
            ENDIF
            BND%ObjMPI(iDomain)%NeighborBackground(iVar,i) = recv_message(iDomain)%Content(counter) 
        ENDDO
      ENDDO
    ENDDO
    
    ! Deallocate all communication variables
    DO iDomain = 1, BND%NoMPIDomains
      DEALLOCATE(send_message(iDomain)%Content)
      DEALLOCATE(recv_message(iDomain)%Content)
    ENDDO
    !
    !
    DEALLOCATE(send_message,       &
               recv_message,       &
               send_imessage,      &
               recv_imessage,      &
               send_request,       &
               recv_request,       &
               send_status_list,   &
               recv_status_list    )

#ifdef VTRACE
!      CALL VTEND(MPI%Trace%hFunctions(3), MPI%iErr)
#endif
#endif

  END SUBROUTINE MPIExchangeBackground


END MODULE MPIExchangeValues_mod
