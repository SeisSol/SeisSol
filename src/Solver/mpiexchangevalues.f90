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
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: MPIExchangeBackground
  !----------------------------------------------------------------------------
CONTAINS
  ! Exchange background values
  SUBROUTINE MPIExchangeBackground(DISC,EQN,BND,MESH,IO,OptionalFields,MPI)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
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
