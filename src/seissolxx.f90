!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2010-2016, SeisSol Group
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
!! Entry point of SeisSol
!!

#include <Initializer/preProcessorMacros.fpp>

module SeisSol
  !----------------------------------------------------------------------------
  use parallel_mpi
  USE ini_SeisSol_mod
  USE calc_SeisSol_mod
  USE close_SeisSol_mod
  USE inioutput_SeisSol_mod
  USE TypesDef
  USE COMMON_operators_mod, ONLY: OpenFile

  use iso_c_binding
use f_ftoc_bind_interoperability
!----------------------------------------------------------------------------
IMPLICIT NONE

contains
subroutine main() bind(c, name='fortran_main')
!----------------------------------------------------------------------------
#ifdef PARALLEL
INCLUDE 'mpif.h'
#endif
!----------------------------------------------------------------------------
REAL                                   :: time
INTEGER                                :: timestep
TYPE (tUnstructDomainDescript), target :: domain
CHARACTER(LEN=600)                     :: name
INTEGER                                :: iTry
LOGICAL                                :: fexist
!----------------------------------------------------------------------------

domain%IO%AbortStatus = 2

CALL CPU_TIME(domain%IO%WallStart)

! Set line size for stdout and stderr (removes auto newlines)
! This does no longer work with C++ main() and produces errors with gfortran
#ifdef __INTEL_COMPILER
open(FORTRAN_STDOUT, recl=FORTRAN_LINE_SIZE)
open(FORTRAN_STDERR, recl=FORTRAN_LINE_SIZE)
#endif

#ifdef PARALLEL
   ! Initialize MPI
   domain%MPI%commWorld = getCommWorld()

   CALL MPI_COMM_RANK(domain%MPI%commWorld, domain%MPI%myrank, domain%MPI%iErr)
   CALL MPI_COMM_SIZE(domain%MPI%commWorld, domain%MPI%nCPU,   domain%MPI%iErr)

   ! Set global variable for the rank, required by the logger
   myrank = domain%MPI%myrank

   WRITE(name,'(a,i5.5,a)') 'StdOut',domain%MPI%myrank,'.txt'

   logInfo0(*) '<--------------------------------------------------------->'
   logInfo0(*) '<                SeisSol MPI initialization               >'
   logInfo0(*) '<--------------------------------------------------------->'
   logInfo(*) ' MPI Communication Initialized. '
   logInfo(*) ' I am processor ', domain%MPI%myrank, ' of ', domain%MPI%nCPU, ' total CPUs.'
   ! Initialize MPI_AUTO_REAL
   domain%MPI%real_kind     = KIND(domain%MPI%real_kindtest)
   domain%MPI%integer_kind  = KIND(domain%MPI%integer_kindtest)
   SELECT CASE(domain%MPI%real_kind)
   CASE(4)
     logInfo0(*) 'Single precision used for real.'
     logInfo0(*) 'Setting MPI_AUTODOUBLE to MPI_REAL'
     domain%MPI%MPI_AUTO_REAL = MPI_REAL
   CASE(8)
     logInfo0(*) ' Double precision used for real.'
     logInfo0(*) ' Setting MPI_AUTODOUBLE to MPI_DOUBLE_PRECISION'
     domain%MPI%MPI_AUTO_REAL = MPI_DOUBLE_PRECISION
   CASE DEFAULT
     logError(*) 'Unknown kind ', domain%MPI%real_kind
     STOP
   END SELECT
   logInfo0(*) ' MPI_AUTO_REAL feature initialized. '
   domain%MPI%MPI_AUTO_INTEGER = 0
   !  
   logInfo0(*) ' MPI initialization done. '
   logInfo0('(a)') ' <--------------------------------------------------------->'
   WRITE(domain%IO%ErrorFile,'(a,i5.5,a)') 'IRREGULARITIES.', domain%MPI%myrank, '.log'    ! Name der Datei in die Fehler geschrieben werden
#else
domain%IO%ErrorFile                    = 'IRREGULARITIES.log'         ! Name der Datei in die Fehler geschrieben werden
myrank = 0
#endif

domain%programTitle                    ='SeisSol'                     ! Name des Programms
domain%IO%Path                         = ''                           ! Standardmaessig liegen alle Dateien im lokalen Verzeichnis
domain%IO%OutInterval%PlotNrGiven      = .FALSE.                      ! PlotNumbers are created by SeisSol

domain%IO%UNIT%FileIn                  = 104
domain%IO%UNIT%FileOut                 = 105
!  domain%IO%UNIT%PMLtest                 = 108
domain%IO%UNIT%other01                 = 110
domain%IO%UNIT%other02                 = 111
domain%IO%UNIT%other03                 = 112
domain%IO%UNIT%other04                 = 113
domain%IO%UNIT%other05                 = 114
domain%IO%UNIT%FileIn_Fault            = 115                          ! Input file for fault parameters
domain%IO%UNIT%FileOut_Tet             = 118
domain%IO%UNIT%FileOut_Hex             = 119
domain%IO%UNIT%receiverStart           = 120                          ! Ab hier koennen Unitnumbers fuer receiver
!                                                                     ! vergeben werden
domain%IO%UNIT%maxThisDom              = 69999                        ! Obere Grenze fuer Unitnumber fuer receiver

  if (domain%MPI%myrank .eq. 0) then
    WRITE(*,*) 'INFORMATION: The assumed unit number is', FORTRAN_STDOUT, 'for stdout and', FORTRAN_STDERR, 'for stderr.'
    WRITE(*,*) '             If no information follows, please change the value.'
  endif

  name = TRIM(domain%IO%ErrorFile)

  DO iTry = 1, 999
      INQUIRE(                                                                  & 
              EXIST = fexist                                                 ,  & 
              FILE  = name                                                      )
      IF(fexist) THEN
#ifdef PARALLEL
         WRITE(domain%IO%ErrorFile,'(a,i3.3,a,i5.5,a)') 'IRREGULARITIES_',iTry,'.', domain%MPI%myrank, '.log'
#else
         WRITE(domain%IO%ErrorFile,'(a,i3.3,a)') 'IRREGULARITIES_',iTry,'.log'
#endif        
         name = TRIM(domain%IO%ErrorFile)
      ELSE
         EXIT
      ENDIF 
  ENDDO

  logInfo0(*) '<--------------------------------------------------------->'  !
  logInfo0(*) '<     Start ini_SeisSol ...                               >'  !
  logInfo0(*) '<--------------------------------------------------------->'  !

  !---- related to fault output.
  !TODO: Move out of main program to inioutput and open it only for dynamic rupture case!
  !  if (domain%MPI%myrank.eq.0) then
  !  	call create_pvd_writer(meta_plotter)
  !  	call open_pvd_writer(meta_plotter)
  !  endif
  !----END  related to fault output.

    ! propagate data structures to c
    call c_interoperability_setDomain( i_domain = c_loc(domain) )

  CALL ini_SeisSol(                                &
       time           =        time              , &
       timestep       =        timestep          , &
       pvar           = domain%pvar              , &
       cvar           = domain%cvar              , &
       EQN            = domain%EQN               , &
       IC             = domain%IC                , &
       MESH           = domain%MESH              , &
       MPI            = domain%MPI               , &
       SOURCE         = domain%SOURCE            , &
       DISC           = domain%DISC              , &
       BND            = domain%BND               , &
       OptionalFields = domain%OptionalFields    , &  
       IO             = domain%IO                , &
       programTitle   = domain%programTitle        )
domain%IO%MPIPickCleaningDone = 0

    logInfo0(*) '<--------------------------------------------------------->'  !
    logInfo0(*) '<     Start inioutput_SeisSol ...                         >'  !
    logInfo0(*) '<--------------------------------------------------------->'  !

    CALL inioutput_SeisSol(                          &
         time           =        time              , &
         timestep       =        timestep          , &
         pvar           = domain%pvar              , &
         cvar           = domain%cvar              , &
         EQN            = domain%EQN               , &
         IC             = domain%IC                , &
         MESH           = domain%MESH              , &
         MPI            = domain%MPI               , &
         SOURCE         = domain%SOURCE            , &
         DISC           = domain%DISC              , & 
         BND            = domain%BND               , &
         OptionalFields = domain%OptionalFields    , &  
         IO             = domain%IO                , &
         programTitle   = domain%programTitle        )

    logInfo0(*) '<--------------------------------------------------------->'  !
    logInfo0(*) '<     Start calc_SeisSol ...                              >'  !
    logInfo0(*) '<--------------------------------------------------------->'  !

    CALL calc_SeisSol(                               &
         time           =        time              , &
         timestep       =        timestep          , &
         pvar           = domain%pvar              , & ! @TODO, breuera: remove not used
         cvar           = domain%cvar              , &
         EQN            = domain%EQN               , &
         MESH           = domain%MESH              , &
         DISC           = domain%DISC              , &
         SOURCE         = domain%SOURCE            , &
         BND            = domain%BND               , &
         IC             = domain%IC                , &
         OptionalFields = domain%OptionalFields    , &  
         IO             = domain%IO                , &
         MPI            = domain%MPI)

#ifdef PARALLEL
   CALL MPI_BARRIER(domain%MPI%commWorld,domain%MPI%iErr)
#endif

  logInfo0(*) '<--------------------------------------------------------->'  !
  logInfo0(*) '<     Start close_SeisSol ...                             >'  !
  logInfo0(*) '<--------------------------------------------------------->'  !

  CALL close_SeisSol(                            &
       pvar           = domain%pvar            , &
       cvar           = domain%cvar            , &
       EQN            = domain%EQN             , &
       IC             = domain%IC              , &
       BND            = domain%BND             , &
       MESH           = domain%MESH            , &
       DISC           = domain%DISC            , &
       SOURCE         = domain%SOURCE          , &
       OptionalFields = domain%OptionalFields  , &  
       IO             = domain%IO              , &
       MPI            = domain%MPI               )

#ifdef PARALLEL
  ! Finalize MPI 

   CALL MPI_BARRIER(domain%MPI%commWorld,domain%MPI%iErr)
#endif
  end subroutine main

END module SeisSol
