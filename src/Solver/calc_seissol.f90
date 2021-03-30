!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
!! Computational entry point.
!!

#include <Initializer/preProcessorMacros.fpp>

MODULE calc_SeisSol_mod
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE calc_SeisSol
     MODULE PROCEDURE calc_SeisSol
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: calc_SeisSol  
  !----------------------------------------------------------------------------
  
CONTAINS

  SUBROUTINE calc_SeisSol(time,timestep,pvar,cvar,EQN,MESH,DISC,SOURCE,BND,IC, &
                         OptionalFields,IO,MPI)
    !--------------------------------------------------------------------------
    USE TypesDef
    USE dg_setup_mod
    USE calc_deltaT_mod
#ifdef HDF
    USE receiver_hdf_mod
#else
    USE receiver_mod
    USE energies_output_mod
#endif
    USE ini_SeisSol_mod
    USE magnitude_output_mod
    USE output_rupturefront_mod
    USE COMMON_operators_mod
#ifdef PARALLEL
    USE MPIExchangeValues_mod
#endif
    use iso_c_binding, only: c_loc
    use monitoring
    use f_ftoc_bind_interoperability

    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tEquations)             :: EQN 
    TYPE (tUnstructMesh)          :: MESH 
    TYPE (tDiscretization)        :: DISC
    TYPE (tSource)                :: SOURCE 
    TYPE (tBoundary)              :: BND 
    TYPE (tInitialCondition)      :: IC
    TYPE (tUnstructOptionalFields):: OptionalFields
    TYPE (tInputOutput)           :: IO
    TYPE (tMPI)                   :: MPI

    REAL                          :: time                 ! current time
    INTEGER                       :: timestep             ! index of time step
    REAL, POINTER                 :: pvar(  : , :  ) 
    REAL, POINTER                 :: cvar(  : , :  )     
    ! local variables declaration
    INTEGER                       :: k_index(4)
    INTEGER                       :: plotNr, minl(1), iter
    INTEGER                       :: iElem, iSide, i, iObject, MPIIndex
    INTEGER(KIND=4)               :: now(3)
    REAL                          :: tol, InitBegin, InitEnd, LoopBegin, LoopEnd
    REAL                          :: MPITime, MPIdt, dt
    REAL                          :: CPUTime
    REAL                          :: AbortStatus, MPIAbortStatus
    INTEGER                       :: iCPU, iErr
    CHARACTER(LEN=60)             :: separator
    LOGICAL                       :: printThisTimeStep
    LOGICAL                       :: ContinueLoop
    INTEGER                       :: nNonZeros(MESH%nElem)
    INTEGER                       :: MaxnNonZeros, LocnNonZeros
    REAL                          :: JT(3,3,10),DIF,x(4),y(4),z(4),JacobiT(3,3)
    INTEGER                       :: counter, counter2,iType,j, countside, UpdType
    REAL(KIND=8)                  :: l_synchronizationPoint;
    INTEGER                       :: iDRupdate
    !--------------------------------------------------------------------------
    INTENT(INOUT)                 :: time,timestep,OptionalFields,EQN,IO,BND,DISC,MESH,SOURCE

    ! register epik/scorep region for time stepping loop
    EPIK_USER_REG(r_time_stepping, "time stepping loop")
    SCOREP_USER_REGION_DEFINE(r_time_stepping)

    ! register epik/scorep region for the communication of the degrees of freedom
    EPIK_USER_REG(r_dofs_communication, "DOFs communication")
    SCOREP_USER_REGION_DEFINE(r_dofs_communication)

    !--------------------------------------------------------------------------
    PARAMETER(separator='<--------------------------------------------------------->')
    !--------------------------------------------------------------------------
    !                                                   !
    !                                                   !
    ! For the first iteration, set printtime to time    !
    DISC%printtime = time+IO%outInterval%TimeInterval   !
    IO%picktime  = time
    IO%picktime_energy = time
    !                                                   !
    logInfo(*) separator                                !
    !                                                   !
    !
    !
    IF (timestep.EQ.0) THEN                             ! Write initial 
       !                                                ! condition       
       logInfo(*) separator                !
       logInfo(*) ' '                      !
       logInfo(*) separator                !
       logInfo(*) '<                   MONITORING TIME STEP                  >'
       logInfo(*) separator                !
       !
    ENDIF

#ifdef USE_MPI
       ! sync all redundant data from the initialization
       call c_interoperability_synchronizeCopyLayerDofs()
#endif

#ifdef PARALLEL
    ! Sync all processes for a concurrent start
    CALL MPI_BARRIER(MPI%commWorld,iErr)
#endif
    DISC%StartCPUTime = dwalltime()

    ! enable dynamic rupture if requested
    if( eqn%dr==1 ) then
      call c_interoperability_enableDynamicRupture()
    endif

    ! check whether the device memory allocated at this point
    ! exceeds the maximum avaliable on a current device
    call c_interoperability_report_device_memory_status()

    ! do the simulation
    call c_interoperability_simulate( i_finalTime = disc%endTime );
    ! End time is currently the only supported abort criteria by GK
    time = disc%endTime
!no generated kernel

    !---------------------------------------------------!
    !---------------- STOP TIME LOOP -------------------!
    !---------------------------------------------------!
    !                                                   !
    ! Get CPU time after integration of the PDE and measure the time that was  
    ! spent in this iteration loop and, if called by KOP, in the ones before.

    ! end epik/scorep region for time stepping loop
    EPIK_USER_END(r_time_stepping)
    SCOREP_USER_REGION_END(r_time_stepping)

#ifdef PARALLEL
    !Ensure that all processes are done with their calculation
    CALL MPI_BARRIER(MPI%commWorld,iErr)
#endif
    DISC%StopCPUTime = dwalltime()

    DISC%LoopCPUTime = DISC%StopCPUTime - DISC%StartCPUTime

    !
    EQN%FinalTime = time !Actual time at which simualation ends, used for time-reversing sources
    EQN%FinalPick = IO%picktime
    !
    logInfo0(*) 'total number of performed time steps: ', timestep
    logInfo0(*) 'final time of the simulation: ', time

    logInfo(*)'<--------------------------------------------------------->'  !
    !
    IF(IO%PGMLocationsFlag.NE.0)THEN
#ifdef HDF
        CALL PGM_output_hdf(IO,MPI)
#else
        CALL PGM_output(IO,MPI)
#endif
    ENDIF

    ! output magnitude for dynamic rupture simulations
    IF (EQN%DR.EQ.1 .AND. DISC%DynRup%magnitude_output_on.EQ.1) CALL magnitude_output(OptionalFields%BackgroundValue,DISC,MESH,MPI,IO)
    ! output GP-wise RF in extra files
    IF (EQN%DR.EQ.1 .AND. DISC%DynRup%RF_output_on.EQ.1) CALL output_rupturefront(DISC,MESH,MPI,IO, BND)

    logInfo(*)'<--------------------------------------------------------->'  !
    logInfo(*)'<     calc_SeisSol successfully finished                  >'  !
    logInfo(*)'<--------------------------------------------------------->'  !
    !
    IO%AbortStatus = 0
    !
  END SUBROUTINE calc_SeisSol                            

END MODULE calc_SeisSol_mod
