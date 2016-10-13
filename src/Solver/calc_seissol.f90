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
                         OptionalFields,IO,MPI,Analyse)
    !--------------------------------------------------------------------------
    USE TypesDef
    USE COMMON_printThisTimeStep_mod
    USE dg_setup_mod
    USE calc_deltaT_mod
#ifdef HDF
    USE receiver_hdf_mod
#else
    USE receiver_mod
    USE energies_output_mod
#endif
    USE data_output_mod
    USE ini_SeisSol_mod
    USE Galerkin3D_solver_mod
    USE magnitude_output_mod
    USE output_rupturefront_mod
    USE COMMON_operators_mod
#ifdef PARALLEL
    USE MPIExchangeValues_mod
#endif
#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc
    use monitoring
    use f_ftoc_bind_interoperability
#endif

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
    TYPE (tAnalyse)               :: Analyse

    REAL                          :: time                 ! current time
    INTEGER                       :: timestep             ! index of time step
    REAL, POINTER                 :: pvar(  : , :  ) 
    REAL, POINTER                 :: cvar(  : , :  )     
    ! local variables declaration
    INTEGER                       :: k_index(4)
    INTEGER                       :: plotNr, minl(1), nIter, iter
    INTEGER                       :: iElem, iSide, i, iObject, MPIIndex
    INTEGER*4                     :: now(3)
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
    INTEGER                       :: counter, counter2,iType,j, countside, cnt, UpdType
#ifdef GENERATEDKERNELS
    real*8 :: l_synchronizationPoint;
    INTEGER                       :: iDRupdate
#endif
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
#ifdef PARALLEL

#ifndef GENERATEDKERNELS
    logInfo0(*) 'communicating star matrices and flux solvers in sparse format' 
    CALL MPIExchangeJacobians_new(DISC,EQN,BND,MESH,IO,OptionalFields,MPI)
#endif
    ! initialize counters for load balance
    cnt   = 0
    nIter = 10

#ifndef GENERATEDKERNELS
    IF (DISC%Galerkin%DGMethod.NE.3) THEN
            CALL MPIExchangeValues_GTS_Init(DISC           = DISC,             &
                                            EQN            = EQN,              &
                                            BND            = BND,              &
                                            MPI            = MPI)

      logInfo0(*) 'GTS communicaition was initialized'
    ENDIF
#endif

#endif

#ifdef GENERATEDKERNELS

#ifdef USE_MPI
       ! sync all redundant data from the initialization
       call c_interoperability_synchronizeCopyLayerDofs()
#endif

#endif

#ifdef PARALLEL
    ! Sync all processes for a concurrent start
    CALL MPI_BARRIER(MPI%commWorld,iErr)
#endif
    DISC%StartCPUTime = dwalltime()

#ifdef GENERATEDKERNELS
    ! enable dynamic rupture if requested
    if( eqn%dr==1 ) then
      call c_interoperability_enableDynamicRupture()
    endif

    ! do the simulation
    call c_interoperability_simulate( i_finalTime = disc%endTime );
    ! End time is currently the only supported abort criteria by GK
    time = disc%endTime
#else
    !
    !                                                   !
    !---------------------------------------------------!
    !----------------START TIME LOOP -------------------!
    !---------------------------------------------------!
    !
    !
    IF((timestep.LT.DISC%MaxIteration).AND.(time.LT.DISC%EndTime)) THEN
        ContinueLoop = .TRUE.
    ELSE
        ContinueLoop = .FALSE.
    ENDIF

    call itime(now)
    logInfo0('(A,I2,A,I2,A,I2)') ' iterating over time loop             system-time: ', now(1), ':', now(2), ':', now(3)

    ! start epik/scorep region for time stepping loop
    EPIK_USER_START(r_time_stepping)
    SCOREP_USER_REGION_BEGIN(r_time_stepping, "time_stepping", SCOREP_USER_REGION_TYPE_COMMON)
    !                                                   !
    DO  WHILE(ContinueLoop)                             ! Time marching GO
     !
     DISC%time = time                                        !

     ! compute the maximum allowed time step width acoording to the CFL-condition for this rank.
     CALL calc_deltaT(                              & ! calc_deltaT
          OptionalFields = OptionalFields         , & ! calc_deltaT
          EQN            = EQN                    , & ! calc_deltaT
          MESH           = MESH                   , & ! calc_deltaT
          DISC           = DISC                   , & ! calc_deltaT
          SOURCE         = SOURCE                 , & ! calc_deltaT
          IO             = IO                     , & ! calc_deltaT
          time           = time                   , & ! calc_deltaT
          printTime      = DISC%printTime           ) ! calc_deltaT   

#ifdef PARALLEL
     IF(MPI%nCPU.GT.1) THEN
        ! get the minimum allowed time step width over all ranks (unstructured grid, so dt may differ in domains)
        CALL MPI_ALLREDUCE(OptionalFields%dt(1),MPIdt,1,MPI%MPI_AUTO_REAL,MPI_MIN,MPI%commWorld,iErr)
        OptionalFields%dt(:) = MPIdt
        !
     ENDIF
#endif

       ! write receiver and energy output
     if(DISC%Galerkin%DGMethod.ne.3) then
#ifdef HDF
          CALL receiver_hdf(                     &
               EQN       = EQN                  ,&
               MESH      = MESH                 ,&
               DISC      = DISC                 ,&
               MPI       = MPI                  ,&
               IO        = IO                   ,&
               time_op   = time                 ,&
               dt_op     = OptionalFields%dt(1)  )
#else
          CALL receiver(                         &
               EQN       = EQN                  ,&
               MESH      = MESH                 ,&
               DISC      = DISC                 ,&
               MPI       = MPI                  ,&
               IO        = IO                   ,&
               time_op   = time                 ,&
               dt_op     = OptionalFields%dt(1)  )

           ! write energy output (time series)
        IF (IO%energy_output_on .EQ. 1) THEN
           CALL energies_output(                 &
               DISC      = DISC                 ,&
               EQN       = EQN                  ,&
               MESH      = MESH                 ,&
               MPI       = MPI                  ,&
               IO        = IO                   ,&
               time_op   = time                 ,&
               dt_op     = OptionalFields%dt(1)  )
        ENDIF
#endif
     endif

    ! 
    IF(DISC%Galerkin%DGMethod.EQ.3) THEN
        CALL CKProcedureForEveryoneLTSUni(                    &  ! Before MPIExchange in the first timestep,                
            time           = time                           , &  ! do CK procedure for everybody in the local timestepping case 
            dt             = OptionalFields%dt(1)           , &  
            iteration      = timestep                       , &  
            MaterialVal    = OptionalFields%BackgroundValue , &  
            OptionalFields = OptionalFields,                  &  
            EQN            = EQN,                             &  
            MESH           = MESH,                            &  
            DISC           = DISC,                            &  
            IC             = IC,                              &  
            SOURCE         = SOURCE,                          &  
            BND            = BND,                             &  
            MPI            = MPI,                             &  
            IO             = IO                               )
    ELSE
        CALL CKProcedureForEveryoneGTSUni(                  &  ! Do the CK procedure 
          time           = time                           , &  ! for all elements in the domain. 
          dt             = OptionalFields%dt(1)           , &  ! Time-integrate the DOF and save 
          iteration      = timestep                       , &  ! the result in dgwork.
          MaterialVal    = OptionalFields%BackgroundValue , &  !  
          OptionalFields = OptionalFields,                  &  !  
          EQN            = EQN,                             &  !  
          MESH           = MESH,                            &  !  
          DISC           = DISC,                            &  !  
          IC             = IC,                              &  !  
          SOURCE         = SOURCE,                          &  !  
          BND            = BND,                             &  !  
          MPI            = MPI,                             &  !  
          IO             = IO                               )  !      
    ENDIF

#ifdef PARALLEL
    ! start epik/scorep region for communication of the degrees of freedom
    EPIK_USER_START( r_dofs_communication )
    SCOREP_USER_REGION_BEGIN(r_dofs_communication, "dofs_communication", SCOREP_USER_REGION_TYPE_COMMON)

    IF(MPI%nCPU.GT.1) THEN
        IF (DISC%Galerkin%DGMethod.EQ.3) THEN
            !
            !   For the parallel version of Seissol:
            !   exchange MPI boundary values before discretization
            !
            CALL MPIExchangeValues_LTS(DISC           = DISC,             &
                                       EQN            = EQN,              &
                                       BND            = BND,              &
                                       MESH           = MESH,             &
                                       IO             = IO,               &
                                       OptionalFields = OptionalFields,   &
                                       MPI            = MPI               )
        ELSE
            CALL MPIExchangeValues_GTS(DISC           = DISC,             &
                                       EQN            = EQN,              &
                                       BND            = BND               )
        ENDIF
       !
    ENDIF

    ! end epik/scorep region for communication of the degrees of freedom
    EPIK_USER_END( r_dofs_communication )
    SCOREP_USER_REGION_END(r_dofs_communication)
#endif

       IF(timestep.EQ.0) THEN
            logInfo(*) '-----------------------------------------------------------'
            logInfo(*) '    Minimum timestep is ', OptionalFields%dt(1)
            logInfo(*) '-----------------------------------------------------------'
       ENDIF

       ! do computations
       SELECT CASE(DISC%Galerkin%DGMethod)
       CASE(1) ! Global time step
            CALL ADERGalerkin3D_GTS(                              &  ! Quadfree ADER-DG
                time           = time                           , &  ! Quadfree ADER-DG
                dt             = OptionalFields%dt(1)           , &  ! Quadfree ADER-DG
                iteration      = timestep                       , &  ! Quadfree ADER-DG
                MaterialVal    = OptionalFields%BackgroundValue , &  ! Quadfree ADER-DG
                OptionalFields = OptionalFields,                  &  ! Quadfree ADER-DG
                EQN            = EQN,                             &  ! Quadfree ADER-DG
                MESH           = MESH,                            &  ! Quadfree ADER-DG
                DISC           = DISC,                            &  ! Quadfree ADER-DG
                IC             = IC,                              &  ! Quadfree ADER-DG
                SOURCE         = SOURCE,                          &  ! Quadfree ADER-DG
                BND            = BND,                             &  ! Quadfree ADER-DG
                MPI            = MPI,                             &  ! Quadfree ADER-DG
                IO             = IO  							  )

       CASE(3) ! Local time step
            CALL ADERGalerkin3D_LTS(                              &  ! Quadfree ADER-DG
                time           = time                           , &  ! Quadfree ADER-DG
                dt             = OptionalFields%dt(1)           , &  ! Quadfree ADER-DG
                iteration      = timestep                       , &  ! Quadfree ADER-DG
                MaterialVal    = OptionalFields%BackgroundValue , &  ! Quadfree ADER-DG
                OptionalFields = OptionalFields,                  &  ! Quadfree ADER-DG
                EQN            = EQN,                             &  ! Quadfree ADER-DG
                MESH           = MESH,                            &  ! Quadfree ADER-DG
                DISC           = DISC,                            &  ! Quadfree ADER-DG
                IC             = IC,                              &  ! Quadfree ADER-DG
                SOURCE         = SOURCE,                          &  ! Quadfree ADER-DG
                BND            = BND,                             &  ! Quadfree ADER-DG
                MPI            = MPI,                             &  ! Quadfree ADER-DG
                IO             = IO                               )  ! Quadfree ADER-DG

       END SELECT

       !
       timestep = timestep + 1                          ! advance time step number 
       DISC%iterationstep = DISC%iterationstep + 1      ! advance iteration number
       !                                                !  
       time      = time + OptionalFields%dt(1)           
       DISC%time = DISC%time + OptionalFields%dt(1) 
       MPITime   = time

       IF(DISC%Galerkin%DGMethod.EQ.3) THEN
         time    = MINVAL( DISC%LocalTime(:) )
         MPITime = time
#ifdef PARALLEL
         IF(MPI%nCPU.GT.1) THEN
            CALL MPI_ALLREDUCE(time,MPItime,1,MPI%MPI_AUTO_REAL,MPI_MIN,MPI%commWorld,iErr)
         ENDIF
#endif
       ENDIF

       !                                            
       CALL CalcPrintThisTimeStep(                    & ! Check if this timestep must be printed
            tol               = tol                 , &
            MaxTolerance      = DISC%MaxTolerance   , &
            timestep          = timestep            , &
            time              = MPITime             , &
            printTime         = DISC%printTime      , &
            printThisTimeStep = printThisTimeStep   , &
            EndTime           = DISC%EndTime        , &
            IO                = IO                    )
       !                                                ! 
       IF ( PrintThisTimeStep ) THEN
          !                                             ! 
          IF (IO%OutInterval%PlotNrGiven) THEN          !
             plotNr = IO%OutInterval%plotNumber         ! If called from Kop
          ELSE                                          !
             plotNr = timestep                          ! If called from seissol
          END IF                                        !
          !                                             !
          CALL data_output(                           & ! data_output
               dt         = OptionalFields%dt(:)    , & ! data_output
               time       = time                    , & ! data_output
               timestep   = plotNr                  , & ! data_output
               EQN        = EQN                     , & ! data_output
               MESH       = MESH                    , & ! data_output
               DISC       = DISC                    , & ! data_output
               SOURCE     = SOURCE                  , & ! data_output
               BND        = BND                     , & ! data_output
               MPI        = MPI                     , & ! data_output
               IO         = IO                      , & ! data_output
               ANALYSE    = ANALYSE                 , & ! data_output
           OptionalFields = OptionalFields            ) ! data_output
       ENDIF                                            !

       IF(modulo(timestep,50).EQ.0) THEN 
           call itime(now)
           logInfo0('(A,I12,A,I2,A,I2,A,I2)') ' done with time step: ', timestep, '    system-time: ', now(1), ':', now(2), ':', now(3)
       ENDIF

       IF((timestep.LT.DISC%MaxIteration).AND.(time.LT.DISC%EndTime)) THEN
            ContinueLoop = .TRUE.
       ELSE
            ContinueLoop = .FALSE.
       ENDIF
       !
#ifdef PARALLEL
       IF(DISC%Galerkin%DGMethod.EQ.3) THEN
           IF(MPI%nCPU.GT.1) THEN
              IF(MPItime.LT.DISC%EndTime) THEN
                 ContinueLoop = .TRUE.
              ELSE
                 ContinueLoop = .FALSE.
              ENDIF
           ENDIF
       ENDIF
#endif
       !  
    ENDDO ! end of time marching    !
#endif
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

    ! Print iteration information after the time loop
    CALL data_output(                           & ! data_output
         dt         = OptionalFields%dt(:)    , & ! data_output
         time       = time                    , & ! data_output
         timestep   = timestep                , & ! data_output
         EQN        = EQN                     , & ! data_output
         MESH       = MESH                    , & ! data_output
         DISC       = DISC                    , & ! data_output
         SOURCE     = SOURCE                  , & ! data_output
         BND        = BND                     , & ! data_output
         MPI        = MPI                     , & ! data_output
         IO         = IO                      , & ! data_output
         ANALYSE    = ANALYSE                 , & ! data_output
     OptionalFields = OptionalFields            ) ! data_output

    ! output magnitude for dynamic rupture simulations
    IF (EQN%DR.EQ.1 .AND. DISC%DynRup%magnitude_output_on.EQ.1) CALL magnitude_output(OptionalFields%BackgroundValue,DISC,MESH,MPI,IO)
    ! output GP-wise RF in extra files
    IF (EQN%DR.EQ.1 .AND. DISC%DynRup%RF_output_on.EQ.1) CALL output_rupturefront(DISC,MESH,MPI,IO)

    logInfo(*)'<--------------------------------------------------------->'  !
    logInfo(*)'<     calc_SeisSol successfully finished                  >'  !
    logInfo(*)'<--------------------------------------------------------->'  !
    !
    IO%AbortStatus = 0
    !
  END SUBROUTINE calc_SeisSol                            

END MODULE calc_SeisSol_mod
