!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) SeisSol Group
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

MODULE COMMON_printThisTimeStep_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE CalcPrintThisTimeStep
     MODULE PROCEDURE CalcPrintThisTimeStep
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  ::  CalcPrintThisTimeStep
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE CalcPrintThisTimeStep(tol,MaxTolerance,timestep,time,printTime,printThisTimeStep,EndTime,IO)    
    !--------------------------------------------------------------------------
    IMPLICIT NONE     
    !--------------------------------------------------------------------------
    TYPE(tInputOutput) :: IO
    INTEGER            :: timestep
    REAL, OPTIONAL     :: tol
    REAL, OPTIONAL     :: MaxTolerance  
    REAL               :: time
    REAL               :: printTime
    REAL               :: EndTime
    LOGICAL            :: printThisTimeStep 
    LOGICAL            :: skip
    !--------------------------------------------------------------------------
    INTENT(IN)         :: tol,MaxTolerance,timestep,time,EndTime,IO
    INTENT(OUT)        :: printThisTimeStep 
    INTENT(INOUT)      :: printTime
    !--------------------------------------------------------------------------

    skip = .FALSE.

    IF (.NOT.skip) THEN

       SELECT CASE(IO%OutInterval%printIntervalCriterion)                     ! Dataoutput 
       CASE(1)                                                                ! Dataoutput 
          IF (MOD(timestep,IO%OutInterval%Interval).EQ.0) THEN                ! Dataoutput 
             printThisTimeStep = .TRUE.                                       ! Dataoutput 
          ELSE                                                                ! Dataoutput 
             printThisTimeStep = .FALSE.                                      ! Dataoutput 
          END IF                                                              ! Dataoutput 
       CASE(2)                                                                ! Dataoutput 
          IF (time.EQ.printtime) THEN                                         ! Dataoutput 
             printThisTimeStep = .TRUE.                                       ! Dataoutput 
             printtime         = printtime + IO%outInterval%TimeInterval      ! Dataoutput 
          ELSE                                                                ! Dataoutput 
             printThisTimeStep = .FALSE.                                      ! Dataoutput 
          END IF                                                              ! Dataoutput 
       CASE(3)                                                                ! Dataoutput 
          IF (     (MOD(timestep,IO%OutInterval%Interval).EQ.0        ) &     ! Dataoutput 
               .OR.(time                                 .EQ.printtime) )THEN ! Dataoutput 
             printThisTimeStep = .TRUE.                                       ! Dataoutput 
             printtime         = printtime + IO%outInterval%TimeInterval      ! Dataoutput 
          ELSE                                                                ! Dataoutput 
             PrintThisTimeStep = .FALSE.                                      ! Dataoutput 
          END IF                                                              ! Dataoutput 
          !                                                                   ! Dataoutput 
       END SELECT                                                             ! Dataoutput      
       !                                                                      ! Dataoutput 
    END IF
    !
    IF(       (.NOT.IO%OutInterval%PlotNrGiven) &                          ! Dataoutput: Falls die Routine nicht vom KOP2d aufgerufen wird
         .AND.(.NOT.printThisTimeStep         ) ) THEN                     ! Dataoutput: und falls dieser Zeitschritt noch nicht geplottet wird  
       IF (time.GE.EndTime ) THEN                                          ! Dataoutput: Falls der Endzeitpunkt erreicht wird, soll dieser auf 
          printThisTimeStep = .TRUE.                                       ! Dataoutput: jeden Fall ausgegeben werden 
       END IF                                                              ! Dataoutput 
    END IF                                                                 ! Dataoutput 

  END SUBROUTINE CalcPrintThisTimeStep

END MODULE COMMON_printThisTimeStep_mod
