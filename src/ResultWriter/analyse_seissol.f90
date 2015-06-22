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

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE analyse_SeisSol_mod
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE analyse_SeisSol
     MODULE PROCEDURE analyse_SeisSol
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: analyse_SeisSol
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE analyse_SeisSol(time,timestep,pvar,EQN,IC,MESH,DISC,BND,     &
       SOURCE,IO,Analyse,OptionalFields,Debug,MPI)
    !--------------------------------------------------------------------------
    USE TypesDef
    USE plot_fields_mod
    USE dg_setup_mod
    USE ini_OptionalFields_mod
    USE COMMON_operators_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tEquations)        :: EQN
    TYPE (tInitialCondition) :: IC
    TYPE (tUnstructMesh)     :: MESH
    TYPE (tDiscretization)   :: DISC
    TYPE (tBoundary)         :: BND
    TYPE (tSource)           :: SOURCE
    TYPE (tInputOutput)      :: IO
    TYPE (tAnalyse)          :: Analyse
    TYPE (tUnstructOptionalFields) :: OptionalFields
    TYPE (tDebug)            :: Debug
    TYPE (tMPI)              :: MPI
    REAL                     :: time
    INTEGER                  :: timestep
    REAL                     :: pvar( MESH%nElem,EQN%nVar )
    ! local variables
    REAL                     :: cvar(      MESH%nElem,EQN%nVar )
    REAL                     :: pvar_init(MESH%nElem,EQN%nVar ) 
    TYPE(tUnstructOptionalFields) :: Dummy
    CHARACTER (LEN=100)      :: mshfile               ! auxiliary output mesh file
    CHARACTER (LEN=100)      :: exactfile
    CHARACTER (LEN=210)      :: pbc_file
    CHARACTER(LEN=10)        :: ExTitleMask(9)        ! Variable names for output of exact solution
    CHARACTER (LEN=4)        :: cmyrank
    INTEGER                  :: ielem, i
    REAL                     :: vector1(3)
    REAL                     :: rhol,rhor,rho,ul,ur,u,pl,pr,p,al,ar,s,sqDiff(1:3)
    REAL                     :: E_FieldVec(3),start_pos(3), tempValue1, tempValue2,tempValue3
    REAL                     :: exactPos(3),numericalPos(3),diffPos(3),NbrOfParticle
    REAL                     :: meanErrorSqx,meanErrorSqy,meanErrorSqz,meanErrorSq
    CHARACTER(LEN=3)         :: varName
    CHARACTER(LEN=610)       :: Name
    CHARACTER (LEN=350)      :: Filename
    !--------------------------------------------------------------------------
    INTENT(IN)               :: time,timestep,Debug
    INTENT(INOUT)            :: IO,pvar,IC,EQN,DISC,MESH,OptionalFields
    !--------------------------------------------------------------------------
    !                                                                          !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) 'CPU-Time:    ', DISC%LoopCPUTime
    logInfo(*) 'CPU-TimeInt: ', DISC%Galerkin%Cpu_TimeInt, ' n: ', DISC%Galerkin%nCpu_TimeInt
    logInfo(*) 'CPU-Volume:  ', DISC%Galerkin%Cpu_Volume, ' n: ', DISC%Galerkin%nCpu_Volume
    logInfo(*) 'CPU-Source:  ', DISC%Galerkin%Cpu_Source, ' n: ', DISC%Galerkin%nCpu_Source
    logInfo(*) 'CPU-Flux:    ', DISC%Galerkin%Cpu_Flux, ' n: ', DISC%Galerkin%nCpu_Flux
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) 'h1 = ', MESH%MaxSQRTVolume                      !
    logInfo(*) 'h2 = ', MESH%MaxCircle                          !
    CALL ini_OptionalFields(OptionalFields = dummy                    , & !
                            SOURCE         = SOURCE                   , & !
                            EQN            = EQN                      , & !
                            MESH           = MESH                     , & !
                            DISC           = DISC                     , & !
                            IO             = IO                         ) !
    !                                                                          !
    logInfo(*) '<--------------------------------------------------------->'
    logInfo(*) '    Average time step of iteration is ',time/timestep
    logInfo(*) '<--------------------------------------------------------->'
    !                                                                          !
    !                                                                          !
    IF(DISC%DiscretizationMethod.EQ.2) THEN                                    ! Analyse DG
            CALL AnalyseGalerkin3D_us_new(                                       & ! Analyse DG
                time           = time                                          , & ! Analyse DG
                ANALYSE        = Analyse                                       , & ! Analyse DG
                EQN            = EQN                                           , & ! Analyse DG
                MESH           = MESH                                          , & ! Analyse DG
                DISC           = DISC                                          , & ! Analyse DG
                BND            = BND                                           , & ! Analyse DG
                SOURCE         = SOURCE                                        , & ! Analyse DG
                IC             = IC                                            , & ! Analyse DG
                IO             = IO                                            , & ! Analyse DG
                OptionalFields = OptionalFields                                , & ! Analyse DG
                MPI            = MPI                                             ) ! Analyse DG
    ENDIF                    
    !                                                                          !
    IF ((timestep.eq.DISC%MaxIteration).OR.(time .eq.DISC%EndTime)) THEN   !
      CALL close_plot_fields(OptionalFields,IO)                                  !
    END IF
    !

    logInfo(*) '<--------------------------------------------------------->'  !
    logInfo(*) '<     analyse_SeisSol successfully finished               >'  !
    logInfo(*) '<--------------------------------------------------------->'  !
                                                                     !
  END SUBROUTINE analyse_SeisSol                                                !

END MODULE analyse_SeisSol_mod
