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

MODULE inioutput_SeisSol_mod
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE inioutput_SeisSol
     MODULE PROCEDURE inioutput_SeisSol
  END INTERFACE

#ifdef PARALLEL
  interface
    function mkdir(path,mode) bind(c,name="mkdir")
      use iso_c_binding
      integer(c_int) :: mkdir
      character(kind=c_char,len=1) :: path(*)
      integer(c_int16_t), value :: mode
    end function mkdir
  end interface
#endif

  !----------------------------------------------------------------------------
  PUBLIC  :: inioutput_SeisSol
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE inioutput_SeisSol(time,timestep,pvar,cvar,EQN,IC,MESH,MPI,      &
       SOURCE,DISC,BND,OptionalFields,IO,Analyse, &
       programTitle) !
    !--------------------------------------------------------------------------
    USE TypesDef
#ifdef HDF
    USE receiver_hdf_mod
#else
    USE receiver_mod
#endif
    USE ini_calcSeisSol_mod
    USE dg_setup_mod
    USE data_output_mod

#ifdef GENERATEDKERNELS
    use iso_c_binding
    use f_ftoc_bind_interoperability
#else
    use WaveFieldWriter
#endif

#ifdef PARALLEL
    use iso_c_binding
#endif
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                                              !
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    REAL                           :: time,x,y,Variable(8),k1,k2               !
    INTEGER                        :: timestep                                 !
    INTEGER                        :: i                                        !
    INTEGER                        :: outputMaskInt(9)                         !
    REAL,POINTER                   :: pvar(:,:)                                ! @TODO, breuera: remove not used
    REAL,POINTER                   :: cvar(:,:)                                !
    TYPE (tEquations)              :: EQN                                      !
    TYPE (tInitialCondition)       :: IC                                       !
    TYPE (tUnstructMesh)           :: MESH                                     !
    TYPE (tMPI), OPTIONAL          :: MPI                                      !
    TYPE (tSource)                 :: SOURCE                                   !
    TYPE (tDiscretization)         :: DISC                                     !
    TYPE (tUnstructOptionalFields) :: OptionalFields                           !
    TYPE (tInputOutput)            :: IO                                       !
    TYPE (tBoundary)               :: BND                                      !
    TYPE (tAnalyse)                :: Analyse                                  !
    CHARACTER(LEN=100)             :: programTitle                             !
    ! local variable declaration                                               !
    CHARACTER(LEN=256)             :: outfile
    CHARACTER(LEN=5)               :: cmyrank
    integer                     :: timestepWavefield
    integer                     :: mkdirRet
    !--------------------------------------------------------------------------
    INTENT(IN)                     :: programTitle                             !
    INTENT(INOUT)                  :: EQN,DISC,IO, OptionalFields, MESH        ! Some values are set in the TypesDef
    INTENT(INOUT)                  :: IC,SOURCE,BND       !
    INTENT(INOUT)                  :: time,timestep, Analyse             !
    !--------------------------------------------------------------------------
    !                                                                          !
    ! register epik/scorep function
    EPIK_FUNC_REG("inioutput")
    SCOREP_USER_FUNC_DEFINE()
    !--------------------------------------------------------------------------
    !                                                                          !
    ! start epik/scorep function
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN("inioutput")

    timestepWavefield = 0

#ifdef HDF
    CALL ini_receiver_hdf(                                &                    ! Initialize receivers
         EQN    = EQN                                   , &                    ! Initialize receivers
         MESH   = MESH                                  , &                    ! Initialize receivers
         DISC   = DISC                                  , &                    ! Initialize receivers
         SOURCE = SOURCE                                , &                    ! Initialize receivers
         IO     = IO                                    , &                    ! Initialize receivers
         MPI    = MPI                                     )                    ! Initialize receivers
    !                                                                          !
#else
    CALL ini_receiver(                                    &                    ! Initialize receivers
         EQN    = EQN                                   , &                    ! Initialize receivers
         MESH   = MESH                                  , &                    ! Initialize receivers
         DISC   = DISC                                  , &                    ! Initialize receivers
         SOURCE = SOURCE                                , &                    ! Initialize receivers
         IO     = IO                                    , &                    ! Initialize receivers
         MPI    = MPI                                     )                    ! Initialize receivers
    !                                                                          !
#endif

#ifdef PARALLEL
    WRITE(cmyrank,'(I5.5)') MPI%myrank     ! myrank -> cmyrank
    outfile = TRIM(IO%OutputFile )//'.'//TRIM(cmyrank)  ! mshfile Name
#else
    outfile = TRIM(IO%OutputFile )                      ! mshfile Name
#endif
    !                                                                          !
    logInfo(*) '<--------------------------------------------------------->'  !
    logInfo(*) '<      Initializing calcSeisSol (OptionalFields)          >'  !
    logInfo(*) '<--------------------------------------------------------->'  !
    !                                                                          !
    CALL ini_calcSeisSol(                                 &                    ! ini_calc
         pvar           = pvar                          , &                    ! ini_calc
         EQN            = EQN                           , &                    ! ini_calc
         MESH           = MESH                          , &                    ! ini_calc
         DISC           = DISC                          , &                    ! ini_calc
         OptionalFields = OptionalFields                , &                    ! ini_calc
         IO             = IO                              )                    ! ini_calc
         !                                                                     !

#ifdef PARALLEL
    if (IO%Format .ne. 5 .and. IO%Format .ne. 6 .and. IO%Format .ne. 10) then
        logInfo0(*) 'Creating output folders'
        WRITE(cmyrank,'(I5.5)') MPI%myrank
        mkdirRet = mkdir(TRIM(IO%OutputFile )//'-'//TRIM(cmyrank)//c_null_char, int(o'777',c_int16_t))
    endif
#endif


#ifdef GENERATEDKERNELS
    do i = 1, 9
        if ( io%OutputMask(3+i) ) then
            outputMaskInt(i) = 1
        else
            outputMaskInt(i) = 0
        end if
    end do
    call c_interoperability_initializeIO(    &
        i_mu        = disc%DynRup%mu,        &
        i_slipRate1 = disc%DynRup%slipRate1, &
        i_slipRate2 = disc%DynRup%slipRate2, &
        i_slip     = disc%DynRup%slip,      &
        i_slip1     = disc%DynRup%slip1,    &
        i_slip2     = disc%DynRup%slip2,    &
        i_state     = disc%DynRup%stateVar,  &
        i_strength  = disc%DynRup%strength,  &
        i_numSides  = mesh%fault%nSide,      &
        i_numBndGP  = disc%galerkin%nBndGP,  &
        i_refinement= io%Refinement,         &
        i_outputMask= outputMaskInt,         &
        i_outputRegionBounds = io%OutputRegionBounds(1:6))
#else
    if (IO%Format .eq. 6) then
        call waveFieldWriterInit(0, disc, eqn, io, mesh, mpi)
    endif

    if (timestep .eq. 0) then
        ! Only write initial data when we don't start from a check point
        logInfo(*) 'Write initial data: '!
        CALL data_output(                                     &                    ! data_output
            dt         = OptionalFields%dt(:)              , &                    ! data_output
            time       = time                              , &                    ! data_output
            timestep   = timestep                          , &                    ! data_output
            EQN        = EQN                               , &                    ! data_output
            MESH       = MESH                              , &                    ! data_output
            DISC       = DISC                              , &                    ! data_output
            SOURCE     = SOURCE                            , &                    ! data_output
            BND        = BND                               , &                    ! data_output
            MPI        = MPI                               , &                    ! data_output
            IO         = IO                                , &                    ! data_output
            ANALYSE    = ANALYSE                           , &                    ! data_output
            OptionalFields = OptionalFields                  )                    ! data_output
    endif
#endif

    ! end epik/scorep function
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()
  END SUBROUTINE inioutput_SeisSol                                                    !

END MODULE inioutput_SeisSol_mod
