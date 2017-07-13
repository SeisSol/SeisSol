!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) 2006-2014, SeisSol Group
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

MODULE ini_calcSeisSol_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE ini_calcSeisSol
     MODULE PROCEDURE ini_calc
  END INTERFACE

  INTERFACE close_calcSeisSol
     MODULE PROCEDURE close_calc
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: ini_calcSeisSol
  PUBLIC  :: close_calcSeisSol
  !----------------------------------------------------------------------------

CONTAINS
  SUBROUTINE ini_calc(pvar,EQN,MESH,DISC,OptionalFields,IO)
    !--------------------------------------------------------------------------
    USE plot_fields_mod,         ONLY : ini_plot_fields
    USE calc_deltaT_mod,         ONLY : ini_calc_deltaT
#ifdef HDF
!    USE hdf_output_utils_mod,         ONLY : ini_hdf_wavefield
#endif

    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)             :: EQN
    TYPE (tUnstructMesh)          :: MESH
    TYPE (tDiscretization)        :: DISC
    TYPE (tUnstructOptionalFields):: OptionalFields
    TYPE (tInputOutput)           :: IO
    REAL                          :: pvar(MESH%nElem,EQN%nVar)
    INTEGER                       :: allocStat
    !--------------------------------------------------------------------------
    INTENT(IN)                    :: pvar,EQN,MESH,DISC
    INTENT(INOUT)                 :: IO
    !--------------------------------------------------------------------------
    !
    ! field pvar_n is used by
    !   - ViscousPart
    !   - face_adjustment
    !
    ALLOCATE(OptionalFields%pvar_n(MESH%nnode,EQN%nVar), &
             STAT = allocStat               )

    IF (allocstat .NE. 0 ) THEN
       logError(*)' ALLOCATE ERROR in ini_calcSeisSol!'
       STOP
    END IF

    CALL ini_plot_fields(     pvar      ,OptionalFields, EQN, DISC, MESH     ,IO)! ok
#ifdef HDF
   ! CALL ini_hdf_wavefield(MESH, IO, DISC)
#endif

  END SUBROUTINE ini_calc

  SUBROUTINE close_calc(OptionalFields,IO,MPI)
    !--------------------------------------------------------------------------
    USE plot_fields_mod,       ONLY : close_plot_fields
    USE calc_deltaT_mod,       ONLY : close_calc_deltaT
#ifdef HDF
    USE hd_output_mod,         ONLY : close_hdf_wavefield
#endif
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tInputOutput)           :: IO
    TYPE (tUnstructOptionalFields):: OptionalFields
    TYPE (tMPI)                   :: MPI
    !--------------------------------------------------------------------------
    !
    ! field pvar_n is used by
    !   - ViscousPart
    !   - face_adjustment
    !
    DEALLOCATE(OptionalFields%pvar_n)
    !
    CALL close_plot_fields(OptionalFields,IO)
    CALL close_calc_DeltaT(OptionalFields)
#ifdef HDF
    if (IO%Format .eq. 5) then
        CALL close_hdf_wavefield(IO,MPI)
    endif
#endif

  END SUBROUTINE close_calc

END MODULE ini_calcSeisSol_mod
