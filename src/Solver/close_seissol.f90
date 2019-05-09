!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) 2006, SeisSol Group
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

MODULE close_SeisSol_mod
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE close_SeisSol
     MODULE PROCEDURE close_SeisSol
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: close_SeisSol
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE close_SeisSol(pvar,cvar,OptionalFields,EQN,IC,BND,MESH,DISC,SOURCE,IO,MPI)
    !--------------------------------------------------------------------------
    USE TypesDef
    USE ini_OptionalFields_mod,  ONLY: close_OptionalFields
    USE dg_setup_mod
    USE allocate_mesh_mod 
    use f_ftoc_bind_interoperability
    USE calc_deltaT_mod,       ONLY : close_calc_deltaT
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    TYPE (tEquations)              :: EQN
    TYPE(tInitialCondition)        :: IC                                                          !
    TYPE(tBoundary)                :: BND                                                         !
    TYPE (tUnstructMesh)           :: MESH
    TYPE (tDiscretization)         :: DISC
    TYPE (tSource)                 :: SOURCE
    TYPE (tInputOutput)            :: IO
    TYPE (tUnstructOptionalFields) :: OptionalFields
    TYPE (tMPI)                    :: MPI
    REAL, POINTER                  :: pvar( :,: )
    REAL, POINTER                  :: cvar( :,: )
    ! local variables
    INTEGER                        :: deallocStat
    !--------------------------------------------------------------------------
    INTENT(IN)                     :: EQN, MPI
    INTENT(INOUT)                  :: MESH,DISC, IO, SOURCE
    !--------------------------------------------------------------------------
    CALL close_calc_DeltaT(OptionalFields)
    !
    ! Call only the necessary deallocate subroutines for 3D unstruct. ADER-DG !!
    ! Only very few of the Finite Volume Mesh information is needed, 
    ! so for 3D ADER-DG only very few mesh information is allocated.
    !
    CALL destruct_mesh_level0_1(MESH)
    CALL destruct_mesh_level0_2(MESH)

    call c_interoperability_finalizeIO()
    !    
    CALL close_OptionalFields(                   &
         OptionalFields = OptionalFields       , &
         SOURCE         = SOURCE               , &
         EQN            = EQN                  , &
         MESH           = MESH                 , &
         DISC           = DISC                 , &
         IO             = IO                     )


    CALL closeGalerkin3D_us(EQN,MESH,DISC,IO)  

    DEALLOCATE (pvar,cvar,STAT=deallocStat)
    IF (deallocStat .NE. 0) THEN                                    ! Falls ein Fehler
       logError(*) 'cannot deallocate all arrays!'                  ! aufgetreten ist
       STOP                                                         ! Programmabbruch
    END IF                                                          !

    logInfo(*) '<--------------------------------------------------------->'  !
    logInfo(*) '<     close_SeisSol successfully finished                 >'  !
    logInfo(*) '<--------------------------------------------------------->'  !
    logInfo(*) ' '
    logInfo(*) ' '
    logInfo(*) '____________  Program finished without errors  ____________'
    logInfo(*) ' '
    logInfo(*) ' '

  END SUBROUTINE close_SeisSol

END MODULE close_SeisSol_mod

