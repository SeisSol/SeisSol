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

MODULE COMMON_InitialField_mod
  !----------------------------------------------------------------------------
  USE TypesDef

  USE iso_c_binding
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE InitialField
     MODULE PROCEDURE InitialField
  SUBROUTINE InitialFieldPlanarWave(time, x, y, z, variables) bind(C, name="initial_field_planarwave")
    use, intrinsic :: iso_c_binding

    real(kind=c_double), value :: time
    real(kind=c_double), value :: x
    real(kind=c_double), value :: y
    real(kind=c_double), value :: z
    real(kind=c_double), dimension(*), intent(inout) :: variables
  END SUBROUTINE

  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC   :: InitialField
  !----------------------------------------------------------------------------
  !
CONTAINS
  !============================================================================
  ! Returns the initial condition in form of a 3D state vector
  ! in the vector Variable for the position x, y, z and time level time
  !============================================================================

  SUBROUTINE InitialField(Variable, Variable_ANE, time, x, y, z, iElem, EQN,IC, SOURCE, IO)
    !--------------------------------------------------------------------------
    USE common_operators_mod
    USE DGBasis_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tEquations)              :: EQN
    TYPE (tInitialCondition)       :: IC
    TYPE (tSource)                 :: SOURCE
    TYPE (tInputOutput)            :: IO
    REAL                           :: Variable(:)
    REAL                           :: Variable_ANE(:)
    REAL                           :: time
    REAL                           :: x, y, z
    INTEGER                        :: i,j
    INTEGER :: iElem
    !
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !INTENT(IN)  :: x, y, z, iElem, EQN, IC, SOURCE, IO
    !INTENT(OUT) :: Variable
    !--------------------------------------------------------------------------
    !
    SELECT CASE(TRIM(IC%cICType))
    CASE('Planarwave')
       Variable_ANE(:) = 0.

       CALL InitialFieldPlanarWave(time,x,y,z,variable)
    CASE('Zero')
      Variable(:) = 0.
      Variable_ANE(:) = 0.
    CASE DEFAULT
       logError(*) 'InitialField: none of the possible initial conditions was chosen'
       logError(*) TRIM(IC%cICType),'|'
       STOP
    END SELECT

  END SUBROUTINE InitialField


END MODULE COMMON_InitialField_mod
