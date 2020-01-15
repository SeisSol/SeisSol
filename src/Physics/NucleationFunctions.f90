!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Thomas Ulrich
!!
!! @section LICENSE
!! Copyright (c) 2007-2020, SeisSol Group
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
!! Smooth step function used to smoothly nucleate a rupture over time Tnuc 
!! similarly with SCEC TPV103 benchmark  

#ifdef BG 
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE NucleationFunctions_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  !PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE Calc_SmoothStepIncrement
     MODULE PROCEDURE Calc_SmoothStepIncrement
  END INTERFACE
  INTERFACE Calc_SmoothStep
     MODULE PROCEDURE Calc_SmoothStep
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Calc_SmoothStepIncrement
  PUBLIC  :: Calc_SmoothStep
  !---------------------------------------------------------------------------!

CONTAINS

  FUNCTION Calc_SmoothStepIncrement(time, Tnuc, dt)result(Gnuc)
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    REAL        :: time, Tnuc, dt, Gnuc, prevtime, Gnuc_prev
    !-------------------------------------------------------------------------!
    IF ((time.GT.0.0D0).AND.(time.LE.Tnuc)) THEN
        Gnuc = Calc_SmoothStep(time, Tnuc)
        prevtime = time - dt
        IF (prevtime.GT.0.0D0) THEN
            Gnuc= Gnuc - Calc_SmoothStep(prevtime, Tnuc)
        ENDIF
    ELSE
        Gnuc=0.0D0
    ENDIF
  END FUNCTION Calc_SmoothStepIncrement

  FUNCTION Calc_SmoothStep(time, Tnuc)result(Gnuc)
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    REAL        :: time, Tnuc, Gnuc
    !-------------------------------------------------------------------------!
    if (time.LE.0) then
        Gnuc=0.0
    else 
        if (time.LT.Tnuc) then
            Gnuc = EXP((time-Tnuc)**2/(time*(time-2.0D0*Tnuc)))
        else
            Gnuc=1d0
        endif
    endif
  END FUNCTION Calc_SmoothStep


END MODULE NucleationFunctions_mod
