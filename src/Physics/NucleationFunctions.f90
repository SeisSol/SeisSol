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
  INTERFACE regularizedYoffe
     MODULE PROCEDURE regularizedYoffe
  END INTERFACE

  !---------------------------------------------------------------------------!
  PUBLIC  :: Calc_SmoothStepIncrement
  PUBLIC  :: Calc_SmoothStep
  PUBLIC  :: regularizedYoffe
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

  FUNCTION C1(t, ts, tr) result(rC1)
    ! C1 to C6 are analytical functions
    ! used for building the regularized Yoffe function
    IMPLICIT NONE
    REAL        :: t, ts, tr, rC1
    rC1 = (0.5 * t + 0.25 * tr) * sqrt(t * (tr - t)) &
        + (t * tr - tr * tr) * asin(sqrt(t / tr)) &
        - 0.75 * tr * tr * atan(sqrt((tr - t) / t))
  END FUNCTION C1

  FUNCTION C2(t, ts, tr) result(rC2)
    IMPLICIT NONE
    REAL        :: t, ts, tr, rC2
    REAL, PARAMETER                :: pi=3.141592653589793
    rC2 = 0.375 * pi * tr * tr
  END FUNCTION C2


  FUNCTION C3(t, ts, tr) result(rC3)
    IMPLICIT NONE
    REAL        :: t, ts, tr, rC3
    rC3 = (ts - t - 0.5 * tr) * sqrt((t - ts) * (tr - t + ts)) &
        + tr * (2 * tr - 2 * t + 2 * ts) * asin(sqrt((t - ts) / tr)) &
        + 1.5 * tr * tr * atan(sqrt((tr - t + ts) / (t - ts)))
  END FUNCTION C3

  FUNCTION C4(t, ts, tr) result(rC4)
    IMPLICIT NONE
    REAL        :: t, ts, tr, rC4
    ! 2 typos fixed in the second term compared with Tinti et al. 2005   
    rC4 = (-ts + 0.5 * t + 0.25 * tr) * sqrt((t - 2.0 * ts) * (tr - t + 2.0 * ts)) &
        - tr * (tr - t + 2.0 * ts) * asin(sqrt((t - 2.0 * ts) / tr)) &
        - 0.75 * tr * tr * atan(sqrt((tr - t + 2.0 * ts) / (t - 2.0 * ts)))
  END FUNCTION C4

  FUNCTION C5(t, ts, tr) result(rC5)
    IMPLICIT NONE
    REAL        :: t, ts, tr, rC5
    REAL, PARAMETER                :: pi=3.141592653589793
    rC5 =  0.5 * pi * tr * (t - tr)
  END FUNCTION C5


  FUNCTION C6(t, ts, tr) result(rC6)
    IMPLICIT NONE
    REAL        :: t, ts, tr, rC6
    REAL, PARAMETER                :: pi=3.141592653589793
    rC6 =  0.5 * pi * tr * (2.0 * ts - t + tr)
  END FUNCTION C6




  FUNCTION regularizedYoffe(t, ts, tr) result(Gnuc)
    ! Implementation of the regularized Yoffe function
    ! defined in Appendix of Tinti et al. (2005)
    IMPLICIT NONE
    REAL        :: t, ts, tr
    REAL        :: K, Gnuc
    REAL, PARAMETER                :: pi=3.141592653589793
    K = 2.0 / (pi * tr * ts * ts)
    if (tr > 2.0 * ts) then
        if (t <= 0) then
            Gnuc = 0
            return
        else if (t <= ts) then
            Gnuc = K * (C1(t, ts, tr) + C2(t, ts, tr))
            return
        else if (t <= 2.0 * ts) then
            Gnuc = K * (C1(t, ts, tr) - C2(t, ts, tr) + C3(t, ts, tr))
            return
        else if (t < tr) then 
            Gnuc = K * (C1(t, ts, tr) + C3(t, ts, tr) + C4(t, ts, tr))
            return
        else if (t < tr + ts) then
            Gnuc = K * (C3(t, ts, tr) + C4(t, ts, tr) + C5(t, ts, tr))
            return
        else if (t < tr + 2.0 * ts) then
            Gnuc = K * (C4(t, ts, tr) + C6(t, ts, tr))
            return
        else
            Gnuc = 0
            return
        endif
    else
        if (t <= 0) then
            Gnuc = 0
            return
        else if (t <= ts) then
            Gnuc = K * (C1(t, ts, tr) + C2(t, ts, tr))
            return
        else if (t < tr) then
            Gnuc = K * (C1(t, ts, tr) - C2(t, ts, tr) + C3(t, ts, tr))
            return
        else if (t <= 2.0 * ts) then
            Gnuc =  K * (C5(t, ts, tr) + C3(t, ts, tr) - C2(t, ts, tr))
            return
        else if (t < tr + ts) then
            Gnuc =  K * (C3(t, ts, tr) + C4(t, ts, tr) + C5(t, ts, tr))
            return
        else if (t < tr + 2.0 * ts) then
            Gnuc = K * (C4(t, ts, tr) + C6(t, ts, tr))
            return
        else
            Gnuc = 0
            return
        endif
     endif
  END FUNCTION regularizedYoffe


END MODULE NucleationFunctions_mod
