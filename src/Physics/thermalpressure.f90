!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Stephanie Wollherr (wollherr AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/wollherr)
!!
!! @section LICENSE
!! Copyright (c) 2007-2019, SeisSol Group
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
!! Pseudo-spectral code for Thermal Pressurization, called for rate-and-state friction

#ifdef BG 
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE Thermalpressure_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  !PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE Calc_ThermalPressure
     MODULE PROCEDURE Calc_ThermalPressure
  END INTERFACE
  INTERFACE heat_source
     MODULE PROCEDURE heat_source
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Calc_ThermalPressure, heat_source

  !---------------------------------------------------------------------------!

CONTAINS

  SUBROUTINE Calc_ThermalPressure(EQN,dt, TP_grid_nz, TP_half_width_shear_zone, alpha_th, alpha_hy, rho_c, &
             Lambda, theta, sigma, Sh, SR, TP_grid, DFinv, temp, pressure)
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)    :: EQN
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: TP_grid_nz
    REAL        :: omega(TP_grid_nz)                                                  ! shear heating source
    REAL        :: Sh, SR, dt, tauV                                           ! shear stress, slip rate
    REAL        :: TP_grid(TP_grid_nz), DFinv(TP_grid_nz)                                         !
    REAL        :: normalized_TP_grid(TP_grid_nz)
    REAL        :: alpha_th, alpha_hy, rho_c                                  ! thermal and hydraulic diffusivities
    REAL        :: theta(TP_grid_nz), sigma(TP_grid_nz)                       ! stored diffusion from previous timestep
    REAL        :: exp_alpha_th(TP_grid_nz), exp_alpha_hy(TP_grid_nz)         ! intermediate variables for optimization
    REAL        :: Lambda, Lambda_prime
    REAL        :: TP_half_width_shear_zone
    REAL        :: temp, pressure                                      ! temperatur, pressure in space domain
    !-------------------------------------------------------------------------!
    INTENT(IN)  :: EQN, dt, TP_grid_nz, TP_half_width_shear_zone, alpha_th, alpha_hy, rho_c, Lambda, Sh, SR, TP_grid, DFinv
    INTENT(INOUT):: theta, sigma
    INTENT(OUT) :: temp, pressure
    !-------------------------------------------------------------------------!

    tauV = Sh*SR !fault strenght*slip rate
    Lambda_prime = Lambda*alpha_th/(alpha_hy-alpha_th)
    !Gaussian shear zone in spectral domain, normalized by w
    normalized_TP_grid(:) = (TP_grid(:)/TP_half_width_shear_zone)**2
    !1. Calculate diffusion of the field at previous timestep

    !temperature
    exp_alpha_th(:) = exp(-alpha_th*dt*normalized_TP_grid(:))
    theta(:) = theta(:) * exp_alpha_th(:)
    !pore pressure + lambda'*temp
    exp_alpha_hy(:) = exp(-alpha_hy*dt*normalized_TP_grid(:))
    sigma(:) = sigma(:) * exp_alpha_hy(:)

    !2. Add current contribution and get new temperature
    omega(:) = heat_source(alpha_th, dt, TP_grid, normalized_TP_grid, exp_alpha_th, TP_grid_nz)
    theta(:) = theta(:) + (tauV / rho_c) * omega(:)
    omega(:) = heat_source(alpha_hy, dt, TP_grid, normalized_TP_grid, exp_alpha_hy, TP_grid_nz)
    sigma(:) = sigma(:) + (((Lambda + Lambda_prime) * tauV) / rho_c) * omega(:)

    !3. Recover temperature and pressure using inverse Fourier
    ! transformation with the calculated fourier coefficients

    temp = SUM(DFinv(:) / TP_half_width_shear_zone * theta(:))
    pressure = SUM(DFinv(:) / TP_half_width_shear_zone * sigma(:))

    !Update pore pressure change (sigma = pore pressure + lambda'*temp)
    !In the BIEM code (Lapusta) they use T without initial value
    pressure = pressure - Lambda_prime * temp

    !Temp and pore pressure change at single GP on the fault + initial values
    temp = temp + EQN%Temp_0
    pressure = -pressure + EQN%Pressure_0
 
  END SUBROUTINE Calc_ThermalPressure

  FUNCTION heat_source(alpha, dt, TP_grid, normalized_TP_grid, exp_alpha, TP_grid_nz) result(omega)

    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: TP_grid_nz                                                 ! number of points at the fault
    REAL        :: alpha, dt                                                  ! difusion parameter
    REAL        :: omega(TP_grid_nz), exp_alpha(:)                            ! exp_alpha is exp(-alpha*dt*normalized_TP_grid(:))
    REAL        :: TP_grid(TP_grid_nz)
    REAL        :: normalized_TP_grid(TP_grid_nz)
    REAL, PARAMETER :: pi=3.141592653589793     ! CONSTANT pi
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: alpha, dt, normalized_TP_grid, TP_grid_nz
    !-------------------------------------------------------------------------!
    !original function in spatial domain
    !omega = 1/(w*sqrt(2*pi))*exp(-0.5*(z/TP_half_width_shear_zone).^2);
    !function in the wavenumber domain *including additional factors in front of the heat source function*
    !omega = 1/(*alpha*TP_grid**2**(sqrt(2.0*pi))*exp(-0.5*(TP_grid*TP_half_width_shear_zone)**2)*(1-exp(-alpha**dt**normalized_TP_grid)) 
    !inserting TP_grid/TP_half_width_shear_zone (scaled) for TP_grid cancels out TP_half_width_shear_zone
    omega = 1.0/(alpha*normalized_TP_grid(:) * sqrt(2.0*pi))*exp(-0.5*(TP_grid(:))**2) * (1.0 - exp_alpha(:))

    RETURN

  END FUNCTION heat_source

END MODULE Thermalpressure_mod
