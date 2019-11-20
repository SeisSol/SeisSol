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
             Lambda, theta, sigma, Sh, SR, Dwn, DFinv, temp, pressure)
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)    :: EQN
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: i
    INTEGER     :: TP_grid_nz
    REAL        :: omega(TP_grid_nz)                                                  ! shear heating source
    REAL        :: Sh, SR, dt, tauV                                           ! shear stress, slip rate
    REAL        :: Dwn(TP_grid_nz), DFinv(TP_grid_nz)                                         !
    REAL        :: alpha_th, alpha_hy, rho_c                                  ! thermal and hydraulic diffusivities
    REAL        :: theta(TP_grid_nz), sigma(TP_grid_nz)                                       ! stored diffusion from previous timestep
    REAL        :: theta_current(TP_grid_nz), sigma_current(TP_grid_nz)                       ! diffusion next timestep
    REAL        :: tmp(TP_grid_nz)
    REAL        :: Lambda, Lambda_prime
    REAL        :: T,p, TP_half_width_shear_zone
    REAL        :: temp, pressure                                      ! temperatur, pressure in space domain
    REAL        :: temp_ini, pressure_ini                                      ! temperatur, pressure in space domain
    !-------------------------------------------------------------------------!
    INTENT(IN)  :: EQN, dt, TP_grid_nz, TP_half_width_shear_zone, alpha_th, alpha_hy, rho_c, Lambda, Sh, SR, Dwn, DFinv
    INTENT(INOUT):: theta, sigma
    INTENT(OUT) :: temp, pressure
    !-------------------------------------------------------------------------!


    tauV = Sh*SR !fault strenght*slip rate
    Lambda_prime = Lambda*alpha_th/(alpha_hy-alpha_th)
    tmp = (Dwn/TP_half_width_shear_zone)**2
    !1. Calculate diffusion of the field at previous timestep

    !temperature
    theta_current = theta*exp(-alpha_th*dt*tmp)
    !pore pressure + lambda'*temp
    sigma_current = sigma*exp(-alpha_hy*dt*tmp)

    !2. Add current contribution and get new temperature
    CALL heat_source(TP_half_width_shear_zone,alpha_th,dt,Dwn,TP_grid_nz,omega)
    theta = theta_current + (tauV/rho_c)*omega
    CALL heat_source(TP_half_width_shear_zone,alpha_hy,dt,Dwn,TP_grid_nz,omega)
    sigma = sigma_current + ((Lambda+Lambda_prime)*tauV)/(rho_c)*omega

    !3. Recover temperature and pressure using inverse Fourier
    ! transformation with the calculated fourier coefficients

    T = 0.0
    p = 0.0

    !new contribution
    DO i=1,TP_grid_nz
       T = T + (DFinv(i)/TP_half_width_shear_zone)*theta(i)
       p = p + (DFinv(i)/TP_half_width_shear_zone)*sigma(i)
    ENDDO

    !Update pore pressure change (sigma = pore pressure + lambda'*temp)
    !In the BIEM code (Lapusta) they use T without initial value
    p = p - Lambda_prime*T

    !Temp and pore pressure change at single GP on the fault + initial values
    temp = T + EQN%Temp_0
    pressure = -p + EQN%Pressure_0
 
  END SUBROUTINE Calc_ThermalPressure

  SUBROUTINE  heat_source(TP_half_width_shear_zone, alpha, dt, Dwn, TP_grid_nz, omega)

    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    REAL        :: TP_half_width_shear_zone
    INTEGER     :: TP_grid_nz                                                         ! number of points at the fault
    REAL        :: Dwn(TP_grid_nz)                                                    ! insert DISC%DynRup%TP_grid
    REAL        :: Sh, SR                                                     ! current shear stress and slip rate
    REAL        :: alpha, dt                                                      ! difussion parameter
    REAL        :: omega(TP_grid_nz)                                                  !
    REAL        :: tmp(TP_grid_nz)
    REAL, PARAMETER :: pi=3.141592653589793     ! CONSTANT pi
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: TP_half_width_shear_zone, alpha, dt, Dwn, TP_grid_nz
    INTENT(OUT)   :: omega
    !-------------------------------------------------------------------------!
    !Gaussian shear zone in spectral domain, normalized by w
    tmp = (Dwn/TP_half_width_shear_zone)**2
    !original function in spatial domain
    !omega = 1/(w*sqrt(2*pi))*exp(-0.5*(z/TP_half_width_shear_zone).^2);
    !function in the wavenumber domain *including additional factors in front of the heat source function*
    !omega = 1/(*alpha*Dwn**2**(sqrt(2.0*pi))*exp(-0.5*(Dwn*TP_half_width_shear_zone)**2)*(1-exp(-alpha**dt**tmp)) 
    !inserting Dwn/TP_half_width_shear_zone (scaled) for Dwn cancels out TP_half_width_shear_zone
    omega = 1.0/(alpha*tmp*(sqrt(2.0*pi)))*exp(-0.5*(Dwn)**2)*(1.0 - exp(-alpha*dt*tmp))

    RETURN

  END SUBROUTINE heat_source

END MODULE Thermalpressure_mod
