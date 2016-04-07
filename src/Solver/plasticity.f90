!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Stephanie Wollherr
!!
!! @section LICENSE
!! Copyright (c) 2007-2014, SeisSol Group
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
!! Plasticity module: checks if an element undergoes plastic yielding or not
#include <Initializer/preProcessorMacros.fpp>

MODULE Plasticity_mod
  !---------------------------------------------------------------------------!
  USE TypesDef

#ifdef GENERATEDKERNELS
    use iso_c_binding, only: c_loc
#endif
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  !PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE Plasticity
     MODULE PROCEDURE Plasticity_3D
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Plasticity_3D
  !---------------------------------------------------------------------------!
  CONTAINS

!yldfac is only caluclated from the first DOF, and all DOF's are adjusted by the same coefficient
  SUBROUTINE Plasticity_3D(DISC, dgvar, DOFStress, nDegFr, nAlignedDegFr, BulkFriction, Tv, dt, mu, lambda, parameters, Energy, pstrain)
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization)    :: DISC
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom       !
    INTEGER     :: nDegFr
    integer     :: nAlignedDegFr

    REAL        :: Stress(1:nDegFr,6)                                         !local stress variable for the yield criterion
    REAL        :: Strain_ini(1:6)                                            !local initial strain variable for elastic strain energy
    REAL        :: Strain_total(1:6)                                          !local strain variable for elastic strain energy
    REAL        :: devStress(1:nDegFr,6)                                      !stress deviator for the yield criterion
    REAL        :: meanStress(1:nDegFr)                                       !mean stress
    REAL        :: angfric, yldfac                                            !Angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !relaxation time for the yield factor  
    REAL        :: mu, lambda                                                 !Lame parameters mu and lambda
    REAL        :: tau,taulim                                                 !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv                                                     !secInv=second Invariant of deviatoric stress
    REAL        :: BulkFriction, Tv
    REAL        :: DOFStress(1:nDegFr,1:6)
    REAL        :: dgvar(1:nAlignedDegFr,1:6)
    REAL        :: dudt_pstrain(1:6)
    REAL        :: pstrain(1:7)
    REAL        :: estrain(1:6), estrain_ini(1:6)                             !total elastic strain
    REAL        :: PlasticEnergy_tmp, EstrainEnergy_tmp
    REAL        :: Energy(1:2)
    REAL        :: parameters(1:3)
    REAL         :: I1,I1_0,I2,I2_0
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, DOFStress, nDegFr, BulkFriction, Tv, dt, mu, lambda, parameters
    INTENT(INOUT) :: dgvar, pstrain, Energy
    !-------------------------------------------------------------------------!

    dudt_pstrain = 0.0


    angfric = ATAN(BulkFriction) !angle of friction
    relaxtime = dt/(Tv) !Tv=dx/V_s with dx=min(dx);  Tv smaller-> stronger plasticity


! ---[ Calculate trial stress tensor ]---
! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz -> need to be specified throughout the whole medium for every element 
! as material values in EQN%IniStress, are mapped to the basis functions in dg_setup

    Stress(:,1:6)= dgvar(1:nDegFr,1:6)  + DOFStress(:,1:6)   !act.Stress + initial stress_xx

    !Calculate the total strain from the elastic stress-strain relation
    Strain_total(1:6) = MATMUL(DISC%Galerkin%Strain_matrix, Stress(1,1:6))

    !Calculate initial strain loading from initial stress loading (elementwise) -> move that outside the routine and calculate beforhand
    Strain_ini(1:6) = MATMUL(DISC%Galerkin%Strain_matrix, DOFStress(1,1:6))


! ---[ Calculate trial yield stress ]---

    ! Mean stress
    meanStress(1:nDegFr) = (Stress(1:nDegFr,1) + Stress(1:nDegFr,2)+ Stress(1:nDegFr,3) )*(1.0D0/3.0D0)

    ! Deviatoric stress
    devStress(1:nDegFr,1) = Stress(1:nDegFr,1) - meanStress(1:nDegFr)
    devStress(1:nDegFr,2) = Stress(1:nDegFr,2) - meanStress(1:nDegFr)
    devStress(1:nDegFr,3) = Stress(1:nDegFr,3) - meanStress(1:nDegFr)

    devStress(1:nDegFr,4:6) = Stress(1:nDegFr,4:6)

    ! Second invariant of stress deviator
    secInv = 0.5*(devStress(1,1)**2+devStress(1,2)**2+devStress(1,3)**2)+devStress(1,4)**2+devStress(1,5)**2+devStress(1,6)**2 

    ! Scalar measure of shear stress
    tau= SQRT(secInv)

    ! Yield stress   
    taulim = parameters(2)*COS(angfric) - meanStress(1)*SIN(angfric)! minus before sinus is for compressional stress=negative.
    taulim = MAX(0.0, taulim)


    ! Stress deviators are adjusted

    IF (tau .GT. taulim) THEN !plastic behaviour, else: elastic and stress tensor=trial stress tensor
       yldfac = 1.0D0- (1.0D0 - taulim/tau)*(1.0D0 - EXP(-relaxtime)) !factor by Duan/Day


       Stress(1:nDegFr,1) = devStress(1:nDegFr,1)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,2) = devStress(1:nDegFr,2)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,3) = devStress(1:nDegFr,3)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,4:6) = devStress(1:nDegFr,4:6)*yldfac


       !----Change the dofs-----
       dgvar(1:nDegFr,1:6) = Stress(1:nDegFr,1:6) - DOFStress(1:nDegFr,1:6)
       !dividing by 2*mu due to tensor convention
       dudt_pstrain(1:6) = ((1-yldfac)/(2*mu))*devStress(1, 1:6) !only the first dof is considered for plastic strain tensor

        
    ENDIF !yield criterion check

    !update plastic strain
    pstrain(1:6) = pstrain(1:6) + dudt_pstrain !plastic strain tensor
    !accumulated plastic strain
    pstrain(7) = pstrain(7)+ dt*sqrt(0.5*(dudt_pstrain(1)**2 + dudt_pstrain(2)**2 &
                                                   + dudt_pstrain(3)**2)+ dudt_pstrain(4)**2 + dudt_pstrain(5)**2 + dudt_pstrain(6)**2)

    !calculate energies
    PlasticEnergy_tmp = Stress(1,1)*dudt_pstrain(1) + Stress(1,2)*dudt_pstrain(2) + Stress(1,3)*dudt_pstrain(3) + 2.0*Stress(1,4)*dudt_pstrain(4) &
                      + 2.0*Stress(1,5)*dudt_pstrain(5) + 2.0*Stress(1,6)*dudt_pstrain(6)

    estrain(1:6) = Strain_total(1:6) - dudt_pstrain(1:6) !total elastic strain, if no plastic yielding -> elastic strain = total strain
    estrain_ini(1:6) = Strain_ini(1:6)

    !first and second invariants of the total elastic strain
    I1 = estrain(1) + estrain(2) + estrain(3)
    I2 = estrain(1)**2 + estrain(2)**2 + estrain(3)**2 + 2.0*estrain(4)**2 + 2.0*estrain(5)**2 + 2.0*estrain(6)**2
    !first and second invariants of the initial strain loading
    I1_0 = estrain_ini(1) + estrain_ini(2) + estrain_ini(3)
    I2_0 = estrain_ini(1)**2 + estrain_ini(2)**2 + estrain_ini(3)**2 + 2.0*estrain_ini(4)**2 + 2.0*estrain_ini(5)**2 + 2.0*estrain_ini(6)**2
    !Elastic strain energy
    !subtracted from the initial elastic strain
    EstrainEnergy_tmp = 0.5*lambda*(I1_0**2-I1**2) + mu*(I2_0-I2)

    Energy(1) = PlasticEnergy_tmp*parameters(1) !multiplied by volume to get integral over element
    Energy(2) = EstrainEnergy_tmp*parameters(1)


 END SUBROUTINE Plasticity_3D

END MODULE Plasticity_mod
