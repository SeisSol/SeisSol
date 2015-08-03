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
  SUBROUTINE Plasticity_3D(dgvar, DOFStress, nDegFr, BulkFriction, Tv, PlastCo, dt, iElem, mu, dudt_plastic, dudt_pstrain)  
    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    !TYPE(tEquations)         :: EQN
    TYPE(tDiscretization)    :: DISC
    !TYPE(tUnstructMesh)      :: MESH
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iElem, LocElemType                                         ! Element number                   !
    INTEGER     :: iSide                                                      ! Local element side number        !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom       !
    INTEGER     :: Bndcase                                                    ! case=1: element on the fault, case=0 element not on the fault
    INTEGER     :: nDegFr

    REAL        :: Stress(1:nDegFr,6)                                         !local stress variable for the yield criterion
    REAL        :: devStress(1:nDegFr,6)                                      !stress deviator for the yield criterion
    REAL        :: meanStress(1:nDegFr)                                       !mean stress
    REAL        :: angfric, yldfac                                            !Angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !relaxation time for the yield factor  
    REAL        :: mu                                                         ! Lame parameter mu
    REAL        :: LocDegFr, LocnVar                                          !
    REAL        :: tau,taulim                                                 !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv                                                     !secInv=second Invariant of deviatoric stress
    REAL        :: BulkFriction, Tv, PlastCo
    REAL        :: DOFStress(1:nDegFr,1:6)
    REAL        :: dgvar(1:nDegFr,1:6)
    REAL        :: dudt_plastic(1:nDegFr,1:6)
    !REAL        :: dudt_pstrain(1:nDegFr,1:9)
    REAL        :: dudt_pstrain(1:6)
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: dgvar, DOFStress, nDegFr, BulkFriction, Tv, PlastCo, dt, iElem, mu
    INTENT(OUT)   :: dudt_plastic, dudt_pstrain
    !-------------------------------------------------------------------------!

    dudt_plastic = 0.0
    dudt_pstrain = 0.0


    LocDegFr = nDegFr
    angfric = ATAN(BulkFriction) !angle of friction
    relaxtime = dt/(Tv) !Tv=dx/V_s with dx=min(dx);  Tv smaller-> stronger plasticity


! ---[ Calculate trial stress tensor ]---
! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz -> need to be specified throughout the whole medium for every element 
! as material values in EQN%IniStress, are mapped to the basis functions in dg_setup

    Stress(:,1:6)= dgvar(:,1:6)  + DOFStress(:,1:6)   !act.Stress + initial stress_xx

! ---[ Calculate trial yield stress ]---

    ! Mean stress
    meanStress(1:LocDegFr) = (Stress(1:LocDegFr,1) + Stress(1:LocDegFr,2)+ Stress(1:LocDegFr,3) )/3

    ! Deviatoric stress
    devStress(1:LocDegFr,1) = Stress(1:LocDegFr,1) - meanStress(1:LocDegFr)
    devStress(1:LocDegFr,2) = Stress(1:LocDegFr,2) - meanStress(1:LocDegFr)
    devStress(1:LocDegFr,3) = Stress(1:LocDegFr,3) - meanStress(1:LocDegFr)

    devStress(1:LocDegFr,4:6) = Stress(1:LocDegFr,4:6)

    ! Second invariant of stress deviator
    secInv = 0.5*(devStress(1,1)**2+devStress(1,2)**2+devStress(1,3)**2)+devStress(1,4)**2+devStress(1,5)**2+devStress(1,6)**2 

    ! Scalar measure of shear stress
    tau= SQRT(secInv)

    ! Yield stress   
    taulim = PlastCo*COS(angfric) - meanStress(1)*SIN(angfric)! minus before sinus is for compressional stress=negative.
    taulim = MAX(0.0, taulim)


    ! Stress deviators are adjusted

    IF (tau .GT. taulim) THEN !plastic behaviour, else: elastic and stress tensor=trial stress tensor
       yldfac = 1.0D0- (1.0D0 - taulim/tau)*(1.0D0 - EXP(-relaxtime)) !factor by Duan/Day

       DO iDegFr=1,LocDegFr

           Stress(iDegFr,1) = devStress(iDegFr,1)*yldfac + meanStress(iDegFr)
           Stress(iDegFr,2) = devStress(iDegFr,2)*yldfac + meanStress(iDegFr)
           Stress(iDegFr,3) = devStress(iDegFr,3)*yldfac + meanStress(iDegFr)
           Stress(iDegFr,4:6) = devStress(iDegFr,4:6)*yldfac


           !----Change of dofs----- 
           dudt_plastic(iDegFr,1:6) = dgvar(iDegFr,1:6) - (Stress(iDegFr,1:6) - DOFStress(iDegFr,1:6))
           !dudt_pstrain(iDegFr,1:6) = ((1-yldfac)/mu)*devStress(iDegFr, 1:6) strain for DOF-wise calculation

       ENDDO
          dudt_pstrain(1:6) = ((1-yldfac)/mu)*devStress(1, 1:6)
        
    ENDIF !yield criterion check

 END SUBROUTINE Plasticity_3D

END MODULE Plasticity_mod
