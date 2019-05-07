!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Stephanie Wollherr
!!
!! @section LICENSE
!! Copyright (c) 2007-2016, SeisSol Group
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

  use iso_c_binding, only: c_loc
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  !PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE Plasticity
     MODULE PROCEDURE Plasticity_3D_high
  END INTERFACE
    INTERFACE Plasticity
     MODULE PROCEDURE Plasticity_3D_avg
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC  :: Plasticity_3D_high
  PUBLIC  :: Plasticity_3D_avg
  !---------------------------------------------------------------------------!
  CONTAINS


  SUBROUTINE Plasticity_3D_high(dgvar, DOFStress, nDegFr, nAlignedDegFr, &
                                Tv, dt, mu, lambda, parameters , Energy, pstrain, &
                                intGaussP, intGaussW, DISC, nVar, nIntGP)


    !-------------------------------------------------------------------------!
  USE DGBasis_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization)    :: DISC
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom
    INTEGER     :: nDegFr
    INTEGER     :: iIntGP, nIntGP
    INTEGER     :: nVar, iPoly
    integer     :: nAlignedDegFr


    REAL        :: stateGP(nVar)                                              !State at GP's
    REAL        :: pstateGP(nVar)                                             !Primitive state at GP's
    REAL        :: Stress_total(nIntGP,6)                                     !local stress variable for yield criterion
    REAL        :: Strain_total(nIntGP,6)                                     !local strain variable for elastic strain energy
    REAL        :: Strain_ini(1:6)                                            !local initial strain variable for elastic strain energy
    REAL        :: devStress(1:nIntGP,6)                                      !deviatoric stress
    REAL        :: meanStress(1:nIntGP)                                       !mean stress
    REAL        :: angfric, yldfac                                            !angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !timestep, relaxation time for the yield factor
    REAL        :: mu, lambda                                                 !Lame parameters
    REAL        :: LocnVar
    REAL        :: tau(1:nIntGP)
    REAL        :: taulim(1:nIntGP)                                           !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv(1:nIntGP)                                           !second Invariant of deviatoric stress
    REAL        :: Tv                                                         !relaxtime
    REAL        :: DOFStress(1:nDegFr,1:6)                                    !initial loading
    REAL        :: dgvar(1:nAlignedDegFr,1:9)                                 !dof's
    REAL        :: dgvar_new(1:nAlignedDegFr,1:6)                             !
    REAL        :: phi                                                        !Value of the base function at GP
    LOGICAL     :: check
    REAL        :: update(1:nIntGP,6)
    REAL        :: newstateGP(1:6)
    REAL        :: dudt_plastic(1:nDegFr,1:6)                                 !stress change due to plasticity
    REAL        :: dudt_pstrain(1:6)                                          !change of plastic strain
    REAL        :: pstrain(1:7)                                               !plastic strain
    REAL        :: estrain(1:6), estrain_ini(1:6)                             !total elastic strain
    REAL        :: PlasticEnergy_tmp, KineticEnergy_tmp, EstrainEnergy_tmp    !temp. energies
    REAL        :: Energy(1:2)                                                !plastic and elastic strain energy
    REAL        :: parameters(1:4)                                            !1=volume of the triangle, 2=plastcohesion, 3=density rho, 4 = bulk friction
    REAL        :: I1,I1_0,I2,I2_0                                            !first and second invariants of strains
    REAL        :: IntGaussP(:,:)
    REAL        :: IntGaussW(:)
    REAL, POINTER :: IntGPBaseFunc(:,:) =>NULL()
    REAL, POINTER :: MassMatrix(:,:)    =>NULL()


    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DOFStress, nDegFr, nAlignedDegFr, Tv, dt, &
                     mu, lambda, parameters, intGaussP, intGaussW, DISC, nVar, nIntGP
    INTENT(INOUT) :: dgvar, Energy, pstrain
    !-------------------------------------------------------------------------!

    !set energies to zero
    dudt_plastic = 0.0D0
    dudt_pstrain = 0.0D0
    Stress_total = 0.0D0
    Energy(1:2)  = 0.0D0

    angfric = ATAN(parameters(4)) !angle of friction = atan(bulk friction)

    IF (Tv .GT. 0) THEN
       relaxtime = 1.0D0 - EXP(-dt/(Tv)) !Tv: direct input via parameter file; Tv smaller-> stronger plasticity
    ELSE
       relaxtime = 1.0
    ENDIF

    iPoly = DISC%Galerkin%nPoly
    ! Basis func values
    IntGPBaseFunc => DISC%Galerkin%IntGPBaseFunc_Tet(1:nDegFr,1:nIntGP,iPoly)
    ! Mass matrix
    MassMatrix    => DISC%Galerkin%MassMatrix_Tet(1:nDegFr,1:nDegFr,iPoly)


    ! ---[ Calculate trial stress tensor ]---
    ! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz
    ! -> need to be specified throughout the whole medium for every element
    ! as material values in EQN%IniStress, are mapped to the basis functions in dg_setup

    ! GP-wise

    DO iIntGP = 1, nIntGP
        !1. get the state at every GP
        stateGP(:) = 0.
        DO iDegFr = 1, nDegFr
           phi = IntGPBaseFunc(iDegFr,iIntGP)
           stateGP(1:nVar) = stateGP(1:nVar) + phi*dgvar(iDegFr,1:nVar)
        ENDDO

        !2. add up all initial loading
        pstateGP(:) = stateGP(:)
        !dofstress are in this case just the elementwise initial stresses
        Stress_total(iIntGP,1:6) = pstateGP(1:6) + DOFStress(1,1:6)

        !Calculate the total strain from the elastic stress-strain relation
        Strain_total(iIntGP, 1:6) = MATMUL(DISC%Galerkin%Strain_matrix, Stress_total(iIntGP,1:6))
    ENDDO

    !Calculate initial strain loading from initial stress loading (elementwise)
    !-> move that outside the routine and calculate before as it is constant over time
    Strain_ini(1:6) = MATMUL(DISC%Galerkin%Strain_matrix,DOFStress(1,1:6))

    ! Mean stress, GP-wise
    meanStress(1:nIntGP) = (Stress_total(:,1) + Stress_total(:,2)+ Stress_total(:,3) )/3

    ! Deviatoric stress, GP-wise
    devStress(1:nIntGP,1) = Stress_total(:,1) - meanStress(:)
    devStress(1:nIntGP,2) = Stress_total(:,2) - meanStress(:)
    devStress(1:nIntGP,3) = Stress_total(:,3) - meanStress(:)

    devStress(1:nIntGP,4:6) = Stress_total(:,4:6)


    ! Second invariant of stress deviator
    secInv(1:nIntGP) = 0.5*(devStress(:,1)**2 + devStress(:,2)**2 + devStress(:,3)**2) + &
                       devStress(:,4)**2 + devStress(:,5)**2 + devStress(:,6)**2

    ! Scalar measure of shear stress
    tau(1:nIntGP)= SQRT(secInv(1:nIntGP))

    ! Yield stress
    ! minus in front of meanstress is for compressional stress=negative.
    taulim(1:nIntGP) = parameters(2)*COS(angfric) - meanStress(1:nIntGP)*SIN(angfric)
    taulim(1:nIntGP) = MAX(0.0, taulim(1:nIntGP))

    check = .false.
    ! Stress deviators are adjusted
    DO iIntGP = 1, nIntGP

       IF (tau(iIntGP) .GT. taulim(iIntGP)) THEN !plastic behaviour, else: elastic and stress tensor=trial stress tensor
           check = .TRUE.
           yldfac = 1.0D0- (1.0D0 - taulim(iIntGP)/tau(iIntGP))*(relaxtime) !factor by Duan/Day 2008

           ! adjustment of stresses, GP-wise for every variable 1-6
           Stress_total(iIntGP,1) = devStress(iIntGP,1)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,2) = devStress(iIntGP,2)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,3) = devStress(iIntGP,3)*yldfac + meanStress(iIntGP)
           Stress_total(iIntGP,4:6) = devStress(iIntGP,4:6)*yldfac
       ENDIF

       !subtract the inital loading; new stress state at every GP
       update(iIntGP,1:6) = Stress_total(iIntGP, 1:6)-DOFStress(1,1:6)
    ENDDO !adjustment over all GP-points


    IF (.NOT. check) THEN !nothing was adjusted
       dudt_plastic(1:nDegFr,1:6) = 0.0
      !
    ELSE !back projection is necessary because at least one GP was adjusted
       dgvar_new = 0.0

       ! back projection to the DOFs
       DO iIntGP = 1, nIntGP
          DO iDegFr = 1, nDegFr
             phi = IntGPBaseFunc(iDegFr,iIntGP) !basis function number idegfr at point iIntGp
             newstateGP = update(iIntGP,1:6) !new value at that GP (=old value if no adjustment)
             !dgvar since it's not zero, only possible for the first iteration
             dgvar_new(iDegFr,1:6) =  dgvar_new(iDegFr,1:6) + IntGaussW(iIntGP)*newstateGP(1:6)*phi
          ENDDO
       ENDDO !nIntGP

       !divide by diagonal mass matrix entries
       DO iDegFr = 1, nDegFr
           dgvar_new(iDegFr,:) = dgvar_new(iDegFr,:) / MassMatrix(iDegFr,iDegFr)
       ENDDO

       !Update
       dudt_plastic(1:nDegFr,1:6) = dgvar(1:nDegFr,1:6)- dgvar_new(1:nDegFr,1:6)
       !just the first dof
       dudt_pstrain(1:6) = (1.0/(2.0*mu))*dudt_plastic(1,1:6)
       !both mu and 2*mu (= tensor formulation) are used in literature


    ENDIF !check = .true.

    dgvar(1:nDegFr,1:6) = dgvar(1:nDegFr,1:6) - dudt_plastic(1:nDegFr,1:6)

    !update plastic strain
    pstrain(1:6) = pstrain(1:6) + dudt_pstrain(1:6)
    !accumulated plastic strain
    pstrain(7) = pstrain(7)+ dt*sqrt(0.5*(dudt_pstrain(1)**2 + dudt_pstrain(2)**2 &
                 + dudt_pstrain(3)**2)+ dudt_pstrain(4)**2 + dudt_pstrain(5)**2 + dudt_pstrain(6)**2)

    !calculate energies
    PlasticEnergy_tmp = Stress_total(1,1)*dudt_pstrain(1) + Stress_total(1,2)*dudt_pstrain(2) + &
                        Stress_total(1,3)*dudt_pstrain(3) + 2.0*Stress_total(1,4)*dudt_pstrain(4) + &
                        2.0*Stress_total(1,5)*dudt_pstrain(5) + 2.0*Stress_total(1,6)*dudt_pstrain(6)

    !total elastic strain, if no plastic yielding -> elastic strain = total strain
    estrain(1:6) = Strain_total(1,1:6) - dudt_pstrain(1:6)
    estrain_ini(1:6) = Strain_ini(1:6)

    !first and second invariants of the total elastic strain
    I1 = estrain(1) + estrain(2) + estrain(3)
    I2 = estrain(1)**2 + estrain(2)**2 + estrain(3)**2 + 2.0*estrain(4)**2 + &
         2.0*estrain(5)**2 + 2.0*estrain(6)**2
    !first and second invariants of the initial strain loading
    I1_0 = estrain_ini(1) + estrain_ini(2) + estrain_ini(3)
    I2_0 = estrain_ini(1)**2 + estrain_ini(2)**2 + estrain_ini(3)**2 + 2.0*estrain_ini(4)**2 + &
           2.0*estrain_ini(5)**2 + 2.0*estrain_ini(6)**2
    !Elastic strain energy
    !subtracted the initial elastic strain
    EstrainEnergy_tmp = 0.5*lambda*(I1**2-I1_0**2) + mu*(I2-I2_0)

    Energy(1) = PlasticEnergy_tmp*parameters(1) !multiplied by volume to get integral over element
    Energy(2) = EstrainEnergy_tmp*parameters(1)

 END SUBROUTINE Plasticity_3D_high


  !yldfac is only caluclated from the first DOF, and all DOF's are adjusted by the same coefficient

  SUBROUTINE Plasticity_3D_avg(DISC, dgvar, DOFStress, nDegFr, nAlignedDegFr, &
                               Tv, dt, mu,lambda, parameters , &
                               Energy, pstrain)

    !-------------------------------------------------------------------------!

    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    TYPE(tDiscretization)    :: DISC
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER     :: iDegFr                                                     ! Index of degree of freedom       !
    INTEGER     :: nDegFr
    INTEGER     :: nAlignedDegFr

    REAL        :: Stress(1:nDegFr,6)                                         !local stress variable for the yield criterion
    REAL        :: Strain_ini(1:6)                                            !local initial strain variable for elastic strain energy
    REAL        :: Strain_total(1:6)                                          !local strain variable for elastic strain energy
    REAL        :: devStress(1:nDegFr,6)                                      !stress deviator for the yield criterion
    REAL        :: meanStress(1:nDegFr)                                       !mean stress
    REAL        :: angfric, yldfac                                            !Angle of friction, yield factor
    REAL        :: dt, relaxtime                                              !timestep, relaxation time for the yield factor
    REAL        :: mu, lambda                                                 !Lame parameters
    REAL        :: tau,taulim                                                 !tau= potential flow function, taulim= drucker-Prager yield stress.
    REAL        :: secInv                                                     !secInv=second Invariant of deviatoric stress
    REAL        :: Tv                                                         !relaxtime
    REAL        :: DOFStress(1:nDegFr,1:6)                                    !initial loading
    REAL        :: dgvar(1:nAlignedDegFr,1:9)                                 !dofs
    REAL        :: dudt_pstrain(1:6)                                          !change in plastic strain
    REAL        :: pstrain(1:7)                                               !plastic strain
    REAL        :: estrain(1:6), estrain_ini(1:6)                             !total elastic strain
    REAL        :: PlasticEnergy_tmp, EstrainEnergy_tmp
    REAL        :: Energy(1:2)
    REAL        :: parameters(1:4)                                            !1=volume of the triangle, 2=plastcohesion, 3=density rho, 4 = bulk friction
    REAL         :: I1,I1_0,I2,I2_0                                           !first and second invariants of strains
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, DOFStress, nDegFr, nAlignedDegFr, Tv, &
                     dt, mu, lambda, parameters
    INTENT(INOUT) :: dgvar, Energy, pstrain
    !-------------------------------------------------------------------------!

    dudt_pstrain = 0.0
    Energy(1:2) = 0.0


    angfric = ATAN(parameters(4)) !angle of friction = atan(bulk friction)

    !Tv: direct input via parameter file; Tv smaller-> stronger plasticity
    IF (Tv .GT. 0) THEN
       relaxtime = 1.0D0 - EXP(-dt/(Tv))
    ELSE
       relaxtime = 1.0
    ENDIF

    ! ---[ Calculate trial stress tensor ]---
    ! Stress tensor components: Initial sxx,syy,szz,sxy,sxz,syz
    ! -> need to be specified throughout the whole medium for every element
    ! as material values in EQN%IniStress, are mapped to the basis functions in dg_setup

    !calculated for every degfr as it is needed later
    Stress(1:nDegFr,1:6)= dgvar(1:nDegFr,1:6)  + DOFStress(1:nDegFr,1:6)   !act.Stress + initial stress_xx

    !Calculate the total strain from the elastic stress-strain relation
    Strain_total(1:6) = MATMUL(DISC%Galerkin%Strain_matrix, Stress(1,1:6))

    !Calculate initial strain loading from initial stress loading (elementwise)
    !-> move that outside the routine and calculate beforhand
    Strain_ini(1:6) = MATMUL(DISC%Galerkin%Strain_matrix, DOFStress(1,1:6))

    ! ---[ Calculate trial yield stress ]---

    ! Mean stress
    meanStress(1:nDegFr) = (Stress(1:nDegFr,1) + Stress(1:nDegFr,2)+ Stress(1:nDegFr,3) )*(1.0D0/3.0D0)

    ! Deviatoric stress
    devStress(1:nDegFr,1) = Stress(1:nDegFr,1) - meanStress(1:nDegFr)
    devStress(1:nDegFr,2) = Stress(1:nDegFr,2) - meanStress(1:nDegFr)
    devStress(1:nDegFr,3) = Stress(1:nDegFr,3) - meanStress(1:nDegFr)

    devStress(1:nDegFr,4:6) = Stress(1:nDegFr,4:6)

    ! Second invariant of stress deviator using only the first degfr = average
    secInv = 0.5*(devStress(1,1)**2 + devStress(1,2)**2 + devStress(1,3)**2) + devStress(1,4)**2 + &
             devStress(1,5)**2 + devStress(1,6)**2

    ! Scalar measure of shear stress
    tau= SQRT(secInv)

    ! Yield stress
    ! minus in front of meanstress is for compressional stress=negative.
    taulim = parameters(2)*COS(angfric) - meanStress(1)*SIN(angfric)
    taulim = MAX(0.0, taulim)


    ! Stress deviators are adjusted

    IF (tau .GT. taulim) THEN !plastic behaviour, else: elastic and stress tensor=trial stress tensor
       yldfac = 1.0D0- (1.0D0 - taulim/tau)*(relaxtime) !factor by Duan/Day


       Stress(1:nDegFr,1) = devStress(1:nDegFr,1)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,2) = devStress(1:nDegFr,2)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,3) = devStress(1:nDegFr,3)*yldfac + meanStress(1:nDegFr)
       Stress(1:nDegFr,4:6) = devStress(1:nDegFr,4:6)*yldfac


       !----update the dofs-----
       dgvar(1:nDegFr,1:6) = Stress(1:nDegFr,1:6) - DOFStress(1:nDegFr,1:6)

       !just the first dof
       dudt_pstrain(1:6) = ((1.0-yldfac)/(2.0*mu))*devStress(1, 1:6)
       !both mu and 2*mu (= tensor formulation) are used in literature

        
    ENDIF !yield criterion check

    !update plastic strain
    pstrain(1:6) = pstrain(1:6) + dudt_pstrain(1:6) !plastic strain tensor

    !accumulated plastic strain
    pstrain(7) = pstrain(7)+ dt*sqrt(0.5*(dudt_pstrain(1)**2 + dudt_pstrain(2)**2 + dudt_pstrain(3)**2) + &
                 dudt_pstrain(4)**2 + dudt_pstrain(5)**2 + dudt_pstrain(6)**2)

    !calculate energies
    PlasticEnergy_tmp = Stress(1,1)*dudt_pstrain(1) + Stress(1,2)*dudt_pstrain(2) + &
                        Stress(1,3)*dudt_pstrain(3) + 2.0*Stress(1,4)*dudt_pstrain(4) + &
                        2.0*Stress(1,5)*dudt_pstrain(5) + 2.0*Stress(1,6)*dudt_pstrain(6)

    estrain(1:6) = Strain_total(1:6) - dudt_pstrain(1:6) !total elastic strain, if no plastic yielding -> elastic strain = total strain
    estrain_ini(1:6) = Strain_ini(1:6)

    !first and second invariants of the total elastic strain
    I1 = estrain(1) + estrain(2) + estrain(3)
    I2 = estrain(1)**2 + estrain(2)**2 + estrain(3)**2 + 2.0*estrain(4)**2 + &
         2.0*estrain(5)**2 + 2.0*estrain(6)**2
    !first and second invariants of the initial strain loading
    I1_0 = estrain_ini(1) + estrain_ini(2) + estrain_ini(3)
    I2_0 = estrain_ini(1)**2 + estrain_ini(2)**2 + estrain_ini(3)**2 + &
           2.0*estrain_ini(4)**2 + 2.0*estrain_ini(5)**2 + 2.0*estrain_ini(6)**2
    !Elastic strain energy
    !subtracted the initial elastic strain
    EstrainEnergy_tmp = 0.5*lambda*(I1**2-I1_0**2) + mu*(I2-I2_0)

    Energy(1) = PlasticEnergy_tmp*parameters(1) !multiplied by volume to get integral over element
    Energy(2) = EstrainEnergy_tmp*parameters(1)


 END SUBROUTINE Plasticity_3D_avg


END MODULE Plasticity_mod
