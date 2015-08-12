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
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE InitialField
     MODULE PROCEDURE InitialField
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
    INTEGER                        :: iElem
    ! Local variable declaration
    INTEGER                        :: i,j
    INTEGER                        :: iVar, iZone, counter
    REAL                           :: cs, cp, rho0
    REAL                           :: RA(EQN%nVar,EQN%nVar), TT(EQN%nVar,EQN%nVar)
    REAL                           :: LambdaA(EQN%nVar)
    REAL                           :: n(3), t1(3), t2(3)
    REAL                           :: kx, ky, kz
    REAL                           :: amplitude,hwidth(3)
    REAL                           :: du
    REAL                           :: dx(3), dxt(3), TrafoMatrix(3,3)              
    REAL                           :: omega, tau
    COMPLEX                        :: IU
    !
    !--------------------------------------------------------------------------
    !--------------------------------------------------------------------------
    !INTENT(IN)  :: x, y, z, iElem, EQN, IC, SOURCE, IO
    !INTENT(OUT) :: Variable
    !--------------------------------------------------------------------------
    !
    SELECT CASE(TRIM(IC%cICType))                                          !
       !                                                                   !
    CASE('Var_Gauss_Puls')
       Variable(:)     = 0.                                                ! Init. with homogeneous background  
       cp = SQRT((EQN%lambda+2.*EQN%mu)/(EQN%rho0))
       cs = SQRT((EQN%mu)/(EQN%rho0))
       amplitude = IC%GP%amplitude                                         ! 
       hwidth    = IC%GP%hwidth                                            ! 
       !
       dx(1)       = x-IC%GP%xc(1)
       dx(2)       = y-IC%GP%xc(2)
       dx(3)       = z-IC%GP%xc(3)
       !
       ! Transpose of Rotation Matrix T => (T^-1)
       TrafoMatrix(1,:) = IC%GP%n(:)
       TrafoMatrix(2,:) = IC%GP%t1(:)
       TrafoMatrix(3,:) = IC%GP%t2(:)           
       !
       ! Compute pointing vector in rotated coordinate system frame
       ! 
       dxt(:)           = MATMUL(TrafoMatrix(:,:),dx(:))
       !
       DO i = 1, IC%GP%SetVar                                              
           Variable(IC%GP%varfield(i)) = Variable(IC%GP%varfield(i))  +   &  
                IC%GP%ampfield(i)*EXP(-0.5*(                              &   
                ((dxt(1))/IC%GP%hwidth(1))**2      + &                       
                ((dxt(2))/IC%GP%hwidth(2))**2      + &                      
                ((dxt(3))/IC%GP%hwidth(3))**2 ) )                           
       ENDDO   
       !
    CASE ('Char_Gauss_Puls','Char_Ricker_Puls')                            ! characteristic pulse
       !                
       Variable(:)     = 0.                                                ! Init. with homogeneous background
       Variable_ANE(:) = 0.                                                ! Init. with homogeneous background
       !
       amplitude = IC%GP%amplitude                                         ! amplitude of pulse
       hwidth    = IC%GP%hwidth                                            ! halfwidth of pulse
       !                                                                   !
       n(:)      = IC%GP%n(:)                                              ! normal vector of plane wave
       t1(:)     = IC%GP%t1(:)                                             ! tangential vector 1 of wave
       t2(:)     = IC%GP%t2(:)                                             ! tangential vector 2 of wave
       !
       cp = SQRT((EQN%lambda+2.*EQN%mu)/(EQN%rho0))                        ! p-wave velocity
       cs = SQRT((EQN%mu)/(EQN%rho0))                                      ! s-wave velocity
       !
       ! Order of Eigenvalues: (-cp,-cs,-cs,0,0,0,+cs,+cs,+cp)
       LambdaA(:) = (/ -cp, -cp, -cs, 0., 0., 0., cs, cs, cp /)            ! Eigenvalues of elastic system
       !
       RA(:,:)  = 0.                                                       ! Initialize with zero
       !
       ! Eigenvectors
       RA(1,1) =  (EQN%lambda + 2.*EQN%mu)/cp
       RA(2,1) =  EQN%lambda/cp
       RA(3,1) =  EQN%lambda/cp
       RA(7,1) =  1.
       RA(4,2) =  EQN%mu/cs
       RA(8,2) =  1.
       RA(6,3) =  EQN%mu/cs
       RA(9,3) =  1.
       RA(5,4) =  1.
       RA(2,5) =  1.
       RA(3,6) =  1.
       RA(6,7) =  EQN%mu/(-cs)
       RA(9,7) =  1.
       RA(4,8) =  EQN%mu/(-cs)
       RA(8,8) =  1.
       RA(1,9) =  (EQN%lambda + 2.*EQN%mu)/(-cp)
       RA(2,9) =  EQN%lambda/(-cp)
       RA(3,9) =  EQN%lambda/(-cp)
       RA(7,9) =  1.
       !
       TT(:,:)  = 0.                                                       ! Initialize with zero
       ! Tensor Rotation Matrix depending on the direction of the wave
       TT(1,1) = n(1)*n(1)      
       TT(1,2) = t1(1)*t1(1)
       TT(1,3) = t2(1)*t2(1)
       TT(1,4) = 2*n(1)*t1(1)
       TT(1,5) = 2*t1(1)*t2(1)
       TT(1,6) = 2*n(1)*t2(1)
       TT(2,1) = n(2)*n(2)
       TT(2,2) = t1(2)*t1(2)
       TT(2,3) = t2(2)*t2(2)
       TT(2,4) = 2*n(2)*t1(2)
       TT(2,5) = 2*t1(2)*t2(2)
       TT(2,6) = 2*n(2)*t2(2)
       TT(3,1) = n(3)*n(3)
       TT(3,2) = t1(3)*t1(3)
       TT(3,3) = t2(3)*t2(3)
       TT(3,4) = 2*n(3)*t1(3)
       TT(3,5) = 2*t1(3)*t2(3)
       TT(3,6) = 2*n(3)*t2(3)
       TT(4,1) = n(2)*n(1)
       TT(4,2) = t1(2)*t1(1)
       TT(4,3) = t2(2)*t2(1)
       TT(4,4) = n(2)*t1(1)+n(1)*t1(2)
       TT(4,5) = t1(2)*t2(1)+t1(1)*t2(2)
       TT(4,6) = n(2)*t2(1)+n(1)*t2(2)
       TT(5,1) = n(3)*n(2)
       TT(5,2) = t1(3)*t1(2)
       TT(5,3) = t2(3)*t2(2)
       TT(5,4) = n(3)*t1(2)+n(2)*t1(3)
       TT(5,5) = t1(3)*t2(2)+t1(2)*t2(3)
       TT(5,6) = n(3)*t2(2)+n(2)*t2(3)
       TT(6,1) = n(3)*n(1)
       TT(6,2) = t1(3)*t1(1)
       TT(6,3) = t2(3)*t2(1)
       TT(6,4) = n(3)*t1(1)+n(1)*t1(3)
       TT(6,5) = t1(3)*t2(1)+t1(1)*t2(3)
       TT(6,6) = n(3)*t2(1)+n(1)*t2(3)
       TT(7,7) = n(1)
       TT(7,8) = t1(1)
       TT(7,9) = t2(1)
       TT(8,7) = n(2)
       TT(8,8) = t1(2)
       TT(8,9) = t2(2)
       TT(9,7) = n(3)
       TT(9,8) = t1(3)
       TT(9,9) = t2(3)
       !
       ! Rotating the eigenvectors into the desired direction
       RA(:,:) = MATMUL(TT,RA)
       !
       ! Compute vector pointing from center of Gausspulse to current location (x,y,z)
       dx(1)       = x-IC%GP%xc(1)
       dx(2)       = y-IC%GP%xc(2)
       dx(3)       = z-IC%GP%xc(3)
       !
       ! Vector Rotation Matrix
       TrafoMatrix(1,:) = IC%GP%n(:)
       TrafoMatrix(2,:) = IC%GP%t1(:)
       TrafoMatrix(3,:) = IC%GP%t2(:)
       !
       ! Compute pointing vector in rotated coordinate system frame
       dxt(:)           = MATMUL(TrafoMatrix(:,:),dx(:))
       ! 
       SELECT CASE(TRIM(IC%cICType))                                           !
       !
       CASE ('Char_Gauss_Puls') 
           DO i = 1, IC%GP%SetVar                                              !
               dxt(1) = dxt(1) - LambdaA(IC%GP%varfield(i))*time
               Variable(:) = Variable(:)  + RA(:,IC%GP%varfield(i))*        &  !
                    IC%GP%ampfield(i)*EXP(-0.5*(                            &  !
                    ((dxt(1))/IC%GP%hwidth(1))**2      + &                     !
                    ((dxt(2))/IC%GP%hwidth(2))**2      + &                     !
                    ((dxt(3))/IC%GP%hwidth(3))**2 ) )                          !
           ENDDO   
       !                                                                       !
       CASE ('Char_Ricker_Puls') 
           DO i = 1, IC%GP%SetVar
               dxt(1) = dxt(1) - LambdaA(IC%GP%varfield(i))*time
               tau =    (dxt(1)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(1) )**2   + &   !
                        (dxt(2)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(2) )**2   + &   !
                        (dxt(3)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(3) )**2         !
                                               !
               Variable(:) = Variable(:)  + RA(:,IC%GP%varfield(i))*                  &   !
                    IC%GP%ampfield(i) * (1 - 2*tau) * EXP(-tau)   
               ! Derivative of Gaussian
               ! Variable(:) = Variable(:)  + RA(:,IC%GP%varfield(i))*        &  !              
               !     IC%GP%ampfield(i) * (- 2* (dxt(1)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(1) )   + &
               !                               (dxt(2)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(2) )   + &
               !                               (dxt(3)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(3) )) * EXP(-tau)
           ENDDO
       
       END SELECT
       !
       !                                                                   !                
       !                                                                   
    CASE ('Gauss_Puls_Rad')                                                    ! Gauss_Puls_Rad for elastic wave equation
 
       Variable(:)     = 0.                                                ! Init. with homogeneous background  
       Variable_ANE(:) = 0.                                                !

 
       amplitude = IC%GP%amplitude                                             
       hwidth    = IC%GP%hwidth                                                 
       !                                                                        
       ! Compute normal vector
       n(1)      =  x-IC%GP%xc(1)                                       
       n(2)      =  y-IC%GP%xc(2)                                       
       n(3)      =  z-IC%GP%xc(3)
       ! Compute length and normalize vector
       du        =  SQRT(SUM(n(:)**2))
       IF(du.GT.0) THEN
          n(:)   =  n(:)/du
       ELSE
          n(:)   = 0.
       ENDIF
       !
       Variable(:) = 0. 
       ! Set the speed vectors (last variables in elastic equations) into the direction of the normal vector
       Variable((EQN%nVar-EQN%Dimension+1):EQN%nVar) = n(:)
       !                                                                    
       Variable((EQN%nVar-EQN%Dimension+1):EQN%nVar) = Variable((EQN%nVar-EQN%Dimension+1):EQN%nVar)*amplitude*  &                  
           EXP(-0.5*(   ((x-IC%GP%xc(1))/IC%GP%hwidth(1))**2      + &                  
                        ((y-IC%GP%xc(2))/IC%GP%hwidth(2))**2      + &                  
                        ((z-IC%GP%xc(3))/IC%GP%hwidth(3))**2 )      )  
       !
       !
       ! Add the homogeneous background
       !
       ! Variable(:) = Variable(:) + IC%GP%Um(:)
       !                                                                   
    CASE('Planarwave')                                                     ! Planarwave
       !
       Variable(:) = IC%PW%Um(:)                                           ! Init. with homogeneous background  
       Variable_ANE(:) = 0.                                                !
       !
       n         = IC%PW%n(:)                                              ! Planarwave 
       !
           cp = SQRT((EQN%lambda+2.*EQN%mu)/(EQN%rho0))
           cs = SQRT((EQN%mu)/(EQN%rho0))

       !
       ! Order of Eigenvalues: (-cp,-cs(y),-cs(z),0,0,0,+cs(y),+cs(z),+cp)
       !
       LambdaA(:) = (/ -cp, -cs, -cs, 0., 0., 0., cs, cs, cp /)
       !
       RA(1,1) = EQN%rho0*(-2*n(2)**2*EQN%mu-2*n(3)**2*EQN%mu+EQN%lambda+2*EQN%mu)      
       RA(1,2) = -2*EQN%mu*n(2)*EQN%rho0*n(1)**2*n(3)
       RA(1,3) = -2*n(3)**2*EQN%mu*EQN%rho0*n(1)**2
       RA(1,4) = -n(3)**2*n(1)
       RA(1,5) = 0
       RA(1,6) = -n(2)*n(1)*n(3)
       RA(1,7) = 2*n(3)**2*EQN%mu*EQN%rho0*n(1)**2
       RA(1,8) = 2*EQN%mu*n(2)*EQN%rho0*n(1)**2*n(3)
       RA(1,9) = -EQN%rho0*(-2*n(2)**2*EQN%mu-2*n(3)**2*EQN%mu+EQN%lambda+2*EQN%mu)
       
       RA(2,1) = EQN%rho0*(2*n(2)**2*EQN%mu+EQN%lambda)
       RA(2,2) = 2*EQN%mu*n(2)*EQN%rho0*n(1)**2*n(3)
       RA(2,3) = 0
       RA(2,4) = 0
       RA(2,5) = -n(3)**2/n(2)*n(1)**2
       RA(2,6) = -1/n(2)*n(1)**3*n(3)
       RA(2,7) = 0
       RA(2,8) = -2*EQN%mu*n(2)*EQN%rho0*n(1)**2*n(3)
       RA(2,9) = -EQN%rho0*(2*n(2)**2*EQN%mu+EQN%lambda)
       
       RA(3,1) = EQN%rho0*(2*n(3)**2*EQN%mu+EQN%lambda)
       RA(3,2) = 0
       RA(3,3) = 2*n(3)**2*EQN%mu*EQN%rho0*n(1)**2
       RA(3,4) = -n(1)**3
       RA(3,5) = -n(2)*n(1)**2
       RA(3,6) = 0
       RA(3,7) = -2*n(3)**2*EQN%mu*EQN%rho0*n(1)**2
       RA(3,8) = 0
       RA(3,9) = -EQN%rho0*(2*n(3)**2*EQN%mu+EQN%lambda)
       
       RA(4,1) = 2*n(2)*EQN%mu*n(1)*EQN%rho0
       RA(4,2) = -EQN%mu*EQN%rho0*n(1)*(2*n(2)**2+n(3)**2-1)*n(3)
       RA(4,3) = -EQN%mu*n(2)*n(3)**2*EQN%rho0*n(1)
       RA(4,4) = 0
       RA(4,5) = 0
       RA(4,6) = n(1)**2*n(3)
       RA(4,7) = EQN%mu*n(2)*n(3)**2*EQN%rho0*n(1)
       RA(4,8) = EQN%mu*EQN%rho0*n(1)*(2*n(2)**2+n(3)**2-1)*n(3)
       RA(4,9) = -2*n(2)*EQN%mu*n(1)*EQN%rho0
       
       RA(5,1) = 2*EQN%mu*n(2)*EQN%rho0*n(3)
       RA(5,2) = n(3)**2*EQN%mu*EQN%rho0*n(1)**2
       RA(5,3) = EQN%mu*n(2)*EQN%rho0*n(1)**2*n(3)
       RA(5,4) = 0
       RA(5,5) = n(1)**2*n(3)
       RA(5,6) = 0
       RA(5,7) = -EQN%mu*n(2)*EQN%rho0*n(1)**2*n(3)
       RA(5,8) = -n(3)**2*EQN%mu*EQN%rho0*n(1)**2
       RA(5,9) = -2*EQN%mu*n(2)*EQN%rho0*n(3)
       
       RA(6,1) = 2*EQN%mu*n(1)*EQN%rho0*n(3)
       RA(6,2) = -EQN%mu*n(2)*n(3)**2*EQN%rho0*n(1)
       RA(6,3) = EQN%mu*EQN%rho0*n(1)*(-2*n(3)**2-n(2)**2+1)*n(3)
       RA(6,4) = n(1)**2*n(3)
       RA(6,5) = 0
       RA(6,6) = 0
       RA(6,7) = -EQN%mu*EQN%rho0*n(1)*(-2*n(3)**2-n(2)**2+1)*n(3)
       RA(6,8) = EQN%mu*n(2)*n(3)**2*EQN%rho0*n(1)
       RA(6,9) = -2*EQN%mu*n(1)*EQN%rho0*n(3)
       
       RA(7,1) = n(1)*sqrt(EQN%rho0*(EQN%lambda+2*EQN%mu))
       RA(7,2) = -n(2)*n(1)*n(3)*sqrt(EQN%rho0*EQN%mu)
       RA(7,3) = -n(3)**2*n(1)*sqrt(EQN%rho0*EQN%mu)
       RA(7,4) = 0
       RA(7,5) = 0
       RA(7,6) = 0
       RA(7,7) = -n(3)**2*n(1)*sqrt(EQN%rho0*EQN%mu)
       RA(7,8) = -n(2)*n(1)*n(3)*sqrt(EQN%rho0*EQN%mu)
       RA(7,9) = n(1)*sqrt(EQN%rho0*(EQN%lambda+2*EQN%mu))
       
       RA(8,1) = n(2)*sqrt(EQN%rho0*(EQN%lambda+2*EQN%mu))
       RA(8,2) = n(1)**2*n(3)*sqrt(EQN%rho0*EQN%mu)
       RA(8,3) = 0
       RA(8,4) = 0
       RA(8,5) = 0
       RA(8,6) = 0
       RA(8,7) = 0
       RA(8,8) = n(1)**2*n(3)*sqrt(EQN%rho0*EQN%mu)
       RA(8,9) = n(2)*sqrt(EQN%rho0*(EQN%lambda+2*EQN%mu))
       
       RA(9,1) = sqrt(EQN%rho0*(EQN%lambda+2*EQN%mu))*n(3)
       RA(9,2) = 0
       RA(9,3) = n(1)**2*n(3)*sqrt(EQN%rho0*EQN%mu)
       RA(9,4) = 0
       RA(9,5) = 0
       RA(9,6) = 0
       RA(9,7) = n(1)**2*n(3)*sqrt(EQN%rho0*EQN%mu)
       RA(9,8) = 0
       RA(9,9) = sqrt(EQN%rho0*(EQN%lambda+2*EQN%mu))*n(3)           
       
       !                                                                    !
       DO i = 1, IC%PW%SetVar                                               !
           omega = LambdaA(IC%PW%varfield(i))*SQRT(SUM(IC%PW%k_vec(:)**2))
           Variable(:) = Variable(:)  + RA(:,IC%PW%varfield(i))*        &   !
                IC%PW%ampfield(i)*SIN(IC%PW%k_vec(1)*x +                &   ! 
                                      IC%PW%k_vec(2)*y +                &   !
                                      IC%PW%k_vec(3)*z - omega*time     )   !
       ENDDO                                                                !
       !
    CASE('Planarwave_Gauss_Puls','Planarwave_Ricker_Puls')          ! CASE Planarwave anelastic
       ! Add the homogeneous background
      
       Variable(:) = IC%GP%Um(:)                                                ! Init. with homogeneous background  
       Variable_ANE(:) = 0. 

       amplitude = IC%GP%amplitude                                         ! 
       hwidth    = IC%GP%hwidth                                            ! 
       !                                                                   !
       !                                                                   ! 
       n(:)      = IC%GP%n(:)                                              ! 

       LambdaA(:)= REAL(IC%PWAN%EigenVal(:))      
       RA(:,:)   = REAL(IC%PWAN%EigenVec(:,:))
       !
       ! Compute vector pointing from center of Gausspulse to current location (x,y,z)
       !
       dx(1)       = x-IC%GP%xc(1)
       dx(2)       = y-IC%GP%xc(2)
       dx(3)       = z-IC%GP%xc(3)
       !
       ! Transpose of Rotation Matrix T => (T^-1)
       TrafoMatrix(1,:) = IC%GP%n(:)
       TrafoMatrix(2,:) = IC%GP%t1(:)
       TrafoMatrix(3,:) = IC%GP%t2(:)           
       !
       ! Compute pointing vector in rotated coordinate system frame
       ! 
       dxt(:)           = MATMUL(TrafoMatrix(:,:),dx(:))
       ! 
       SELECT CASE(TRIM(IC%cICType))                                           !
       !
       CASE ('Planarwave_Gauss_Puls') 
           DO i = 1, IC%GP%SetVar                                              !
               dxt(1) = dxt(1) - LambdaA(IC%GP%varfield(i))*time
               Variable(:) = Variable(:)  + RA(:,IC%GP%varfield(i))*        &  !
                    IC%GP%ampfield(i)*EXP(-0.5*(                            &  ! 
                    ((dxt(1))/IC%GP%hwidth(1))**2      + &                     !  
                    ((dxt(2))/IC%GP%hwidth(2))**2      + &                     ! 
                    ((dxt(3))/IC%GP%hwidth(3))**2 ) )                          ! 
           ENDDO   
       !                                                                       !
       CASE ('Planarwave_Ricker_Puls') 
           DO i = 1, IC%GP%SetVar
               dxt(1) = dxt(1) - LambdaA(IC%GP%varfield(i))*time
               tau =    (dxt(1)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(1) )**2   + &   !  
                        (dxt(2)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(2) )**2   + &   ! 
                        (dxt(3)/LambdaA(IC%GP%varfield(i)) * IC%GP%hwidth(3) )**2         ! 
                                               !
               Variable(:) = Variable(:)  + RA(:,IC%GP%varfield(i))*        &  !              
                    IC%GP%ampfield(i) * (1 - 2*tau) * EXP(-tau)   
           ENDDO
       !
       END SELECT
       !   
    CASE('PlanarwaveAnel','PlanarwaveAn')                                       ! CASE Planarwave anelastic
       ! Add the homogeneous background
      
       Variable(:) = IC%PW%Um(:)                                                ! Init. with homogeneous background  
       Variable_ANE(:) = 0. 
      
       ! Wavenumbers      
            kx    = IC%PWAN%Wavenumbers(1)
            ky    = IC%PWAN%Wavenumbers(2)
            kz    = IC%PWAN%Wavenumbers(3)

       ! Imaginary Unit
       IU = (0.,1.)

         DO i = 1, IC%PW%SetVar                                                   !
           Variable(:) = Variable(:)  +                                     &
                         real(IC%PW%ampfield(i) * IC%PWAN%EigenVec(1:EQN%nVar,IC%PW%varfield(i)) * &
                                   exp( IU * ( - kx*x - ky*y - kz*z )) ) 
           Variable_ANE(:) = Variable_ANE(:)  +                             &
                         real(IC%PW%ampfield(i) * IC%PWAN%EigenVec(EQN%nVar+1:EQN%nVarTotal,IC%PW%varfield(i)) * &
                                   exp( IU * ( - kx*x - ky*y - kz*z )) ) 
       ENDDO
       !   
    CASE('PlanarwaveAniso')                                                     ! CASE Planarwave anelastic
       ! Add the homogeneous background
      
       Variable(:) = IC%PW%Um(:)                                                ! Init. with homogeneous background  
       Variable_ANE(:) = 0. 
       ! Imaginary Unit
       IU = (0.,1.)
 
       DO j = 1,3
           ! Wavenumbers      
              kx    = IC%PWANISO(j)%Wavenumbers(1)
              ky    = IC%PWANISO(j)%Wavenumbers(2)
              kz    = IC%PWANISO(j)%Wavenumbers(3)

           DO i = 1, IC%PWANISO(j)%SetVar                                                   !
             Variable(:) = Variable(:)  +                                     &
                           real(IC%PWANISO(j)%ampfield(i) * IC%PWANISO(j)% &
				EigenVec(1:EQN%nVar,IC%PWANISO(j)%varfield(i)) * &
                                   exp( IU * ( - kx*x - ky*y - kz*z )) ) 
             Variable_ANE(:) = Variable_ANE(:)  +                             &
                           real(IC%PWANISO(j)%ampfield(i) * IC%PWANISO(j)% &
				EigenVec(EQN%nVar+1:EQN%nVarTotal,IC%PWANISO(j)%varfield(i)) * &
                                   exp( IU * ( - kx*x - ky*y - kz*z )) ) 
           ENDDO
       ENDDO
       !   
    CASE DEFAULT                                                                ! DEFAULT  
       logError(*) 'InitialField: none of the possible initial conditions was chosen'
       logError(*) TRIM(IC%cICType),'|'                         ! DEFAULT
       STOP 
    END SELECT                                                                  !

  END SUBROUTINE InitialField


END MODULE COMMON_InitialField_mod
