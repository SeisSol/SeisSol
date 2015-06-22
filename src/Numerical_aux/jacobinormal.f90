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

  MODULE JacobiNormal_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE JacobiNormal3D
     MODULE PROCEDURE JacobiNormal3D
  END INTERFACE
  INTERFACE RotationMatrix3D
     MODULE PROCEDURE RotationMatrix3D
  END INTERFACE
  INTERFACE BuildJacobians3D_Ela
     MODULE PROCEDURE BuildJacobians3D_Ela
  END INTERFACE
  INTERFACE BuildJacobians3D_Ani
     MODULE PROCEDURE BuildJacobians3D_Ani
  END INTERFACE
  INTERFACE BuildJacobians3D_Poro
     MODULE PROCEDURE BuildJacobians3D_Poro
  END INTERFACE
  INTERFACE AnelasticAuxMatrix3D_Iso
     MODULE PROCEDURE AnelasticAuxMatrix3D_Iso
  END INTERFACE
  INTERFACE AnelasticAuxMatrix3D_Ani
     MODULE PROCEDURE AnelasticAuxMatrix3D_Ani
  END INTERFACE

  !---------------------------------------------------------------------------!
  PUBLIC  :: JacobiNormal3D 
  PUBLIC  :: RotationMatrix3D
  PUBLIC  :: BuildJacobians3D_Ela
  PUBLIC  :: BuildJacobians3D_Ani
  PUBLIC  :: BuildJacobians3D_Poro
  PUBLIC  :: AnelasticAuxMatrix3D_Iso
  PUBLIC  :: AnelasticAuxMatrix3D_Ani
  !---------------------------------------------------------------------------!

CONTAINS

  ! looking in the normal direction (nx,ny,nz) and returns the matrix and its
  ! absolute value (needed for Riemann-solver)

  SUBROUTINE JacobiNormal3D(A,absA,W0,n0,s0,t0,EQN,IO,w_speed,AniVec)
    !-------------------------------------------------------------------------!
    
    USE COMMON_operators_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations)   :: EQN                                                 ! 
    TYPE(tInputOutput) :: IO
    REAL               :: A(EQN%nVar,EQN%nVar)                                ! Jacobian in local (nx,ny) system
    REAL               :: absA(EQN%nVar,EQN%nVar)                             ! Abs. value of A in (nx,ny) sys. 
    REAL               :: W0(EQN%nBackgroundVar)                              ! State for eval A in glob. system 
    REAL               :: n0(3), s0(3), t0(3)                                 ! Normal and tangent vectors
    REAL               :: w_speed(EQN%nNonZeroEV)                             ! Wave speeds normal to the surface
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER         :: i,j,rus,PoroFlux                                       ! 
    REAL            :: rho, u, v, w, p, c0, cp, cs                            ! 
    REAL            :: mu, lambda
    REAL            :: LA(EQN%nVar,EQN%nVar)                                  ! Eigenvalue matrix
    REAL            :: RA(EQN%nVar,EQN%nVar)                                  ! Right eigenvector matrix
    REAL            :: iRA(EQN%nVar,EQN%nVar)                                 ! Inverse right eigenvector matrix
    REAL            :: Test(EQN%nVar,EQN%nVar)
    REAL            :: Unity(EQN%nVar,EQN%nVar)
    REAL            :: T(6,6), iT(6,6), TT(6,6)
    REAL            :: nx, ny, nz, sx, sy, sz, tx, ty, tz, RES(6)
    REAL            :: checksum
    REAL            :: C(6,6)                                                 ! Material constants (Voigt matrix)
    REAL            :: Voigt_rot(6,6)                                         ! Voigt matrix (for rotation)
    REAL            :: amax                                                   ! max. wavespeeds
    REAL            :: E(9), GG(3), HH(3), II(3), NORMA, XX(3,3)              ! elements of the anisotropic eigenvectors
    REAL            :: AniVec(3,3),VV(3)    
    REAL            :: temp(3,3), TOL
    INTEGER         :: ZEROS(3), J2, J3
    REAL            :: coefficients(4)
    COMPLEX         :: solution(3)
    REAL            :: Re_solution(3), Im_solution(3)
    REAL            :: K_F, K_S, K_mean, MM, Alpha(6), rho_S                                    ! Porous parameters
    REAL            :: rho_F, nu, Poro, Kappa(3), Tor(3), rho1(3), rho2(3), beta1(3), beta2(3)  ! Porous parameters
    REAL            :: X(2),Y(2),Z(2),r1(2),r7(2),r10(2),r11(2),r2(2),NORMA1,NORMA2
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: EQN, IO, W0, n0, s0, t0, w_speed
    INTENT(OUT)     :: A, absA, AniVec
    !-------------------------------------------------------------------------!

        ! Elastic wave equation
        !
        ! Filling of isotropic (G,H,I) values
        !
        GG = (/ 1., 0., 0./) / SQRT(W0(1))
        HH = (/ 0., 1., 0./) / SQRT(W0(1))
        II = (/ 0., 0., 1./) / SQRT(W0(1))
        !
        PoroFlux = 0

        IF(EQN%Poroelasticity.NE.0)THEN ! If material is non-porous, use elastic jacobians
          IF(W0(27).NE.0.)THEN
             PoroFlux = 1
          ENDIF
        ENDIF

        IF(EQN%Anisotropy.EQ.0.AND.EQN%Poroelasticity.EQ.0)THEN
            !
            ! Initialization
            !
            A(:,:)      = 0.
            LA(:,:)     = 0.
            RA(:,:)     = 0.
            iRA(:,:)    = 0.
            ! 
            rho    = W0(1)
            mu     = W0(2)
            lambda = W0(3)
            !
            cp = SQRT((lambda+2*mu)/rho)
            cs = SQRT(mu/rho)
            !
            ! Jacobian in x-direction (rotated!)
            !
            A(1,7) = -lambda-2*mu
            A(2,7) = -lambda
            A(3,7) = -lambda
            A(4,8) = -mu
            A(6,9) = -mu
            A(7,1) = -1./rho
            A(8,4) = -1./rho
            A(9,6) = -1./rho
            !
            ! Eigenvalue Matrix
            !
            LA(1,1) = -cp
            LA(2,2) = -cs
            LA(3,3) = -cs
            LA(4,4) =  0.
            LA(5,5) =  0.
            LA(6,6) =  0.
            LA(7,7) =  cs
            LA(8,8) =  cs
            LA(9,9) =  cp
            !
            IF(mu.EQ.0.) THEN
                absA(:,:) = 0.
                absA(1,1) = w_speed(1)
                absA(2,1) = 1./sqrt(lambda)*lambda/sqrt(rho)
                absA(3,1) = 1./sqrt(lambda)*lambda/sqrt(rho)
                absA(7,7) = sqrt(lambda)/sqrt(rho)
                RETURN
            ENDIF
            !
            RA(1,1) =  lambda + 2.*mu
            RA(2,1) =  lambda
            RA(3,1) =  lambda
            RA(7,1) = +cp      
            RA(4,2) =  mu
            RA(8,2) =  cs
            RA(6,3) =  mu
            RA(9,3) =  cs
            RA(5,4) =  1.      
            RA(2,5) =  1. 
            RA(3,6) =  1. 
            RA(6,7) =  mu
            RA(9,7) = -cs
            RA(4,8) =  mu
            RA(8,8) = -cs
            RA(1,9) =  lambda + 2.*mu
            RA(2,9) =  lambda
            RA(3,9) =  lambda
            RA(7,9) = -cp
            !
            iRA(1,1) = 1./(lambda+2*mu)/2.      
            iRA(1,7) = 1./cp/2.
            iRA(2,4) = 1./mu/2.
            iRA(2,8) = 1./cs/2.
            iRA(3,6) = 1./mu/2.
            iRA(3,9) = 1./cs/2.
            iRA(4,5) = 1
            iRA(5,1) = -lambda/(lambda+2*mu)
            iRA(5,2) = 1.
            iRA(6,1) = -lambda/(lambda+2*mu)
            iRA(6,3) = 1.
            iRA(7,6) = 1./mu/2.
            iRA(7,9) = -1./cs/2.
            iRA(8,4) = 1./mu/2.
            iRA(8,8) = -1./cs/2.
            iRA(9,1) = 1./(lambda+2.*mu)/2.
            iRA(9,7) = -1./cp/2.

         ELSE   !Anisotropic fluxes
            SELECT CASE(PoroFlux)
            CASE(0)
            !
            ! Initialization
            !
            A(:,:)          = 0.
            T(:,:)          = 0.
            iT(:,:)         = 0.
            absA(:,:)       = 0.
            rho             = W0(1)
            !
            nx = n0(1)
            ny = n0(2)
            nz = n0(3)
            sx = s0(1)
            sy = s0(2)
            sz = s0(3)
            tx = t0(1)
            ty = t0(2)
            tz = t0(3)
            ! 
            ! Rotation matrices for the transformation of the Voigt matrix
            ! C' = T*C*T^T
            !
            ! Transformation matrix TO the rotated system              
            T(1,1) = nx**2
            T(1,2) = ny**2
            T(1,3) = nz**2
            T(1,4) = 2*nz*ny
            T(1,5) = 2*nz*nx
            T(1,6) = 2*ny*nx
            T(2,1) = sx**2
            T(2,2) = sy**2
            T(2,3) = sz**2
            T(2,4) = 2*sz*sy
            T(2,5) = 2*sz*sx
            T(2,6) = 2*sy*sx
            T(3,1) = tx**2
            T(3,2) = ty**2
            T(3,3) = tz**2
            T(3,4) = 2*tz*ty
            T(3,5) = 2*tz*tx
            T(3,6) = 2*ty*tx
            T(4,1) = sx*tx
            T(4,2) = sy*ty
            T(4,3) = sz*tz
            T(4,4) = sz*ty+sy*tz
            T(4,5) = sz*tx+sx*tz
            T(4,6) = sy*tx+sx*ty
            T(5,1) = nx*tx
            T(5,2) = ny*ty
            T(5,3) = nz*tz
            T(5,4) = nz*ty+ny*tz
            T(5,5) = nz*tx+nx*tz
            T(5,6) = ny*tx+nx*ty
            T(6,1) = nx*sx
            T(6,2) = ny*sy
            T(6,3) = nz*sz
            T(6,4) = nz*sy+ny*sz
            T(6,5) = nz*sx+nx*sz
            T(6,6) = ny*sx+nx*sy
            !
            ! Transpose of Transformation matrix
            !
            DO i = 1, 6 
              DO j = 1, 6
                TT(i,j) = T(j,i)
              ENDDO
            ENDDO
            !           
            ! Material constants in the global xyz system
            !            
            rho    = W0( 1)
            c(1,1) = W0( 2)
            c(1,2) = W0( 3)
            c(1,3) = W0( 4)
            c(1,4) = W0( 5)
            c(1,5) = W0( 6)
            c(1,6) = W0( 7)
            c(2,2) = W0( 8)
            c(2,3) = W0( 9)
            c(2,4) = W0(10)
            c(2,5) = W0(11)
            c(2,6) = W0(12)
            c(3,3) = W0(13)
            c(3,4) = W0(14)
            c(3,5) = W0(15)
            c(3,6) = W0(16)
            c(4,4) = W0(17)
            c(4,5) = W0(18)
            c(4,6) = W0(19)
            c(5,5) = W0(20)
            c(5,6) = W0(21)
            c(6,6) = W0(22)
            !
            ! Initialize lower half of Voigt matrix using the symmetry property
            !
            DO i = 1, 6 
              DO j = 1, i-1
                c(i,j) = c(j,i)
              ENDDO
            ENDDO
            !
            ! Initialize and rotate Voigt matrix to get local material properties
            !
            Voigt_rot(:,:) = MATMUL( T(:,:), MATMUL(c(:,:),TT(:,:)) )
            !
            c(:,:) = Voigt_rot(:,:)
            ! 
            ! Jacobian in normal direction
            !
            A(1,7:9) = (/ -c(1,1), -c(1,6), -c(1,5) /)
            A(2,7:9) = (/ -c(1,2), -c(2,6), -c(2,5) /)
            A(3,7:9) = (/ -c(1,3), -c(3,6), -c(3,5) /)
            A(4,7:9) = (/ -c(1,6), -c(6,6), -c(5,6) /)
            A(5,7:9) = (/ -c(1,4), -c(4,6), -c(4,5) /)
            A(6,7:9) = (/ -c(1,5), -c(5,6), -c(5,5) /)            
            A(7,1)   = -1./rho
            A(8,4)   = -1./rho
            A(9,6)   = -1./rho
            !
            rus=0  ! TOGGLE: rus=0 (GODUNOV) OR rus=1 (RUSANOV)
            !
            IF(rus.EQ.0) THEN
             ! Godunov's flux anisotropic
             !
             XX(1,:) = (c(1,1)-w_speed(1:3)**2*rho)
             XX(2,:) = (c(6,6)-w_speed(1:3)**2*rho)
             XX(3,:) = (c(5,5)-w_speed(1:3)**2*rho)
             !
             ! To avoid zero vectors (G,H,I) we shuffle through 3 different G-J solutions
             TOL = (c(1,1)**2) / 1.0e10 !(maximum tolerance)
             !
             J=1
             !
             DO WHILE (J<=3) !Loops through 3 different results obtained by Gauss-Jordan Elimination
               GG(J) = XX(2,J)*XX(3,J)-c(5,6)**2
               HH(J) = c(5,6)*c(1,5)-c(1,6)*XX(3,J)
               II(J) = c(1,6)*c(5,6)-c(1,5)*XX(2,J)

               IF((ABS(GG(J))+ABS(HH(J))+ABS(II(J))).LE.TOL) THEN
                 GG(J) = c(1,5)*c(5,6)-XX(3,J)*c(1,6)
                 HH(J) = XX(3,J)*XX(1,J)-c(1,5)**2
                 II(J) = c(1,6)*c(1,5)-XX(1,J)*c(5,6)
               ENDIF
               IF((ABS(GG(J))+ABS(HH(J))+ABS(II(J))).LE.TOL) THEN
                 GG(J) = c(5,6)*c(1,6)-XX(2,J)*c(1,5)
                 HH(J) = c(1,5)*c(1,6)-XX(1,J)*c(5,6)
                 II(J) = XX(1,J)*XX(2,J)-c(1,6)**2
               ENDIF
               IF((ABS(GG(J))+ABS(HH(J))+ABS(II(J))).LE.TOL) THEN !In this case we assume an isotropic-like solution
                  GG(:)=0.d0
                  HH(:)=0.d0
                  II(:)=0.d0
                  GG(1)=1.d0 
                  HH(2)=1.d0 
                  II(3)=1.d0 
                  J=4  !Terminate DO WHILE
               ENDIF
               J=J+1
             ENDDO
             !
             ! Normalization
             DO J=1,3
               NORMA = SQRT(GG(J)**2+HH(J)**2+II(J)**2)
               GG(J) = GG(J) / NORMA
               HH(J) = HH(J) / NORMA
               II(J) = II(J) / NORMA
             ENDDO
             !
             ! Looking for linear-independence of the solutions
             J=1
             TOL=1.0e-6 
             DO WHILE (J<=3)
               J2=MOD(J,3)+1
               NORMA=GG(J)*GG(MOD(J,3)+1)+HH(J)*HH(MOD(J,3)+1)+II(J)*II(MOD(J,3)+1)
               IF(ABS(NORMA).GE.TOL) THEN
                 J3=MOD(J2,3)+1
                 GG(J) = HH(J2)*II(J3)-HH(J3)*II(J2)
                 HH(J) = II(J2)*GG(J3)-II(J3)*GG(J2)
                 II(J) = GG(J2)*HH(J3)-GG(J3)*HH(J2)
               ENDIF
               J=J+1
             ENDDO
             !
             GG(:) = GG(:) / SQRT(rho)
             HH(:) = HH(:) / SQRT(rho)
             II(:) = II(:) / SQRT(rho)
             !
             DO I=1,3
               E(7) = GG(I)
               E(8) = HH(I)
               E(9) = II(I)
               E(1) = (c(1,1)*E(7) + c(1,6)*E(8) + c(1,5)*E(9)) / w_speed(I)
               E(2) = (c(1,2)*E(7) + c(2,6)*E(8) + c(2,5)*E(9)) / w_speed(I)
               E(3) = (c(1,3)*E(7) + c(3,6)*E(8) + c(3,5)*E(9)) / w_speed(I)
               E(4) = (c(1,6)*E(7) + c(6,6)*E(8) + c(5,6)*E(9)) / w_speed(I)
               E(5) = (c(1,4)*E(7) + c(4,6)*E(8) + c(4,5)*E(9)) / w_speed(I)
               E(6) = (c(1,5)*E(7) + c(5,6)*E(8) + c(5,5)*E(9)) / w_speed(I)
               !
               RA(:,I)=E(:)
               !
               absA(1,1) = absA(1,1) + E(1)*E(7)
               absA(2,1) = absA(2,1) + E(2)*E(7)
               absA(3,1) = absA(3,1) + E(3)*E(7)
               absA(4,1) = absA(4,1) + E(4)*E(7)
               absA(5,1) = absA(5,1) + E(5)*E(7)
               absA(6,1) = absA(6,1) + E(6)*E(7)
               absA(2,4) = absA(2,4) + E(2)*E(8)
               absA(3,4) = absA(3,4) + E(3)*E(8)
               absA(4,4) = absA(4,4) + E(4)*E(8)
               absA(5,4) = absA(5,4) + E(5)*E(8)
               absA(6,4) = absA(6,4) + E(6)*E(8)
               absA(2,6) = absA(2,6) + E(2)*E(9)
               absA(3,6) = absA(3,6) + E(3)*E(9)
               absA(5,6) = absA(5,6) + E(5)*E(9)
               absA(6,6) = absA(6,6) + E(6)*E(9)
               !
             ENDDO
             !
             ! Complete the right eigenvector matrix by symmetry
             DO I=1,3
               RA(1:6,10-I) =-RA(1:6,I)
               RA(7:9,10-I) = RA(7:9,I)
             ENDDO
             RA(2,4) = 1.0
             RA(3,5) = 1.0
             RA(6,6) = 1.0
             !
             absA(1,4) = absA(4,1)
             absA(1,6) = absA(6,1)
             absA(4,6) = absA(6,4)
             absA(7,7) = absA(1,1)
             absA(7,8) = absA(4,1)
             absA(7,9) = absA(6,1)
             absA(8,7) = absA(4,1)
             absA(8,8) = absA(4,4)
             absA(8,9) = absA(4,6)
             absA(9,7) = absA(6,1)
             absA(9,8) = absA(4,6)
             absA(9,9) = absA(6,6)
             !
            ELSEIF (rus.eq.1) THEN
             ! Matrix of numerical viscosity for the Rusanov flux
             ! 
             absa(:,:) = 0. 
             amax = w_speed(1)
             Unity(:,:) = 0.
             DO i = 1, EQN%nVar
                 Unity(i,i) = 1.
             ENDDO
             absA(:,:) = amax*Unity(:,:)
             !
            ENDIF !finishes "rus" if
            !
            AniVec(1,:)=GG(:)
            AniVec(2,:)=HH(:)
            AniVec(3,:)=II(:)
            !
            RETURN
            ! 
         CASE(1)  !Porous media for the moment just Rusanov fluxes
            !
            ! Initialization
            !
            A(:,:)          = 0.
            T(:,:)          = 0.
            iT(:,:)         = 0.
            absA(:,:)       = 0.
            rho             = W0(1)
            !
            nx = n0(1)
            ny = n0(2)
            nz = n0(3)
            sx = s0(1)
            sy = s0(2)
            sz = s0(3)
            tx = t0(1)
            ty = t0(2)
            tz = t0(3)
            ! 
            ! Rotation matrices for the transformation of the Voigt matrix
            ! C' = T*C*T^T
            !
            ! Transformation matrix TO the rotated system              
            T(1,1) = nx**2
            T(1,2) = ny**2
            T(1,3) = nz**2
            T(1,4) = 2*nz*ny
            T(1,5) = 2*nz*nx
            T(1,6) = 2*ny*nx
            T(2,1) = sx**2
            T(2,2) = sy**2
            T(2,3) = sz**2
            T(2,4) = 2*sz*sy
            T(2,5) = 2*sz*sx
            T(2,6) = 2*sy*sx
            T(3,1) = tx**2
            T(3,2) = ty**2
            T(3,3) = tz**2
            T(3,4) = 2*tz*ty
            T(3,5) = 2*tz*tx
            T(3,6) = 2*ty*tx
            T(4,1) = sx*tx
            T(4,2) = sy*ty
            T(4,3) = sz*tz
            T(4,4) = sz*ty+sy*tz
            T(4,5) = sz*tx+sx*tz
            T(4,6) = sy*tx+sx*ty
            T(5,1) = nx*tx
            T(5,2) = ny*ty
            T(5,3) = nz*tz
            T(5,4) = nz*ty+ny*tz
            T(5,5) = nz*tx+nx*tz
            T(5,6) = ny*tx+nx*ty
            T(6,1) = nx*sx
            T(6,2) = ny*sy
            T(6,3) = nz*sz
            T(6,4) = nz*sy+ny*sz
            T(6,5) = nz*sx+nx*sz
            T(6,6) = ny*sx+nx*sy
            !
            ! Transpose of Transformation matrix
            !
            DO i = 1, 6 
              DO j = 1, 6
                TT(i,j) = T(j,i)
              ENDDO
            ENDDO
             !           
             ! Material constants in the global xyz system
             !            
             rho_S  = W0( 1)
             c(1,1) = W0( 2)
             c(1,2) = W0( 3)
             c(1,3) = W0( 4)
             c(1,4) = W0( 5)
             c(1,5) = W0( 6)
             c(1,6) = W0( 7)
             c(2,2) = W0( 8)
             c(2,3) = W0( 9)
             c(2,4) = W0(10)
             c(2,5) = W0(11)
             c(2,6) = W0(12)
             c(3,3) = W0(13)
             c(3,4) = W0(14)
             c(3,5) = W0(15)
             c(3,6) = W0(16)
             c(4,4) = W0(17)
             c(4,5) = W0(18)
             c(4,6) = W0(19)
             c(5,5) = W0(20)
             c(5,6) = W0(21)
             c(6,6) = W0(22)
             !
             ! Initialize lower half of Voigt matrix using the symmetry property
             !
             DO i = 1, 6 
               DO j = 1, i-1
                 c(i,j) = c(j,i)
               ENDDO
             ENDDO
             !
             ! Initialize and rotate Voigt matrix to get local material properties
             !
             Voigt_rot(:,:) = MATMUL( T(:,:), MATMUL(c(:,:),TT(:,:)) )
             !
             c(:,:) = Voigt_rot(:,:)
             ! 
             rho_F    = W0(23)
             K_F      = W0(24)
             nu       = W0(25)
             K_S      = W0(26)
             Poro     = W0(27)

              ! Permeabilities and tortuosities are rotated assuming that they represent an ellipsoid in the reference system!
             ! The ellipsoid is aligned with the mesh reference system (no rotation in .par file)
             Kappa(1) = SQRT((W0(28)*nx)**2+(W0(29)*ny)**2+(W0(30)*nz)**2)
             Kappa(2) = SQRT((W0(28)*sx)**2+(W0(29)*sy)**2+(W0(30)*sz)**2)
             Kappa(3) = SQRT((W0(28)*tx)**2+(W0(29)*ty)**2+(W0(30)*tz)**2)
             Tor(1)   = SQRT((W0(31)*nx)**2+(W0(32)*ny)**2+(W0(33)*nz)**2)
             Tor(2)   = SQRT((W0(31)*sx)**2+(W0(32)*sy)**2+(W0(33)*sz)**2)
             Tor(3)   = SQRT((W0(31)*tx)**2+(W0(32)*ty)**2+(W0(33)*tz)**2)

             rho      = rho_S * (1 - Poro) + Poro * rho_F
             K_Mean   = 1./9.*(c(1,1)+c(2,2)+c(3,3)+2*(c(1,2)+c(1,3)+c(2,3)))
             MM       = K_S / ((1 - K_Mean/K_S) - Poro*(1 - K_S/K_F))
             Alpha(1) = 1 - (c(1,1)+c(1,2)+c(1,3)) / (3.*K_S)
             Alpha(2) = 1 - (c(1,2)+c(2,2)+c(2,3)) / (3.*K_S)
             Alpha(3) = 1 - (c(1,3)+c(2,3)+c(3,3)) / (3.*K_S)
             Alpha(4) = - (c(1,4)+c(2,4)+c(3,4)) / (3.*K_S)
             Alpha(5) = - (c(1,5)+c(2,5)+c(3,5)) / (3.*K_S)
             Alpha(6) = - (c(1,6)+c(2,6)+c(3,6)) / (3.*K_S)
             Rho1(:)  = rho - (rho_F**2 / (rho_F * Tor(:) / Poro))
             Rho2(:)  = rho_F - (rho_F * Tor(:) / Poro) * rho /rho_F 
             Beta1(:) = rho_F / (rho_F * Tor(:) / Poro)
             Beta2(:)  = rho / rho_F
             !
             ! Computation of undrained ce(i,j) coeffs
             DO i=1,6
               DO j=1,6
                 c(i,j) = c(i,j) + MM * Alpha(i)*Alpha(j)
               ENDDO
             ENDDO
             ! 
             ! Jacobian in x-direction
             !
             A(1,7:9) = (/ -c(1,1), -c(1,6), -c(1,5) /)
             A(2,7:9) = (/ -c(1,2), -c(2,6), -c(2,5) /)
             A(3,7:9) = (/ -c(1,3), -c(3,6), -c(3,5) /)
             A(4,7:9) = (/ -c(1,6), -c(6,6), -c(5,6) /)
             A(5,7:9) = (/ -c(1,4), -c(4,6), -c(4,5) /)
             A(6,7:9) = (/ -c(1,5), -c(5,6), -c(5,5) /)
             A(7,1) = -1./rho1(1)
             A(8,4) = -1./rho1(1)
             A(9,6) = -1./rho1(1)
             A(10,7:9) = (/ Alpha(1)*MM, Alpha(6)*MM, Alpha(5)*MM /)
             A(11,1) = -1./rho2(1)
             A(12,4) = -1./rho2(1)
             A(13,6) = -1./rho2(1)
             A(7,10) = -beta1(1)/rho1(1)
             A(11,10)= -beta2(1)/rho2(1)
             A(1,11) = -Alpha(1)*MM
             A(2,11) = -Alpha(2)*MM
             A(3,11) = -Alpha(3)*MM
             A(4,11) = -Alpha(6)*MM
             A(5,11) = -Alpha(4)*MM
             A(6,11) = -Alpha(5)*MM
             A(10,11)= MM
             !             
             ! Matrix of numerical viscosity for the Rusanov flux
             ! 
             absa(:,:) = 0. 
             !
             rus=1
             !
             IF(rus.EQ.1) THEN
               amax = w_speed(1)
               Unity(:,:) = 0.
               DO i = 1, EQN%nVar
                   Unity(i,i) = 1.
               ENDDO
               absA(:,:) = amax*Unity(:,:)
               !
             ELSEIF (rus.eq.0) THEN !Poroelastic Godunov fluxes (CAUTION: only work for the isotropic case!)
               !
               IF((c(6,6)/c(1,1)).GE.1e-10) THEN !Poroelastic case
                 !
                 X(1) = c(1,1)-rho1(1)*(w_speed(1))**2*(1.d0+Beta1(1)/(Beta2(1)-Beta1(1)))
                 X(2) = c(1,1)-rho1(1)*(w_speed(4))**2*(1.d0+Beta1(1)/(Beta2(1)-Beta1(1)))
                 Z(1) = Alpha(1)*MM-rho1(1)*(w_speed(1))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 Z(2) = Alpha(1)*MM-rho1(1)*(w_speed(4))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 Y(1) = MM+rho2(1)*(w_speed(1))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 Y(2) = MM+rho2(1)*(w_speed(4))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 !
                 r7(1)  = 1.d0
                 r7(2)  = 1.d0
                 r11(:) = -X(:)/Z(:)
                 r10(1) = -(Alpha(1)*MM*r7(1)+MM*r11(1))/w_speed(1)
                 r10(2) = -(Alpha(1)*MM*r7(2)+MM*r11(2))/w_speed(4)
                 r1(1)  = (c(1,1)*r7(1)+Alpha(1)*MM*r11(1))/w_speed(1)
                 r1(2)  = (c(1,1)*r7(2)+Alpha(1)*MM*r11(2))/w_speed(4)
                 r2(1)  = r1(1)-2*c(6,6)*r7(1)/w_speed(1)
                 r2(2)  = r1(2)-2*c(6,6)*r7(2)/w_speed(4)
                 !
                 NORMA1 = (w_speed(1))/(2.d0*r1(1)-2.d0*r10(1)*r11(1))
                 NORMA2 = (w_speed(4))/(2.d0*r1(2)-2.d0*r10(2)*r11(2))
                 !
                 absa(1,1)   =  2*r1(1)*r7(1)*NORMA1    + 2*r1(2)*r7(2)*NORMA2
                 absa(2,1)   =  2*r2(1)*r7(1)*NORMA1    + 2*r2(2)*r7(2)*NORMA2
                 absa(3,1)   =  2*r2(1)*r7(1)*NORMA1    + 2*r2(2)*r7(2)*NORMA2
                 absa(10,1)  =  2*r10(1)*r7(1)*NORMA1   + 2*r10(2)*r7(2)*NORMA2
                 absa(1,10)  = -2*r1(1)*r11(1)*NORMA1   - 2*r1(2)*r11(2)*NORMA2
                 absa(10,10) = -2*r10(1)*r11(1)*NORMA1  - 2*r10(2)*r11(2)*NORMA2
                 absa(7,7)   =  absa(1,1)
                 absa(11,11) =  absa(10,10)
                 absa(11,7)  = -absa(1,10)
                 absa(7,11)  = -absa(10,1)
                 !
                 absa(4,4)   =  w_speed(2)
                 absa(6,6)   =  w_speed(2)
                 absa(8,8)   =  w_speed(2)
                 absa(9,9)   =  w_speed(2)
                 !
               ELSEIF((c(6,6)/c(1,1)).LT.1e-10) THEN !Poroacoustic case
                 !
                 X(1) = c(1,1)-rho1(1)*(w_speed(1))**2*(1.d0+Beta1(1)/(Beta2(1)-Beta1(1)))
                 X(2) = c(1,1)-rho1(1)*(w_speed(2))**2*(1.d0+Beta1(1)/(Beta2(1)-Beta1(1)))
                 Z(1) = Alpha(1)*MM-rho1(1)*(w_speed(1))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 Z(2) = Alpha(1)*MM-rho1(1)*(w_speed(2))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 Y(1) = MM+rho2(1)*(w_speed(1))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 Y(2) = MM+rho2(1)*(w_speed(2))**2*(1.d0/(Beta2(1)-Beta1(1)))
                 !
                 r7(1)  = 1.d0
                 r7(2)  = 1.d0
                 r11(:) = -X(:)/Z(:)
                 r10(1) = -(Alpha(1)*MM*r7(1)+MM*r11(1))/w_speed(1)
                 r10(2) = -(Alpha(1)*MM*r7(2)+MM*r11(2))/w_speed(2)
                 r1(1)  = (c(1,1)*r7(1)+Alpha(1)*MM*r11(1))/w_speed(1)
                 r1(2)  = (c(1,1)*r7(2)+Alpha(1)*MM*r11(2))/w_speed(2)
                 r2(1)  = r1(1)
                 r2(2)  = r1(2)
                 !
                 NORMA1 = (w_speed(1))/(2.d0*r1(1)-2.d0*r10(1)*r11(1))
                 NORMA2 = (w_speed(2))/(2.d0*r1(2)-2.d0*r10(2)*r11(2))
                 !
                 absa(1,1)   =  2*r1(1)*r7(1)*NORMA1    + 2*r1(2)*r7(2)*NORMA2
                 absa(2,1)   =  2*r2(1)*r7(1)*NORMA1    + 2*r2(2)*r7(2)*NORMA2
                 absa(3,1)   =  2*r2(1)*r7(1)*NORMA1    + 2*r2(2)*r7(2)*NORMA2
                 absa(10,1)  =  2*r10(1)*r7(1)*NORMA1   + 2*r10(2)*r7(2)*NORMA2
                 absa(1,10)  = -2*r1(1)*r11(1)*NORMA1   - 2*r1(2)*r11(2)*NORMA2
                 absa(10,10) = -2*r10(1)*r11(1)*NORMA1  - 2*r10(2)*r11(2)*NORMA2
                 absa(7,7)   =  absa(1,1)
                 absa(11,11) =  absa(10,10)
                 absa(11,7)  = -absa(1,10)
                 absa(7,11)  = -absa(10,1)
                 !
               ENDIF
               !
             ENDIF

             RETURN
             !
             END SELECT
             !
          ENDIF

    !TEST(:,:) = MATMUL(iRA,RA)
    !Unity(:,:) = 0.
    !DO i = 1, EQN%nVar
    !  Unity(i,i) = 1.
    !ENDDO
    !checksum  = SUM((TEST(:,:)-Unity(:,:))**2)
    !IF(checksum.GT.1e-6) THEN
    !  print *, 'fatal error: iR*R .NE. 1 but ', checksum, W0(:)
    !  stop
    !ENDIF
    !TEST(:,:) = MATMUL(MATMUL(iRA,A),RA)
    !checksum  = SUM((TEST(:,:)-LA(:,:))**2)
    !IF(checksum.GT.1e-6) THEN
    !  print *, 'fatal error: iR*A*R .NE. Lambda'
    !  stop
    !ENDIF

    !absA  = MATMUL(MATMUL(RA,ABS(LA)),iRA)

    absA(:,:) = 0.
    Test(:,:) = RA(:,:)
    DO i = 1, EQN%nVar
      Test(:,i) = Test(:,i)*ABS(LA(i,i))
    ENDDO

    DO j = 1, EQN%nVar
      DO i = 1, EQN%nVar
         absA(:,j) = absA(:,j) + Test(:,i)*iRA(i,j)
      ENDDO
    ENDDO
    !
    AniVec(1,:)=GG(:)
    AniVec(2,:)=HH(:)
    AniVec(3,:)=II(:)


  END SUBROUTINE JacobiNormal3D

  !-------------------------------------------------------------------------!


  ! Compute the rotation matrix on 3D
  SUBROUTINE RotationMatrix3D(n1,n2,n3,T,iT,EQN)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations):: EQN                                                    ! 
    REAL            :: n1(3),n2(3),n3(3)                                      ! Component of normal vector
    REAL            :: T(EQN%nVar,EQN%nVar), iT(EQN%nVar,EQN%nVar)            ! Rotation matrix
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: n1, n2, n3, EQN                                        !
    INTENT(OUT)     :: T, iT                                                  !
    !-------------------------------------------------------------------------!
    REAL            :: nx,ny,nz
    REAL            :: sx,sy,sz
    REAL            :: tx,ty,tz
    INTEGER         :: i,j

    nx=n1(1); ny=n1(2); nz=n1(3);
    sx=n2(1); sy=n2(2); sz=n2(3);
    tx=n3(1); ty=n3(2); tz=n3(3);

    T(:,:)  = 0.
    iT(:,:) = 0.             
    ! u' = iT * u 
    ! u   = T * u'
    ! Inverse Transformation matrix
    T(1,1) = nx * nx
    T(1,2) = sx * sx
    T(1,3) = tx * tx
    T(1,4) = 2 * nx * sx
    T(1,5) = 2 * sx * tx
    T(1,6) = 2 * nx * tx
    T(2,1) = ny * ny
    T(2,2) = sy * sy
    T(2,3) = ty * ty
    T(2,4) = 2 * ny * sy
    T(2,5) = 2 * sy * ty
    T(2,6) = 2 * ny * ty
    T(3,1) = nz * nz
    T(3,2) = sz * sz
    T(3,3) = tz * tz
    T(3,4) = 2 * nz * sz
    T(3,5) = 2 * sz * tz
    T(3,6) = 2 * nz * tz
    T(4,1) = ny * nx
    T(4,2) = sy * sx
    T(4,3) = ty * tx
    T(4,4) = ny * sx + nx * sy
    T(4,5) = sy * tx + sx * ty
    T(4,6) = ny * tx + nx * ty
    T(5,1) = nz * ny
    T(5,2) = sz * sy
    T(5,3) = tz * ty
    T(5,4) = nz * sy + ny * sz
    T(5,5) = sz * ty + sy * tz
    T(5,6) = nz * ty + ny * tz
    T(6,1) = nz * nx
    T(6,2) = sz * sx
    T(6,3) = tz * tx
    T(6,4) = nz * sx + nx * sz
    T(6,5) = sz * tx + sx * tz
    T(6,6) = nz * tx + nx * tz
    T(7,7) = nx
    T(7,8) = sx
    T(7,9) = tx
    T(8,7) = ny
    T(8,8) = sy
    T(8,9) = ty
    T(9,7) = nz
    T(9,8) = sz
    T(9,9) = tz
       !
    ! Transformation matrix
    iT(1,1) = nx * nx
    iT(1,2) = ny * ny
    iT(1,3) = nz * nz
    iT(1,4) = 2 * ny * nx
    iT(1,5) = 2 * nz * ny
    iT(1,6) = 2 * nz * nx
    iT(2,1) = sx * sx
    iT(2,2) = sy * sy
    iT(2,3) = sz * sz
    iT(2,4) = 2 * sy * sx
    iT(2,5) = 2 * sz * sy
    iT(2,6) = 2 * sz * sx
    iT(3,1) = tx * tx
    iT(3,2) = ty * ty
    iT(3,3) = tz * tz
    iT(3,4) = 2 * ty * tx
    iT(3,5) = 2 * tz * ty
    iT(3,6) = 2 * tz * tx
    iT(4,1) = nx * sx
    iT(4,2) = ny * sy
    iT(4,3) = nz * sz
    iT(4,4) = ny * sx + nx * sy
    iT(4,5) = nz * sy + ny * sz
    iT(4,6) = nz * sx + nx * sz
    iT(5,1) = sx * tx
    iT(5,2) = sy * ty
    iT(5,3) = sz * tz
    iT(5,4) = sy * tx + sx * ty
    iT(5,5) = sz * ty + sy * tz
    iT(5,6) = sz * tx + sx * tz
    iT(6,1) = nx * tx
    iT(6,2) = ny * ty
    iT(6,3) = nz * tz
    iT(6,4) = ny * tx + nx * ty
    iT(6,5) = nz * ty + ny * tz
    iT(6,6) = nz * tx + nx * tz
    iT(7,7) = nx
    iT(7,8) = ny
    iT(7,9) = nz
    iT(8,7) = sx
    iT(8,8) = sy
    iT(8,9) = sz
    iT(9,7) = tx
    iT(9,8) = ty
    iT(9,9) = tz

    IF(EQN%Poroelasticity.NE.0) THEN
      T(10,10)  = 1.d0
      iT(10,10) = 1.d0
      DO i=1,3
        DO j=1,3
          T(10+i,10+j)  = T(6+i,6+j)
          iT(10+i,10+j) = iT(6+i,6+j)
        ENDDO
      ENDDO
    ENDIF


  END SUBROUTINE RotationMatrix3D


  ! Compute the elastic Jacobians in 3D
  SUBROUTINE BuildJacobians3D_Ela(Mat,A,B,C,EQN)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations):: EQN                                                    !      
    REAL            :: Mat(EQN%nBackgroundVar) 
    REAL            :: A(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: B(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: C(EQN%nVarTotal,EQN%nVarTotal)                         ! Jacobian matrices
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: Mat,EQN                                                !
    INTENT(OUT)     :: A, B, C                                                !
    !-------------------------------------------------------------------------!
    REAL            :: rhoinv, mu, lambda
    !-------------------------------------------------------------------------!
    rhoinv = 1.0/Mat(1)  ! rho = 1, mu = 2, lambda = 3
    mu     = Mat(2)
    lambda = Mat(3)
    !
    A(:,:) = 0.
    B(:,:) = 0.
    C(:,:) = 0.
    ! 
    ! Jacobian in x-direction
    !  
    A(1,7) = -lambda-2*mu
    A(2,7) = -lambda
    A(3,7) = -lambda
    A(4,8) = -mu
    A(6,9) = -mu
    A(7,1) = -rhoinv
    A(8,4) = -rhoinv
    A(9,6) = -rhoinv
    !
    ! Jacobian in y-direction
    !
    B(1,8) = -lambda
    B(2,8) = -lambda-2*mu
    B(3,8) = -lambda
    B(4,7) = -mu
    B(5,9) = -mu
    B(7,4) = -rhoinv
    B(8,2) = -rhoinv
    B(9,5) = -rhoinv
    !
    ! Jacobian in z-direction
    !
    C(1,9) = -lambda      
    C(2,9) = -lambda
    C(3,9) = -lambda-2*mu
    C(5,8) = -mu
    C(6,7) = -mu
    C(7,6) = -rhoinv
    C(8,5) = -rhoinv
    C(9,3) = -rhoinv


  END SUBROUTINE BuildJacobians3D_Ela


  ! Compute the anisotropic Jacobians in 3D
  SUBROUTINE BuildJacobians3D_Ani(Mat,A,B,C,EQN)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations):: EQN                                                    ! 
    REAL            :: Mat(EQN%nBackgroundVar) 
    REAL            :: A(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: B(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: C(EQN%nVarTotal,EQN%nVarTotal)                         ! Jacobian matrices
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: Mat,EQN                                                !
    INTENT(OUT)     :: A, B, C                                                !
    !-------------------------------------------------------------------------!
    REAL            :: rhoinv, ce(6,6)
    !-------------------------------------------------------------------------!
    !
    A(:,:) = 0.
    B(:,:) = 0.
    C(:,:) = 0.
    !
    rhoinv  = 1.0/Mat( 1)
    ce(1,1) = Mat( 2)
    ce(1,2) = Mat( 3)
    ce(1,3) = Mat( 4)
    ce(1,4) = Mat( 5)
    ce(1,5) = Mat( 6)
    ce(1,6) = Mat( 7)
    ce(2,2) = Mat( 8)
    ce(2,3) = Mat( 9)
    ce(2,4) = Mat(10)
    ce(2,5) = Mat(11)
    ce(2,6) = Mat(12)
    ce(3,3) = Mat(13)
    ce(3,4) = Mat(14)
    ce(3,5) = Mat(15)
    ce(3,6) = Mat(16)
    ce(4,4) = Mat(17)
    ce(4,5) = Mat(18)
    ce(4,6) = Mat(19)
    ce(5,5) = Mat(20)
    ce(5,6) = Mat(21)
    ce(6,6) = Mat(22)
    ! 
    ! Jacobian in x-direction
    !
    A(1,7:9) = (/ -ce(1,1), -ce(1,6), -ce(1,5) /)
    A(2,7:9) = (/ -ce(1,2), -ce(2,6), -ce(2,5) /)
    A(3,7:9) = (/ -ce(1,3), -ce(3,6), -ce(3,5) /)
    A(4,7:9) = (/ -ce(1,6), -ce(6,6), -ce(5,6) /)
    A(5,7:9) = (/ -ce(1,4), -ce(4,6), -ce(4,5) /)
    A(6,7:9) = (/ -ce(1,5), -ce(5,6), -ce(5,5) /)
    A(7,1) = -rhoinv
    A(8,4) = -rhoinv
    A(9,6) = -rhoinv
    !
    ! Jacobian in y-direction
    !
    B(1,7:9) = (/ -ce(1,6), -ce(1,2), -ce(1,4) /)
    B(2,7:9) = (/ -ce(2,6), -ce(2,2), -ce(2,4) /)
    B(3,7:9) = (/ -ce(3,6), -ce(2,3), -ce(3,4) /)
    B(4,7:9) = (/ -ce(6,6), -ce(2,6), -ce(4,6) /)
    B(5,7:9) = (/ -ce(4,6), -ce(2,4), -ce(4,4) /)
    B(6,7:9) = (/ -ce(5,6), -ce(2,5), -ce(4,5) /)
    B(7,4) = -rhoinv
    B(8,2) = -rhoinv
    B(9,5) = -rhoinv
    !
    ! Jacobian in z-direction
    !
    C(1,7:9) = (/ -ce(1,5), -ce(1,4), -ce(1,3) /)
    C(2,7:9) = (/ -ce(2,5), -ce(2,4), -ce(2,3) /)
    C(3,7:9) = (/ -ce(3,5), -ce(3,4), -ce(3,3) /)
    C(4,7:9) = (/ -ce(5,6), -ce(4,6), -ce(3,6) /)
    C(5,7:9) = (/ -ce(4,5), -ce(4,4), -ce(3,4) /)
    C(6,7:9) = (/ -ce(5,5), -ce(4,5), -ce(3,5) /)
    C(7,6) = -rhoinv
    C(8,5) = -rhoinv
    C(9,3) = -rhoinv


  END SUBROUTINE BuildJacobians3D_Ani

  ! Compute the poroelastic Jacobians in 3D
  SUBROUTINE BuildJacobians3D_Poro(Mat,A,B,C,EC,EQN)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations):: EQN                                                    ! 
    REAL            :: Mat(EQN%nBackgroundVar) 
    REAL            :: A(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: B(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: C(EQN%nVarTotal,EQN%nVarTotal)                         ! Jacobian matrices
    REAL            :: EC(2,3) 
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: Mat,EQN                                                !
    INTENT(OUT)     :: A, B, C, EC                                            !
    !-------------------------------------------------------------------------!
    REAL            :: rho_S, ce(6,6), rho_F, K_F, K_S, rho, Poro, Kappa(3)
    REAL            :: K_Mean, MM, Alpha(6), Rho1(3), Rho2(3), Beta1(3), nu
    REAL            :: Tor(3), Beta2(3), T(EQN%nVarTotal,EQN%nVarTotal), iT(EQN%nVarTotal,EQN%nVarTotal)
    INTEGER         :: i, j
    !-------------------------------------------------------------------------!
    !
    A(:,:) = 0.
    B(:,:) = 0.
    C(:,:) = 0.
    EC(:,:)= 0.
    !
    rho_S   = Mat( 1)
    ce(1,1) = Mat( 2)
    ce(1,2) = Mat( 3)
    ce(1,3) = Mat( 4)
    ce(1,4) = Mat( 5)
    ce(1,5) = Mat( 6)
    ce(1,6) = Mat( 7)
    ce(2,2) = Mat( 8)
    ce(2,3) = Mat( 9)
    ce(2,4) = Mat(10)
    ce(2,5) = Mat(11)
    ce(2,6) = Mat(12)
    ce(3,3) = Mat(13)
    ce(3,4) = Mat(14)
    ce(3,5) = Mat(15)
    ce(3,6) = Mat(16)
    ce(4,4) = Mat(17)
    ce(4,5) = Mat(18)
    ce(4,6) = Mat(19)
    ce(5,5) = Mat(20)
    ce(5,6) = Mat(21)
    ce(6,6) = Mat(22)
    rho_F   = Mat(23)
    K_F     = Mat(24)
    nu      = Mat(25)
    K_S     = Mat(26)
    Poro    = Mat(27)
    Kappa(1) = Mat(28) 
    Kappa(2) = Mat(29) 
    Kappa(3) = Mat(30) 
    Tor(1)  = Mat(31) 
    Tor(2)  = Mat(32) 
    Tor(3)  = Mat(33) 
    rho     = rho_S * (1 - Poro) + Poro * rho_F
    K_Mean  = 1./9.*(ce(1,1)+ce(2,2)+ce(3,3)+2*(ce(1,2)+ce(1,3)+ce(2,3)))
    MM      = K_S / ((1 - K_Mean/K_S) - Poro*(1 - K_S/K_F))
    Alpha(1) = 1 - (ce(1,1)+ce(1,2)+ce(1,3)) / (3.*K_S)

    Alpha(2) = 1 - (ce(1,2)+ce(2,2)+ce(2,3)) / (3.*K_S)
    Alpha(3) = 1 - (ce(1,3)+ce(2,3)+ce(3,3)) / (3.*K_S)
    Alpha(4) = - (ce(1,4)+ce(2,4)+ce(3,4)) / (3.*K_S)
    Alpha(5) = - (ce(1,5)+ce(2,5)+ce(3,5)) / (3.*K_S)
    Alpha(6) = - (ce(1,6)+ce(2,6)+ce(3,6)) / (3.*K_S)
    Rho1(:)  = rho - (rho_F**2 / (rho_F * Tor(:) / Poro))
    Rho2(:)  = rho_F - (rho_F * Tor(:) / Poro) * rho /rho_F 
    Beta1(:) = rho_F / (rho_F * Tor(:) / Poro)
    Beta2(:) = rho / rho_F
    !
    ! Computation of undrained ce(i,j) coeffs
    DO i=1,6
      DO j=1,6
         ce(i,j) = ce(i,j) + MM * Alpha(i)*Alpha(j)
      ENDDO
    ENDDO
    ! 
    ! Jacobian in x-direction
    !
    A(1,7:9) = (/ -ce(1,1), -ce(1,6), -ce(1,5) /)
    A(2,7:9) = (/ -ce(1,2), -ce(2,6), -ce(2,5) /)
    A(3,7:9) = (/ -ce(1,3), -ce(3,6), -ce(3,5) /)
    A(4,7:9) = (/ -ce(1,6), -ce(6,6), -ce(5,6) /)
    A(5,7:9) = (/ -ce(1,4), -ce(4,6), -ce(4,5) /)
    A(6,7:9) = (/ -ce(1,5), -ce(5,6), -ce(5,5) /)
    A(7,1) = -1./rho1(1)
    A(8,4) = -1./rho1(1)
    A(9,6) = -1./rho1(1)
    A(10,7:9) = (/ Alpha(1)*MM, Alpha(6)*MM, Alpha(5)*MM /)
    A(11,1) = -1./rho2(1)
    A(12,4) = -1./rho2(1)
    A(13,6) = -1./rho2(1)
    A(7,10) = -beta1(1)/rho1(1)
    A(11,10)= -beta2(1)/rho2(1)
    A(1,11) = -Alpha(1)*MM
    A(2,11) = -Alpha(2)*MM
    A(3,11) = -Alpha(3)*MM
    A(4,11) = -Alpha(6)*MM
    A(5,11) = -Alpha(4)*MM
    A(6,11) = -Alpha(5)*MM
    A(10,11)= MM
    !
    ! Jacobian in y-direction
    !
    B(1,7:9) = (/ -ce(1,6), -ce(1,2), -ce(1,4) /)
    B(2,7:9) = (/ -ce(2,6), -ce(2,2), -ce(2,4) /)
    B(3,7:9) = (/ -ce(3,6), -ce(2,3), -ce(3,4) /)
    B(4,7:9) = (/ -ce(6,6), -ce(2,6), -ce(4,6) /)
    B(5,7:9) = (/ -ce(4,6), -ce(2,4), -ce(4,4) /)
    B(6,7:9) = (/ -ce(5,6), -ce(2,5), -ce(4,5) /)
    B(7,4) = -1./rho1(2)
    B(8,2) = -1./rho1(2)
    B(9,5) = -1./rho1(2)
    B(10,7:9) = (/ Alpha(6)*MM, Alpha(2)*MM, Alpha(4)*MM /)
    B(11,4) = -1./rho2(2)
    B(12,2) = -1./rho2(2)
    B(13,5) = -1./rho2(2)
    B(8,10) = -beta1(2)/rho1(2)
    B(12,10)= -beta2(2)/rho2(2)
    B(1,12) = -Alpha(1)*MM
    B(2,12) = -Alpha(2)*MM
    B(3,12) = -Alpha(3)*MM
    B(4,12) = -Alpha(6)*MM
    B(5,12) = -Alpha(4)*MM
    B(6,12) = -Alpha(5)*MM
    B(10,12)= MM
    !
    ! Jacobian in z-direction
    !
    C(1,7:9) = (/ -ce(1,5), -ce(1,4), -ce(1,3) /)
    C(2,7:9) = (/ -ce(2,5), -ce(2,4), -ce(2,3) /)
    C(3,7:9) = (/ -ce(3,5), -ce(3,4), -ce(3,3) /)
    C(4,7:9) = (/ -ce(5,6), -ce(4,6), -ce(3,6) /)
    C(5,7:9) = (/ -ce(4,5), -ce(4,4), -ce(3,4) /)
    C(6,7:9) = (/ -ce(5,5), -ce(4,5), -ce(3,5) /)
    C(7,6) = -1./rho1(3)
    C(8,5) = -1./rho1(3)
    C(9,3) = -1./rho1(3)
    C(10,7:9) = (/ Alpha(5)*MM, Alpha(4)*MM, Alpha(3)*MM /)
    C(11,6) = -1./rho2(3)
    C(12,5) = -1./rho2(3)
    C(13,3) = -1./rho2(3)
    C(9,10) = -beta1(3)/rho1(3)
    C(13,10)= -beta2(3)/rho2(3)
    C(1,13) = -Alpha(1)*MM
    C(2,13) = -Alpha(2)*MM
    C(3,13) = -Alpha(3)*MM
    C(4,13) = -Alpha(6)*MM
    C(5,13) = -Alpha(4)*MM
    C(6,13) = -Alpha(5)*MM
    C(10,13)= MM
    !
    EC(1,:) = beta1(:)/rho1(:)*nu/Kappa(:) !Entries of the reaction matrix E
    EC(2,:) = beta2(:)/rho2(:)*nu/Kappa(:)

  END SUBROUTINE BuildJacobians3D_Poro


  ! Compute the anelastic auxiliary matrices on 3D
  SUBROUTINE AnelasticAuxMatrix3D_Iso(AuxMatrix1,n1,n2,n3,Mat,EQN)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations):: EQN                                                    ! 
    REAL            :: n1(3),n2(3),n3(3)                                      ! Component of normal vector
    REAL            :: AuxMatrix1(6,9)
    REAL            :: Mat(EQN%nBackgroundVar) 
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: Mat, n1, n2, n3,EQN                                    ! 
    INTENT(OUT)     :: AuxMatrix1                                             !
    !-------------------------------------------------------------------------!
    REAL            :: nx,ny,nz
    REAL            :: sx,sy,sz
    REAL            :: tx,ty,tz
    REAL            :: den1, den2, den3
    !-------------------------------------------------------------------------!

    nx=n1(1); ny=n1(2); nz=n1(3);
    sx=n2(1); sy=n2(2); sz=n2(3);
    tx=n3(1); ty=n3(2); tz=n3(3);

    den1=1.d0/(SQRT((Mat(3)+2*Mat(2))/Mat(1))*Mat(1))
    den2=1.d0/(SQRT(Mat(2)/Mat(1))*Mat(1))
    den3=1.d0/(2.0*SQRT(Mat(2)/Mat(1))*Mat(1))
    !
    auxMatrix1(:,:) = 0.
    auxMatrix1(1,1) = (nx*nx*nx*nx)   *den1 + (nx*nx*sx*sx)                 *den2 + (nx*nx*tx*tx)                 *den2
    auxMatrix1(1,2) = (nx*nx*ny*ny)   *den1 + (nx*ny*sx*sy)                 *den2 + (nx*ny*tx*ty)                 *den2
    auxMatrix1(1,3) = (nx*nx*nz*nz)   *den1 + (nx*nz*sz*sx)                 *den2 + (nx*nz*tx*tz)                 *den2
    auxMatrix1(1,4) = (2*nx*nx*nx*ny) *den1 + (nx*sx*(ny*sx+nx*sy))         *den2 + (nx*tx*(ny*tx+nx*ty))         *den2
    auxMatrix1(1,5) = (2*nx*nx*ny*nz) *den1 + (nx*sx*(nz*sy+ny*sz))         *den2 + (nx*tx*(nz*ty+ny*tz))         *den2
    auxMatrix1(1,6) = (2*nx*nx*nx*nz) *den1 + (nx*sx*(nz*sx+nx*sz))         *den2 + (nx*tx*(nz*tx+nx*tz))         *den2
    auxMatrix1(2,1) = (nx*nx*ny*ny)   *den1 + (nx*ny*sx*sy)                 *den2 + (nx*ny*tx*ty)                 *den2
    auxMatrix1(2,2) = (ny*ny*ny*ny)   *den1 + (ny*ny*sy*sy)                 *den2 + (ny*ny*ty*ty)                 *den2
    auxMatrix1(2,3) = (ny*ny*nz*nz)   *den1 + (ny*nz*sz*sy)                 *den2 + (ny*nz*ty*tz)                 *den2
    auxMatrix1(2,4) = (2*nx*ny*ny*ny) *den1 + (ny*sy*(ny*sx+nx*sy))         *den2 + (ny*ty*(ny*tx+nx*ty))         *den2
    auxMatrix1(2,5) = (2*ny*ny*ny*nz) *den1 + (ny*sy*(nz*sy+ny*sz))         *den2 + (ny*ty*(nz*ty+ny*tz))         *den2
    auxMatrix1(2,6) = (2*nx*ny*ny*nz) *den1 + (ny*sy*(nz*sx+nx*sz))         *den2 + (ny*ty*(nz*tx+nx*tz))         *den2
    auxMatrix1(3,1) = (nx*nx*nz*nz)   *den1 + (nx*nz*sx*sz)                 *den2 + (nx*nz*tx*tz)                 *den2
    auxMatrix1(3,2) = (ny*ny*nz*nz)   *den1 + (ny*nz*sy*sz)                 *den2 + (ny*nz*ty*tz)                 *den2
    auxMatrix1(3,3) = (nz*nz*nz*nz)   *den1 + (nz*nz*sz*sz)                 *den2 + (nz*nz*tz*tz)                 *den2
    auxMatrix1(3,4) = (2*nx*ny*nz*nz) *den1 + (nz*sz*(ny*sx+nx*sy))         *den2 + (nz*tz*(ny*tx+nx*ty))         *den2
    auxMatrix1(3,5) = (2*ny*nz*nz*nz) *den1 + (nz*sz*(nz*sy+ny*sz))         *den2 + (nz*tz*(nz*ty+ny*tz))         *den2
    auxMatrix1(3,6) = (2*nx*nz*nz*nz) *den1 + (nz*sz*(nz*sx+nx*sz))         *den2 + (nz*tz*(nz*tx+nx*tz))         *den2
    auxMatrix1(4,1) = (nx*nx*nx*ny)   *den1 + (nx*sx*(ny*sx+nx*sy))         *den3 + (nx*tx*(ny*tx+nx*ty))         *den3
    auxMatrix1(4,2) = (nx*ny*ny*ny)   *den1 + (ny*sy*(ny*sx+nx*sy))         *den3 + (ny*ty*(ny*tx+nx*ty))         *den3
    auxMatrix1(4,3) = (nx*ny*nz*nz)   *den1 + (nz*sz*(ny*sx+nx*sy))         *den3 + (nz*tz*(ny*tx+nx*ty))         *den3
    auxMatrix1(4,4) = (2*nx*nx*ny*ny) *den1 + ((ny*sx+nx*sy)*(ny*sx+nx*sy)) *den3 + ((ny*tx+nx*ty)*(ny*tx+nx*ty)) *den3
    auxMatrix1(4,5) = (2*nx*ny*ny*nz) *den1 + (ny*sx+nx*sy)*(nz*sy+ny*sz)   *den3 + (ny*tx+nx*ty)*(nz*ty+ny*tz)   *den3
    auxMatrix1(4,6) = (2*nx*nx*ny*nz) *den1 + (ny*sx+nx*sy)*(nz*sx+nx*sz)   *den3 + (ny*tx+nx*ty)*(nz*tx+nx*tz)   *den3
    auxMatrix1(5,1) = (nx*nx*ny*nz)   *den1 + (nx*sx*(nz*sy+ny*sz))         *den3 + (nx*tx*(nz*ty+ny*tz))         *den3
    auxMatrix1(5,2) = (nz*ny*ny*ny)   *den1 + (ny*sy*(nz*sy+ny*sz))         *den3 + (ny*ty*(nz*ty+ny*tz))         *den3
    auxMatrix1(5,3) = (ny*nz*nz*nz)   *den1 + (nz*sz*(nz*sy+ny*sz))         *den3 + (nz*tz*(nz*ty+ny*tz))         *den3
    auxMatrix1(5,4) = (2*ny*ny*nx*nz) *den1 + (ny*sx+nx*sy)*(nz*sy+ny*sz)   *den3 + (ny*tx+nx*ty)*(nz*ty+ny*tz)   *den3
    auxMatrix1(5,5) = (2*ny*ny*nz*nz) *den1 + ((nz*sy+ny*sz)*(nz*sy+ny*sz)) *den3 + ((nz*ty+ny*tz)*(nz*ty+ny*tz)) *den3
    auxMatrix1(5,6) = (2*nz*nz*ny*nx) *den1 + (nz*sx+nx*sz)*(nz*sy+ny*sz)   *den3 + (nz*tx+nx*tz)*(nz*ty+ny*tz)   *den3
    auxMatrix1(6,1) = (nx*nx*nx*nz)   *den1 + (nx*sx*(nz*sx+nx*sz))         *den3 + (nx*tx*(nz*tx+nx*tz))         *den3
    auxMatrix1(6,2) = (nx*ny*ny*nz)   *den1 + (ny*sy*(nz*sx+nx*sz))         *den3 + (ny*ty*(nz*tx+nx*tz))         *den3
    auxMatrix1(6,3) = (nx*nz*nz*nz)   *den1 + (nz*sz*(nx*sz+nz*sx))         *den3 + (nz*tz*(nz*tx+nx*tz))         *den3
    auxMatrix1(6,4) = (2*nx*nx*ny*nz) *den1 + (ny*sx+nx*sy)*(nz*sx+nx*sz)   *den3 + (ny*tx+nx*ty)*(nz*tx+nx*tz)   *den3
    auxMatrix1(6,5) = (2*nz*nz*nx*ny) *den1 + (nz*sx+nx*sz)*(nz*sy+ny*sz)   *den3 + (nz*tx+nx*tz)*(nz*ty+ny*tz)   *den3
    auxMatrix1(6,6) = (2*nz*nz*nx*nx) *den1 + ((nz*sx+nx*sz)*(nz*sx+nx*sz)) *den3 + ((nz*tx+nx*tz)*(nz*tx+nx*tz)) *den3
    !
    auxMatrix1(1,7) = -nx*nx*nx   - nx*sx*sx - nx*tx*tx
    auxMatrix1(1,8) = -nx*nx*ny   - nx*sx*sy - nx*tx*ty
    auxMatrix1(1,9) = -nx*nx*nz   - nx*sx*sz - nx*tx*tz
    auxMatrix1(2,7) = -nx*ny*ny   - ny*sx*sy - ny*tx*ty
    auxMatrix1(2,8) = -ny*ny*ny   - ny*sy*sy - ny*ty*ty
    auxMatrix1(2,9) = -ny*ny*nz   - ny*sy*sz - ny*ty*tz
    auxMatrix1(3,7) = -nx*nz*nz   - nz*sz*sx - nz*tx*tz
    auxMatrix1(3,8) = -ny*nz*nz   - nz*sz*sy - nz*ty*tz
    auxMatrix1(3,9) = -nz*nz*nz   - nz*sz*sz - nz*tz*tz
    auxMatrix1(4,7) = -nx*nx*ny - 0.5*sx*(ny*sx+nx*sy) - 0.5*tx*(ny*tx+nx*ty)
    auxMatrix1(4,8) = -nx*ny*ny - 0.5*sy*(ny*sx+nx*sy) - 0.5*ty*(ny*tx+nx*ty)
    auxMatrix1(4,9) = -nx*ny*nz - 0.5*sz*(ny*sx+nx*sy) - 0.5*tz*(ny*tx+nx*ty)
    auxMatrix1(5,7) = -nx*ny*nz - 0.5*sx*(nz*sy+ny*sz) - 0.5*tx*(nz*ty+ny*tz)
    auxMatrix1(5,8) = -nz*ny*ny - 0.5*sy*(nz*sy+ny*sz) - 0.5*ty*(nz*ty+ny*tz)
    auxMatrix1(5,9) = -ny*nz*nz - 0.5*sz*(nz*sy+ny*sz) - 0.5*tz*(nz*ty+ny*tz)
    auxMatrix1(6,7) = -nx*nx*nz - 0.5*sx*(nz*sx+nx*sz) - 0.5*tx*(nz*tx+nx*tz)
    auxMatrix1(6,8) = -nz*ny*nx - 0.5*sy*(nz*sx+nx*sz) - 0.5*ty*(nz*tx+nx*tz)
    auxMatrix1(6,9) = -nx*nz*nz - 0.5*sz*(nz*sx+nx*sz) - 0.5*tz*(nz*tx+nx*tz)

  END SUBROUTINE AnelasticAuxMatrix3D_Iso

  ! Compute the anelastic auxiliary matrices on 3D
  SUBROUTINE AnelasticAuxMatrix3D_Ani(AuxMatrix1,n1,n2,n3,w_speed,AniVec,EQN)
    !-------------------------------------------------------------------------!
    
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tEquations):: EQN                                                    ! 
    REAL            :: n1(3),n2(3),n3(3)                                      ! Component of normal vector
    REAL            :: AuxMatrix1(6,9)
    REAL            :: w_speed(EQN%nNonZeroEV),AniVec(3,3)
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: w_speed, n1, n2, n3, AniVec,EQN                        ! 
    INTENT(OUT)     :: AuxMatrix1                                             !
    !-------------------------------------------------------------------------!
    REAL            :: AuxMatrix3(EQN%nVarTotal,EQN%nVarTotal), T(EQN%nVarTotal,EQN%nVarTotal), iT(EQN%nVarTotal,EQN%nVarTotal)
    REAL            :: nx,ny,nz
    REAL            :: sx,sy,sz
    REAL            :: tx,ty,tz
    REAL            :: inv_speed
    INTEGER         :: i,j
    !-------------------------------------------------------------------------!

    nx=n1(1); ny=n1(2); nz=n1(3);
    sx=n2(1); sy=n2(2); sz=n2(3);
    tx=n3(1); ty=n3(2); tz=n3(3);

    !
    auxMatrix3(:,:)=0.
    !
    T(:,:)  = 0.
    iT(:,:) = 0.             
    ! u' = iT * u 
    ! u   = T * u'
    ! Inverse Transformation matrix
    T(1,1) = nx * nx
    T(1,2) = sx * sx
    T(1,3) = tx * tx
    T(1,4) = 2 * nx * sx
    T(1,5) = 2 * sx * tx
    T(1,6) = 2 * nx * tx
    T(2,1) = ny * ny
    T(2,2) = sy * sy
    T(2,3) = ty * ty
    T(2,4) = 2 * ny * sy
    T(2,5) = 2 * sy * ty
    T(2,6) = 2 * ny * ty
    T(3,1) = nz * nz
    T(3,2) = sz * sz
    T(3,3) = tz * tz
    T(3,4) = 2 * nz * sz
    T(3,5) = 2 * sz * tz
    T(3,6) = 2 * nz * tz
    T(4,1) = ny * nx
    T(4,2) = sy * sx
    T(4,3) = ty * tx
    T(4,4) = ny * sx + nx * sy
    T(4,5) = sy * tx + sx * ty
    T(4,6) = ny * tx + nx * ty
    T(5,1) = nz * ny
    T(5,2) = sz * sy
    T(5,3) = tz * ty
    T(5,4) = nz * sy + ny * sz
    T(5,5) = sz * ty + sy * tz
    T(5,6) = nz * ty + ny * tz
    T(6,1) = nz * nx
    T(6,2) = sz * sx
    T(6,3) = tz * tx
    T(6,4) = nz * sx + nx * sz
    T(6,5) = sz * tx + sx * tz
    T(6,6) = nz * tx + nx * tz
    T(7,7) = nx
    T(7,8) = sx
    T(7,9) = tx
    T(8,7) = ny
    T(8,8) = sy
    T(8,9) = ty
    T(9,7) = nz
    T(9,8) = sz
    T(9,9) = tz
    !
    ! Transformation matrix
    iT(1,1) = nx * nx
    iT(1,2) = ny * ny
    iT(1,3) = nz * nz
    iT(1,4) = 2 * ny * nx
    iT(1,5) = 2 * nz * ny
    iT(1,6) = 2 * nz * nx
    iT(2,1) = sx * sx
    iT(2,2) = sy * sy
    iT(2,3) = sz * sz
    iT(2,4) = 2 * sy * sx
    iT(2,5) = 2 * sz * sy
    iT(2,6) = 2 * sz * sx
    iT(3,1) = tx * tx
    iT(3,2) = ty * ty
    iT(3,3) = tz * tz
    iT(3,4) = 2 * ty * tx
    iT(3,5) = 2 * tz * ty
    iT(3,6) = 2 * tz * tx
    iT(4,1) = nx * sx
    iT(4,2) = ny * sy
    iT(4,3) = nz * sz
    iT(4,4) = ny * sx + nx * sy
    iT(4,5) = nz * sy + ny * sz
    iT(4,6) = nz * sx + nx * sz
    iT(5,1) = sx * tx
    iT(5,2) = sy * ty
    iT(5,3) = sz * tz
    iT(5,4) = sy * tx + sx * ty
    iT(5,5) = sz * ty + sy * tz
    iT(5,6) = sz * tx + sx * tz
    iT(6,1) = nx * tx
    iT(6,2) = ny * ty
    iT(6,3) = nz * tz
    iT(6,4) = ny * tx + nx * ty
    iT(6,5) = nz * ty + ny * tz
    iT(6,6) = nz * tx + nx * tz
    iT(7,7) = nx
    iT(7,8) = ny
    iT(7,9) = nz
    iT(8,7) = sx
    iT(8,8) = sy
    iT(8,9) = sz
    iT(9,7) = tx
    iT(9,8) = ty
    iT(9,9) = tz

    DO j=1,3
      inv_speed = 1.d0 / w_speed(j)
      auxMatrix3(1,1) = auxMatrix3(1,1) + AniVec(1,j) * AniVec(1,j) * inv_speed
      auxMatrix3(1,4) = auxMatrix3(1,4) + AniVec(1,j) * AniVec(2,j) * inv_speed
      auxMatrix3(1,6) = auxMatrix3(1,6) + AniVec(1,j) * AniVec(3,j) * inv_speed
      auxMatrix3(4,1) = auxMatrix3(4,1) + AniVec(2,j) * AniVec(1,j) * inv_speed * 0.5
      auxMatrix3(4,4) = auxMatrix3(4,4) + AniVec(2,j) * AniVec(2,j) * inv_speed * 0.5
      auxMatrix3(4,6) = auxMatrix3(4,6) + AniVec(2,j) * AniVec(3,j) * inv_speed * 0.5
      auxMatrix3(6,1) = auxMatrix3(6,1) + AniVec(3,j) * AniVec(1,j) * inv_speed * 0.5
      auxMatrix3(6,4) = auxMatrix3(6,4) + AniVec(3,j) * AniVec(2,j) * inv_speed * 0.5
      auxMatrix3(6,6) = auxMatrix3(6,6) + AniVec(3,j) * AniVec(3,j) * inv_speed * 0.5
    ENDDO
    !
    auxMatrix3(1,7) = -1.d0
    auxMatrix3(4,8) = -0.5d0
    auxMatrix3(6,9) = -0.5d0
    !
    auxMatrix3 = MATMUL( T,MATMUL(auxMatrix3,iT) )
    !
    auxMatrix1(:,:) = auxMatrix3(1:6,1:9)


  END SUBROUTINE AnelasticAuxMatrix3D_Ani


END MODULE JacobiNormal_mod
