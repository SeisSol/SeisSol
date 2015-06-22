!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) 2006, SeisSol Group
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

MODULE ConvertXiEta2XY_mod
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE ConvertXiEta2XY
     MODULE PROCEDURE ConvertXiEta2XY
  END INTERFACE
  INTERFACE ConvertXiEtaZeta2XYZ
     MODULE PROCEDURE ConvertXiEtaZeta2XYZ
  END INTERFACE
  INTERFACE ConvertSuperXiEta2XY
     MODULE PROCEDURE ConvertSuperXiEta2XY
  END INTERFACE
  INTERFACE MirrorXYWallPolynomial
     MODULE PROCEDURE MirrorXYWallPolynomial
  END INTERFACE
  INTERFACE XiEtaMonomialBasis
     MODULE PROCEDURE XiEtaMonomialBasis
  END INTERFACE
  INTERFACE MovePolynomial
     MODULE PROCEDURE MovePolynomial
  END INTERFACE
  INTERFACE RotatePolynomial
     MODULE PROCEDURE RotatePolynomial
  END INTERFACE
  INTERFACE XYBoundaryPolynomial
     MODULE PROCEDURE XYBoundaryPolynomial
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: ConvertXiEta2XY
  PUBLIC  :: ConvertXiEtaZeta2XYZ
  PUBLIC  :: ConvertSuperXiEta2XY
  PUBLIC  :: MirrorXYWallPolynomial
  PUBLIC  :: XiEtaMonomialBasis
  PUBLIC  :: MovePolynomial
  PUBLIC  :: RotatePolynomial
  PUBLIC  :: XYBoundaryPolynomial
  !----------------------------------------------------------------------------

CONTAINS

  !
  ! Converts the DG polynomial in xi-eta space into a polynomial in x-y space
  ! Needed by nonlinear ADER-DG method for the calculation of x-y derivatives.
  !
  SUBROUTINE ConvertXiEta2XY(XYcPoly,u_hat,cpoly,x,y,nVar,nDegFr,    &                           
       nOrdPoly,etype,IO_errUnit)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: nVar                                      ! Number of variables
    INTEGER   :: nDegFr                                    ! No. of degrees of freedom for DG
    INTEGER   :: nOrdPoly                                  ! Order of the basis polynomial
    INTEGER   :: etype                                     ! Elementtype (3=tri, 4=qua)
    INTEGER   :: IO_errUnit                                ! Unit number for writing errors
    REAL      :: cPoly(0:nOrdPoly,0:nOrdPoly,0:nDegFr-1)   ! DG DOF Poly coefficients
    REAL      :: XYcPoly(nVar,0:nOrdPoly,0:nOrdPoly)       ! Converted xy polynomial
    REAL      :: u_hat(0:nDegFr-1,nVar)                    ! Values of the DOF in DG xi-eta base
    REAL      :: x(eType), y(eType)                        ! Element vertices
    ! Local variable declaration
    INTEGER   :: i,j,k,l,iVar                              ! Loop counters 
    REAL      :: xi0,eta0,xi_x,xi_y,eta_x,eta_y,JacobiDet  ! Metrics
    REAL      :: XiEta2XY(0:nOrdPoly,0:nOrdPoly,0:nOrdPoly,0:nOrdPoly) ! Conversion Tensor x^k,y^l,xi^i,eta^j
    REAL      :: a(nVar,0:nOrdPoly,0:nOrdPoly)             ! xi-eta coefficients 
    !--------------------------------------------------------------------------
    INTENT(IN)      :: nVar, nDegFr, nOrdPoly, etype, IO_errUnit
    INTENT(IN)      :: u_hat, cpoly, x, y
    INTENT(OUT)     :: XYcpoly
    !--------------------------------------------------------------------------

    ! Determinant of the Jacobian of the x-y to xi-eta transformation
    JacobiDet = (x(2)-x(1))*(y(eType)-y(1)) - (x(eType)-x(1))*(y(2)-y(1))
    ! xi = xi0 + xi_x*x + xi_y*y
    xi0   = ( x(eType)*y(1) - x(1)*y(eType) )/JacobiDet
    xi_x  = ( y(eType)- y(1) ) / JacobiDet
    xi_y  = ( x(1)- x(eType) ) / JacobiDet
    ! eta = eta0 + eta_x*x + eta_y*y
    eta0  = ( x(1)*y(2) - x(2)*y(1) ) / JacobiDet
    eta_x = ( y(1) - y(2) ) / JacobiDet
    eta_y = ( x(2) - x(1) ) / JacobiDet

    ! Calculate xi-eta polynomial coefficients a_ij for each variable
    DO i = 0, nOrdPoly
       DO j = 0, nOrdPoly-i
          DO iVar = 1, nVar
             a(iVar,i,j) = DOT_PRODUCT( cPoly(i,j,:), u_hat(:,iVar) )
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion tensor
    XiEta2XY(:,:,:,:) = 0.
    ! First element for both exponents equal to 0 is 1
    XiEta2XY(0,0,0,0) = 1.

    ! For the polynomial, we have the following identity: P(xi(x,y),eta(x,y)) = P(x,y) =>
    ! sum_i sum_j a_ij * (xi0 + xi_x*x + xi_y*y)^i * (eta0 + eta_x*x + eta_y*y)^j = 
    ! sum_k sum_l c_kl * x^k * y^l

    ! Initialize conversion for pure xi^i
    j = 0
    DO i = 1, nOrdPoly 
       DO l = 0, i+j-1
          DO k = 0,i+j-1-l
             XiEta2XY(k  ,l  ,i,j)   = XiEta2XY(k  ,l  ,i,j)  +  XiEta2XY(k,l,i-1,j) * xi0
             XiEta2XY(k+1,l  ,i,j)   = XiEta2XY(k+1,l  ,i,j)  +  XiEta2XY(k,l,i-1,j) * xi_x
             XiEta2XY(k  ,l+1,i,j)   = XiEta2XY(k  ,l+1,i,j)  +  XiEta2XY(k,l,i-1,j) * xi_y
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion for mixed xi^i * eta^j
    DO j = 1, nOrdPoly
       DO i = 0, nOrdPoly-j 
          DO l = 0, i+j-1
             DO k = 0, i+j-1-l
                XiEta2XY(k  ,l  ,i,j)   = XiEta2XY(k  ,l  ,i,j)  +  XiEta2XY(k,l,i,j-1) * eta0
                XiEta2XY(k+1,l  ,i,j)   = XiEta2XY(k+1,l  ,i,j)  +  XiEta2XY(k,l,i,j-1) * eta_x
                XiEta2XY(k  ,l+1,i,j)   = XiEta2XY(k  ,l+1,i,j)  +  XiEta2XY(k,l,i,j-1) * eta_y
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from xi, eta exponents to x, y exponents.
    XYcPoly(:,:,:) = 0.
    DO l = 0, nOrdPoly
       DO k = 0, nOrdPoly-l ! If total polynomial degree is reached, then end loop
          DO j = 0, nOrdPoly
             ! Start the next loop so that xi and eta exponents can reach the 
             ! xy exponent but such they do not become negative
             DO i = MAX(0,l+k-j), nOrdPoly-j 
                XYcPoly(:,k,l) = XYcPoly(:,k,l) + a(:,i,j)*XiEta2XY(k,l,i,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE ConvertXiEta2XY

  !
  ! Converts the DG polynomial in xi-eta-zeta space into a polynomial in x-y-z space
  ! Needed by nonlinear ADER-DG method for the calculation of x-y-z derivatives.
  !
  SUBROUTINE ConvertXiEtaZeta2XYZ(XYZcPoly,u_hat,cpoly,x,y,z,nVar,nDegFr,    &                           
       nOrdPoly)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: nVar                                      ! Number of variables
    INTEGER   :: nDegFr                                    ! No. of degrees of freedom for DG
    INTEGER   :: nOrdPoly                                  ! Order of the basis polynomial
    INTEGER   :: etype                                     ! Elementtype (3=tri, 4=qua)
    REAL      :: cPoly(0:nOrdPoly,0:nOrdPoly,0:nOrdPoly,0:nDegFr-1)   ! DG DOF Poly coefficients
    REAL      :: XYZcPoly(nVar,0:nOrdPoly,0:nOrdPoly,0:nOrdPoly)       ! Converted xyz polynomial
    REAL      :: u_hat(0:nDegFr-1,nVar)                    ! Values of the DOF in DG xi-eta base
    REAL      :: x(4), y(4), z(4)                          ! Element vertices
    ! Local variable declaration
    INTEGER   :: i,j,k,l,m,n,iVar                          ! Loop counters 
    REAL      :: xi0,eta0,zeta0                            ! Metrics
    REAL      :: xi_x,xi_y,xi_z                            ! Metrics
    REAL      :: eta_x,eta_y,eta_z                         ! Metrics
    REAL      :: zeta_x,zeta_y,zeta_z                      ! Metrics
    REAL      :: JacobiDet                                 ! Metrics
    REAL      :: XiEtaZeta2XYZ(0:nOrdPoly,0:nOrdPoly, &    !
                               0:nOrdPoly,0:nOrdPoly, &    !
                               0:nOrdPoly,0:nOrdPoly  )    ! Conversion Tensor x^l,y^m,z^n, xi^i,eta^j,zeta^k
    REAL      :: a(nVar,0:nOrdPoly,0:nOrdPoly,0:nOrdPoly)  ! xi-eta-zeta coefficients 
    REAL      :: point(3)
    !--------------------------------------------------------------------------
    INTENT(IN)      :: nVar, nDegFr, nOrdPoly
    INTENT(IN)      :: u_hat, cpoly, x, y, z
    INTENT(OUT)     :: XYZcPoly
    !--------------------------------------------------------------------------
    !
    ! Determinant of the Jacobian of the x-y-z to xi-eta-zeta transformation
    !
    JacobiDet = y(1)*x(4)*z(2)+y(1)*x(3)*z(4)-y(1)*x(4)*z(3)+x(2)*y(3)*z(4)- &
                y(1)*x(2)*z(4)+x(1)*y(2)*z(4)+x(3)*z(2)*y(4)-x(3)*z(1)*y(4)- &
                x(3)*y(2)*z(4)-x(1)*y(4)*z(2)+x(1)*y(4)*z(3)+y(1)*x(2)*z(3)+ &
                x(4)*z(1)*y(3)-x(1)*y(2)*z(3)-x(4)*z(1)*y(2)+x(1)*y(3)*z(2)+ &
                x(3)*z(1)*y(2)-x(2)*y(4)*z(3)-x(1)*y(3)*z(4)-y(1)*x(3)*z(2)+ &
                x(2)*z(1)*y(4)-x(2)*z(1)*y(3)+x(4)*y(2)*z(3)-x(4)*z(2)*y(3)
    !
    ! xi = xi0 + xi_x*x + xi_y*y + xi_z*z
    !
    xi0   = -(-z(4)*y(1)*x(3)-z(1)*x(4)*y(3)-x(1)*z(3)*y(4)+x(1)*z(4)*y(3) + &
               z(3)*y(1)*x(4)+z(1)*x(3)*y(4))/JacobiDet
    xi_x  = -(z(1)*y(3)-z(1)*y(4)-z(4)*y(3)-z(3)*y(1)+z(3)*y(4)+z(4)*y(1))/JacobiDet
    xi_y  = -(z(1)*x(4)-z(1)*x(3)+x(1)*z(3)-z(4)*x(1)-z(3)*x(4)+z(4)*x(3))/JacobiDet
    xi_z  = -(y(1)*x(3)-y(3)*x(1)-y(1)*x(4)-y(4)*x(3)+y(3)*x(4)+y(4)*x(1))/JacobiDet
    !
    ! eta = eta0 + eta_x*x + eta_y*y + eta_z*z
    !
    eta0  = ( y(2)*x(1)*z(4)-y(1)*z(4)*x(2)-y(4)*x(1)*z(2)+z(1)*y(4)*x(2) - &
              z(1)*y(2)*x(4)+y(1)*z(2)*x(4))/JacobiDet
    eta_x = ( y(1)*z(4)-y(2)*z(4)+y(2)*z(1)-y(1)*z(2)+y(4)*z(2)-z(1)*y(4))/JacobiDet
    eta_y = (-x(2)*z(1)+x(4)*z(1)+x(1)*z(2)+x(2)*z(4)-x(4)*z(2)-x(1)*z(4))/JacobiDet
    eta_z = (-x(4)*y(1)-y(2)*x(1)+x(1)*y(4)+y(2)*x(4)+y(1)*x(2)-y(4)*x(2))/JacobiDet
    !
    ! zeta = zeta0 + zeta_x*x + zeta_y*y + zeta_z*z
    !
    zeta0  = -(y(2)*x(1)*z(4)+y(4)*x(1)*z(3)-y(3)*x(1)*z(4)-y(4)*x(1)*z(2)+ &
               y(1)*z(4)*x(3)+y(4)*z(2)*x(3)-y(4)*x(2)*z(3)+y(3)*x(2)*z(4)- &
               z(1)*y(4)*x(3)-y(2)*z(4)*x(3)-z(1)*y(2)*x(4)-y(3)*z(2)*x(4)+ &
               z(1)*y(4)*x(2)+z(1)*y(3)*x(4)+y(2)*z(3)*x(4)+y(1)*z(2)*x(4)- &
               y(1)*z(3)*x(4)-y(1)*z(4)*x(2)-JacobiDet)/JacobiDet
    zeta_x = -( y(2)*z(1)-y(1)*z(2)-y(3)*z(1)+y(1)*z(3)+y(3)*z(2)-y(2)*z(3))/JacobiDet
    zeta_y = -(-x(2)*z(1)+x(2)*z(3)-x(1)*z(3)+x(3)*z(1)-x(3)*z(2)+x(1)*z(2))/JacobiDet
    zeta_z = -( x(1)*y(3)+y(1)*x(2)+y(2)*x(3)-x(3)*y(1)-y(3)*x(2)-y(2)*x(1))/JacobiDet  
    !
    ! Calculate xi-eta polynomial coefficients a_ijk for each variable
    !
    DO i = 0, nOrdPoly
       DO j = 0, nOrdPoly-i
          DO k = 0, nOrdPoly-i-j
              DO iVar = 1, nVar
                 a(iVar,i,j,k) = DOT_PRODUCT( cPoly(i,j,k,:), u_hat(:,iVar) )
              ENDDO
          ENDDO
       ENDDO
    ENDDO
    !
    ! Initialize conversion tensor
    !
    XiEtaZeta2XYZ(:,:,:,:,:,:) = 0.
    !
    ! First element for both exponents equal to 0 is 1
    ! 
    XiEtaZeta2XYZ(0,0,0,0,0,0) = 1.
    !
    ! For the polynomial, we have the following identity: P(xi(x,y,z),eta(x,y,z),zeta(x,y,z)) = P(x,y,z) =>
    ! sum_i sum_j sum_k a_ijk * (xi0 + xi_x*x + xi_y*y)^i * (eta0 + eta_x*x + eta_y*y)^j * (zeta0 + zeta_x*x + zeta_y*y)^k= 
    ! sum_l sum_m sum_n c_lmn * x^l * y^m * z^n

    ! Initialize conversion for pure xi^i
    j = 0
    k = 0
    DO i = 1, nOrdPoly 
     DO n = 0, i + j + k - 1
      DO m = 0, i + j + k - 1 - n
       DO l = 0, i + j + k - 1 - n - m
         XiEtaZeta2XYZ(l,m,n,i,j,k)   = XiEtaZeta2XYZ(l,m,n,i,j,k)   + XiEtaZeta2XYZ(l,m,n,i-1,j,k) * xi0
         XiEtaZeta2XYZ(l+1,m,n,i,j,k) = XiEtaZeta2XYZ(l+1,m,n,i,j,k) + XiEtaZeta2XYZ(l,m,n,i-1,j,k) * xi_x
         XiEtaZeta2XYZ(l,m+1,n,i,j,k) = XiEtaZeta2XYZ(l,m+1,n,i,j,k) + XiEtaZeta2XYZ(l,m,n,i-1,j,k) * xi_y
         XiEtaZeta2XYZ(l,m,n+1,i,j,k) = XiEtaZeta2XYZ(l,m,n+1,i,j,k) + XiEtaZeta2XYZ(l,m,n,i-1,j,k) * xi_z
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    k = 0
    ! Initialize conversion for mixed xi^i * eta^j but NOT zeta^k
    DO j = 1, nOrdPoly
     DO i = 0, nOrdPoly - k - j 
      DO n = 0, i + j + k - 1
       DO m = 0, i + j + k - 1 - n
        DO l = 0, i + j + k - 1 - n - m
          XiEtaZeta2XYZ(l,m,n,i,j,k)   = XiEtaZeta2XYZ(l,m,n,i,j,k)   + XiEtaZeta2XYZ(l,m,n,i,j-1,k) * eta0
          XiEtaZeta2XYZ(l+1,m,n,i,j,k) = XiEtaZeta2XYZ(l+1,m,n,i,j,k) + XiEtaZeta2XYZ(l,m,n,i,j-1,k) * eta_x
          XiEtaZeta2XYZ(l,m+1,n,i,j,k) = XiEtaZeta2XYZ(l,m+1,n,i,j,k) + XiEtaZeta2XYZ(l,m,n,i,j-1,k) * eta_y
          XiEtaZeta2XYZ(l,m,n+1,i,j,k) = XiEtaZeta2XYZ(l,m,n+1,i,j,k) + XiEtaZeta2XYZ(l,m,n,i,j-1,k) * eta_z
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    ! Initialize conversion for mixed xi^i * eta^j * zeta^k
    DO k = 1, nOrdPoly
     DO j = 0, nOrdPoly - k
      DO i = 0, nOrdPoly - k - j 
       DO n = 0, i + j + k - 1
        DO m = 0, i + j + k - 1 - n
         DO l = 0, i + j + k - 1 - n - m
           XiEtaZeta2XYZ(l,m,n,i,j,k)   = XiEtaZeta2XYZ(l,m,n,i,j,k)   + XiEtaZeta2XYZ(l,m,n,i,j,k-1) * zeta0
           XiEtaZeta2XYZ(l+1,m,n,i,j,k) = XiEtaZeta2XYZ(l+1,m,n,i,j,k) + XiEtaZeta2XYZ(l,m,n,i,j,k-1) * zeta_x
           XiEtaZeta2XYZ(l,m+1,n,i,j,k) = XiEtaZeta2XYZ(l,m+1,n,i,j,k) + XiEtaZeta2XYZ(l,m,n,i,j,k-1) * zeta_y
           XiEtaZeta2XYZ(l,m,n+1,i,j,k) = XiEtaZeta2XYZ(l,m,n+1,i,j,k) + XiEtaZeta2XYZ(l,m,n,i,j,k-1) * zeta_z
         ENDDO
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from xi, eta, zeta exponents to x, y, z exponents.
    !
    XYZcPoly(:,:,:,:) = 0.
    !
    DO n = 0, nOrdPoly
     DO m = 0, nOrdPoly - n
      DO l = 0, nOrdPoly - n - m    
         !
         DO k = 0, nOrdPoly 
          DO j = 0, nOrdPoly - k 
           DO i = 0, nOrdPoly - k - j 
                XYZcPoly(:,l,m,n) = XYZcPoly(:,l,m,n) + a(:,i,j,k)*XiEtaZeta2XYZ(l,m,n,i,j,k)
           ENDDO
          ENDDO
         ENDDO
         !
      ENDDO
     ENDDO
    ENDDO
    !
  END SUBROUTINE ConvertXiEtaZeta2XYZ

  !
  ! Converts the DG polynomial for superparametric elements 
  ! from the xi-eta space into a polynomial in x-y space 
  ! Needed by nonlinear ADER-DG method for the calculation of x-y derivatives.
  !
  SUBROUTINE ConvertSuperXiEta2XY(XYcPoly,u_hat,cpoly,nVar,nDegFr,    &                           
       nOrdPoly,alpha,beta,nOrdBnd,IO_errUnit)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: nVar                                                                   ! Number of variables
    INTEGER   :: nDegFr                                                                 ! No. of degrees of freedom for DG
    INTEGER   :: nOrdPoly                                                               ! Order of the basis polynomial
    INTEGER   :: nOrdBnd                                                                ! Order of boundary approximation
    INTEGER   :: IO_errUnit                                                             ! Unit number for writing errors
    REAL      :: XYcPoly(nVar,0:nOrdPoly*nOrdBnd,0:nOrdPoly*nOrdBnd)                    ! Converted xy polynomial   
    REAL      :: u_hat(0:nDegFr-1,nVar)                                                 ! Values of the DOF in DG xi-eta base
    REAL      :: cPoly(0:nOrdPoly,0:nOrdPoly,0:nDegFr-1)                                ! DG DOF Poly coefficients
    REAL      :: alpha(0:nOrdBnd,0:nOrdBnd)                                             ! Inverse mapping coefficients for xi
    REAL      :: beta(0:nOrdBnd,0:nOrdBnd)                                              ! Inverse mapping coefficients for eta
    ! Local variable declaration
    INTEGER   :: i,j,k,l,p,q,iVar                                                       ! Loop counters 
    REAL      :: XiEta2XY(0:nOrdPoly*nOrdBnd,0:nOrdPoly*nOrdBnd,0:nOrdPoly,0:nOrdPoly)  ! Conversion Tensor x^k,y^l,xi^i,eta^j
    REAL      :: a(nVar,0:nOrdPoly,0:nOrdPoly)                                          ! xi-eta coefficients 
    !--------------------------------------------------------------------------
    INTENT(IN)     :: nVar, nDegFr, nOrdPoly, nOrdBnd, IO_errUnit
!    INTENT(IN)     :: u_hat, cpoly, alpha, beta
    INTENT(OUT)    :: XYcpoly
    !--------------------------------------------------------------------------

    ! Calculate xi-eta polynomial coefficients a_ij for each variable
    DO i = 0, nOrdPoly
       DO j = 0, nOrdPoly-i
          DO iVar = 1, nVar
             a(iVar,i,j) = DOT_PRODUCT( cPoly(i,j,:), u_hat(:,iVar) )
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion tensor
    XiEta2XY(:,:,:,:) = 0.
    ! First element for both exponents equal to 0 is 1
    XiEta2XY(0,0,0,0) = 1.

    ! For the polynomial, we have the following identity: P(xi(x,y),eta(x,y)) = P(x,y) =>
    ! sum_i sum_j a_ij * ( sum_p sum_q alpha_pq x^p y^q )^i * ( sum_p sum_q beta_pq x^p y^q )^j = 
    ! sum_k sum_l c_kl * x^k * y^l

    ! Initialize conversion for pure xi^i
    j = 0
    DO i = 1, nOrdPoly 
       DO l = 0, (nOrdPoly-1)*nOrdBnd
          DO k = 0, (nOrdPoly-1)*nOrdBnd
             !
             DO p = 0, nOrdBnd
                DO q = 0, nOrdBnd-p
                   XiEta2XY(k+p,l+q,i,j) = XiEta2XY(k+p,l+q,i,j)  +  XiEta2XY(k,l,i-1,j) * alpha(p,q)
                ENDDO
             ENDDO
             !
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion for mixed xi^i * eta^j
    DO j = 1, nOrdPoly
       DO i = 0, nOrdPoly-j 
          DO l = 0, (nOrdPoly-1)*nOrdBnd
             DO k = 0, (nOrdPoly-1)*nOrdBnd
                ! 
                DO p = 0, nOrdBnd
                   DO q = 0, nOrdBnd-p
                      XiEta2XY(k+p,l+q,i,j) = XiEta2XY(k+p,l+q,i,j) + XiEta2XY(k,l,i,j-1)*beta(p,q)
                   ENDDO
                ENDDO
                !
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from xi, eta exponents to x, y exponents.
    XYcPoly(:,:,:) = 0.
    DO l = 0, nOrdPoly*nOrdBnd
       DO k = 0, nOrdPoly*nOrdBnd-l ! If total polynomial degree is reached, then end loop
          DO j = 0, nOrdPoly
             DO i = 0, nOrdPoly-j 
                XYcPoly(:,k,l) = XYcPoly(:,k,l) + a(:,i,j)*XiEta2XY(k,l,i,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO


  END SUBROUTINE ConvertSuperXiEta2XY

  !
  ! Mirrors the Taylor DG polynomial in x-y space, given by the Taylor series 
  ! at the wall, with respect to the outward unit normal vector n(:)
  ! We suppose that the polynomial contains the state vector with velocities
  ! already transformed to the edge-local system, i.e. UWall(2,:) contains the
  ! wall-normal velocity component (and its derivatives) and UWall(3,;) the
  ! tangential components.
  !
  SUBROUTINE MirrorXYWallPolynomial(UMirror,UWall,n,nVar,nDegFr,    &                           
       nOrdPoly,etype,Faculty,IO_errUnit)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: nVar                                      ! Number of variables
    INTEGER   :: nDegFr                                    ! No. of degrees of freedom for DG
    INTEGER   :: nOrdPoly                                  ! Order of the basis polynomial
    INTEGER   :: etype                                     ! Elementtype (3=tri, 4=qua)
    INTEGER   :: IO_errUnit                                ! Unit number for writing errors
    REAL      :: UMirror(1:nVar,0:nOrdPoly,0:nOrdPoly)     ! Mirror series poly coefficients
    REAL      :: Uwall(1:nVar,0:nOrdPoly,0:nOrdPoly)       ! Taylor series poly coefficients
    REAL      :: n(2)                                      ! Outward unit normal vector
    REAL      :: Faculty(0:nOrdPoly)                       ! Precalculated faculties
    ! Local variable declaration
    INTEGER   :: i,j,k,l,iVar                              ! Loop counters 
    REAL      :: XiEta2XY(0:nOrdPoly,0:nOrdPoly,0:nOrdPoly,0:nOrdPoly) ! Conversion Tensor x^k,y^l,xi^i,eta^j
    REAL      :: a(1:nVar,0:nOrdPoly,0:nOrdPoly)           ! Original xi-eta coefficients 
    REAL      :: b(1:nVar,0:nOrdPoly,0:nOrdPoly)           ! Mirrored xi-eta coefficients 
    !--------------------------------------------------------------------------
    INTENT(IN)      :: nVar, nDegFr, nOrdPoly, etype, IO_errUnit
    INTENT(IN)      :: UWall, n
    INTENT(OUT)     :: UMirror
    !--------------------------------------------------------------------------

    ! The original polynomial is given by the Taylor series at the wall
    ! The transformation into the local edge coordinate system is given by
    ! 
    ! x = xi*n(1) - eta*n(2)
    ! y = xi*n(2) + eta*n(1)
    ! 
    ! xi  =  x*n(1) + y*n(2) 
    ! eta = -x*n(2) + y*n(1)

    ! Initialize conversion tensor
    XiEta2XY(:,:,:,:) = 0.
    ! First element for both exponents equal to 0 is 1
    XiEta2XY(0,0,0,0) = 1.

    ! For the polynomial, we have the following identity: P(xi,eta) = P(x(xi,eta),y(xi,eta)) =>
    ! sum_k sum_l a_kl * xi^k * eta^l = 
    ! sum_i sum_j W_ij/(i! * j!) * [n(1)*xi - n(2)*eta]^i * [n(2)*xi + n(1)*eta]^j

    ! Initialize conversion for pure x^i
    j = 0
    DO i = 1, nOrdPoly 
       DO l = 0, i+j-1
          DO k = 0,i+j-1-l
             XiEta2XY(k+1,l  ,i,j) = XiEta2XY(k+1,l  ,i,j) + XiEta2XY(k,l,i-1,j)*n(1) 
             XiEta2XY(k  ,l+1,i,j) = XiEta2XY(k  ,l+1,i,j) - XiEta2XY(k,l,i-1,j)*n(2) 
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion for mixed x^i * y^j
    DO j = 1, nOrdPoly
       DO i = 0, nOrdPoly-j 
          DO l = 0, i+j-1
             DO k = 0, i+j-1-l
                XiEta2XY(k+1,l  ,i,j)   = XiEta2XY(k+1,l  ,i,j) + XiEta2XY(k,l,i,j-1)*n(2) 
                XiEta2XY(k  ,l+1,i,j)   = XiEta2XY(k  ,l+1,i,j) + XiEta2XY(k,l,i,j-1)*n(1) 
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from x, y exponents to xi, eta exponents.
    a(:,:,:) = 0.
    DO l = 0, nOrdPoly
       DO k = 0, nOrdPoly-l ! If total polynomial degree is reached, then end loop
          DO j = 0, nOrdPoly
             ! Start the next loop so that xi and eta exponents can reach the 
             ! xy exponent but such they do not become negative
             DO i = 0, nOrdPoly-j 
                a(:,k,l) = a(:,k,l) + UWall(:,i,j)/(Faculty(i)*Faculty(j))*XiEta2XY(k,l,i,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Now perform the mirroring:
    ! rho,v,p: multiply all even xi derivatives with +1, all odd derivatives with -1
    ! u      : multiply all even xi derivatives with -1, all odd derivatives with +1
    DO k = 0, nOrdPoly
     IF(MOD(k,2).EQ.0) THEN
       ! Even derivatives of xi
       DO l = 0, nOrdPoly     
         b(1,k,l) =  a(1,k,l)
         b(2,k,l) = -a(2,k,l)
         b(3,k,l) =  a(3,k,l)
         b(4,k,l) =  a(4,k,l)
       ENDDO
     ELSE
       ! Odd derivatives of xi
       DO l = 0, nOrdPoly     
         b(1,k,l) = -a(1,k,l)
         b(2,k,l) =  a(2,k,l)
         b(3,k,l) = -a(3,k,l)
         b(4,k,l) = -a(4,k,l)
       ENDDO
     ENDIF
    ENDDO

    ! ...and finally transform back to the x-y system

    ! Initialize conversion tensor
    XiEta2XY(:,:,:,:) = 0.
    ! First element for both exponents equal to 0 is 1
    XiEta2XY(0,0,0,0) = 1.

    ! For the polynomial, we have the following identity: P(xi(x,y),eta(x,y)) = P(x,y) =>
    ! sum_i sum_j b_ij * [n(1)*x + n(2)*y]^i * [-n(2)*x + n(1)*y]^j = 
    ! sum_k sum_l W_kl/(k! * l!) * x^k * y^l

    ! Initialize conversion for pure xi^i
    j = 0
    DO i = 1, nOrdPoly 
       DO l = 0, i+j-1
          DO k = 0,i+j-1-l
             XiEta2XY(k+1,l  ,i,j) = XiEta2XY(k+1,l  ,i,j) + XiEta2XY(k,l,i-1,j) * n(1) 
             XiEta2XY(k  ,l+1,i,j) = XiEta2XY(k  ,l+1,i,j) + XiEta2XY(k,l,i-1,j) * n(2) 
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion for mixed xi^i * eta^j
    DO j = 1, nOrdPoly
       DO i = 0, nOrdPoly-j 
          DO l = 0, i+j-1
             DO k = 0, i+j-1-l
                XiEta2XY(k+1,l  ,i,j)   = XiEta2XY(k+1,l  ,i,j) - XiEta2XY(k,l,i,j-1) * n(2) 
                XiEta2XY(k  ,l+1,i,j)   = XiEta2XY(k  ,l+1,i,j) + XiEta2XY(k,l,i,j-1) * n(1) 
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from xi, eta exponents to x, y exponents.
    UMirror(:,:,:) = 0.
    DO l = 0, nOrdPoly
       DO k = 0, nOrdPoly-l ! If total polynomial degree is reached, then end loop
          DO j = 0, nOrdPoly
             ! Start the next loop so that xi and eta exponents can reach the 
             ! xy exponent but such they do not become negative
             DO i = 0, nOrdPoly-j 
                UMirror(:,k,l) = UMirror(:,k,l) + b(:,i,j)*XiEta2XY(k,l,i,j)*(Faculty(k)*Faculty(l))
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE MirrorXYWallPolynomial


  !
  ! Converts the DG polynomial in xi-eta space into the monomial basis in 
  ! xi-eta space
  !
  SUBROUTINE XiEtaMonomialBasis(MonomialBasis,u_hat,cpoly,nVar,nDegFr,    &                           
       nOrdPoly,IO_errUnit)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: nVar                                      ! Number of variables
    INTEGER   :: nDegFr                                    ! No. of degrees of freedom for DG
    INTEGER   :: nOrdPoly                                  ! Order of the basis polynomial
    INTEGER   :: IO_errUnit                                ! Unit number for writing errors
    REAL      :: u_hat(0:nDegFr-1,nVar)                    ! Values of the DOF in DG xi-eta base
    REAL      :: cPoly(0:nOrdPoly,0:nOrdPoly,0:nDegFr-1)   ! DG DOF Poly coefficients
    REAL      :: MonomialBasis(nVar,0:nOrdPoly,0:nOrdPoly) ! xi-eta coefficients in monomial basis:
    !                                                      ! P(xi,eta) = a_ij * xi^i * eta^j
    ! Local variable declaration
    INTEGER   :: i,j,iVar                                  ! Loop counters 
    !--------------------------------------------------------------------------
    INTENT(IN)      :: nVar, nDegFr, nOrdPoly, IO_errUnit
    INTENT(IN)      :: u_hat, cpoly
    INTENT(OUT)     :: MonomialBasis
    !--------------------------------------------------------------------------

    ! Calculate xi-eta polynomial coefficients a_ij for each variable
    ! (monomial basis)
    DO i = 0, nOrdPoly
       DO j = 0, nOrdPoly-i
          DO iVar = 1, nVar
             MonomialBasis(iVar,i,j) = DOT_PRODUCT( cPoly(i,j,:), u_hat(:,iVar) )
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE XiEtaMonomialBasis


  !
  ! Moves a polynomial (defined by monomials) by the transformation
  ! xi   = x  - delta_xi
  ! eta  = y  - delta_eta
  ! The input polynomial is supposed to be given in the xi-eta system,
  ! the output polynomial is then given in the x-y system
  !
  SUBROUTINE MovePolynomial(Poly_out,Poly_in,delta_xi,delta_eta,nVar, &
       nOrdPoly,IO_errUnit,XiEta2XY_out)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: nVar                                      ! Number of variables
    INTEGER   :: nOrdPoly                                  ! Order of the basis polynomial
    INTEGER   :: IO_errUnit                                ! Unit number for writing errors
    REAL      :: Poly_out(nVar,0:nOrdPoly,0:nOrdPoly)      ! xi-eta coefficients 
    REAL      :: Poly_in(nVar,0:nOrdPoly,0:nOrdPoly)       ! xi-eta coefficients 
    REAL      :: delta_xi, delta_eta                       ! Deplacement of the polynomial
    REAL, OPTIONAL :: XiEta2XY_out(0:nOrdPoly,0:nOrdPoly,0:nOrdPoly,0:nOrdPoly) ! Conversion Tensor x^k,y^l,xi^i,eta^j
    ! Local variable declaration
    INTEGER   :: i,j,k,l                                   ! Loop counters 
    REAL      :: XiEta2XY(0:nOrdPoly,0:nOrdPoly,0:nOrdPoly,0:nOrdPoly) ! Conversion Tensor x^k,y^l,xi^i,eta^j
    !--------------------------------------------------------------------------
    INTENT(IN)      :: nVar, nOrdPoly, IO_errUnit
    INTENT(IN)      :: delta_xi, delta_eta
    INTENT(IN)      :: Poly_in
    INTENT(OUT)     :: Poly_out
    !--------------------------------------------------------------------------

    ! Initialize conversion tensor
    XiEta2XY(:,:,:,:) = 0.
    ! First element for both exponents equal to 0 is 1
    XiEta2XY(0,0,0,0) = 1.

    ! For the polynomial, we have the following identity: P(xi(x,y),eta(x,y)) = P(x,y) =>
    ! sum_i sum_j a_ij * (x-delta_xi)^i * (y-delta_eta)^j = 
    ! sum_k sum_l c_kl * x^k * y^l

    ! Initialize conversion for pure xi^i
    j = 0
    DO i = 1, nOrdPoly 
       DO l = 0, i+j-1
          DO k = 0,i+j-1-l
             XiEta2XY(k  ,l  ,i,j)      = XiEta2XY(k  ,l  ,i,j)  -  XiEta2XY(k,l,i-1,j) * delta_xi
             XiEta2XY(k+1,l  ,i,j)      = XiEta2XY(k+1,l  ,i,j)  +  XiEta2XY(k,l,i-1,j) * 1.          ! here: xi_x = 1. 
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion for mixed xi^i * eta^j
    DO j = 1, nOrdPoly
       DO i = 0, nOrdPoly-j 
          DO l = 0, i+j-1
             DO k = 0, i+j-1-l
                XiEta2XY(k  ,l  ,i,j)   = XiEta2XY(k  ,l  ,i,j)  -  XiEta2XY(k,l,i,j-1) * delta_eta
                XiEta2XY(k  ,l+1,i,j)   = XiEta2XY(k  ,l+1,i,j)  +  XiEta2XY(k,l,i,j-1) * 1.         ! here: eta_y = 1.
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    IF(PRESENT(XiEta2XY_out)) THEN
      XiEta2XY_out(:,:,:,:) = XiEta2XY(:,:,:,:)
      RETURN
    ENDIF

    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from xi, eta exponents to x, y exponents.
    Poly_out(:,:,:) = 0.
    DO l = 0, nOrdPoly
       DO k = 0, nOrdPoly-l ! If total polynomial degree is reached, then end loop
          DO j = 0, nOrdPoly
             DO i = 0, nOrdPoly-j 
                Poly_out(:,k,l) = Poly_out(:,k,l) + Poly_in(:,i,j)*XiEta2XY(k,l,i,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE MovePolynomial

  !
  ! Rotates a polynomial (defined by monomials) by the transformation
  !
  ! x  = nx*xi - ny*eta
  ! y  = ny*xi + nx*eta
  ! 
  ! xi   =  nx*x + ny*y
  ! eta  = -ny*x + nx*y
  !
  ! The input polynomial is supposed to be given in the xi-eta system,
  ! the output polynomial is then given in the x-y system
  !
  SUBROUTINE RotatePolynomial(Poly_out,Poly_in,nx,ny,nVar,nOrdPoly,IO_errUnit)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: nVar                                      ! Number of variables
    INTEGER   :: nOrdPoly                                  ! Order of the basis polynomial
    INTEGER   :: IO_errUnit                                ! Unit number for writing errors
    REAL      :: Poly_out(nVar,0:nOrdPoly,0:nOrdPoly)      ! xi-eta coefficients 
    REAL      :: Poly_in(nVar,0:nOrdPoly,0:nOrdPoly)       ! xi-eta coefficients 
    REAL      :: nx, ny                                    ! Deplacement of the polynomial
    ! Local variable declaration
    INTEGER   :: i,j,k,l                                   ! Loop counters 
    REAL      :: XiEta2XY(0:nOrdPoly,0:nOrdPoly,0:nOrdPoly,0:nOrdPoly) ! Conversion Tensor x^k,y^l,xi^i,eta^j
    !--------------------------------------------------------------------------
    INTENT(IN)      :: nVar, nOrdPoly, IO_errUnit
    INTENT(IN)      :: nx, ny
    INTENT(IN)      :: Poly_in
    INTENT(OUT)     :: Poly_out
    !--------------------------------------------------------------------------

    ! Initialize conversion tensor
    XiEta2XY(:,:,:,:) = 0.
    ! First element for both exponents equal to 0 is 1
    XiEta2XY(0,0,0,0) = 1.

    ! For the polynomial, we have the following identity: P(xi(x,y),eta(x,y)) = P(x,y) =>
    ! sum_i sum_j a_ij * (nx*x + ny*y)^i * (-ny*x + nx*y)^j = 
    ! sum_k sum_l c_kl * x^k * y^l

    ! Initialize conversion for pure xi^i
    j = 0
    DO i = 1, nOrdPoly 
       DO l = 0, i+j-1
          DO k = 0,i+j-1-l
             XiEta2XY(k+1,l  ,i,j)   = XiEta2XY(k+1,l  ,i,j)  +  XiEta2XY(k,l,i-1,j) * nx
             XiEta2XY(k  ,l+1,i,j)   = XiEta2XY(k  ,l+1,i,j)  +  XiEta2XY(k,l,i-1,j) * ny
          ENDDO
       ENDDO
    ENDDO

    ! Initialize conversion for mixed xi^i * eta^j
    DO j = 1, nOrdPoly
       DO i = 0, nOrdPoly-j 
          DO l = 0, i+j-1
             DO k = 0, i+j-1-l
                XiEta2XY(k+1,l  ,i,j)   = XiEta2XY(k+1,l  ,i,j)  -  XiEta2XY(k,l,i,j-1) * ny
                XiEta2XY(k  ,l+1,i,j)   = XiEta2XY(k  ,l+1,i,j)  +  XiEta2XY(k,l,i,j-1) * nx
             ENDDO
          ENDDO
       ENDDO
    ENDDO

    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from xi, eta exponents to x, y exponents.
    Poly_out(:,:,:) = 0.
    DO l = 0, nOrdPoly
       DO k = 0, nOrdPoly-l ! If total polynomial degree is reached, then end loop
          DO j = 0, nOrdPoly
             DO i = 0, nOrdPoly-j 
                Poly_out(:,k,l) = Poly_out(:,k,l) + Poly_in(:,i,j)*XiEta2XY(k,l,i,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE RotatePolynomial


  !
  ! Returns the coefficients of the boundary basis polynomials for all DOF
  ! in physical coordinates. The polynomial is expanded about the side midpoint 
  ! and is defined in the interval [-L/2 ; +L/2]. 
  !
  ! [ xi (chi) ] = [xi0  ] + chi/L * [dir_xi ]
  ! [ eta(chi) ] = [eta0 ] + chi/L * [dir_eta]
  !
  SUBROUTINE XYBoundaryPolynomial(Boundary_Basis, Length, cPoly, iSide, eType, &
     nVar, nDegFr, nOrdPoly, IO_errUnit)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER   :: iSide                                     ! Local side number
    INTEGER   :: eType                                     ! Element type
    INTEGER   :: nVar                                      ! Number of variables
    INTEGER   :: iDegFr                                    ! Number of DOF
    INTEGER   :: nDegFr                                    ! Number of DOF
    INTEGER   :: nOrdPoly                                  ! Order of the basis polynomial
    INTEGER   :: IO_errUnit                                ! Unit number for writing errors
    REAL      :: Length                                    ! Physical side length
    REAL      :: cPoly(0:nOrdPoly,0:nOrdPoly,0:nDegFr-1)   ! DG DOF Poly coefficients
    REAL      :: Boundary_Basis(0:nOrdPoly,1:nDegFr)             ! chi coefficients in x-y space
    REAL      :: nx, ny                                    ! Deplacement of the polynomial
    ! Local variable declaration
    INTEGER   :: i,j,l                                     ! Loop counters 
    REAL      :: xi0,eta0,dir_xi,dir_eta
    REAL      :: XiEta2Chi(0:nOrdPoly,0:nOrdPoly,0:nOrdPoly) ! Conversion Tensor chi^k, xi^i,eta^j
    !--------------------------------------------------------------------------
    INTENT(IN)      :: Length, nVar, nOrdPoly, nDegFr, cPoly, IO_errUnit
    INTENT(IN)      :: iSide, eType
    INTENT(OUT)     :: Boundary_Basis
    !--------------------------------------------------------------------------

    SELECT CASE(eType)
    CASE(3)
      SELECT CASE(iSide)
      CASE(1)
        xi0     =  0.5
        eta0    =  0.0
        dir_xi  =  1.0 
        dir_eta =  0.0
      CASE(2)
        xi0     =  0.5
        eta0    =  0.5
        dir_xi  = -1.0 
        dir_eta =  1.0
      CASE(3)
        xi0     =  0.0
        eta0    =  0.5
        dir_xi  =  0.0 
        dir_eta = -1.0
      END SELECT
    CASE(4)
      SELECT CASE(iSide)
      CASE(1)
        xi0     =  0.5
        eta0    =  0.0
        dir_xi  =  1.0 
        dir_eta =  0.0
      CASE(2)
        xi0     =  1.0
        eta0    =  0.5
        dir_xi  =  0.0 
        dir_eta =  1.0
      CASE(3)
        xi0     =  0.5
        eta0    =  1.0
        dir_xi  = -1.0
        dir_eta =  0.0
      CASE(4)
        xi0     =  0.0
        eta0    =  0.5
        dir_xi  =  0.0 
        dir_eta = -1.0
      END SELECT
    END SELECT

    ! Initialize conversion tensor
    XiEta2Chi(:,:,:) = 0.
    ! First element for both exponents equal to 0 is 1
    XiEta2Chi(0,0,0) = 1.

    ! For the polynomial, we have the following identity: P(xi(chi),eta(chi)) = P(chi) =>
    ! sum_i sum_j a_ij * (xi0 + chi/L * dir_xi)^i * (eta0 + chi/L * dir_eta)^j = 
    ! sum_k sum_l c_l * chi^l

    ! Initialize conversion for pure xi^i
    j = 0
    DO i = 1, nOrdPoly 
       DO l = 0, i+j-1
          XiEta2Chi(l  ,i,j) = XiEta2Chi(l  ,i,j)  +  XiEta2Chi(l,i-1,j) * xi0
          XiEta2Chi(l+1,i,j) = XiEta2Chi(l+1,i,j)  +  XiEta2Chi(l,i-1,j) * dir_xi/Length
       ENDDO
    ENDDO

    ! Initialize conversion for mixed xi^i * eta^j
    DO j = 1, nOrdPoly
      DO i = 0, nOrdPoly-j 
         DO l = 0, i+j-1
            XiEta2Chi(l  ,i,j)  = XiEta2Chi(l  ,i,j)  +  XiEta2Chi(l,i,j-1) * eta0
            XiEta2Chi(l+1,i,j)  = XiEta2Chi(l+1,i,j)  +  XiEta2Chi(l,i,j-1) * dir_eta/Length
         ENDDO
      ENDDO
    ENDDO

    ! Convert Polynomial by using the conversion tensor i.e. by using 
    ! the mapping from xi, eta exponents to x, y exponents.
    Boundary_Basis(:,:) = 0.

    DO iDegFr = 1, nDegFr
       DO l = 0, nOrdPoly
          DO j = 0, nOrdPoly
             DO i = 0, nOrdPoly-j 
                Boundary_Basis(l,iDegFr) = Boundary_Basis(l,iDegFr) + cPoly(i,j,iDegFr-1)*XiEta2Chi(l,i,j)
             ENDDO
          ENDDO
       ENDDO
    ENDDO

  END SUBROUTINE XYBoundaryPolynomial



END MODULE ConvertXiEta2XY_mod
