!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Josep de la Puente Alvarez (josep.delapuente AT bsc.es, http://www.geophysik.uni-muenchen.de/Members/jdelapuente)
!!
!! @section LICENSE
!! Copyright (c) 2008, SeisSol Group
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

MODULE QuadPoints_mod
  !-------------------------------------------------------------------------!
  USE TypesDef
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE TriangleQuadraturePoints
     MODULE PROCEDURE TriangleQuadraturePoints
  END INTERFACE

  INTERFACE TetrahedronQuadraturePoints
     MODULE PROCEDURE TetrahedronQuadraturePoints
  END INTERFACE

  INTERFACE HypercubeQuadraturePoints
     MODULE PROCEDURE HypercubeQuadraturePoints
  END INTERFACE

  INTERFACE TriangleSubPoints
     MODULE PROCEDURE TriangleSubPoints
  END INTERFACE
  
  interface CellCentresOfSubdivision
    module procedure CellCentresOfSubdivision
  end interface
  !----------------------------------------------------------------------------
  PUBLIC  :: TetrahedronQuadraturePoints
  PUBLIC  :: TriangleQuadraturePoints
  PUBLIC  :: HypercubeQuadraturePoints
  PUBLIC  :: TriangleSubPoints
  public  :: CellCentresOfSubdivision
  !----------------------------------------------------------------------------

CONTAINS

  ! Quadrature formula of arbitrary accuracy on the reference triangle
  ! consisting of the nodes (0,0), (1,0), (0,1)
  SUBROUTINE TriangleQuadraturePoints(nIntGP,IntGaussP,IntGaussW,M,IO,quiet)
    !-------------------------------------------------------------------------!
    USE gauss_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER               :: nIntGP         ! Number of 2D integration points !
    REAL, POINTER         :: IntGaussP(:,:) ! Positions of 2D int. points     !
    REAL, POINTER         :: IntGaussW(:)   ! Weights of 2D int. points       !
    INTEGER               :: M              ! Number of 1D Gausspoints        !
    TYPE(tInputOutput)    :: IO             ! IO structure                    !
    LOGICAL, OPTIONAL     :: quiet          ! if quiet then less warnings     !
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER         :: i,j                  ! Loop counters                   !
    INTEGER         :: iIntGP               ! Loop counter                    !
    INTEGER         :: allocstat            ! Allocation status return int.   !
    REAL            :: tol                  ! Tolerance of 0.0                !
    REAL, POINTER   :: mu1(:)               ! 1D quadrature positions in y1   !
    REAL, POINTER   :: mu2(:)               ! 1D quadrature positions in y2   !
    REAL, POINTER   :: A1(:)                ! 1D quadrature weights for y1    !
    REAL, POINTER   :: A2(:)                ! 1D quadrature weights for y2    !
    LOGICAL         :: quiet_internal       ! map of OPTIONAL quiet           !
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: M, IO, quiet
    INTENT(OUT)     :: nIntGP
    !-------------------------------------------------------------------------!
    tol = 1.0 / (10.0**(PRECISION(1.0)-2) )                                   !
    !-------------------------------------------------------------------------!
    !                                                                         !
    IF (.NOT.PRESENT(quiet)) THEN
       quiet_Internal=.FALSE.
    ELSE
       quiet_Internal=quiet
    END IF
    ! Quadrature points are defined by the conical product of 1D Gauss-Jacobi 
    ! formulas with M quadrature points. See Stround, p. 28ff for details.
    nIntGP = M*M

    IF(.NOT.ASSOCIATED(IntGaussP)) THEN
       IF (.NOT.quiet_Internal) THEN
          logWarning(*) 'Warning: IntGaussP not allocated before call of'
          logWarning(*) '         TriangleQuadraturePoints.'
          logWarning(*) '         Allocating now.'
       END IF
       ALLOCATE(IntGaussP(2,nIntGP),STAT=allocstat)
       IF (allocStat .NE. 0) THEN
          logError(*) 'TriangleQuadraturePoints: could not allocate all variables!'
          STOP
       END IF
    ENDIF
    IF(.NOT.ASSOCIATED(IntGaussW)) THEN
       IF (.NOT.quiet_Internal) THEN     
          logWarning(*) 'Warning: IntGaussW not allocated before call of'
          logWarning(*) '         TriangleQuadraturePoints.'
          logWarning(*) '         Allocating now.'
       END IF
       ALLOCATE(IntGaussW(nIntGP),STAT=allocstat)
       IF (allocStat .NE. 0) THEN
          logError(*) 'TriangleQuadraturePoints: could not allocate all variables!'
          STOP
       END IF
    ENDIF

    ALLOCATE(mu1(M), mu2(M), A1(M), A2(M), STAT = allocstat)
    IF (allocStat .NE. 0) THEN
       logError(*) 'TriangleQuadraturePoints: could not allocate all variables!'
       STOP
    END IF

    !CALL gaujac(mu1,A1,M,1.,0.)     ! Get the Gauss-Jacobi positions and weights
    !CALL gaujac(mu2,A2,M,0.,0.)     ! Get the Gauss-Jacobi positions and weights
    CALL gauss_jacobi(1.,0.,mu1,A1)     ! Get the Gauss-Jacobi positions and weights
    CALL gauss_jacobi(0.,0.,mu2,A2)     ! Get the Gauss-Jacobi positions and weights

    mu1(:) = 0.5*mu1(:) + 0.5       ! Shift and rescale positions, because Stroud
    A1(:)  = 0.5**2*A1(:)           ! integrates over the interval [0,1] and
    mu2(:) = 0.5*mu2(:) + 0.5       ! the function gaujac of the num. recipes
    A2(:)  = 0.5**1*A2(:)           ! integrates over the interval [-1,1].

    iIntGP = 1

    DO i = 1, M
       DO j = 1, M
          intGaussP(1,iIntGP) = mu1(i)
          intGaussP(2,iIntGP) = mu2(j)*(1.-mu1(i))
          intGaussW(iIntGP)   = A1(i)*A2(j)      
          iIntGP              = iIntGP + 1
       ENDDO
    ENDDO

    DEALLOCATE(mu1, mu2, A1, A2)

    IF (     ((ABS(SUM(intGaussW(:)))-0.5).GT.tol) &
         .OR.(.NOT.quiet_internal)               ) THEN
       logInfo(*) 'Integration points in TRIANGLE calculated with conical product.'
       logInfo(*) 'Number of integration points is ', nIntGP
       logInfo(*) 'SUM of Weights is ', SUM(intGaussW(:)), ' and must be 0.5! '
    END IF

  END SUBROUTINE TriangleQuadraturePoints

  ! Quadrature formula of arbitrary accuracy on the reference tetrahedron
  ! consisting of the nodes (0,0,0), (1,0,0), (0,1,0), (0,0,1)

  SUBROUTINE TetrahedronQuadraturePoints(nIntGP,IntGaussP,IntGaussW,M,IO,quiet)
    !-------------------------------------------------------------------------!
    USE gauss_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER               :: nIntGP         ! Number of 2D integration points !
    REAL, POINTER         :: IntGaussP(:,:) ! Positions of 2D int. points     !
    REAL, POINTER         :: IntGaussW(:)   ! Weights of 2D int. points       !
    INTEGER               :: M              ! Number of 1D Gausspoints        !
    TYPE(tInputOutput)    :: IO             ! IO structure                    !
    LOGICAL, OPTIONAL     :: quiet          ! if quiet then less warnings     !
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER         :: i,j,k                ! Loop counters                   !
    INTEGER         :: iIntGP               ! Loop counter                    !
    INTEGER         :: allocstat            ! Allocation status return int.   !
    REAL            :: tol                  ! Tolerance of 0.0                !
    REAL, POINTER   :: mu1(:)               ! 1D quadrature positions in y1   !
    REAL, POINTER   :: mu2(:)               ! 1D quadrature positions in y2   !
    REAL, POINTER   :: mu3(:)               ! 1D quadrature positions in y3   !
    REAL, POINTER   :: A1(:)                ! 1D quadrature weights for y1    !
    REAL, POINTER   :: A2(:)                ! 1D quadrature weights for y2    !
    REAL, POINTER   :: A3(:)                ! 1D quadrature weights for y3    !
    LOGICAL         :: quiet_internal       ! map of OPTIONAL quiet           !
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: M, IO, quiet
    INTENT(OUT)     :: nIntGP
    !-------------------------------------------------------------------------!
    tol = 1.0 / (10.0**(PRECISION(1.0)-2) )                                   !
    !-------------------------------------------------------------------------!
    !                                                                         !
    IF (.NOT.PRESENT(quiet)) THEN
       quiet_Internal=.FALSE.
    ELSE
       quiet_Internal=quiet
    END IF
    ! Quadrature points are defined by the conical product of 1D Gauss-Jacobi 
    ! formulas with M quadrature points. See Stround, p. 28ff for details.
    nIntGP = M*M*M

    IF(.NOT.ASSOCIATED(IntGaussP)) THEN
       IF (.NOT.quiet_Internal) THEN
          logWarning(*) 'IntGaussP not allocated before call of TriangleQuadraturePoints. Allocating now.'
       END IF
       ALLOCATE(IntGaussP(3,nIntGP),STAT=allocstat)
       IF (allocStat .NE. 0) THEN
          logError(*) 'TriangleQuadraturePoints: could not allocate all variables!'
          STOP
       END IF
    ENDIF
    IF(.NOT.ASSOCIATED(IntGaussW)) THEN
       IF (.NOT.quiet_Internal) THEN     
          logWarning(*) 'IntGaussW not allocated before call of TriangleQuadraturePoints. Allocating now.'
       END IF
       ALLOCATE(IntGaussW(nIntGP),STAT=allocstat)
       IF (allocStat .NE. 0) THEN
          logError(*) 'TriangleQuadraturePoints: could not allocate all variables!'
          STOP
       END IF
    ENDIF

    ALLOCATE(mu1(M), mu2(M), mu3(M), A1(M), A2(M), A3(M), STAT = allocstat)
    IF (allocStat .NE. 0) THEN
       logError(*) 'TriangleQuadraturePoints: could not allocate all variables!'
       STOP
    END IF

    !CALL gaujac(mu1,A1,M,2.,0.)     ! Get the Gauss-Jacobi positions and weights
    !CALL gaujac(mu2,A2,M,1.,0.)     ! Get the Gauss-Jacobi positions and weights
    !CALL gaujac(mu3,A3,M,0.,0.)     ! Get the Gauss-Jacobi positions and weights
    CALL gauss_jacobi(2.,0.,mu1,A1)     ! Get the Gauss-Jacobi positions and weights
    CALL gauss_jacobi(1.,0.,mu2,A2)     ! Get the Gauss-Jacobi positions and weights
    CALL gauss_jacobi(0.,0.,mu3,A3)     ! Get the Gauss-Jacobi positions and weights

    mu1(:) = 0.5*mu1(:) + 0.5       ! Shift and rescale positions, because Stroud
    A1(:)  = 0.5**3*A1(:)           ! integrates over the interval [0,1] and
    mu2(:) = 0.5*mu2(:) + 0.5       ! the function gaujac of the num. recipes
    A2(:)  = 0.5**2*A2(:)           ! integrates over the interval [-1,1].
    mu3(:) = 0.5*mu3(:) + 0.5       ! the function gaujac of the num. recipes
    A3(:)  = 0.5**1*A3(:)           ! integrates over the interval [-1,1].

    iIntGP = 1

    DO i = 1, M
       DO j = 1, M
         DO k = 1, M
              intGaussP(1,iIntGP) = mu1(i)
              intGaussP(2,iIntGP) = mu2(j)*(1.-mu1(i))
              intGaussP(3,iIntGP) = mu3(k)*(1.-mu2(j))*(1.-mu1(i))
              intGaussW(iIntGP)   = A1(i)*A2(j)*A3(k)      
              iIntGP              = iIntGP + 1
         ENDDO
       ENDDO
    ENDDO

    DEALLOCATE(mu1, mu2, mu3, A1, A2, A3)

    IF (     ((ABS(SUM(intGaussW(:)))-1./6.).GT.tol) &
         .OR.(.NOT.quiet_internal)               ) THEN
       logInfo(*) 'Integration points in TETRAHEDRON calculated with conical product.'
       logInfo(*) 'Number of integration points is  ', nIntGP
       logInfo(*) 'SUM of Weights is ', SUM(intGaussW(:)), ' and must be 0.1666666...! '
    END IF

  END SUBROUTINE TetrahedronQuadraturePoints

  ! Quadrature formula of arbitrary accuracy on the physical hypercube 
  ! [x1_1, x2_1] x [x1_2, x2_2] x ... x [x1_n, x2_n]
  SUBROUTINE HypercubeQuadraturePoints(nIntGP,IntGaussP,IntGaussW, M,nDim,X1,X2,IO,quiet )
    !-------------------------------------------------------------------------!
    use gauss_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER               :: nIntGP         ! Number of nD integration points !
    REAL, POINTER         :: IntGaussP(:,:) ! Positions of nD int. points     !
    REAL, POINTER         :: IntGaussW(:)   ! Weights of nD int. points       !
    INTEGER               :: M              ! Number of 1D Gausspoints        !
    INTEGER               :: nDim           ! Number of dimensions            !
    REAL                  :: X1(nDim)       ! Physical hypercube corner 1     !
    REAL                  :: X2(nDim)       ! Physical hypercube corner 2     !
    REAL                  :: tol            ! Tolerance of 0.0                !
    TYPE(tInputOutput)    :: IO             ! IO structure                    !
    LOGICAL, OPTIONAL     :: quiet          ! if quiet then less warnings     !
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER         :: iDim                 ! Loop counters                   !
    INTEGER         :: iIntGP               ! Loop counter                    !
    INTEGER         :: allocstat            ! Allocation status return int.   !
    INTEGER         :: remain               ! Remainder (aux. variable)       !
    INTEGER         :: multiindex(1:nDim)   ! n-dimensional multi-index       !
    REAL, POINTER   :: mu(:,:)              ! 1D quadrature positions         !
    REAL, POINTER   :: A(:,:)               ! 1D quadrature weights           !
    LOGICAL         :: quiet_internal       ! map of OPTIONAL quiet           !
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: M, X1, X2, IO, quiet
    INTENT(OUT)     :: nIntGP
    !-------------------------------------------------------------------------!
    tol = 1.0 / (10.0**(PRECISION(1.0)-2) )                                   !
    !-------------------------------------------------------------------------!
    !                                                                         !
    IF (.NOT.PRESENT(quiet)) THEN                                             !
       quiet_Internal=.FALSE.                                                 !
    ELSE                                                                      !
       quiet_Internal=quiet                                                   !
    END IF                                                                    !
    !                                                                         !
    ! Quadrature points are defined by the tensor product of 1D Gauss-Legendre! 
    ! formulas with M quadrature points on a hypercube in physical space      !
    !                                                                         !
    nIntGP = M**nDim                                                          !
    !                                                                         !    
    IF(.NOT.ASSOCIATED(IntGaussP)) THEN                                       !
       IF (.NOT.quiet_Internal) THEN                                          !
          logWarning(*)                                 &          !
               'IntGaussP not allocated before the ',     &          !
               'call of CubeQuadraturePoints. Allocating now.'                !
       END IF                                                                 !
       ALLOCATE(IntGaussP(nDim,nIntGP),STAT=allocstat)                        !
       IF (allocStat .NE. 0) THEN                                             !
          logError(*)                             &               !
               'Error in HypercubeQuadraturePoints: ',        &               !
               'could not allocate all variables!'                            !
          STOP                                                                !
       END IF                                                                 !
    ENDIF                                                                     !
    IF(.NOT.ASSOCIATED(IntGaussW)) THEN                                       !
       IF (.NOT.quiet_Internal) THEN                                          !
          logWarning(*)                               &            !
               'IntGaussW not allocated before the ',   &            !
               'call of CubeQuadraturePoints. Allocating now.'                !
       END IF                                                                 !
       ALLOCATE(IntGaussW(nIntGP),STAT=allocstat)                             !
       IF (allocStat .NE. 0) THEN                                             !
          logError(*)                             &               !
               'Error in HypercubeQuadraturePoints: ',        &               !
               'could not allocate all variables!'                            !
          STOP                                                                !
       END IF                                                                 !
    ENDIF                                                                     !
    !                                                                         !
    ALLOCATE(mu(M,nDim), A(M,nDim), STAT = allocstat)                         !
    !                                                                         !
    IF (allocStat .NE. 0) THEN                                                !
       logError(*)                               &                !
            'Error in HypercubeQuadraturePoints: ',          &                !
            'could not allocate all variables!'                               !
       STOP                                                                   !
    END IF                                                                    !
    !                                                                         !
    DO iDim = 1, nDim                                                         !
       CALL gausslegendre(X1(iDim), X2(iDim),mu(:,iDim),A(:,iDim))            !
    ENDDO                                                                     !
    !                                                                         !
    iIntGP = 1                                                                !
    !                                                                         !
    DO iIntGP = 1, nIntGP                                                     !
       ! Compute the multiindex inside the hypercube                          !
       remain = iIntGP                                                        !
       DO iDim = nDim, 1, -1                                                  !
          multiindex(iDim) = (remain-1)/M**(iDim-1)+1                         !
          remain           = MOD((remain-1),M**(iDim-1))+1                    !
       ENDDO                                                                  !
       IntGaussW(iIntGP) = 1.                                                 !
       DO iDim = 1, nDim                                                      !
          IntGaussP(iDim,iIntGP) = mu(multiindex(iDim),iDim)                  !
          IntGaussW(iIntGP) = IntGaussW(iIntGP)*A(multiindex(iDim),iDim)      !
       ENDDO                                                                  !
    ENDDO                                                                     !
    !                                                                         !
    DEALLOCATE(mu, A)                                                         !
    !                                                                         !
    IF (ABS((SUM(IntGaussW(:))-PRODUCT(X2(:)-X1(:)))/SUM(IntGaussW(:))).GT.tol &
         .OR. .NOT.quiet_Internal                                            ) THEN
       logInfo(*) ' Quadrature formula using tensor product '   !
       logInfo(*) ' for the hypercube of dimension  ', nDim     !
       logInfo(*) ' Number of integration points is  ', nIntGP  !
       logInfo(*) ' SUM of Weights is ', SUM(IntGaussW(:))  , & !
            ' and must be ', PRODUCT( X2(:)-X1(:) )                           !
    END IF                                                                    !
    !                                                                         !
  END SUBROUTINE HypercubeQuadraturePoints



  ! Subroutine that subdivides a triangle with vertex coordinates X(:) and Y(:)
  ! regularly. Resulting number of points is 4^M
  SUBROUTINE TriangleSubPoints(nPoints,Points,M,IO)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    TYPE(tInputOutput)    :: IO             ! IO structure                    !
    INTEGER               :: nPoints        ! Number of 2D points             !
    REAL, POINTER         :: Points(:,:)    ! Positions of 2D points          !
    INTEGER               :: M              ! Number of recursions            !
    !-------------------------------------------------------------------------!
    INTEGER         :: counter              ! Loop counter                    !
    INTEGER         :: recdepth             ! Recursion counter
    INTEGER         :: allocstat            ! Allocation status return int.   !
    REAL            :: XY(3,2)
    !-------------------------------------------------------------------------!
    INTENT(IN)      :: M, IO
    INTENT(OUT)     :: nPoints
    !-------------------------------------------------------------------------!
    !                                                                         !
    nPoints = 4**M

    IF(.NOT.ASSOCIATED(Points)) THEN
       ALLOCATE(Points(nPoints,2),STAT=allocstat)
       IF (allocStat .NE. 0) THEN
          logError(*) 'TriangleQuadraturePoints: could not allocate all variables!'
          STOP
       END IF
    ENDIF

    XY(1,:) = (/0.,0./)
    XY(2,:) = (/1.,0./)
    XY(3,:) = (/0.,1./)
    !
    recdepth = 0
    counter  = 0
    !
    CALL Divide(Points, counter, recdepth, M, nPoints, XY)
    !
  END SUBROUTINE TriangleSubPoints

  RECURSIVE SUBROUTINE Divide(Points, counter, recdepth, M, nPoints, XY)
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    INTEGER               :: counter,recdepth, M, nPoints
    REAL                  :: XY(3,2) 
    REAL                  :: Points(nPoints,2)
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER         :: i
    REAL            :: MidSide(3,2) 
    REAL            :: NewXY(4,3,2)
    !-------------------------------------------------------------------------!
    IF(recdepth.EQ.M) THEN
       counter = counter  + 1
       IF(counter.GT.nPoints) THEN
          logError(*) 'Error in Divide. counter > nPoints', counter
          STOP
       ENDIF
       Points(counter,:) = 1./3.*(XY(1,:)+XY(2,:)+XY(3,:))
       RETURN
    ENDIF
    MidSide(1,:) = 0.5*( XY(1,:) + XY(2,:) )
    MidSide(2,:) = 0.5*( XY(2,:) + XY(3,:) )
    MidSide(3,:) = 0.5*( XY(3,:) + XY(1,:) )
    
    ! New triangle I
    NewXY(1,1,:) = XY(1,:)
    NewXY(1,2,:) = MidSide(1,:)
    NewXY(1,3,:) = MidSide(3,:)

    ! New triangle II
    NewXY(2,1,:) = MidSide(1,:)
    NewXY(2,2,:) = XY(2,:)
    NewXY(2,3,:) = MidSide(2,:)

    ! New triangle III
    NewXY(3,1,:) = MidSide(2,:)
    NewXY(3,2,:) = XY(3,:)
    NewXY(3,3,:) = MidSide(3,:)

    ! New triangle IV
    NewXY(4,1,:) = MidSide(1,:)
    NewXY(4,2,:) = MidSide(2,:)
    NewXY(4,3,:) = MidSide(3,:)

    DO i = 1, 4
       CALL Divide(Points,counter,recdepth+1,M,nPoints,NewXY(i,:,:))
    ENDDO
        
  END SUBROUTINE Divide
  
  subroutine CellCentresOfSubdivision(order,nodes)
    implicit none
    integer               :: order
    real, pointer         :: nodes(:,:)
    intent(in)            :: order
    intent(out)           :: nodes
    
    select case(order)
      case (2)
        nodes = reshape((/0.16666666666666665741,0.16666666666666665741,0.66666666666666662966,0.16666666666666665741,0.16666666666666665741,0.66666666666666662966,0.33333333333333331483,0.33333333333333331483/), (/ 2, size(nodes,2) /))
      case (3)
        nodes = reshape((/0.08333333333333332871,0.08333333333333332871,0.33333333333333331483,0.08333333333333332871,0.08333333333333332871,0.33333333333333331483,0.16666666666666665741,0.16666666666666665741,0.58333333333333337034,0.08333333333333332871,0.83333333333333337034,0.08333333333333332871,0.58333333333333337034,0.33333333333333331483,0.66666666666666662966,0.16666666666666665741,0.08333333333333332871,0.58333333333333337034,0.33333333333333331483,0.58333333333333337034,0.08333333333333332871,0.83333333333333337034,0.16666666666666665741,0.66666666666666662966,0.41666666666666668517,0.41666666666666668517,0.16666666666666665741,0.41666666666666668517,0.41666666666666668517,0.16666666666666665741,0.33333333333333331483,0.33333333333333331483/), (/ 2, size(nodes,2) /))
      case (4)
        nodes = reshape((/0.08333333333333332871,0.08333333333333332871,0.33333333333333331483,0.08333333333333332871,0.08333333333333332871,0.33333333333333331483,0.16666666666666665741,0.16666666666666665741,0.58333333333333337034,0.08333333333333332871,0.83333333333333337034,0.08333333333333332871,0.58333333333333337034,0.33333333333333331483,0.66666666666666662966,0.16666666666666665741,0.08333333333333332871,0.58333333333333337034,0.33333333333333331483,0.58333333333333337034,0.08333333333333332871,0.83333333333333337034,0.16666666666666665741,0.66666666666666662966,0.41666666666666668517,0.41666666666666668517,0.16666666666666665741,0.41666666666666668517,0.41666666666666668517,0.16666666666666665741,0.33333333333333331483,0.33333333333333331483/), (/ 2, size(nodes,2) /))
      case (5)
        nodes = reshape((/0.08333333333333332871,0.08333333333333332871,0.33333333333333331483,0.08333333333333332871,0.08333333333333332871,0.33333333333333331483,0.16666666666666665741,0.16666666666666665741,0.58333333333333337034,0.08333333333333332871,0.83333333333333337034,0.08333333333333332871,0.58333333333333337034,0.33333333333333331483,0.66666666666666662966,0.16666666666666665741,0.08333333333333332871,0.58333333333333337034,0.33333333333333331483,0.58333333333333337034,0.08333333333333332871,0.83333333333333337034,0.16666666666666665741,0.66666666666666662966,0.41666666666666668517,0.41666666666666668517,0.16666666666666665741,0.41666666666666668517,0.41666666666666668517,0.16666666666666665741,0.33333333333333331483,0.33333333333333331483/), (/ 2, size(nodes,2) /))
      case (6)
        nodes = reshape((/0.04166666666666666435,0.04166666666666666435,0.16666666666666665741,0.04166666666666666435,0.04166666666666666435,0.16666666666666665741,0.08333333333333332871,0.08333333333333332871,0.29166666666666668517,0.04166666666666666435,0.41666666666666668517,0.04166666666666666435,0.29166666666666668517,0.16666666666666665741,0.33333333333333331483,0.08333333333333332871,0.04166666666666666435,0.29166666666666668517,0.16666666666666665741,0.29166666666666668517,0.04166666666666666435,0.41666666666666668517,0.08333333333333332871,0.33333333333333331483,0.20833333333333334259,0.20833333333333334259,0.08333333333333332871,0.20833333333333334259,0.20833333333333334259,0.08333333333333332871,0.16666666666666665741,0.16666666666666665741,0.54166666666666662966,0.04166666666666666435,0.66666666666666662966,0.04166666666666666435,0.54166666666666662966,0.16666666666666665741,0.58333333333333337034,0.08333333333333332871,0.79166666666666662966,0.04166666666666666435,0.91666666666666662966,0.04166666666666666435,0.79166666666666662966,0.16666666666666665741,0.83333333333333337034,0.08333333333333332871,0.54166666666666662966,0.29166666666666668517,0.66666666666666662966,0.29166666666666668517,0.54166666666666662966,0.41666666666666668517,0.58333333333333337034,0.33333333333333331483,0.70833333333333337034,0.20833333333333334259,0.58333333333333337034,0.20833333333333334259,0.70833333333333337034,0.08333333333333332871,0.66666666666666662966,0.16666666666666665741,0.04166666666666666435,0.54166666666666662966,0.16666666666666665741,0.54166666666666662966,0.04166666666666666435,0.66666666666666662966,0.08333333333333332871,0.58333333333333337034,0.29166666666666668517,0.54166666666666662966,0.41666666666666668517,0.54166666666666662966,0.29166666666666668517,0.66666666666666662966,0.33333333333333331483,0.58333333333333337034,0.04166666666666666435,0.79166666666666662966,0.16666666666666665741,0.79166666666666662966,0.04166666666666666435,0.91666666666666662966,0.08333333333333332871,0.83333333333333337034,0.20833333333333334259,0.70833333333333337034,0.08333333333333332871,0.70833333333333337034,0.20833333333333334259,0.58333333333333337034,0.16666666666666665741,0.66666666666666662966,0.45833333333333331483,0.45833333333333331483,0.33333333333333331483,0.45833333333333331483,0.45833333333333331483,0.33333333333333331483,0.41666666666666668517,0.41666666666666668517,0.20833333333333334259,0.45833333333333331483,0.08333333333333332871,0.45833333333333331483,0.20833333333333334259,0.33333333333333331483,0.16666666666666665741,0.41666666666666668517,0.45833333333333331483,0.20833333333333334259,0.33333333333333331483,0.20833333333333334259,0.45833333333333331483,0.08333333333333332871,0.41666666666666668517,0.16666666666666665741,0.29166666666666668517,0.29166666666666668517,0.41666666666666668517,0.29166666666666668517,0.29166666666666668517,0.41666666666666668517,0.33333333333333331483,0.33333333333333331483/), (/ 2, size(nodes,2) /))
      case (7)
        nodes = reshape((/0.04166666666666666435,0.04166666666666666435,0.16666666666666665741,0.04166666666666666435,0.04166666666666666435,0.16666666666666665741,0.08333333333333332871,0.08333333333333332871,0.29166666666666668517,0.04166666666666666435,0.41666666666666668517,0.04166666666666666435,0.29166666666666668517,0.16666666666666665741,0.33333333333333331483,0.08333333333333332871,0.04166666666666666435,0.29166666666666668517,0.16666666666666665741,0.29166666666666668517,0.04166666666666666435,0.41666666666666668517,0.08333333333333332871,0.33333333333333331483,0.20833333333333334259,0.20833333333333334259,0.08333333333333332871,0.20833333333333334259,0.20833333333333334259,0.08333333333333332871,0.16666666666666665741,0.16666666666666665741,0.54166666666666662966,0.04166666666666666435,0.66666666666666662966,0.04166666666666666435,0.54166666666666662966,0.16666666666666665741,0.58333333333333337034,0.08333333333333332871,0.79166666666666662966,0.04166666666666666435,0.91666666666666662966,0.04166666666666666435,0.79166666666666662966,0.16666666666666665741,0.83333333333333337034,0.08333333333333332871,0.54166666666666662966,0.29166666666666668517,0.66666666666666662966,0.29166666666666668517,0.54166666666666662966,0.41666666666666668517,0.58333333333333337034,0.33333333333333331483,0.70833333333333337034,0.20833333333333334259,0.58333333333333337034,0.20833333333333334259,0.70833333333333337034,0.08333333333333332871,0.66666666666666662966,0.16666666666666665741,0.04166666666666666435,0.54166666666666662966,0.16666666666666665741,0.54166666666666662966,0.04166666666666666435,0.66666666666666662966,0.08333333333333332871,0.58333333333333337034,0.29166666666666668517,0.54166666666666662966,0.41666666666666668517,0.54166666666666662966,0.29166666666666668517,0.66666666666666662966,0.33333333333333331483,0.58333333333333337034,0.04166666666666666435,0.79166666666666662966,0.16666666666666665741,0.79166666666666662966,0.04166666666666666435,0.91666666666666662966,0.08333333333333332871,0.83333333333333337034,0.20833333333333334259,0.70833333333333337034,0.08333333333333332871,0.70833333333333337034,0.20833333333333334259,0.58333333333333337034,0.16666666666666665741,0.66666666666666662966,0.45833333333333331483,0.45833333333333331483,0.33333333333333331483,0.45833333333333331483,0.45833333333333331483,0.33333333333333331483,0.41666666666666668517,0.41666666666666668517,0.20833333333333334259,0.45833333333333331483,0.08333333333333332871,0.45833333333333331483,0.20833333333333334259,0.33333333333333331483,0.16666666666666665741,0.41666666666666668517,0.45833333333333331483,0.20833333333333334259,0.33333333333333331483,0.20833333333333334259,0.45833333333333331483,0.08333333333333332871,0.41666666666666668517,0.16666666666666665741,0.29166666666666668517,0.29166666666666668517,0.41666666666666668517,0.29166666666666668517,0.29166666666666668517,0.41666666666666668517,0.33333333333333331483,0.33333333333333331483/), (/ 2, size(nodes,2) /))
      case default
        logError(*) 'Unsupported order in CellCentresOfSubdivision'
    end select
  end subroutine CellCentresOfSubdivision

END MODULE QuadPoints_mod


