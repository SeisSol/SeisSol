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

MODULE TrilinearInterpolation_mod
    !************************************************************************
    IMPLICIT NONE
    !************************************************************************
    INTERFACE TrilinearFromRaster
        MODULE PROCEDURE TrilinearFromRaster
    END INTERFACE

    INTERFACE EvalTrilinearPolynomia
        MODULE PROCEDURE EvalTrilinearPolynomia
    END INTERFACE
    INTERFACE TrilinearInterpolationPolynomia
        MODULE PROCEDURE TrilinearInterpolationPolynomia
    END INTERFACE

    !************************************************************************

    PUBLIC :: TrilinearFromRaster
    
    PRIVATE :: EvalTrilinearPolynomia
    PRIVATE :: TrilinearInterpolationPolynomia

CONTAINS

!****************************************************************************
! Evaluate the Trilinear interpolation for a set of points from a raster array
! P(iPoint)     iPoint from SetPoints
!
! This subroutine assume that the raster contain the set of points
! If is unable to evaluate the plynomia assign 0
! 
! Raster is a 3d array regular spaced in x (dx), y (dy) and z (dz)
!  with first index for x direction, second index for y and third index for z direction
! The lower left coordinate is (x0,y0,z0)
!
! If present, dP has spatial derivatives
! x- dP(:,1)   y- dP(:,2)   z- dP(:,3)
    SUBROUTINE TrilinearFromRaster(P,Set_x,Set_y,Set_z,Raster,dx,dy,dz,x0,y0,z0,dP)
        !-------------------------------------------------------------------!
        REAL                        :: P(:)                                 ! Interpolated value
        REAL                        :: Set_x(:)                             ! Interpolated position x
        REAL                        :: Set_y(:)                             ! Interpolated position y
        REAL                        :: Set_z(:)                             ! Interpolated position z
        REAL                        :: Raster(:,:,:)                        ! Raster 3 rank matrix
        REAL                        :: dx, dy, dz                           ! Delta in the raster
        REAL                        :: x0,y0, z0                            ! Origin of the raster
        REAL                        :: dP(:,:)                              ! Derivatives of interpolated value
        !-------------------------------------------------------------------!
        INTENT(OUT)                 :: P
        INTENT(IN)                  :: Set_x, Set_y, Set_z
        INTENT(IN)                  :: Raster
        INTENT(IN)                  :: dx,dy, dz
        INTENT(IN)                  :: x0,y0,z0
        !-------------------------------------------------------------------!
        OPTIONAL                    :: dP
        !-------------------------------------------------------------------!
        ! Local variables
        REAL                        :: LocalRaster(1:2,1:2,1:2)             ! Sub set of the raster conaining 2x2x2 elements
                                                                            ! is used to evaluate derivatives on vertexes
        REAL                        :: A(0:1,0:1,0:1)                       ! Coeficients of the trilinear polynomia
        
        REAL                        :: x,y,z                                ! Data point to evaluate
        REAL                        :: x1,y1,z1                             ! Actual point to evaluate

        REAL                        :: F   (0:1,0:1,0:1)                    ! Values of function on vertexes

        INTEGER                     :: ix, iy, iz                           ! Index in the raster array
        INTEGER                     :: ix_old, iy_old, iz_old               ! previous index in the raster array
        INTEGER                     :: nx, ny, nz                           ! Number of cells in the raster
        INTEGER                     :: nPoints                              ! Number of points to evaluate
        INTEGER                     :: iPoint                               ! Actual point to evaluate
        INTEGER                     :: fi, fj, fk, li, lj, lk               !
        !-------------------------------------------------------------------!

        ! Sub set of the raster
        LocalRaster = 0.0d0

        ! Number of points where to evaluate P
        nPoints = SIZE(Set_x,DIM=1)

        ! Dimension of the raster
        nx = SIZE(Raster,DIM=1)
        ny = SIZE(Raster,DIM=2)
        nz = SIZE(Raster,DIM=3)

        ix_old = huge(ix_old); iy_old = huge(iy_old); iz_old = huge(iz_old);

        DO iPoint=1,nPoints
            x = Set_x(iPoint);
            y = Set_y(iPoint);
            z = Set_z(iPoint);

            ! Find the cell
            ! I am between raster point [ix:ix+1,iy:iy+1,iz:iz+1]
            ix = INT((x-x0)/dx)+1
            iy = INT((y-y0)/dy)+1
            iz = INT((z-z0)/dz)+1

            ! Evaluation point for the trilinear polynomia
            ! Note is scaled to the unitary square
            x1 = (x-(x0+(ix-1)*dx))/dx
            y1 = (y-(y0+(iy-1)*dy))/dy
            z1 = (z-(z0+(iz-1)*dz))/dz

            IF (ix.EQ.ix_old .AND. iy.EQ.iy_old .AND. iz.EQ.iz_old) THEN
            
                ! If the new points belong to the same raster cell
                ! Eval polynomia
                IF (PRESENT(dP)) THEN
                    CALL EvalTrilinearPolynomia(P(iPoint),x1,y1,z1,A,dpx=dP(iPoint,1),dpy=dP(iPoint,2),dpz=dP(iPoint,3))
                    dP(iPoint,1) = dP(iPoint,1)/dx
                    dP(iPoint,2) = dP(iPoint,2)/dy
                    dP(iPoint,3) = dP(iPoint,3)/dz
                ELSE
                    CALL EvalTrilinearPolynomia(P(iPoint),x1,y1,z1,A)
                ENDIF

            ELSE

                do fi=0,1
                  do fj=0,1
                    do fk=0,1
                      li = min(max(ix+fi,1),nx)
                      lj = min(max(iy+fj,1),ny)
                      lk = min(max(iz+fk,1),nz)
                      F(fi,fj,fk) = Raster(li,lj,lk)
                    enddo
                  enddo
                enddo
                
                ! Eval coeficients
                CALL TrilinearInterpolationPolynomia(A,F)
                ! Eval polynomia
                IF (PRESENT(dP)) THEN
                    CALL EvalTrilinearPolynomia(P(iPoint),x1,y1,z1,A,dpx=dP(iPoint,1),dpy=dP(iPoint,2),dpz=dP(iPoint,3))
                    dP(iPoint,1) = dP(iPoint,1)/dx
                    dP(iPoint,2) = dP(iPoint,2)/dy
                    dP(iPoint,3) = dP(iPoint,3)/dz
                ELSE
                    CALL EvalTrilinearPolynomia(P(iPoint),x1,y1,z1,A)
                ENDIF

                ix_old = ix; iy_old = iy; iz_old = iz;

            ENDIF

        ENDDO ! iPoint

END SUBROUTINE TrilinearFromRaster

!****************************************************************************
! Eval a Trilinear polynomia from coefficients A
SUBROUTINE EvalTrilinearPolynomia(p,x,y,z,A,dpx,dpy,dpz)
        !-------------------------------------------------------------------!
        REAL                        :: p                                    ! P(x,y,z)
        REAL                        :: x,y,z                                ! Evaluation coordinates
        REAL                        :: A(0:1,0:1,0:1)                       ! Coeficients of the polynomia
        REAL                        :: dpx, dpy, dpz                        ! If present, spatial derivatives
        !-------------------------------------------------------------------!
        INTENT(OUT)                 :: p
        INTENT(IN)                  :: x,y,z
        INTENT(IN)                  :: A
        !-------------------------------------------------------------------!
        OPTIONAL                    :: dpx,dpy,dpz
        !-------------------------------------------------------------------!
        ! Local variables
        INTEGER                     :: i
        INTEGER                     :: j
        INTEGER                     :: k
        !-------------------------------------------------------------------!

        p = 0.0d0
        DO i=0,1
            DO j=0,1
                DO k=0,1
                    p = p + A(i,j,k)*(x**i)*(y**j)*(z**k)
                ENDDO ! k
            ENDDO ! j
        ENDDO ! i

        IF (PRESENT(dpx)) THEN
            ! x-derivative
            dpx = 0.0d0
            DO i=1,1
                DO j=0,1
                    DO k=0,1
                        dpx = dpx + A(i,j,k)*i*(x**(i-1))*(y**(j))*(z**(k))
                    ENDDO ! k
                ENDDO ! j
            ENDDO ! i
        ENDIF
        IF (PRESENT(dpy)) THEN
            ! y-derivative
            dpy = 0.0d0
            DO i=0,1
                DO j=1,1
                    DO k=0,1
                        dpy = dpy + A(i,j,k)*(x**(i))*j*(y**(j-1))*(z**(k))
                    ENDDO ! k
                ENDDO ! j
            ENDDO ! i
        ENDIF
        IF (PRESENT(dpz)) THEN
            ! z-derivative
            dpz = 0.0d0
            DO i=0,1
                DO j=0,1
                    DO k=1,1
                        dpz = dpz + A(i,j,k)*(x**(i))*(y**(j))*k*(z**(k-1))
                    ENDDO ! k
                ENDDO ! j
            ENDDO ! i
        ENDIF
        CONTINUE

END SUBROUTINE EvalTriLinearPolynomia
!****************************************************************************

!****************************************************************************
! Compute the Trilinear interpolation polynomia
! P(x,y,z) = SUM_i=0..1  SUM_j=0..1  SUM_k=0..1  A[i,j,k] x^i y^j z^k
! 
! Vertexes of unitary square [0,1]x[0,1]x[0,1]
! x=0,y=0,z=0 -> [0,0,0]
! x=1,y=0,z=0 -> [1,0,0]
! x=1,y=1,z=0 -> [1,1,0]
! x=0,y=1,z=0 -> [0,1,0]
! 
! x=0,y=0,z=1 -> [0,0,1]
! x=1,y=0,z=1 -> [1,0,1]
! x=1,y=1,z=1 -> [1,1,1]
! x=0,y=1,z=1 -> [0,1,1]
!
! F     values of the function P at vertexes
!
SUBROUTINE TrilinearInterpolationPolynomia(A,F)
        !-------------------------------------------------------------------!
        IMPLICIT NONE
        !-------------------------------------------------------------------!
        REAL                        :: A   (0:1,0:1,0:1)
        REAL                        :: F   (0:1,0:1,0:1)
        !-------------------------------------------------------------------!
        INTENT(OUT)                 :: A
        INTENT(IN)                  :: F
        !-------------------------------------------------------------------!
        ! Local variables
        !-------------------------------------------------------------------!

        A(0,0,0) =  F(0,0,0)
        A(1,0,0) = -F(0,0,0)+F(1,0,0)
        A(0,1,0) = -F(0,0,0)+F(0,1,0)
        A(0,0,1) = -F(0,0,0)+F(0,0,1)
        A(1,1,0) =  F(0,0,0)-F(0,1,0)-F(1,0,0)+F(1,1,0)
        A(1,0,1) =  F(0,0,0)-F(0,0,1)-F(1,0,0)+F(1,0,1)
        A(0,1,1) =  F(0,0,0)-F(0,0,1)-F(0,1,0)+F(0,1,1)
        A(1,1,1) = -F(0,0,0)+F(1,0,0)+F(0,1,0)+F(0,0,1)-F(1,1,0)-F(1,0,1)-F(0,1,1)+F(1,1,1)

END SUBROUTINE TrilinearInterpolationPolynomia
!****************************************************************************

END MODULE TrilinearInterpolation_mod
