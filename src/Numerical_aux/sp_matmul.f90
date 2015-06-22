!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Michael Dumbser (michael.dumbser AT unitn.it, https://www5.unitn.it/People/it/Web/Persona/PER0029602#INFO)
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
!!
!! @section DESCRIPTION
!! Sparse matrix-matrix multiplication routines


MODULE SP_MATMUL_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  INTERFACE SPL_MATMUL
     MODULE PROCEDURE SPL_MATMUL
  END INTERFACE
  INTERFACE SPR_MATMUL
     MODULE PROCEDURE SPR_MATMUL
  END INTERFACE
  INTERFACE SPT_MATMUL_P
     MODULE PROCEDURE SPT_MATMUL_P
  END INTERFACE
  INTERFACE SPT_MATMUL_M
     MODULE PROCEDURE SPT_MATMUL_M
  END INTERFACE
  INTERFACE SP_ADD3D_MULT
     MODULE PROCEDURE SP_ADD3D_MULT
  END INTERFACE
  INTERFACE SPT_T3MUL_M
     MODULE PROCEDURE SPT_T3MUL_M
  END INTERFACE
  INTERFACE SPT_T3MUL_P
     MODULE PROCEDURE SPT_T3MUL_P
  END INTERFACE
  INTERFACE SPL_T3MUL
     MODULE PROCEDURE SPL_T3MUL
  END INTERFACE
  INTERFACE SPL_T3MULB
     MODULE PROCEDURE SPL_T3MULB
  END INTERFACE

  ! 
  !---------------------------------------------------------------------------!
  PUBLIC  :: SPL_MATMUL
  PUBLIC  :: SPR_MATMUL
  PUBLIC  :: SPT_MATMUL_P
  PUBLIC  :: SPT_MATMUL_M
  PUBLIC  :: SP_ADD3D_MULT
  PUBLIC  :: SPT_T3MUL_M
  PUBLIC  :: SPT_T3MUL_P
  PUBLIC  :: SPL_T3MUL
  PUBLIC  :: SPL_T3MULB

  !---------------------------------------------------------------------------!

CONTAINS
    !
    !> The subroutine SPL_MATMUL multiplies a sparse matrix A from the left
    !! to a full matrix B and adds the result to the matrix C:
    !! 
    !! C_ik = C_ik + A_ij * B_jk
    !<
      SUBROUTINE SPL_MATMUL(C, leadc, A, B, leadb, n, o)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      INTEGER                   :: n,o, leadb, leadc
      TYPE(tSparseMatrix)       :: A                                            ! (n,n)
      REAL                      :: C(leadc,o)                                   ! (n,o)
      REAL                      :: B(leadb,o)                                   ! (n,o)
!    ! Local variable declaration
      INTEGER                   :: iNonZero,i,j
      INTEGER                   :: kk
      REAL                      :: BT(o,n)
      REAL                      :: CT(o,n)
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: n,o,A,B, leadb, leadc
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      CT = 0.
      BT = TRANSPOSE(B(1:n,:))
      DO iNonZero = 1, A%nNonZero(n)
        i = A%NonZeroIndex1(iNonZero)
        j = A%NonZeroIndex2(iNonZero)
        CT(:,i) = CT(:,i) + BT(:,j)*A%NonZero(iNonZero) 
      ENDDO
      C(1:n,:) = C(1:n,:) + TRANSPOSE(CT)
!    !
      END SUBROUTINE SPL_MATMUL

    !
    !> The subroutine SPR_MATMUL multiplies a sparse matrix B from the right
    !! to a full matrix A and adds the result to the matrix C
    !< C_ik = C_ik + A_ij * B_jk
    !
      PURE SUBROUTINE SPR_MATMUL(C, A, B, n)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseMatrix)       :: B                                                ! (n,n)
      INTEGER                   :: n
      REAL                      :: C(n,n)                                           ! (n,n)
      REAL                      :: A(n,n)                                           ! (n,n)
!    ! Local variable declaration
      INTEGER                   :: iNonZero,j,k
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: n,A,B
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, B%nNonZero(n) 
        j = B%NonZeroIndex1(iNonZero)
        k = B%NonZeroIndex2(iNonZero)
        C(:,k) = C(:,k) + A(:,j)*B%NonZero(iNonZero)
      ENDDO
!    !
      END SUBROUTINE SPR_MATMUL

    !
    !> The subroutine SPT_MATMUL multiplies a sparse matrix A from the left
    !! to the transpose of a full matrix B and adds the result to the 
    !! transposed matrix C
    !! 
    !< C_ki = C_ki + A_ij * B_kj
    !
      PURE SUBROUTINE SPT_MATMUL_P(C, A, B, n, o)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseMatrix)       :: A                                                  ! Dimension (n,n)
      INTEGER                   :: n,o
      REAL                      :: C(o,n)
      REAL                      :: B(o,n)
!    ! Local variable declaration
      INTEGER                   :: i,j,iNonZero
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: n,o,A,B
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, A%nNonZero(n) 
        i = A%NonZeroIndex1(iNonZero)
        j = A%NonZeroIndex2(iNonZero)
        C(:,i) = C(:,i) + B(:,j)*A%NonZero(iNonZero)
      ENDDO
!    !
      END SUBROUTINE SPT_MATMUL_P
    !
    !> The subroutine SPT_MATMUL multiplies a sparse matrix A from the left
    !! to the transpose of a full matrix B and subtracts the result from the 
    !! transposed matrix C
    !! 
    !< C_ki = C_ki - A_ij * B_kj
    !
      PURE SUBROUTINE SPT_MATMUL_M(C, leadc, A, B, leadb, n, o)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseMatrix)       :: A                                                  ! Dimension (n,n)
      INTEGER                   :: n,o, leadb, leadc
      REAL                      :: C(leadc,n)
      REAL                      :: B(leadb,n)
!    ! Local variable declaration
      INTEGER                   :: i,j,iNonZero
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: n,o,A,B, leadb, leadc
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, A%nNonZero(n) 
        i = A%NonZeroIndex1(iNonZero)
        j = A%NonZeroIndex2(iNonZero)
        C(1:o,i) = C(1:o,i) - B(1:o,j)*A%NonZero(iNonZero)
      ENDDO
!    !
      END SUBROUTINE SPT_MATMUL_M

    !
    !> The subroutine SP_ADD3D adds three sparse matrices A,B,C weighted with
    !! some linear weights x,y,z, adds a scalar s to the diagonal, multiplies with a
    !! full matrix E and adds the result to a full matrix D.
    !! 
    !< D_ki = D_ki + (x*A_ij + y*B_ij + z*C_ij + s*delta_ij) * E_kj
    !
      PURE SUBROUTINE SP_ADD3D_MULT(D, A, B, C, E, s, x, y, z, n, o)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseMatrix)       :: A                                                  ! Dimension (n,n)
      TYPE(tSparseMatrix)       :: B                                                  ! Dimension (n,n)
      TYPE(tSparseMatrix)       :: C                                                  ! Dimension (n,n)      
      REAL                      :: s
      INTEGER                   :: n,o
      REAL                      :: D(o,n)
      REAL                      :: E(o,n)
      REAL                      :: x, y, z
!    ! Local variable declaration
      INTEGER                   :: i,j,iNonZero
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: n,o,A,B,C,E,s,x,y,z
      INTENT(INOUT)             :: D
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, A%nNonZero(n)   
        i = A%NonZeroIndex1(iNonZero)
        j = A%NonZeroIndex2(iNonZero)
        D(:,i) = D(:,i) + x*E(:,j)*A%NonZero(iNonZero)
      ENDDO
      DO iNonZero = 1, B%nNonZero(n) 
        i = B%NonZeroIndex1(iNonZero)
        j = B%NonZeroIndex2(iNonZero)
        D(:,i) = D(:,i) + y*E(:,j)*B%NonZero(iNonZero)
      ENDDO
      DO iNonZero = 1, C%nNonZero(n) 
        i = C%NonZeroIndex1(iNonZero)
        j = C%NonZeroIndex2(iNonZero)
        D(:,i) = D(:,i) + z*E(:,j)*C%NonZero(iNonZero)
      ENDDO
    
      DO i = 1, n
        D(:,i) = D(:,i) + s*E(:,i)
      ENDDO
!    !
      END SUBROUTINE SP_ADD3D_MULT

    !
    !> The subroutine SPT_T3MUL multiplies a sparse 3D tensor A from the left
    !! to the full matrix B and subtracts the result from the 
    !! 3D tensor C
    !! 
    !< C_lpm = C_lpm - A_pqm * B_lq 
    !
      PURE SUBROUTINE SPT_T3MUL_M(C, A, B, p, q, l, m)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseTensor3)      :: A                                                  ! Dimension (p,q,l)
      INTEGER                   :: p,q,l,m
      REAL                      :: C(l,p,m)
      REAL                      :: B(l,q)
!    ! Local variable declaration
      INTEGER                   :: ii,i,j,k,iNonZero
      REAL                      :: Temp
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: p,q,l,m,A,B
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, A%nIndex3
          i = A%Index3(iNonZero)
          C(:,:,i) = 0.
      ENDDO
      DO iNonZero = 1, A%nNonZero 
        i = A%NonZeroIndex1(iNonZero)
        j = A%NonZeroIndex2(iNonZero)
        k = A%NonZeroIndex3(iNonZero)
        C(:,i,k) = C(:,i,k) - B(:,j)*A%NonZero(iNonZero)
      ENDDO
!    !
      END SUBROUTINE SPT_T3MUL_M

    !
    !> The subroutine SPT_T3MUL multiplies a sparse 3D tensor A from the left
    !! to the full matrix B and adds the result to the 
    !! 3D tensor C
    !! 
    !< C_lpm = C_lpm + A_pqm * B_lq 
    !
      PURE SUBROUTINE SPT_T3MUL_P(C, A, B, p, q, l, m)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseTensor3)      :: A                                                  ! Dimension (p,q,l)
      INTEGER                   :: p,q,l,m
      REAL                      :: C(l,p,m)
      REAL                      :: B(l,q)
      REAL                      :: Temp
!    ! Local variable declaration
      INTEGER                   :: i,j,k,iNonZero
      INTEGER                   :: ii
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: p,q,l,m,A,B
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, A%nIndex3
          i = A%Index3(iNonZero)
          C(:,:,i) = 0.
      ENDDO
      DO iNonZero = 1, A%nNonZero 
        i = A%NonZeroIndex1(iNonZero)
        j = A%NonZeroIndex2(iNonZero)
        k = A%NonZeroIndex3(iNonZero)
        C(:,i,k) = C(:,i,k) + B(:,j)*A%NonZero(iNonZero)
      ENDDO
!    !
      END SUBROUTINE SPT_T3MUL_P

    !
    !> The subroutine SPL_T3MUL multiplies a sparse 3D tensor A from the left
    !! to the full 3D tensor B and subtracts the result from the 
    !! 3D tensor C
    !! 
    !< C_kp = C_kp + A_klm * B_lpm 
    !
      PURE SUBROUTINE SPL_T3MUL(C, A, B, p, k, l, m)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseTensor3)      :: A                                                  ! Dimension (p,q,l)
      INTEGER                   :: p,k,l,m
      REAL                      :: C(k,p)
      REAL                      :: B(l,p,m)
!    ! Local variable declaration
      INTEGER                   :: ii,jj,kk,iNonZero
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: p,k,l,m,A,B
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, A%nNonZero 
        ii = A%NonZeroIndex1(iNonZero)
        jj = A%NonZeroIndex2(iNonZero)
        kk = A%NonZeroIndex3(iNonZero)
        C(ii,:) = C(ii,:) + B(jj,:,kk)*A%NonZero(iNonZero)
      ENDDO
!    !
      END SUBROUTINE SPL_T3MUL

    !
    !> The subroutine SPL_T3MULB multiplies a sparse 3D tensor A from the left
    !! to the full 3D tensor B and subtracts the result from the 
    !! 3D tensor C
    !! 
    !< C_kp = C_kp + A_klm * B_lpm 
    !
      SUBROUTINE SPL_T3MULB(C, A, B, p, k, l, m,nIndex3,Index3)
!    !--------------------------------------------------------------------------
      
      IMPLICIT NONE
!    !--------------------------------------------------------------------------
!    ! Argument list declaration
      TYPE(tSparseTensor3b)     :: A                                                  ! Dimension (p,q,l)
      INTEGER                   :: p,k,l,m,i,j
      INTEGER                   :: nIndex3 
      INTEGER                   :: Index3(nIndex3)
      REAL                      :: C(k,p)
      REAL                      :: CT(p,k)
      REAL                      :: B(l,p,m)
      REAL                      :: BT(p,l)
!    ! Local variable declaration
      INTEGER                   :: kk,iNonZero, iNonZero2
!    !--------------------------------------------------------------------------
      INTENT(IN)                :: p,k,l,m,A,B,nIndex3,Index3
      INTENT(INOUT)             :: C
!    !--------------------------------------------------------------------------
!    !
      DO iNonZero = 1, nIndex3
        kk = Index3(iNonZero)
        !CALL SPL_MATMUL(C, k, A%SpSubMatrix(kk), B(1,1,kk), l, k, p)    
        CT = 0.
        BT = TRANSPOSE(B(1:k,:,kk))
        DO iNonZero2 = 1, A%SpSubMatrix(kk)%nNonZero(k)
          i = A%SpSubMatrix(kk)%NonZeroIndex1(iNonZero2)
          j = A%SpSubMatrix(kk)%NonZeroIndex2(iNonZero2)
        CT(:,i) = CT(:,i) + BT(:,j)*A%SpSubMatrix(kk)%NonZero(iNonZero2) 
        ENDDO
        C(1:k,:) = C(1:k,:) + TRANSPOSE(CT)

      ENDDO
!    !
      END SUBROUTINE SPL_T3MULB

END MODULE SP_MATMUL_mod
