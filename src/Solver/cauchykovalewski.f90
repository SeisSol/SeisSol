!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Verena Hermann (hermann AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/hermann)
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

MODULE CauchyKovalewski_mod
  !-------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  INTERFACE CauchyKovalewskiLinear
     MODULE PROCEDURE CauchyKovalewskiLinear3D,  & 
                      CauchyKovalewskiHexa3D
  END INTERFACE
  ! 
  !---------------------------------------------------------------------------!
  PUBLIC  :: CauchyKovalewskiLinear
  PUBLIC  :: CauchyKovalewski3D
  PUBLIC  :: TimeIntTaylor
  !---------------------------------------------------------------------------!

CONTAINS

    SUBROUTINE CauchyKovalewskiLinear3D(TimeIntDOF,TaylorDOF,DOF,tref,t1,t2,              &
                                           LocPoly,LocDegFr,Anelasticity,A,B,C,E,EQN,DISC,IO )
        !-------------------------------------------------------------------------!
        
        USE SP_MATMUL_mod
        !-------------------------------------------------------------------------!
        IMPLICIT NONE
        !-------------------------------------------------------------------------!
        ! Argument list declaration                                               !
        TYPE(tEquations)      :: EQN
        TYPE(tDiscretization) :: DISC
        TYPE(tInputOutput)    :: IO 
        INTEGER               :: LocPoly,LocDegFr,Anelasticity
        TYPE(tSparseMatrix)   :: A
        TYPE(tSparseMatrix)   :: B
        TYPE(tSparseMatrix)   :: C
        TYPE(tSparseMatrix)   :: E
        REAL, OPTIONAL        :: TimeIntDOF(LocDegFr,EQN%nVarTotal)               ! timeintegrated DOF
        REAL                  :: TaylorDOF(LocDegFr,EQN%nVarTotal,0:LocPoly)      ! time - taylorseries for DOF
        REAL                  :: DOF(LocDegFr,EQN%nVarTotal)                      ! Incoming degrees of freedom
        REAL                  :: tref, t1, t2                                     ! Reference, start and end time
        !-------------------------------------------------------------------------!
        ! Local variable declaration                                              !
        INTEGER     :: iOrd,i,r1,r2,r3,j,k                                        ! Loop variables              
        INTEGER     :: iVar,iDegFr                                                ! Index of degree of freedom  
        REAL        :: DOFT(EQN%nVarTotal,DISC%Galerkin%nDegFrRec)                ! Transposed degrees of freedom 
        REAL        :: Temp(LocDegFr,EQN%nVarTotal)                     
        REAL        :: VarTemp(EQN%nVarTotal)
        REAL        :: t1k(0:LocPoly), t2k(0:LocPoly)
        REAL        :: dtk(0:LocPoly)
        !-------------------------------------------------------------------------!
        INTENT(IN)    :: DOF, tref, t1, t2, LocPoly, LocDegFr
        INTENT(IN)    :: EQN,DISC,IO
        INTENT(OUT)   :: TimeIntDOF, TaylorDOF
        !-------------------------------------------------------------------------!
        !
        TaylorDOF(:,:,0) = DOF(:,:)
        !
        SELECT CASE(Anelasticity)
        CASE(0)
            DO iOrd = 0, LocPoly - 1
                ! Initialization
                TaylorDOF(:,:,iOrd + 1) = 0.            
                ! A matrix 
                Temp = 0.
                CALL SPT_MATMUL_M(Temp, LocDegFr, A, TaylorDOF(:,:,iOrd), LocDegFr, EQN%nVarTotal, LocDegFr)            
                CALL SPL_MATMUL(TaylorDOF(:,:,iOrd+1), LocDegFr, DISC%Galerkin%ADG_xi,   Temp, LocDegFr, LocDegFr, EQN%nVarTotal)
                ! B matrix 
                Temp = 0.
                CALL SPT_MATMUL_M(Temp, LocDegFr, B, TaylorDOF(:,:,iOrd), LocDegFr, EQN%nVarTotal, LocDegFr)            
                CALL SPL_MATMUL(TaylorDOF(:,:,iOrd+1), LocDegFr, DISC%Galerkin%ADG_eta,  Temp, LocDegFr, LocDegFr, EQN%nVarTotal)
                ! C matrix 
                Temp = 0.
                CALL SPT_MATMUL_M(Temp, LocDegFr, C, TaylorDOF(:,:,iOrd), LocDegFr, EQN%nVarTotal, LocDegFr)            
                CALL SPL_MATMUL(TaylorDOF(:,:,iOrd+1), LocDegFr, DISC%Galerkin%ADG_zeta, Temp, LocDegFr, LocDegFr, EQN%nVarTotal)
            ENDDO
        CASE(1)
            DO iOrd = 0, LocPoly - 1
                ! Source term
                Temp = 0.
                CALL SPT_MATMUL_M(Temp, LocDegFr, E, TaylorDOF(:,:,iOrd), LocDegFr, EQN%nVarTotal, LocDegFr)            
                TaylorDOF(:,:,iOrd + 1) = Temp(:,:)            
                ! A matrix 
                Temp = 0.
                CALL SPT_MATMUL_M(Temp, LocDegFr, A, TaylorDOF(:,:,iOrd), LocDegFr, EQN%nVarTotal, LocDegFr)            
                CALL SPL_MATMUL(TaylorDOF(:,:,iOrd+1), LocDegFr, DISC%Galerkin%ADG_xi,   Temp, LocDegFr, LocDegFr, EQN%nVarTotal)
                ! B matrix 
                Temp = 0.
                CALL SPT_MATMUL_M(Temp, LocDegFr, B, TaylorDOF(:,:,iOrd), LocDegFr, EQN%nVarTotal, LocDegFr)            
                CALL SPL_MATMUL(TaylorDOF(:,:,iOrd+1), LocDegFr, DISC%Galerkin%ADG_eta,  Temp, LocDegFr, LocDegFr, EQN%nVarTotal)
                ! C matrix 
                Temp = 0.
                CALL SPT_MATMUL_M(Temp, LocDegFr, C, TaylorDOF(:,:,iOrd), LocDegFr, EQN%nVarTotal, LocDegFr)            
                CALL SPL_MATMUL(TaylorDOF(:,:,iOrd+1), LocDegFr, DISC%Galerkin%ADG_zeta, Temp, LocDegFr, LocDegFr, EQN%nVarTotal)
            ENDDO
        END SELECT
        !
        IF(PRESENT(TimeIntDOF)) THEN
            t1k(0) = (t1-tref)
            t2k(0) = (t2-tref)
            dtk(0) = t2k(0)-t1k(0)
            DO i = 1, LocPoly
                t1k(i) = t1k(i-1)*(t1-tref)/REAL(i+1)
                t2k(i) = t2k(i-1)*(t2-tref)/REAL(i+1)
                dtk(i) = t2k(i)-t1k(i)
            ENDDO
            TimeIntDOF(:,:) = 0.
            DO k = 0, LocPoly
                 TimeIntDOF(:,:) = TimeIntDOF(:,:) + TaylorDOF(:,:,k)*dtk(k)
            ENDDO
        ENDIF
        !
    END SUBROUTINE CauchyKovalewskiLinear3D

    SUBROUTINE CauchyKovalewskiHexa3D(TimeIntDOF,TaylorDOF,DOF,tref,t1,t2,              &
                                      LocPoly,LocDegFr,Anelasticity,A,B,C,E,src_deriv,  &
                                      EQN,DISC,IO                                       )
        !-------------------------------------------------------------------------!
        
        USE SP_MATMUL_mod
        !-------------------------------------------------------------------------!
        IMPLICIT NONE
        !-------------------------------------------------------------------------!
        ! Argument list declaration                                               !
        TYPE(tEquations)      :: EQN
        TYPE(tDiscretization) :: DISC
        TYPE(tInputOutput)    :: IO 
        INTEGER               :: LocPoly,LocDegFr,Anelasticity
        TYPE(tSparseTensor3)  :: A
        TYPE(tSparseTensor3)  :: B
        TYPE(tSparseTensor3)  :: C
        TYPE(tSparseTensor3)  :: E
        REAL, OPTIONAL        :: src_deriv(LocDegFr,EQN%nVarTotal,0:LocPoly)      ! source term derivatives 
        REAL, OPTIONAL        :: TimeIntDOF(LocDegFr,EQN%nVarTotal)               ! timeintegrated DOF
        REAL                  :: TaylorDOF(LocDegFr,EQN%nVarTotal,0:LocPoly)      ! time - taylorseries for DOF
        REAL                  :: DOF(LocDegFr,EQN%nVarTotal)                      ! Incoming degrees of freedom
        REAL                  :: tref, t1, t2                                     ! Reference, start and end time
        !-------------------------------------------------------------------------!
        ! Local variable declaration                                              !
        INTEGER     :: iOrd,i,r1,r2,r3,j,k                                        ! Loop variables              
        INTEGER     :: iVar,iDegFr                                                ! Index of degree of freedom  
        INTEGER     :: iNonZero
        REAL        :: DOFT(EQN%nVarTotal,DISC%Galerkin%nDegFrRec)                ! Transposed degrees of freedom 
        REAL        :: Temp(LocDegFr,EQN%nVarTotal,DISC%Galerkin%nDegFrMat)                     
        REAL        :: VarTemp(EQN%nVarTotal)
        REAL        :: t1k(0:LocPoly), t2k(0:LocPoly)
        REAL        :: dtk(0:LocPoly)
        REAL        :: E0(EQN%nVarTotal,EQN%nVarTotal), E1(EQN%nVarTotal,EQN%nVarTotal)
        !-------------------------------------------------------------------------!
        INTENT(IN)    :: DOF, tref, t1, t2, LocPoly, LocDegFr
        INTENT(IN)    :: EQN,DISC,IO
        INTENT(OUT)   :: TimeIntDOF, TaylorDOF
        !-------------------------------------------------------------------------!
        !
        TaylorDOF(:,:,0) = DOF(:,:)
        Temp = 0.d0
        !
        SELECT CASE(Anelasticity)
        CASE(0)
          !
          IF(PRESENT(src_deriv)) THEN
            DO iOrd = 0, LocPoly - 1
                ! Initialization
                TaylorDOF(:,:,iOrd + 1) = 0.            
                ! A matrix 
                CALL SPT_T3MUL_M(Temp, A, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGxiHexa_Sp,   Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,A%nIndex3,A%Index3)
                ! B matrix 
                CALL SPT_T3MUL_M(Temp, B, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGetaHexa_Sp,  Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,B%nIndex3,B%Index3)
                ! C matrix 
                CALL SPT_T3MUL_M(Temp, C, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGzetaHexa_Sp, Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,C%nIndex3,C%Index3)
                ! source derivatives
                TaylorDOF(:,:,iOrd+1) = TaylorDOF(:,:,iOrd+1) + src_deriv(:,:,iOrd)
            ENDDO
          ELSE
            DO iOrd = 0, LocPoly - 1
                ! Initialization
                TaylorDOF(:,:,iOrd + 1) = 0.            
                ! A matrix 
                CALL SPT_T3MUL_M(Temp, A, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGxiHexa_Sp,   Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,A%nIndex3,A%Index3)
                ! B matrix 
                CALL SPT_T3MUL_M(Temp, B, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGetaHexa_Sp,  Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,B%nIndex3,B%Index3)
                ! C matrix 
                CALL SPT_T3MUL_M(Temp, C, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGzetaHexa_Sp, Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,C%nIndex3,C%Index3)
            ENDDO
          ENDIF
          !
        CASE(1)
          !
          IF(PRESENT(src_deriv)) THEN
            DO iOrd = 0, LocPoly - 1
                ! Initialization            
                ! E matrix 
                CALL SPT_T3MUL_M(Temp, E, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                TaylorDOF(:,:,iOrd+1) = Temp(:,:,1) !Still only for piecewise constant E matrices
                ! A matrix 
                CALL SPT_T3MUL_M(Temp, A, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGxiHexa_Sp,   Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,A%nIndex3,A%Index3)
                ! B matrix 
                CALL SPT_T3MUL_M(Temp, B, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGetaHexa_Sp,  Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,B%nIndex3,B%Index3)
                ! C matrix 
                CALL SPT_T3MUL_M(Temp, C, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGzetaHexa_Sp, Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,C%nIndex3,C%Index3)
                ! source derivatives
                TaylorDOF(:,:,iOrd+1) = TaylorDOF(:,:,iOrd+1) + src_deriv(:,:,iOrd)
            ENDDO
          ELSE
            DO iOrd = 0, LocPoly - 1
                ! Initialization
                ! E matrix 
                CALL SPT_T3MUL_M(Temp, E, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                TaylorDOF(:,:,iOrd+1) = Temp(:,:,1) !Still only for piecewise constant E matrices
                ! A matrix 
                CALL SPT_T3MUL_M(Temp, A, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGxiHexa_Sp,   Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,A%nIndex3,A%Index3)
                ! B matrix 
                CALL SPT_T3MUL_M(Temp, B, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGetaHexa_Sp,  Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,B%nIndex3,B%Index3)
                ! C matrix 
                CALL SPT_T3MUL_M(Temp, C, TaylorDOF(:,:,iOrd), EQN%nVarTotal, EQN%nVarTotal, LocDegFr, DISC%Galerkin%nDegFrMat )            
                CALL SPL_T3MULB(TaylorDOF(:,:,iOrd+1), DISC%Galerkin%ADGzetaHexa_Sp, Temp, EQN%nVarTotal, LocDegFr, LocDegFr, DISC%Galerkin%nDegFrMat,C%nIndex3,C%Index3)
            ENDDO
          ENDIF
          !
        END SELECT
        !
        IF(PRESENT(TimeIntDOF)) THEN
            !SELECT CASE(E%nNonZero)
            !CASE(0)
                t1k(0) = (t1-tref)
                t2k(0) = (t2-tref)
                dtk(0) = t2k(0)-t1k(0)
                DO i = 1, LocPoly
                    t1k(i) = t1k(i-1)*(t1-tref)/REAL(i+1)
                    t2k(i) = t2k(i-1)*(t2-tref)/REAL(i+1)
                    dtk(i) = t2k(i)-t1k(i)
                ENDDO
                TimeIntDOF(:,:) = 0.
                DO k = 0, LocPoly
                     TimeIntDOF(:,:) = TimeIntDOF(:,:) + TaylorDOF(:,:,k)*dtk(k)
                ENDDO
            !CASE DEFAULT

            !END SELECT
        ENDIF
        !
    END SUBROUTINE CauchyKovalewskiHexa3D


  !===========================================================================!
  !!                                                                         !!
  !!  CauchyKovalewski3D performs the time integral of the                   !!
  !!  degrees of freedom.                                                    !!
  !!                                                                         !!
  !===========================================================================!
    
    SUBROUTINE CauchyKovalewski3D(TimeIntDof,TimeDerDof,Dof,Dt,             & !
                                  A_Sp,B_Sp,C_Sp,E_Sp,                      & !
                                  ADGxi_Sp, ADGeta_Sp, ADGzeta_Sp, Mklm_Sp, & !
                                  ReactionTerm,LocDegFr,LocDegFrMat,        & !
                                  LocPoly,nVar,Src_Deriv                      )
    !-------------------------------------------------------------------------!
    USE SP_MATMUL_mod
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    REAL                           :: TimeIntDof(:,:)                         ! Time integrated dof
    REAL                           :: TimeDerDof(:,:,0:)                      ! Time derivatives of the dof
    REAL                           :: Dof(:,:)                                ! Actual dof
    REAL                           :: Dt                                      ! Time interval
    TYPE(tSparseTensor3)           :: A_Sp                                    ! Star matrix A (Sparse version)
    TYPE(tSparseTensor3)           :: B_Sp                                    ! Star matrix B (Sparse version)
    TYPE(tSparseTensor3)           :: C_Sp                                    ! Star matrix C (Sparse version)
    TYPE(tSparseTensor3)           :: E_Sp                                    ! Star matrix E (Sparse version)
    TYPE(tSparseTensor3b)          :: ADGxi_Sp                                ! Stiff matrix in xi and l basis function   / Mass matrix (Sparse version)
    TYPE(tSparseTensor3b)          :: ADGeta_Sp                               ! Stiff matrix in eta and l basis function  / Mass matrix (Sparse version)
    TYPE(tSparseTensor3b)          :: ADGzeta_Sp                              ! Stiff matrix in zeta and l basis function / Mass matrix (Sparse version)
    TYPE(tSparseTensor3b)          :: Mklm_Sp                                 ! Int phi_k phi_l phi_m / Mass matrix (Sparse version)
    INTEGER                        :: ReactionTerm                            ! If 1 use the E matrix in CK procedure
    INTEGER                        :: LocDegFr                                ! Number of degrees of freedoms for unknown variable
    INTEGER                        :: LocDegFrMat                             ! Number of degrees of freedoms for material variable
    INTEGER                        :: LocPoly                                 ! Number of time derivatives
    INTEGER                        :: nVar                                    ! Number of variables
    REAL                           :: Src_Deriv(:,:,0:)                       ! Time derivatives of space-time dependent source function
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: Dof, Dt                                                  !
    INTENT(IN)    :: A_Sp, B_Sp, C_Sp, E_Sp                                   !
    INTENT(IN)    :: ADGxi_Sp, ADGeta_Sp, ADGzeta_Sp, Mklm_Sp                 !
    INTENT(IN)    :: ReactionTerm                                             !
    INTENT(IN)    :: LocDegFr                                                 !
    INTENT(IN)    :: LocDegFrMat                                              !
    INTENT(IN)    :: LocPoly                                                  !
    INTENT(IN)    :: nVar                                                     !
    INTENT(OUT)   :: TimeIntDof                                               !
    !-------------------------------------------------------------------------!
    OPTIONAL      :: TimeIntDof, TimeDerDof                                   !
    OPTIONAL      :: Src_Deriv                                                !
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    !-------------------------------------------------------------------------!
    REAL                           :: Aux_TimeDerDof(LocDegFr,nVar,0:LocPoly) ! Time derivatives of fegrees of freedom
    REAL                           :: Temp(LocDegFr,nVar,LocDegFrMat)         ! Temporal tensor
    REAL                           :: dtk(0:LocPoly)                          ! Dt^(k+1)/(k+1)!
    !
    INTEGER                        :: iOrd                                    ! Counter
    INTEGER                        :: i                                       ! Counter
    INTEGER                        :: k                                       ! Counter
    !-------------------------------------------------------------------------!

    Aux_TimeDerDof(:,:,0) = DOF(:,:)
    !
    ! Compute time derivatives of the degrees of freedom    
    DO iOrd = 0, LocPoly - 1
        Aux_TimeDerDof(:,:,iOrd + 1) = 0.0d0                              ! Initialize to zero (iOrd+1)-time derivative
        IF (ReactionTerm.EQ.1) THEN
            ! Sparse E star multiplication
            Temp = 0.
            CALL SPT_T3MUL_P(Temp, E_Sp, Aux_TimeDerDof(:,:,iOrd),      & ! Temp_lpm = Temp_lpm + E_pqm * u_lq
                             nVar, nVar, LocDegFr, LocDegFrMat )          !            
            CALL SPL_T3MULB(Aux_TimeDerDof(:,:,iOrd+1), Mklm_Sp, Temp,  & ! Taylor_kp = Taylor_lp + MassMatrix_klm * Temp_lpm
                             nVar, LocDegFr, LocDegFr, LocDegFrMat,     & !
                             E_Sp%nIndex3, E_Sp%Index3 )                  !
        ENDIF
        !
        ! Sparse A star multiplication
        Temp = 0.
        CALL SPT_T3MUL_M(Temp, A_Sp, Aux_TimeDerDof(:,:,iOrd),          & ! Temp_lpm = Temp_lpm - A_pqm * u_lq
                         nVar, nVar, LocDegFr, LocDegFrMat )              !            
        CALL SPL_T3MULB(Aux_TimeDerDof(:,:,iOrd+1), ADGxi_Sp, Temp,     & ! Taylor_kp = Taylor_lp + ADGxi_klm * Temp_lpm
                         nVar, LocDegFr, LocDegFr, LocDegFrMat,         & !
                         A_Sp%nIndex3, A_Sp%Index3 )                      !
        !
        ! Sparse B star multiplication
        Temp = 0.
        CALL SPT_T3MUL_M(Temp, B_Sp, Aux_TimeDerDof(:,:,iOrd),          & ! Temp_lpm = Temp_lpm - B_pqm * u_lq
                         nVar, nVar, LocDegFr, LocDegFrMat )              !            
        CALL SPL_T3MULB(Aux_TimeDerDof(:,:,iOrd+1), ADGeta_Sp, Temp,    & ! Taylor_kp = Taylor_lp + ADGeta_klm * Temp_lpm
                         nVar, LocDegFr, LocDegFr, LocDegFrMat,         & !
                         B_Sp%nIndex3, B_Sp%Index3 )                      !
        !                 
        ! Sparse C star multiplication
        Temp = 0.
        CALL SPT_T3MUL_M(Temp, C_Sp, Aux_TimeDerDof(:,:,iOrd),          & ! Temp_lpm = Temp_lpm - C_pqm * u_lq
                         nVar, nVar, LocDegFr, LocDegFrMat )              !            
        CALL SPL_T3MULB(Aux_TimeDerDof(:,:,iOrd+1), ADGzeta_Sp, Temp,   & ! Taylor_kp = Taylor_lp + ADGzeta_klm * Temp_lpm
                         nVar, LocDegFr, LocDegFr, LocDegFrMat,         & !
                         C_Sp%nIndex3, C_Sp%Index3 )                      !
        !                 
        ! Source derivatives
        IF(PRESENT(src_deriv)) THEN
            Aux_TimeDerDof(:,:,iOrd+1) = Aux_TimeDerDof(:,:,iOrd+1) +   & ! Add dof time derivatives of the source term
                                         src_deriv(:,:,iOrd)              !
        ENDIF
        !
    ENDDO ! iOrd
    !
    IF (PRESENT(TimeDerDof)) THEN
        ! Compute the time derivatives of the degrees of freedom
        TimeDerDof = Aux_TimeDerDof
    ENDIF
    !
    IF (PRESENT(TimeIntDOF)) THEN
        ! Compute the time integrated degrees of freedom
        dtk(0) = Dt                                                           ! Dt^1/(1)
        DO i = 1, LocPoly                                                     !
            dtk(i) = dtk(i-1)*Dt/REAL(i+1)                                    ! Dt^{k+1}/(k+1)!
        ENDDO
        !
        TimeIntDOF(:,:) = 0.0d0                                               ! Initialization
        DO k = 0, LocPoly                                                     ! 
             TimeIntDOF(:,:) = TimeIntDOF(:,:) + Aux_TimeDerDof(:,:,k)*dtk(k) ! Integral over [0,Dt] of dof time derivatives
        ENDDO                                                                 !
    ENDIF
    !
    CONTINUE
    !
    END SUBROUTINE CauchyKovalewski3D

    ! Integrate a Taylor series constructed at tref
    ! in the time interval [t1,t2]
    SUBROUTINE TimeIntTaylor(TimeIntDof,TimeDerDof,t0,t1,t2,       &
                                        LocnPoly,LocnVar,LocnDegFr            )
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration                                               !
    REAL                           :: TimeIntDof(LocnDegFr,LocnVar)           ! Time integrated dof
    REAL                           :: TimeDerDof(LocnDegFr,LocnVar,0:LocnPoly)! Time derivatives of the dof
    REAL                           :: t0                                      ! Time lever of Taylor series construction
    REAL                           :: t1                                      ! Lower integration time
    REAL                           :: t2                                      ! Higer integration time
    INTEGER                        :: LocnPoly
    INTEGER                        :: LocnVar
    INTEGER                        :: LocnDegFr
    !-------------------------------------------------------------------------!
    INTENT(IN)                     :: TimeDerDof,t0,t1,t2
    INTENT(IN)                     :: LocnPoly,LocnVar,LocnDegFr
    INTENT(OUT)                    :: TimeIntDof
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    !-------------------------------------------------------------------------!
    REAL                           :: t1k(0:LocnPoly)
    REAL                           :: t2k(0:LocnPoly)
    REAL                           :: dtk(0:LocnPoly)
    INTEGER                        :: iDegPol
    !-------------------------------------------------------------------------!
    !
    TimeIntDof = 0.0d0
    ! Time-integration of the element itself, using Taylor expansion
    t1k(0) = (t1-t0)
    t2k(0) = (t2-t0)
    dtk(0) = t2k(0)-t1k(0)
    DO iDegPol = 1, LocnPoly
        t1k(iDegPol) = t1k(iDegPol-1)*(t1-t0)/REAL(iDegPol+1)
        t2k(iDegPol) = t2k(iDegPol-1)*(t2-t0)/REAL(iDegPol+1)
        dtk(iDegPol) = t2k(iDegPol)-t1k(iDegPol)
    ENDDO
    DO iDegPol = 0, LocnPoly
        TimeIntDof = TimeIntDof + TimeDerDof(1:LocnDegFr,1:LocnVar,iDegPol) * dtk(iDegPol)
    ENDDO
    !
    CONTINUE
    !
    END SUBROUTINE TimeIntTaylor

END MODULE CauchyKovalewski_mod
