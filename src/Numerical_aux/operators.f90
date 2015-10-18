!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2009, SeisSol Group
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

MODULE COMMON_operators_mod
  !----------------------------------------------------------------------------
  USE typesDef
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------

  INTERFACE OPERATOR (.CA.)
     MODULE PROCEDURE Circa
  END INTERFACE

  INTERFACE OPERATOR (.NE.)
     MODULE PROCEDURE notEqualLogicals
  END INTERFACE

  INTERFACE OPERATOR (.EQ.)
     MODULE PROCEDURE EqualLogicals
  END INTERFACE

  INTERFACE OPERATOR (.x.)
     MODULE PROCEDURE Kreuzprodukt
  END INTERFACE

  INTERFACE OPERATOR (.cp.)
     MODULE PROCEDURE cross_product
  END INTERFACE

  INTERFACE OPERATOR (.vn.)
     MODULE PROCEDURE vector_norm
  END INTERFACE

  INTERFACE OPERATOR (.vp.)
     MODULE PROCEDURE scalar_product
  END INTERFACE

  INTERFACE OPERATOR (.skalar.)
     MODULE PROCEDURE skalarProdukt
  END INTERFACE

  INTERFACE OPERATOR (.dyadic.)
     MODULE PROCEDURE dyadicProduct
  END INTERFACE

  INTERFACE OPERATOR (.IN.)
     MODULE PROCEDURE isValueInList
  END INTERFACE

  INTERFACE OPERATOR(.im.)
      MODULE PROCEDURE MatrixMatrixProductInv
  END INTERFACE

  INTERFACE OPERATOR(.m.)
      MODULE PROCEDURE MatrixMatrixProduct
  END INTERFACE

  INTERFACE Heaviside
     MODULE PROCEDURE Heaviside
  END INTERFACE

  INTERFACE int_to_logical
     MODULE PROCEDURE int_to_logical
  END INTERFACE

  INTERFACE ASSIGNMENT(=)
     MODULE PROCEDURE logical_to_integer     , &
                      integer_to_logical     
  END INTERFACE

  INTERFACE STRING
     MODULE PROCEDURE STRING_INT, STRING_REAL
  END INTERFACE

  INTERFACE OpenFile
     MODULE PROCEDURE OpenFile
  END INTERFACE

  INTERFACE winkel
     MODULE PROCEDURE winkel
  END INTERFACE

  INTERFACE determinant
     MODULE PROCEDURE determinant2d, determinant3d, determinant5d
  END INTERFACE

  INTERFACE AddNodeToList
     MODULE PROCEDURE AddNodeToList
  END INTERFACE

  INTERFACE FindNodeInList
     MODULE PROCEDURE FindNodeInList
  END INTERFACE

  INTERFACE ZerosPolyO2
     MODULE PROCEDURE ZerosPolyO2
  END INTERFACE

  INTERFACE ZerosPolyO3
     MODULE PROCEDURE ZerosPolyO3
  END INTERFACE

  INTERFACE ZerosPolyO4
     MODULE PROCEDURE ZerosPolyO4
  END INTERFACE

  INTERFACE AngelBetweenVectors
     MODULE PROCEDURE AngelBetweenVectors
  END INTERFACE

  INTERFACE IniSparseMatrix
     MODULE PROCEDURE IniSparseMatrix
  END INTERFACE

  INTERFACE UnpackSparseMatrix
     MODULE PROCEDURE UnpackSparseMatrix
  END INTERFACE

  INTERFACE CloseSparseMatrix
     MODULE PROCEDURE CloseSparseMatrix
  END INTERFACE

  INTERFACE IniSparseMatrix2
     MODULE PROCEDURE IniSparseMatrix2
  END INTERFACE

  INTERFACE CloseSparseMatrix2
     MODULE PROCEDURE CloseSparseMatrix2
  END INTERFACE

  INTERFACE IniSparseVector
     MODULE PROCEDURE IniSparseVector
  END INTERFACE

  INTERFACE CloseSparseVector
     MODULE PROCEDURE CloseSparseVector
  END INTERFACE

  INTERFACE IniSparseTensor3
     MODULE PROCEDURE IniSparseTensor3
  END INTERFACE

  INTERFACE UnpackSparseTensor3
     MODULE PROCEDURE UnpackSparseTensor3
  END INTERFACE

  INTERFACE CloseSparseTensor3
     MODULE PROCEDURE CloseSparseTensor3
  END INTERFACE

  INTERFACE IniSparseTensor3b
     MODULE PROCEDURE IniSparseTensor3b
  END INTERFACE

  INTERFACE CloseSparseTensor3b
     MODULE PROCEDURE CloseSparseTensor3b
  END INTERFACE

  INTERFACE IniSparseTensor4
     MODULE PROCEDURE IniSparseTensor4
  END INTERFACE

  INTERFACE CloseSparseTensor4
     MODULE PROCEDURE CloseSparseTensor4
  END INTERFACE

  INTERFACE MatrixInverse3x3
     MODULE PROCEDURE MatrixInverse3x3
  END INTERFACE
  
  INTERFACE MatrixInverse2x2
     MODULE PROCEDURE MatrixInverse2x2
  END INTERFACE
  
  INTERFACE TF_MISFITS
     MODULE PROCEDURE TF_MISFITS
  END INTERFACE

  INTERFACE FCOOLR
     MODULE PROCEDURE FCOOLR
  END INTERFACE

  INTERFACE CWT
     MODULE PROCEDURE CWT
  END INTERFACE

  INTERFACE MORLET
     MODULE PROCEDURE MORLET
  END INTERFACE
  
  INTERFACE dwalltime
     MODULE PROCEDURE dwalltime
  END INTERFACE

  INTERFACE XYinTriangle
     MODULE PROCEDURE XYinTriangle
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC :: OPERATOR(.CA.)
  PUBLIC :: OPERATOR(.NE.)
  PUBLIC :: OPERATOR(.EQ.)
  PUBLIC :: OPERATOR(.x.)
  PUBLIC :: OPERATOR(.cp.)
  PUBLIC :: OPERATOR(.vn.)
  PUBLIC :: OPERATOR(.vp.)
  PUBLIC :: OPERATOR(.skalar.)
  PUBLIC :: OPERATOR(.dyadic.)
  PUBLIC :: OPERATOR(.IN.)
  PUBLIC :: OPERATOR(.im.)
  PUBLIC :: OPERATOR(.m.)
  PUBLIC :: ASSIGNMENT(=)
  PUBLIC :: winkel
  PUBLIC :: int_to_logical
  PUBLIC :: STRING
  PUBLIC :: OpenFile
  PUBLIC :: determinant
  PUBLIC :: AddNodeToList
  PUBLIC :: FindNodeInList
  PUBLIC :: ZerosPolyO2
  PUBLIC :: ZerosPolyO3
  PUBLIC :: ZerosPolyO4
  PUBLIC :: AngelBetweenVectors
  PUBLIC :: Heaviside
  PUBLIC :: IniSparseVector
  PUBLIC :: CloseSparseVector
  PUBLIC :: IniSparseMatrix
  PUBLIC :: UnpackSparseMatrix
  PUBLIC :: CloseSparseMatrix
  PUBLIC :: IniSparseMatrix2
  PUBLIC :: CloseSparseMatrix2
  PUBLIC :: IniSparseTensor3
  PUBLIC :: UnpackSparseTensor3
  PUBLIC :: CloseSparseTensor3
  PUBLIC :: IniSparseTensor3b
  PUBLIC :: CloseSparseTensor3b
  PUBLIC :: IniSparseTensor4
  PUBLIC :: CloseSparseTensor4
  PUBLIC :: MatrixInverse3x3
  PUBLIC :: MatrixInverse2x2
  PUBLIC :: TF_MISFITS
  PUBLIC :: dwalltime
  PUBLIC :: XYinTriangle
  !----------------------------------------------------------------------------
  
  INTEGER, PARAMETER :: ik = selected_int_kind(6)  ! needed for dwalltime
  LOGICAL, SAVE      :: first = .true.             ! needed for dwalltime
  INTEGER(ik), SAVE  :: count_rate, count_max      ! needed for dwalltime
  REAL, SAVE         :: conversion = 0.0d0         ! needed for dwalltime

CONTAINS

  FUNCTION Circa(x,y)
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    REAL    :: x,y
    LOGICAL :: Circa
    REAL    :: tol
    !--------------------------------------------------------------------------
    INTENT(IN) :: x,y
    !--------------------------------------------------------------------------
    
    tol = 1e-4
    IF(ABS(x-y).LT.tol) THEN
      Circa = .TRUE.
    ELSE
      Circa = .FALSE.
    ENDIF

  END FUNCTION Circa
  
  FUNCTION notEqualLogicals(log1,log2)
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    LOGICAL :: log1
    LOGICAL :: log2
    LOGICAL :: notEqualLogicals
    !--------------------------------------------------------------------------
    INTENT(IN) :: log1, log2
    !--------------------------------------------------------------------------

    notEqualLogicals = log1.NEQV.log2

  END FUNCTION notEqualLogicals

  FUNCTION EqualLogicals(log1,log2)
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    LOGICAL :: log1
    LOGICAL :: log2
    LOGICAL :: EqualLogicals
    !--------------------------------------------------------------------------
    INTENT(IN) :: log1, log2
    !--------------------------------------------------------------------------

    EqualLogicals = .NOT.(log1.NEQV.log2)

  END FUNCTION EqualLogicals

  FUNCTION Kreuzprodukt(vec_1,vec_2)
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    REAL :: Kreuzprodukt(3)
    REAL :: vec_1(3)
    REAL :: vec_2(3)
    !--------------------------------------------------------------------------
    INTENT(IN) :: vec_1,vec_2
    !--------------------------------------------------------------------------
    Kreuzprodukt(1) = vec_1(2)*vec_2(3) - vec_1(3)*vec_2(2)
    Kreuzprodukt(2) = vec_1(3)*vec_2(1) - vec_1(1)*vec_2(3)
    Kreuzprodukt(3) = vec_1(1)*vec_2(2) - vec_1(2)*vec_2(1)
  END FUNCTION Kreuzprodukt

  FUNCTION cross_product(vec_1,vec_2) 
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    TYPE(tVector)    :: cross_product,help
    TYPE(tVector)    :: vec_1, vec_2
    !--------------------------------------------------------------------------
    INTENT(IN) :: vec_1,vec_2
    !--------------------------------------------------------------------------
    help%x(1:vec_1%VecLength) = vec_1%y*vec_2%z - vec_1%z*vec_2%y
    help%y(1:vec_1%VecLength) = vec_1%z*vec_2%x - vec_1%x*vec_2%z
    help%z(1:vec_1%VecLength) = vec_1%x*vec_2%y - vec_1%y*vec_2%x
    cross_product = help
  END FUNCTION cross_product

  FUNCTION vector_norm(vec_1,dummy)
     
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    TYPE(tVector)    :: vec_1
    REAL             :: vector_norm(vec_1%VecLength)
    INTEGER          :: dummy
    !--------------------------------------------------------------------------
    INTENT(IN)       :: vec_1,dummy
    !--------------------------------------------------------------------------
    vector_norm(1) = dummy                                                           ! Only to avoid Info Message
    vector_norm(1:vec_1%VecLength) = vec_1%x*vec_1%x
    vector_norm(1:vec_1%VecLength) = vector_norm(1:vec_1%VecLength) + vec_1%y*vec_1%y
    vector_norm(1:vec_1%VecLength) = vector_norm(1:vec_1%VecLength) + vec_1%z*vec_1%z
    vector_norm(1:vec_1%VecLength) = SQRT(vector_norm(1:vec_1%VecLength))
  END FUNCTION vector_norm

  FUNCTION scalar_product(vec_1,vec_2)
     
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    TYPE(tVector)    :: vec_1
    TYPE(tVector)    :: vec_2
    REAL             :: scalar_product(vec_1%VecLength)
    !--------------------------------------------------------------------------
    INTENT(IN)       :: vec_1,vec_2
    !--------------------------------------------------------------------------
    scalar_product(1:vec_1%VecLength) = vec_1%x*vec_2%x
    scalar_product(1:vec_1%VecLength) = scalar_product(1:vec_1%VecLength) + vec_1%y*vec_2%y
    scalar_product(1:vec_1%VecLength) = scalar_product(1:vec_1%VecLength) + vec_1%z*vec_1%z
  END FUNCTION scalar_product

  PURE FUNCTION dyadicproduct(vec_1,vec_2)
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    REAL :: dyadicproduct(3,3)
    REAL :: vec_1(3)
    REAL :: vec_2(3)
    INTEGER::i,j
    !--------------------------------------------------------------------------
    INTENT(IN) :: vec_1,vec_2
    !--------------------------------------------------------------------------
    !                                                  !
    DO i=1,3 
      DO j=1,3
        dyadicproduct(i,j)=vec_1(i)*vec_2(j)
      ENDDO
    ENDDO
    !                                                  !
  END FUNCTION dyadicproduct

  PURE FUNCTION skalarprodukt(vec_1,vec_2)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    REAL :: skalarprodukt
    REAL :: vec_1(:)
    REAL :: vec_2(:)
    INTENT(IN) :: vec_1,vec_2
    !--------------------------------------------------------------------------
    skalarprodukt=SUM(vec_1(:)*vec_2(:))
    !                                                  !
  END FUNCTION skalarprodukt

  FUNCTION isValueInList(value,list)
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    LOGICAL :: isValueInList
    INTEGER :: value
    INTEGER :: list(:)
    !local Variables
    INTEGER :: fieldSize
    INTEGER :: I
    !--------------------------------------------------------------------------
    INTENT(IN) :: value, list
    !--------------------------------------------------------------------------

    fieldSize     = SIZE(list)
    isValueInList = .FALSE.   

    DO I=1,fieldSize
       IF (list(I).EQ.value) THEN
          isValueInList = .TRUE.
       END IF
    END DO

  END FUNCTION isValueInList
    
  SUBROUTINE logical_to_integer(integerVar,logicalVar)
    !-----------------------------------------------------------------------
    LOGICAL     :: logicalVar
    INTEGER     :: integerVar
    !-----------------------------------------------------------------------
    INTENT(IN)  :: logicalVar
    INTENT(OUT) :: integerVar
    !-----------------------------------------------------------------------
    
    IF (logicalVar) THEN
       integerVar = 1
    ELSE
       integerVar = 0
    END IF
    
  END SUBROUTINE logical_to_integer
  
  PURE ELEMENTAL SUBROUTINE integer_to_logical(logicalVar,integerVar)
    !-----------------------------------------------------------------------
    LOGICAL     :: logicalVar
    INTEGER     :: integerVar
    !-----------------------------------------------------------------------
    INTENT(IN)  :: integerVar
    INTENT(OUT) :: logicalVar
    !-----------------------------------------------------------------------
    
    IF (integerVar .EQ. 0) THEN
       logicalVar = .FALSE.
    ELSE
       logicalVar = .TRUE.
    END IF
    
  END SUBROUTINE integer_to_logical
  
  PURE ELEMENTAL FUNCTION Int_To_Logical(int,TrueValue,FalseValue,Default) 
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER     :: int
    INTEGER     :: TrueValue
    INTEGER     :: FalseValue
    LOGICAL     :: DEFAULT
    LOGICAL     :: Int_To_Logical
    !------------------------------------------------------------------------!
    INTENT(IN)  :: int,TrueValue,FalseValue,Default
    !------------------------------------------------------------------------!
    
    IF (int.EQ.TrueValue) THEN
       Int_To_Logical = .TRUE.
    ELSE IF (int.EQ. FalseValue) THEN
       Int_To_Logical = .FALSE.
    ELSE
       Int_To_Logical = Default
    END IF
    
  END FUNCTION Int_To_Logical
  
  PURE FUNCTION STRING_INT(int,digits)
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER           :: int
    INTEGER, OPTIONAL :: digits
    CHARACTER(LEN=10) :: STRING_int
    CHARACTER(LEN=15) :: Format
    CHARACTER(LEN=5)  :: cDigits
    !------------------------------------------------------------------------!
    INTENT(IN)        :: int
    INTENT(IN)        :: digits
    !------------------------------------------------------------------------!

    IF (PRESENT(digits)) THEN

       WRITE(cDigits,'(I5.2)') digits

       Format='(I'//TRIM(cDigits)//'.'//TRIM(cDigits)//')'

    ELSE

       Format='(I9.9)'

    END IF

    WRITE(STRING_int,Format) int                                                       

  END FUNCTION STRING_INT

  PURE FUNCTION STRING_real(number,digits)
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL              :: number
    INTEGER, OPTIONAL :: digits
    CHARACTER(LEN=30) :: STRING_real
    CHARACTER(LEN=18) :: Format
    CHARACTER(LEN=5)  :: cDigits
    CHARACTER(LEN=5)  :: cDigits2
    !------------------------------------------------------------------------!
    INTENT(IN)        :: number
    INTENT(IN)        :: digits
    !------------------------------------------------------------------------!

    IF (PRESENT(digits)) THEN

       WRITE(cDigits,'(I5.2)') digits
       WRITE(cDigits2,'(I5.2)') digits+7

       Format='(SP,E'//TRIM(cDigits2)//'.'//TRIM(cDigits)//')'

    ELSE

       Format='(SP,E16.9)'

    END IF

    WRITE(STRING_real,Format) number

  END FUNCTION STRING_REAL

  SUBROUTINE OpenFile(UnitNr,Name,create)
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    INTEGER            :: UnitNr
    CHARACTER(LEN=*)   :: Name
    LOGICAL            :: create
    ! local variables
    INTEGER            :: status
    LOGICAL            :: fexist
    CHARACTER(LEN=100) :: AllreadyOpenFileName
    CHARACTER(LEN=10)  :: FileStatus
    !------------------------------------------------------------------------!
    INTENT(IN)         :: UnitNr, Name, create
    !------------------------------------------------------------------------!
    !                                                                        !        
    INQUIRE(                                                               & !
         UNIT  = UnitNr                                                 ,  & !
         EXIST = fexist                                                 ,  & !
         NAME  = AllreadyOpenFileName                                      ) !
    !                                                                        !        
    IF(.NOT.fexist) THEN                                                     !
       logError(*) 'Error opening unit    : ',UnitNr              !
       logError(*) 'Unit is already open under the name:',TRIM(AllreadyOpenFileName)
       STOP                                                                  !
    ENDIF                                                                    !
    !                                                                        !        
    IF (create) THEN                                                         !
       FileStatus = 'UNKNOWN'                                                !
    ELSE                                                                     !
       !                                                                     !            
       FileStatus = 'OLD'                                                    !
       !                                                                     !            
       INQUIRE(                                                            & !
            FILE  = TRIM(Name)                                          ,  & !
            EXIST = fexist                                                 ) !
       !                                                                     !        
       IF(.NOT.fexist) THEN                                                  !
          logError(*) 'Error opening file    : ', Name            !
          logError(*) 'File does not exist.'                       !
          STOP                                                               !
       ENDIF                                                                 !
       !                                                                     !
    END IF                                                                   !
    !                                                                        !        
    OPEN(UNIT   = UnitNr                                                 , & !
         FILE   = TRIM(Name)                                             , & !
         STATUS = TRIM(FileStatus)                                       , & !
         RECL   = 700                                                    , & !
         IOSTAT = status                                                   ) !
    !                                                                        !        
    IF (status .NE. 0) THEN                                                  !
       logError(*) 'could not open ',Name                   !
       logError(*) 'File does exist but cannot be opened.'         !
       logError(*) 'IOSTAT:',status                                !
       STOP                                                                  !
    END IF                                                                   !
    !                                                                        !        
  END SUBROUTINE OpenFile

    PURE FUNCTION MatrixMatrixProduct(A,B)
      !------------------------------------------------------------------------!
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      ! Argument list declaration      
      REAL             :: A(1:5,1:5), B(1:5,1:5)
      ! Return variable declaration      
      REAL             :: MatrixMatrixProduct(1:5,1:5)
      ! Local variable declaration      
      INTEGER          :: i
      !------------------------------------------------------------------------!
      INTENT(IN)       :: A, B
      !------------------------------------------------------------------------!
      
      DO i=1,5
        MatrixMatrixProduct(:,i) =                     &                       
                                    A(:,1)*B(1,i) +    &
                                    A(:,2)*B(2,i) +    &
                                    A(:,3)*B(3,i) +    &
                                    A(:,4)*B(4,i) +    &
                                    A(:,5)*B(5,i)
      ENDDO

    END FUNCTION MatrixMatrixProduct

    PURE FUNCTION Heaviside(t)
    !--------------------------------------------------------------------------
    IMPLICIT NONE                                            
    !--------------------------------------------------------------------------
    ! Argument list declaration  
    REAL             :: t
    !--------------------------------------------------------------------------
    ! Return variable declaration    
    INTEGER          :: Heaviside
    !--------------------------------------------------------------------------
    INTENT(IN)       :: t
    !--------------------------------------------------------------------------
     IF(t.LT.0.)THEN
       Heaviside = 0
     ELSE
       Heaviside = 1
     ENDIF 
   
    END FUNCTION Heaviside

    PURE FUNCTION MatrixMatrixProductInv(B,A)
      !------------------------------------------------------------------------!
      IMPLICIT NONE
      !------------------------------------------------------------------------!
      ! Argument list declaration      
      REAL             :: A(1:5,1:5), B(1:5,1:5)
      ! Return variable declaration      
      REAL             :: MatrixMatrixProductInv(1:5,1:5)
      ! Local variable declaration      
      INTEGER          :: i
      !------------------------------------------------------------------------!
      INTENT(IN)       :: A, B
      !------------------------------------------------------------------------!
      
      DO i=1,5
        MatrixMatrixProductInv(:,i) =                     &                       
                                    A(:,1)*B(1,i) +    &
                                    A(:,2)*B(2,i) +    &
                                    A(:,3)*B(3,i) +    &
                                    A(:,4)*B(4,i) +    &
                                    A(:,5)*B(5,i)
      ENDDO

    END FUNCTION MatrixMatrixProductInv

    SUBROUTINE winkel(p1,p2,phi,warnOut,quiet)
      !--------------------------------------------------------------------------
      IMPLICIT NONE                                            
      !--------------------------------------------------------------------------
      INTEGER        ,OPTIONAL :: warnOut
      REAL                     :: p1(2)
      REAL                     :: p2(2)
      REAL                     :: phi
      REAL                     :: dx,dy,r,phi1,phi2
      LOGICAL, OPTIONAL        :: quiet
      ! local Variables
      INTEGER                  :: UnitNr
      REAL                     :: arg1, arg2
      REAL                     :: pi
      LOGICAL                  :: quiet_internal
      !--------------------------------------------------------------------------
      INTENT(IN)               :: p1,p2,warnOut,quiet
      INTENT(OUT)              :: phi
      !--------------------------------------------------------------------------
      !                                                   !
      pi = ACOS(-1.0)                                     !
      !                                                   !   
      IF (PRESENT(quiet)) THEN                            ! 
         quiet_internal = quiet                           ! 
      ELSE                                                ! 
         quiet_internal = .FALSE.                         ! 
      END IF                                              ! 
      !                                                   ! 
      IF (PRESENT(warnOut)) THEN                          !   
         UnitNr = warnOut                                 !   
      ELSE                                                ! Falls die Winkel Routine von   
         UnitNr = 6                                       ! den Testroutinen aus aufgerufen wird
      END IF                                              ! soll die Ausgabe auf den Bildschirm
      !                                                   ! erfolgen. Ansonsten in die Warnings Datei
      !                                                   ! (Da in der Regel h�ufig kleinere Abweichungen
      !                                                   ! auftreten  
      dx   = p2(1)-p1(1)                                  ! 
      dy   = p2(2)-p1(2)                                  ! 
      r    = SQRT(dx*dx+dy*dy)                            ! 
      !                                                   !   
      IF (r.EQ.0.0) THEN                                  ! 
         phi = 0.0                                        ! 
         IF (.NOT.quiet_internal) THEN                    !
            WRITE(UnitNr,*)'WARNING: winkel'              ! 
         END IF                                           !
      ELSE                                                ! 
         !                                                !
         arg1 = dy/r                                      !
         arg2 = dx/r                                      !
         !                                                !
         IF (arg1.GT.1.0) THEN                            !
            arg1 = 1.0                                    !
         END IF                                           !
         IF (arg1.LT.-1.0) THEN                           !
            arg1 = -1.0                                   !
         END IF                                           !
         IF (arg2.GT.1.0) THEN                            !
            arg2 = 1.0                                    !
         END IF                                           !
         IF (arg2.LT.-1.0) THEN                           !
            arg2 = -1.0                                   !
         END IF                                           !
         !                                                !
         phi1 = ASIN(arg1)                                ! 
         phi2 = ACOS(arg2)                                ! 
         !                                                !
         IF (phi1.GE.0.0 ) THEN                           ! x+y+ oder x-y+
            !                                             !      Quadrant 
            IF (phi2.LE.pi/2.0) THEN                      ! x+y+ Quadrant
               !                                          !
               phi = phi2                                 !
               !                                          !
!!$               IF (ABS(phi1-phi).GT.tol) THEN             ! Fehlerkontrolle
!!$                  IF (.NOT.quiet_internal) THEN           ! Fehlerkontrolle
!!$                     WRITE(UnitNr,*)'Warning: winkel'     ! Fehlerkontrolle
!!$                     WRITE(UnitNr,*)'     phi1         =',phi1 ! Fehlerkontrolle
!!$                     WRITE(UnitNr,*)'     phi2         =',phi2 ! Fehlerkontrolle
!!$                     WRITE(UnitNr,*)'     tol          =',tol  ! Fehlerkontrolle
!!$                     WRITE(UnitNr,*)'     ABS(phi1-phi)=',ABS(phi1-phi)
!!$                  END IF                                  ! Fehlerkontrolle
!!$               END IF                                     ! Fehlerkontrolle
               !                                          !
            ELSE                                          ! x-y+ Quadrant
               phi = phi2                                 !
            END IF                                        !
            !                                             !
         ELSE                                             ! x-y- oder x+y- 
            !                                             !      Quadrant
            IF (phi2.LE.pi/2.0) THEN                      ! x+y- Quadrant
               phi = 2*pi + phi1                          ! Achtung: phi1 ist 
               !                                          !          negativ
            ELSE                                          ! x-y- Quadrant
               phi =   pi - phi1                          ! Achtung: phi1 ist
               !                                          !          negativ
            END IF                                        !
            !                                             !
         END IF                                           !
         !                                                ! 
      END IF                                              !
      !                                                   !   
    END SUBROUTINE winkel                                 !

  PURE FUNCTION determinant2d(A,B)
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL              :: determinant2d
    REAL              :: A(2),B(2)
    !------------------------------------------------------------------------!
    INTENT(IN)        :: A,B
    !------------------------------------------------------------------------!

    determinant2d=A(1)*B(2)-A(2)*B(1)

  END FUNCTION determinant2d

  PURE FUNCTION determinant3d(A,B,C)
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL              :: determinant3d
    REAL              :: A(3),B(3),C(3)
    !------------------------------------------------------------------------!
    INTENT(IN)        :: A,B,C
    !------------------------------------------------------------------------!
    determinant3d=A(1)*B(2)*C(3)-A(1)*B(3)*C(2)+A(2)*B(3)*C(1)-A(2)*B(1)*C(3) &
                  +A(3)*B(1)*C(2)-A(3)*B(2)*C(1)
  END FUNCTION determinant3d

  PURE FUNCTION determinant4d(A,B,C,D)
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL              :: determinant4d
    REAL              :: A(4),B(4),C(4),D(4)
    !------------------------------------------------------------------------!
    INTENT(IN)        :: A,B,C,D
    !------------------------------------------------------------------------!
    determinant4d = A(1)*determinant3d(B(2:4),C(2:4),D(2:4)) - &
                    B(1)*determinant3d(A(2:4),C(2:4),D(2:4)) + &
                    C(1)*determinant3d(A(2:4),B(2:4),D(2:4)) - &
                    D(1)*determinant3d(A(2:4),B(2:4),C(2:4))
  END FUNCTION determinant4d

  PURE FUNCTION determinant5d(A,B,C,D,E)
    !------------------------------------------------------------------------!
    IMPLICIT NONE
    !------------------------------------------------------------------------!
    REAL              :: determinant5d                                       !
    REAL              :: A(5),B(5),C(5),D(5),E(5)                            !
    !------------------------------------------------------------------------!
    INTENT(IN)        :: A,B,C,D,E
    !------------------------------------------------------------------------!
    determinant5d = A(1)*determinant4d(B(2:5),C(2:5),D(2:5),E(2:5)) - &
                    B(1)*determinant4d(A(2:5),C(2:5),D(2:5),E(2:5)) + &
                    C(1)*determinant4d(A(2:5),B(2:5),D(2:5),E(2:5)) - &
                    D(1)*determinant4d(A(2:5),B(2:5),C(2:5),E(2:5)) + &
                    E(1)*determinant4d(A(2:5),B(2:5),C(2:5),D(2:5)) 
  END FUNCTION determinant5d


  SUBROUTINE AddNodeToList(Nodelist,VertexNr,counter,size)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER                   :: size, VertexNr, counter
    INTEGER                   :: Nodelist(size)
    ! Local variable declaration
    INTEGER                   :: iNode
    !--------------------------------------------------------------------------
    INTENT(IN)                :: VertexNr,size
    INTENT(INOUT)             :: counter,Nodelist
    !--------------------------------------------------------------------------
    DO iNode = 1, counter
       IF(Nodelist(iNode).EQ.VertexNr) THEN
         ! Node already exists in list, so return without doing anything
         RETURN
       ENDIF
    ENDDO
    ! Obviously the node has not been found, so add it
    counter = counter + 1
    Nodelist(counter) = VertexNr
  END SUBROUTINE AddNodeToList

  
  SUBROUTINE FindNodeInList(index,Nodelist,VertexNr,counter,size)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    INTEGER                   :: index,size, VertexNr, counter
    INTEGER                   :: Nodelist(size)
    ! Local variable declaration
    INTEGER                   :: iNode
    !--------------------------------------------------------------------------
    INTENT(IN)                :: Nodelist,VertexNr,counter,size
    INTENT(OUT)               :: index
    !--------------------------------------------------------------------------
    DO iNode = 1, counter
       IF(Nodelist(iNode).EQ.VertexNr) THEN
         ! Node found in list, so return index
         index = iNode
         RETURN
       ENDIF
    ENDDO
    ! Obviously the node has not been found, so return -1
    index = -1
  END SUBROUTINE FindNodeInList


  ! 
  ! Analytically compute all the roots of a non-degenerate second order polynomial
  ! a=coefficient(1) must be non-zero.
  ! 
  SUBROUTINE ZerosPolyO2(solution,coefficients)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    REAL                      :: coefficients(3)
    REAL                      :: solution(2)
    ! Local variable declaration
    REAL                      :: a,b,c
    REAL                      :: d
    !--------------------------------------------------------------------------
    INTENT(IN)                :: coefficients
    INTENT(OUT)               :: solution
    !--------------------------------------------------------------------------
    !
    a = coefficients(1)
    b = coefficients(2)
    c = coefficients(3)
    !
    d = b**2-4*a*c
    !
    IF(d.LT.0.) THEN
        PRINT *, ' d negative in ZerosPolyO2. No real roots found! '
        PRINT *, d
        STOP
    ENDIF
    !
    solution(1) = (-b-SQRT(d))/(2.*a)
    solution(2) = (-b+SQRT(d))/(2.*a)
    !
  END SUBROUTINE ZerosPolyO2

 
  ! 
  ! Analytically compute all the roots of a non-degenerate third order polynomial
  ! a=coefficient(1) must be non-zero.
  ! 
  SUBROUTINE ZerosPolyO3(solution,coefficients)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    REAL                      :: coefficients(4)
    COMPLEX                   :: solution(3)
    ! Local variable declaration
    COMPLEX                   :: a,b,c,d
    !--------------------------------------------------------------------------
    INTENT(IN)                :: coefficients
    INTENT(OUT)               :: solution
    !--------------------------------------------------------------------------
    !
    a = CMPLX( coefficients(1), 0. )
    b = CMPLX( coefficients(2), 0. )
    c = CMPLX( coefficients(3), 0. )
    d = CMPLX( coefficients(4), 0. )
    !
    solution(1) = 1./a*(36*c*b*a-108*d*a**2-8*b**3+12*SQRT(3.)*SQRT(4*c**3*a      &
                  -c**2*b**2-18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)**(1./3.)/6.     & 
                  -2./3.*(3*c*a-b**2)/a/(36*c*b*a-108*d*a**2-8*b**3+12*SQRT(3.)   &
                  *SQRT(4*c**3*a-c**2*b**2-18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)   &
                  **(1./3.)-b/a/3.

    solution(2) = -1/a*(36*c*b*a-108*d*a**2-8*b**3+12*SQRT(3.)*SQRT(4*c**3*a-     &
                   c**2*b**2-18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)**(1./3.)/12+    &
                   (3*c*a-b**2)/a/(36*c*b*a-108*d*a**2-8*b**3+12*SQRT(3.)*SQRT(4* &
                   c**3*a-c**2*b**2-18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)**(1./3.) &
                   /3.-b/a/3+CMPLX(0.,1./2.)*SQRT(3.)*(1/a*(36*c*b*a-108*d*a**2-  &
                   8*b**3+12*SQRT(3.)*SQRT(4*c**3*a-c**2*b**2-18*c*b*a*d+27*d**2* &
                   a**2+4*d*b**3)*a)**(1./3.)/6+2./3.*(3*c*a-b**2)/a/(36*c*b*a-   &
                   108*d*a**2-8*b**3+12*SQRT(3.)*SQRT(4*c**3*a-c**2*b**2-18*c*b*a &
                   *d+27*d**2*a**2+4*d*b**3)*a)**(1./3.))

    solution(3) = -1/a*(36*c*b*a-108*d*a**2-8*b**3+12*SQRT(3.)*SQRT(4*c**3*a-     &
                   c**2*b**2-18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)**(1./3.)/12.    &
                   +(3*c*a-b**2)/a/(36*c*b*a-108*d*a**2-8*b**3+12*SQRT(3.)*       &
                   SQRT(4*c**3*a-c**2*b**2-18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)** &
                   (1./3.)/3-b/a/3.+CMPLX(0.,-1./2.)*SQRT(3.)*(1/a*(36*c*b*a-     &
                   108*d*a**2-8*b**3+12*SQRT(3.)*SQRT(4*c**3*a-c**2*b**2-         &
                   18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)**(1./3.)/6+2./3.*(3*c*a-  &
                   b**2)/a/(36*c*b*a-108*d*a**2-8*b**3+12*SQRT(3.)*SQRT(4*c**3*a- &
                   c**2*b**2-18*c*b*a*d+27*d**2*a**2+4*d*b**3)*a)**(1./3.))
  END SUBROUTINE ZerosPolyO3


  ! 
  ! Analytically compute all the roots of a non-degenerate fourth order polynomial
  ! a=coefficient(1) must be non-zero.
  ! 
  SUBROUTINE ZerosPolyO4(solution,coefficients)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    REAL                      :: coefficients(5)
    REAL                      :: solution(4)
    ! Local variable declaration
    REAL                      :: b,c,d,e
    COMPLEX                   :: a3,b3,c3,d3
    REAL                      :: a2,b2,c2,d2
    REAL                      :: y,Ap,Am
    COMPLEX                   :: s, sqrarg
    !--------------------------------------------------------------------------
    INTENT(IN)                :: coefficients
    INTENT(OUT)               :: solution
    !--------------------------------------------------------------------------
    !
    ! Normalized coefficients of the fourth order polynomial
    !
    b = coefficients(2)/coefficients(1)
    c = coefficients(3)/coefficients(1)
    d = coefficients(4)/coefficients(1)
    e = coefficients(5)/coefficients(1)
    !
    ! Coefficients of the auxilliary third order polynomial
    !
    a3 = CMPLX(  8.             , 0. )
    b3 = CMPLX( -4.*c           , 0. )
    c3 = CMPLX( 2.*b*d-8*e      , 0. )
    d3 = CMPLX( e*(4*c-b*b)-d*d , 0. )
    !
    ! (Any) real root of the aux. third order polynomial
    !
    s = 1./a3*(36*c3*b3*a3-108*d3*a3**2-8*b3**3+12*SQRT(3.)*SQRT(4*c3**3*a3      &
        -c3**2*b3**2-18*c3*b3*a3*d3+27*d3**2*a3**2+4*d3*b3**3)*a3)**(1./3.)/6.     & 
        -2./3.*(3*c3*a3-b3**2)/a3/(36*c3*b3*a3-108*d3*a3**2-8*b3**3+12*SQRT(3.)   &
        *SQRT(4*c3**3*a3-c3**2*b3**2-18*c3*b3*a3*d3+27*d3**2*a3**2+4*d3*b3**3)*a3)   &
        **(1./3.)-b3/a3/3.
    !
    y = REAL(s)
    IF(ABS(AIMAG(s)).GT.1e-6*ABS(y)) THEN
        PRINT *, ' Real root of aux. poly not real in ZerosPolyO4! '
        PRINT *, s
        STOP
    ENDIF
    ! Aux. variables A+/-
    Ap = SQRT(8.*y+b*b-4*c) 
    Am = -Ap
    !
    ! Coefficients of an aux. second order polynomial which yields the final solution 
    !
    a2 = 1.
    b2 = 0.5*(b+Ap)
    c2 = y+(b*y-d)/Ap
    !
    d2 = b2**2-4*a2*c2
    !
    IF(d2.LT.0.) THEN
        PRINT *, ' d2 (using A+) negative in ZerosPolyO4. No real roots found! '
        PRINT *, d2
        STOP
    ENDIF
    !
    solution(1) = (-b2-SQRT(d2))/(2.*a2)
    solution(2) = (-b2+SQRT(d2))/(2.*a2)
    !
    a2 = 1.
    b2 = 0.5*(b+Am)
    c2 = y+(b*y-d)/Am
    !
    d2 = b2**2-4*a2*c2
    !
    IF(d2.LT.0.) THEN
        PRINT *, ' d2 (using A-) negative in ZerosPolyO4. No real roots found! '
        PRINT *, d2
        STOP
    ENDIF
    !
    solution(3) = (-b2-SQRT(d2))/(2.*a2)
    solution(4) = (-b2+SQRT(d2))/(2.*a2)
    !
  END SUBROUTINE ZerosPolyO4


SUBROUTINE AngelBetweenVectors(cosinus,a,b)                                                        !
!--------------------------------------------------------------------------------------------------!
                                                                                       !
!--------------------------------------------------------------------------------------------------!
   IMPLICIT NONE                                                                                   !
!--------------------------------------------------------------------------------------------------!
! argument list declaration                                                                        !
   REAL                             :: a(3),b(3), cosinus                                          !
! Local variable declaration                                                                       !
   REAL                             :: c,d,e                                                       !
!--------------------------------------------------------------------------------------------------!
   INTENT(IN)                       :: a,b                                                         !
   INTENT(OUT)                      :: cosinus                                                     !
!--------------------------------------------------------------------------------------------------!
                                                                                                   !
   c = a.skalar.b                                                                                  !
   d = a.skalar.a                                                                                  !
   e = b.skalar.b                                                                                  !
   cosinus = c/SQRT(d*e)                                                                           !
 RETURN                                                                                            !
 END SUBROUTINE AngelBetweenVectors                                                                !
                                                                                                
  !
  ! The subroutine IniSparseVector prepares the storage 
  ! and administration of a sparse vector
  !
  SUBROUTINE IniSparseVector(SparseVector, FullVector, m)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseVector)       :: SparseVector
    INTEGER                   :: m
    REAL                      :: FullVector(m)
    ! Local variable declaration
    INTEGER                   :: i,counter,iNonZero
    INTEGER, POINTER          :: TempIndexArray(:)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: FullVector,m
    INTENT(OUT)               :: SparseVector
    !--------------------------------------------------------------------------
    !
    SparseVector%m = m
    !
    ALLOCATE( TempIndexArray(SparseVector%m) )
    !
    ! Get the non-zero entries
    !
    counter = 0    
    DO i = 1, SparseVector%m
         IF(FullVector(i).NE.0.) THEN
            counter = counter + 1
            TempIndexArray(counter) = i
         ENDIF
    ENDDO
    !
    SparseVector%nNonZero = counter
    ALLOCATE( SparseVector%NonZero(     SparseVector%nNonZero) )
    ALLOCATE( SparseVector%NonZeroIndex(SparseVector%nNonZero) )
    !
    DO iNonZero = 1, SparseVector%nNonZero
      SparseVector%NonZeroIndex(iNonZero)   = TempIndexArray(iNonZero)
      SparseVector%NonZero(iNonZero)        = FullVector( SparseVector%NonZeroIndex(iNonZero) ) 
    ENDDO
    !
    DEALLOCATE(TempIndexArray)
    !
  END SUBROUTINE IniSparseVector

  !
  ! The subroutine CloseSparseVector deallocates the storage 
  ! of a sparse vector
  !
  SUBROUTINE CloseSparseVector(SparseVector)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseVector)       :: SparseVector
    !--------------------------------------------------------------------------
    INTENT(INOUT)             :: SparseVector
    !--------------------------------------------------------------------------
    !
    DEALLOCATE( SparseVector%NonZero )
    DEALLOCATE( SparseVector%NonZeroIndex )
    !
  END SUBROUTINE CloseSparseVector

  !
  ! The subroutine IniSparseMatrix prepares the storage 
  ! and administration of a sparse matrix
  !
  SUBROUTINE IniSparseMatrix(SparseMatrix, FullMatrix, n)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseMatrix)       :: SparseMatrix
    INTEGER                   :: n
    REAL                      :: FullMatrix(n,n)
    ! Local variable declaration
    INTEGER                   :: i,j,k,counter,iNonZero
    INTEGER, POINTER          :: TempIndexArray(:,:)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: n
    INTENT(OUT)               :: SparseMatrix
    !--------------------------------------------------------------------------
    !
    SparseMatrix%m = n
    SparseMatrix%n = n
    !
    ALLOCATE( TempIndexArray(2,SparseMatrix%n*SparseMatrix%n) )
    ALLOCATE( SparseMatrix%nNonZero(n) )
    SparseMatrix%nNonZero(:) = 0
    !
    ! Get the non-zero entries
    !
    counter = 0
    DO k = 1, n    
        DO j = 1, k
          DO i = 1, k
             IF((i.EQ.k).OR.(j.EQ.k)) THEN
                 IF(FullMatrix(i,j).NE.0.) THEN
                    counter = counter + 1
                    TempIndexArray(:,counter) = (/i,j/)
                 ENDIF
             ENDIF
          ENDDO
        ENDDO
        SparseMatrix%nNonZero(k) = counter
    ENDDO
    !
    ALLOCATE( SparseMatrix%NonZero(      SparseMatrix%nNonZero(n))  )
    ALLOCATE( SparseMatrix%NonZeroIndex1(SparseMatrix%nNonZero(n))  )
    ALLOCATE( SparseMatrix%NonZeroIndex2(SparseMatrix%nNonZero(n))  )
    !
    DO iNonZero = 1, SparseMatrix%nNonZero(n)
      SparseMatrix%NonZeroIndex1(iNonZero) = TempIndexArray(1,iNonZero)
      SparseMatrix%NonZeroIndex2(iNonZero) = TempIndexArray(2,iNonZero)
      SparseMatrix%NonZero(iNonZero) = FullMatrix(SparseMatrix%NonZeroIndex1(iNonZero), &
                                                  SparseMatrix%NonZeroIndex2(iNonZero)  ) 
    ENDDO
    !
    DEALLOCATE(TempIndexArray)
    !
  END SUBROUTINE IniSparseMatrix

  ! The subroutine UnpackSparseMatrix unpacks a sparse 
  ! matrix of dimension 2 into a full matrix
  !
  SUBROUTINE UnpackSparseMatrix(FullMatrix, SparseMatrix, n)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseMatrix)       :: SparseMatrix
    INTEGER                   :: n
    REAL                      :: FullMatrix(n,n)
    ! Local variable declaration
    INTEGER                   :: i,j,counter,iNonZero
    !--------------------------------------------------------------------------
    INTENT(IN)                :: n,SparseMatrix
    INTENT(OUT)               :: FullMatrix
    !--------------------------------------------------------------------------
    !
    FullMatrix = 0.
    DO iNonZero = 1, SparseMatrix%nNonZero(n)
      i = SparseMatrix%NonZeroIndex1(iNonZero) 
      j = SparseMatrix%NonZeroIndex2(iNonZero) 
      FullMatrix(i,j) = SparseMatrix%NonZero(iNonZero)
    ENDDO
    !
  END SUBROUTINE UnpackSparseMatrix  
  !
  ! The subroutine CloseSparseMatrix deallocates the storage 
  ! of a sparse matrix
  !
  SUBROUTINE CloseSparseMatrix(SparseMatrix)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseMatrix)       :: SparseMatrix
    !--------------------------------------------------------------------------
    INTENT(INOUT)               :: SparseMatrix
    !--------------------------------------------------------------------------
    !
    DEALLOCATE( SparseMatrix%NonZero       )
    DEALLOCATE( SparseMatrix%NonZeroIndex1 )
    DEALLOCATE( SparseMatrix%NonZeroIndex2 )
    DEALLOCATE( SparseMatrix%nNonZero      )
    !
  END SUBROUTINE CloseSparseMatrix

  !
  ! The subroutine IniSparseMatrix prepares the storage 
  ! and administration of a sparse matrix
  !
  SUBROUTINE IniSparseMatrix2(SparseMatrix, FullMatrix, m, n)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseMatrix2)      :: SparseMatrix
    INTEGER                   :: m,n
    REAL                      :: FullMatrix(m,n)
    ! Local variable declaration
    INTEGER                   :: i
    !--------------------------------------------------------------------------
    INTENT(IN)                :: m,n,FullMatrix
    INTENT(OUT)               :: SparseMatrix
    !--------------------------------------------------------------------------
    !
    SparseMatrix%m = m
    SparseMatrix%n = n
    !
    ALLOCATE( SparseMatrix%RowVector(m) )
    !
    DO i = 1, m
        CALL IniSparseVector(SparseMatrix%RowVector(i), FullMatrix(i,:), n)
    ENDDO
    !
  END SUBROUTINE IniSparseMatrix2

  !
  ! The subroutine CloseSparseMatrix deallocates the storage 
  ! of a sparse matrix
  !
  SUBROUTINE CloseSparseMatrix2(SparseMatrix)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseMatrix2)        :: SparseMatrix
    INTEGER                     :: i
    !--------------------------------------------------------------------------
    INTENT(INOUT)               :: SparseMatrix
    !--------------------------------------------------------------------------
    !
    DO i = 1, SparseMatrix%m
     CALL  CloseSparseVector( SparseMatrix%RowVector(i) )
    ENDDO
    !
  END SUBROUTINE CloseSparseMatrix2

  !
  ! The subroutine IniSparseTensor3 prepares the storage 
  ! and administration of a sparse tensor of dimension 3
  !
  SUBROUTINE IniSparseTensor3(SparseTensor, FullTensor, p, q, l)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseTensor3)      :: SparseTensor
    INTEGER                   :: p,q,l
    REAL                      :: FullTensor(p,q,l)
    ! Local variable declaration
    INTEGER                   :: i,j,k,counter,iNonZero
    INTEGER, POINTER          :: TempIndexArray(:,:)
    INTEGER, POINTER          :: TmpIdx3Struct(:)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: p,q,l
    INTENT(OUT)               :: SparseTensor
    !--------------------------------------------------------------------------
    !
    SparseTensor%p = p
    SparseTensor%q = q
    SparseTensor%l = l
    !
    ALLOCATE( TempIndexArray(3,SparseTensor%p*SparseTensor%q*SparseTensor%l) )
    ALLOCATE( TmpIdx3Struct(SparseTensor%l) )
    TmpIdx3Struct(:) = 0. 
    SparseTensor%nNonZero = 0
    !
    ! Get the non-zero entries
    !
    counter = 0
    DO k = 1, l
     DO j = 1, q
      DO i = 1, p
         IF(FullTensor(i,j,k).NE.0.) THEN
            counter = counter + 1
            TempIndexArray(:,counter) = (/i,j,k/)
            TmpIdx3Struct(k) = 1
         ENDIF        
      ENDDO
     ENDDO
    ENDDO
    SparseTensor%nNonZero = counter   
    SparseTensor%nIndex3  = SUM( TmpIdx3Struct(:) ) 
    !
    ALLOCATE( SparseTensor%NonZero(      SparseTensor%nNonZero)  )
    ALLOCATE( SparseTensor%NonZeroIndex1(SparseTensor%nNonZero)  )
    ALLOCATE( SparseTensor%NonZeroIndex2(SparseTensor%nNonZero)  )
    ALLOCATE( SparseTensor%NonZeroIndex3(SparseTensor%nNonZero)  )
    ALLOCATE( SparseTensor%Index3(SparseTensor%nIndex3)          )
    !
    iNonZero = 0 
    DO k = 1, l
       IF(TmpIdx3Struct(k).EQ.1) THEN
         iNonZero = iNonZero + 1
         SparseTensor%Index3(iNonZero) = k 
       ENDIF
    ENDDO
    !
    DO iNonZero = 1, SparseTensor%nNonZero
      SparseTensor%NonZeroIndex1(iNonZero) = TempIndexArray(1,iNonZero)
      SparseTensor%NonZeroIndex2(iNonZero) = TempIndexArray(2,iNonZero)
      SparseTensor%NonZeroIndex3(iNonZero) = TempIndexArray(3,iNonZero)
      SparseTensor%NonZero(iNonZero)       = FullTensor(SparseTensor%NonZeroIndex1(iNonZero), &
                                                        SparseTensor%NonZeroIndex2(iNonZero), &
                                                        SparseTensor%NonZeroIndex3(iNonZero)  ) 
    ENDDO
    !
    DEALLOCATE(TempIndexArray)
    DEALLOCATE(TmpIdx3Struct) 
    !
  END SUBROUTINE IniSparseTensor3
  !
  !
  ! The subroutine UnpackSparseTensor3 unpacks a sparse 
  ! tensor of dimension 3 into a full tensor
  !
  SUBROUTINE UnpackSparseTensor3(FullTensor, SparseTensor, p, q, l)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseTensor3)      :: SparseTensor
    INTEGER                   :: p,q,l
    REAL                      :: FullTensor(p,q,l)
    ! Local variable declaration
    INTEGER                   :: i,j,k,counter,iNonZero
    !--------------------------------------------------------------------------
    INTENT(IN)                :: p,q,l,SparseTensor
    INTENT(OUT)               :: FullTensor
    !--------------------------------------------------------------------------
    !
    FullTensor = 0.
    DO iNonZero = 1, SparseTensor%nNonZero
      i = SparseTensor%NonZeroIndex1(iNonZero) 
      j = SparseTensor%NonZeroIndex2(iNonZero) 
      k = SparseTensor%NonZeroIndex3(iNonZero) 
      FullTensor(i,j,k) = SparseTensor%NonZero(iNonZero)
    ENDDO
    !
  END SUBROUTINE UnpackSparseTensor3  
  !
  ! The subroutine CloseSparseTensor3 deallocates the storage 
  ! of a sparse tensor of dimension 3
  !
  SUBROUTINE CloseSparseTensor3(SparseTensor)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseTensor3)        :: SparseTensor
    !--------------------------------------------------------------------------
    INTENT(INOUT)               :: SparseTensor
    !--------------------------------------------------------------------------
    !
    DEALLOCATE( SparseTensor%NonZero       )
    DEALLOCATE( SparseTensor%NonZeroIndex1 )
    DEALLOCATE( SparseTensor%NonZeroIndex2 )
    DEALLOCATE( SparseTensor%NonZeroIndex3 )
    DEALLOCATE( SparseTensor%Index3        ) 
    !
  END SUBROUTINE CloseSparseTensor3
  !
  ! The subroutine IniSparseTensor3b prepares the storage 
  ! and administration of a sparse tensor of dimension 3 (b)
  !
  SUBROUTINE IniSparseTensor3b(SparseTensor, FullTensor, p, q, l)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseTensor3b)     :: SparseTensor
    INTEGER                   :: p,q,l
    REAL                      :: FullTensor(p,q,l)
    ! Local variable declaration
    INTEGER                   :: i,j,k,counter,iNonZero
    INTEGER, POINTER          :: TempIndexArray(:,:)
    INTEGER, POINTER          :: TmpIdx3Struct(:)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: p,q,l,FullTensor
    INTENT(OUT)               :: SparseTensor
    !--------------------------------------------------------------------------
    !
    SparseTensor%p = p
    SparseTensor%q = q
    SparseTensor%l = l
    !
    ALLOCATE( SparseTensor%SpSubMatrix(l) )
    DO k = 1, l
       CALL IniSparseMatrix(SparseTensor%SpSubMatrix(k), FullTensor(:,:,k), p)
    ENDDO
    !
  END SUBROUTINE IniSparseTensor3b
  !
  ! The subroutine IniSparseTensor3b prepares the storage 
  ! and administration of a sparse tensor of dimension 3 (b)
  !
  SUBROUTINE CloseSparseTensor3b(SparseTensor)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseTensor3b)     :: SparseTensor
    ! Local variable declaration
    INTEGER                   :: k
    !--------------------------------------------------------------------------
    INTENT(INOUT)             :: SparseTensor
    !--------------------------------------------------------------------------
    !
    DO k = 1, SparseTensor%l
       CALL CloseSparseMatrix(SparseTensor%SpSubMatrix(k))
    ENDDO
    DEALLOCATE( SparseTensor%SpSubMatrix )
    !
  END SUBROUTINE CloseSparseTensor3b


  !
  ! The subroutine IniSparseTensor4 prepares the storage 
  ! and administration of a sparse tensor
  !
  SUBROUTINE IniSparseTensor4(SparseTensor4,FullTensor4,p,q,l,m)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseTensor4)      :: SparseTensor4
    INTEGER                   :: p,q,l,m
    REAL                      :: FullTensor4(p,q,l,m)
    ! Local variable declaration
    INTEGER                   :: a,b,i,j,counter,iNonZero
    INTEGER, POINTER          :: TempIndexArray(:,:)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: p,q,l,m,FullTensor4
    INTENT(OUT)               :: SparseTensor4
    !--------------------------------------------------------------------------
    !
    SparseTensor4%p = p
    SparseTensor4%q = q
    SparseTensor4%l = l
    SparseTensor4%m = m
    !
    ALLOCATE( TempIndexArray(4,p*q*l*m) )
    !
    ! Get the non-zero entries
    !
    counter = 0    
    DO b = 1, SparseTensor4%m
     DO a = 1, SparseTensor4%l
      DO j = 1, SparseTensor4%q
       DO i = 1, SparseTensor4%p
          IF(FullTensor4(i,j,a,b).NE.0.) THEN
             counter = counter + 1
             TempIndexArray(:,counter) = (/i,j,a,b/)
          ENDIF
       ENDDO
      ENDDO
     ENDDO
    ENDDO
    !
    SparseTensor4%nNonZero = counter
    ALLOCATE( SparseTensor4%NonZero(       SparseTensor4%nNonZero) )
    ALLOCATE( SparseTensor4%NonZeroIndex(4,SparseTensor4%nNonZero) )
    !
    DO iNonZero = 1, SparseTensor4%nNonZero
      SparseTensor4%NonZeroIndex(:,iNonZero) = TempIndexArray(:,iNonZero)
      SparseTensor4%NonZero(iNonZero) = FullTensor4(SparseTensor4%NonZeroIndex(1,iNonZero), &
                                                    SparseTensor4%NonZeroIndex(2,iNonZero), &
                                                    SparseTensor4%NonZeroIndex(3,iNonZero), &
                                                    SparseTensor4%NonZeroIndex(4,iNonZero)  )
    ENDDO
    !
    DEALLOCATE(TempIndexArray)
    !
  END SUBROUTINE IniSparseTensor4

  !
  ! The subroutine CloseSparseMatrix deallocates the storage 
  ! of a sparse matrix
  !
  SUBROUTINE CloseSparseTensor4(SparseTensor4)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE(tSparseTensor4)        :: SparseTensor4
    !--------------------------------------------------------------------------
    INTENT(INOUT)               :: SparseTensor4
    !--------------------------------------------------------------------------
    !
    DEALLOCATE( SparseTensor4%NonZero )
    DEALLOCATE( SparseTensor4%NonZeroIndex )
    !
  END SUBROUTINE CloseSparseTensor4


  !
  ! The subroutine MatrixInverse3x3 computes the
  ! inverse of a 3 x 3 matrix analytically 
  !
  SUBROUTINE MatrixInverse3x3(iA,A)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    REAL                        :: iA(3,3)
    REAL                        :: A(3,3)
    ! Local variable declaration
    REAL                        :: t4,t6,t8,t10,t12,t14,t17
    !--------------------------------------------------------------------------
    INTENT(IN)                  :: A
    INTENT(OUT)                 :: iA
    !--------------------------------------------------------------------------
    !
    t4 = A(1,1)*A(2,2)      
    t6 = A(1,1)*A(2,3)
    t8 = A(1,2)*A(2,1)
    t10 = A(1,3)*A(2,1)
    t12 = A(1,2)*A(3,1)
    t14 = A(1,3)*A(3,1)
    t17 = 1/(t4*A(3,3)-t6*A(3,2)-t8*A(3,3)+t10*A(3,2)+t12*A(2,3)-t14*A(2,2))
    iA(1,1) = (A(2,2)*A(3,3)-A(2,3)*A(3,2))*t17
    iA(1,2) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))*t17
    iA(1,3) = (A(1,2)*A(2,3)-A(1,3)*A(2,2))*t17
    iA(2,1) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))*t17
    iA(2,2) = (A(1,1)*A(3,3)-t14)*t17
    iA(2,3) = -(t6-t10)*t17
    iA(3,1) = (A(2,1)*A(3,2)-A(2,2)*A(3,1))*t17
    iA(3,2) = -(A(1,1)*A(3,2)-t12)*t17
    iA(3,3) = (t4-t8)*t17
    !
  END SUBROUTINE MatrixInverse3x3


  !
  ! The subroutine MatrixInverse2x2 computes the
  ! inverse of a 2 x 2 matrix analytically 
  !
  SUBROUTINE MatrixInverse2x2(iA,A)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    REAL                        :: iA(2,2)
    REAL                        :: A(2,2)
    ! Local variable declaration
    REAL                        :: det
    !--------------------------------------------------------------------------
    INTENT(IN)                  :: A
    INTENT(OUT)                 :: iA
    !--------------------------------------------------------------------------
    !
    det = A(1,1)*A(2,2)-A(2,1)*A(1,2) 
    
    iA(1,1) =  A(2,2) / det;        iA(1,2) = -A(1,2) / det;
    iA(2,1) = -A(2,1) / det;        iA(2,2) =  A(1,1) / det;
    !
  END SUBROUTINE MatrixInverse2x2

  !
  ! Time-frequency misfits computed following Kristekova et al. 2006,
  ! taken from the TF-MISFITS software, available at www.nuquake.eu.
  !
  SUBROUTINE TF_MISFITS(S,S_REF,dt,mt,nf_tf,fmin,fmax,TFEM,TFPM,CWT_REF,TEM,TPM,TIEM,TIPM,EM,PM)
    !--------------------------------------------------------------------------
    
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    INTEGER , PARAMETER     :: KN = 16
    INTEGER , PARAMETER     :: NF = 2**KN
    REAL,     PARAMETER     :: PI = 3.14159265359
    REAL,     PARAMETER     :: W0 = 6. 
    REAL                :: dt, fmin, fmax        ! Timestep and frequency range of signals
    INTEGER             :: mt, nf_tf             ! Number of time samples and frequency samples
    REAL                :: S(mt), S_REF(mt)      ! Input signal and reference signal
    !--------------------------------------------------------------------------
    INTEGER             :: I,J, L,M
    REAL                :: NOMSE, NOMSP, DENOMS, DF, FF, MAXDENOM, MAXTF 
    REAL                :: DE(1:MT,1:NF_TF)  , DP(1:MT,1:NF_TF) 
    REAL                :: DET(1:MAX(MT,NF_TF)), DPT(1:MAX(MT,NF_TF)), DEF(1:MAX(MT,NF_TF)), DPF(1:MAX(MT,NF_TF))                                        
    COMPLEX             :: WV(1:MT,1:NF_TF), WV_REF(1:MT,1:NF_TF)
    COMPLEX             :: SS(NF), SS_REF(NF)
    REAL                :: SS_tmp(2*NF), SS_REF_tmp(2*NF)
    REAL                :: EM,   PM                ! Single-valued envelope and phase misfits
    REAL, OPTIONAL      :: TEM(1:MT),  TPM(1:MT)   ! Time dependent envelope and phase misfits
    REAL, OPTIONAL      :: TIEM, TIPM              ! Time integrated envelope and phase misfits
    REAL                :: TFEM(1:MT,1:NF_TF), TFPM(1:MT,1:NF_TF), CWT_REF(1:MT,1:NF_TF) 
    !--------------------------------------------------------------------------
    INTENT(IN)          :: S, S_REF, dt,mt,nf_tf,fmin,fmax
    INTENT(INOUT)       :: TFEM,TFPM,CWT_REF,TEM,TPM,TIEM,TIPM,EM,PM
    !--------------------------------------------------------------------------

    ! COMPUTATION OF CWT
    M=NF
    DF  = 1./(DT*REAL(NF))
    FF = EXP( LOG(FMAX/FMIN) / REAL(NF_TF - 1) )
   
    SS     = CMPLX(0.,0.)
    SS_REF = CMPLX(0.,0.)
    DO I = 1, MT
      SS     (I) = CMPLX( S   (I), 0. )
      SS_REF (I) = CMPLX( S_REF(I), 0. )
    END DO

    !Modified input for FTT (must be real!)
    DO I=1,NF
      SS_tmp(2*I-1)     = REAL(SS(I))
      SS_tmp(2*I)       = 0.0
      SS_REF_tmp(2*I-1) = REAL(SS_REF(I))
      SS_REF_tmp(2*I)   = 0.0
    ENDDO
  
    CALL FCOOLR (KN, SS_tmp    , -1.)
    CALL FCOOLR (KN, SS_REF_tmp, -1.)

    DO I=1,NF
      SS(I)=CMPLX(SS_tmp(2*I-1), SS_tmp(2*I))
      SS_REF(I)=CMPLX(SS_REF_tmp(2*I-1), SS_REF_tmp(2*I))
    ENDDO

    !WV     (1:MT,1:NF_TF) = CWT(SS    , MT ,NF_TF, DF, FF, FMIN)
    !WV_REF (1:MT,1:NF_TF) = CWT(SS_REF, MT, NF_TF, DF, FF, FMIN)
    CALL CWT(WV    ,  SS    , MT ,NF_TF, DF, FF, FMIN)
    CALL CWT(WV_REF,  SS_REF, MT ,NF_TF, DF, FF, FMIN)
    
    CWT_REF = ABS(WV_REF)
  
    ! COMPUTATION OF TFEM AND TFPM

    MAXTF = MAXVAL(ABS(WV_REF(1:MT,1:NF_TF)))
        
    DE = ( ABS(WV) - ABS(WV_REF) )
    !WV_REF=WV_REF+CMPLX(1.0e-20,1.0e-20)
    DP = ABS(WV_REF) * (ATAN2(IMAG(WV/WV_REF),REAL(WV/WV_REF)))/PI
    TFEM = DE / MAXTF
    TFPM = DP / MAXTF
        
    ! COMPUTATION OF AUXILIARY VARIABLES

    DO I = 1, MT
      DET  (I) = 0.
      DPT  (I) = 0.
      DO L = 1, NF_TF
        DET (I) = DET (I) +     DE(I,L)
        DPT (I) = DPT (I) +     DP(I,L)
      END DO
      DET (I) = DET (I)/REAL(NF_TF)
      DPT (I) = DPT (I)/REAL(NF_TF)
    END DO

    DO L = 1, NF_TF
      DEF   (L) = 0.
      DPF   (L) = 0.
      DO I = 1, MT
        DEF (L) = DEF (L) +     DE(I,L)
        DPF (L) = DPF (L) +     DP(I,L)
      END DO
      DEF (L) = DEF (L)/REAL(MT)
      DPF (L) = DPF (L)/REAL(MT)
    END DO  
  
    MAXDENOM = 0.
    DO I = 1, MT
      DENOMS = 0.
      DO L = 1, NF_TF
        DENOMS = DENOMS + ABS(WV_REF(I,L))
      END DO
      DENOMS = DENOMS / REAL(NF_TF)
      IF ( DENOMS > MAXDENOM ) MAXDENOM = DENOMS
    END DO

    ! COMPUTATION OF TEM
    IF(PRESENT(TEM)) THEN
      DO I = 1, MT
        TEM(I) = DET (I) / MAXDENOM
      END DO
    ENDIF

    ! COMPUTATION OF TPM
    IF(PRESENT(TPM)) THEN
      DO I = 1, MT
        TPM(I) = DPT (I) / MAXDENOM
      END DO
    ENDIF

    ! COMPUTATION OF TIPM
    IF(PRESENT(TIPM)) THEN
      TIPM = 0.0
      DO I = 1, MT
        TIPM = TIPM + DT*TPM(I)
      END DO
    ENDIF

    ! COMPUTATION OF TIEM
    IF(PRESENT(TIPM)) THEN
      TIEM = 0.0
      DO I = 1, MT
        TIEM = TIEM + DT*TEM(I)
      END DO
    ENDIF
  
    ! COMPUTATION OF EM AND PM
 
    NOMSE  = 0.
    NOMSP  = 0.
    DENOMS = 0.

    DO I = 1, MT
      DO L = 1, NF_TF

        NOMSE  = NOMSE  + ABS(DE    (I,L))*ABS(DE    (I,L))
        NOMSP  = NOMSP  + ABS(DP    (I,L))*ABS(DP    (I,L))
        DENOMS = DENOMS + ABS(WV_REF(I,L))*ABS(WV_REF(I,L))

      END DO
    END DO
  
    EM = SQRT ( NOMSE / DENOMS )
    PM = SQRT ( NOMSP / DENOMS )

  END SUBROUTINE TF_MISFITS
  !
  ! Fast Fourier transform algorithm
  !
  SUBROUTINE FCOOLR(K,D,SN)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    INTEGER                  :: LX, K, IL, I, NKK, LA, NCK, LCK, L2K, NW, ICK, LS
    INTEGER                  :: J1, J, JH, JH1, ID, JJ
    REAL,          PARAMETER :: PI  = 3.141592654, PI2 = 6.283185308
    REAL                     :: SN, SH, Q1, Q2, FNW, AA, W1, W2
    INTEGER , DIMENSION (32) :: INU
    REAL    , DIMENSION ( *) :: D
    !--------------------------------------------------------------------------
    LX = 2**K
    Q1 = LX
    IL = LX
    SH = SN * PI2 / Q1
    DO I = 1, K
      IL     = IL / 2
      INU(I) = IL
    END DO
    NKK = 1
    DO LA = 1, K
      NCK = NKK
      NKK = NKK + NKK
      LCK = LX  / NCK
      L2K = LCK + LCK
      NW = 0
      DO ICK = 1, NCK
        FNW = NW
        AA  = SH * FNW
        W1  = COS(AA)
        W2  = SIN(AA)
        LS  = L2K * (ICK-1)
        DO I = 2, LCK, 2
          J1  = I  + LS
          J   = J1 - 1
          JH  = J  + LCK
          JH1 = JH + 1
          Q1  = D(JH)*W1 - D(JH1)*W2
          Q2  = D(JH)*W2 + D(JH1)*W1
          D(JH ) = D(J ) - Q1
          D(JH1) = D(J1) - Q2
          D(J  ) = D(J ) + Q1
          D(J1 ) = D(J1) + Q2
        END DO
        DO I = 2, K
          ID = INU(I)
          IL = ID + ID
          IF ( (NW-ID-IL*(NW/IL)) < 0 ) EXIT
          NW = NW - ID
        END DO
        NW = NW + ID
      END DO
    END DO
    NW = 0
    DO J = 1, LX
      IF ( (NW-J) >= 0 ) THEN 
        JJ  = NW  + NW + 1
        J1  = JJ  + 1
        JH1 = J   + J
        JH  = JH1 - 1
        Q1  = D(JJ)
        D(JJ ) = D(JH )
        D(JH ) = Q1
        Q1     = D(J1 )
        D(J1 ) = D(JH1)
        D(JH1) = Q1
      END IF
      DO I = 1, K
        ID = INU(I)
        IL = ID + ID
        IF ( (NW-ID-IL*(NW/IL)) < 0 ) EXIT
        NW = NW - ID
      END DO
      NW = NW + ID
    END DO
       
    RETURN
 
  END SUBROUTINE FCOOLR
  !
  ! Continuous Wavelet Transform, for wavelet analysis
  !
  SUBROUTINE CWT(WV,F_V, MT, NF_TF, DFA, FF, FMIN)
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    !
    INTEGER , PARAMETER     :: KN = 16
    INTEGER , PARAMETER     :: NF = 2**KN
    REAL,     PARAMETER     :: PI = 3.14159265359
    REAL,     PARAMETER     :: W0 = 6. 
    INTEGER , INTENT(IN)    :: MT, NF_TF
    REAL,     INTENT(IN)    :: DFA, FF, FMIN
    COMPLEX , INTENT(INOUT), DIMENSION(NF) :: F_V
    INTEGER    :: I,J, NF21
    REAL       :: S, F
    COMPLEX, DIMENSION (1:MT, 1:NF_TF) :: WV
    COMPLEX, DIMENSION (1:NF)          :: FPWV, FWV
    REAL, DIMENSION(1:2*NF)            :: FWV_tmp
    !--------------------------------------------------------------------------

    NF21 = NF/2 + 1
  
    F = FMIN/FF

    DO I = 1, NF_TF

      F = F*FF
      S = W0 / (2.*PI*F)

      DO J = 1, NF21
        FPWV(J) = MORLET( S, REAL(J-1)*DFA )
      END DO

      DO  J = NF21+1, NF
        FPWV(J) = (0.0,0.0)
      END DO

      DO J = 1, NF
        FWV (J) = F_V (J)*CONJG( FPWV(J) )
      END DO

      !Modified input for FTT (must be real!)
      DO J=1,NF
        FWV_tmp(2*J-1) = REAL(FWV(J))
        FWV_tmp(2*J)   = 0.0
      ENDDO

      CALL FCOOLR (KN,FWV_tmp ,1.0)

      DO J=1,NF
        FWV(J)=CMPLX(FWV_tmp(2*J-1), FWV_tmp(2*J))
      ENDDO

      WV(1:MT,I) = FWV(1:MT)*DFA

    END DO
  
  END SUBROUTINE CWT
  !
  ! Morlet wavelet, used as analyzing wavelet in CWT
  !
  FUNCTION MORLET ( S, FA )
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    REAL,     PARAMETER  :: PI = 3.14159265359
    REAL,     PARAMETER  :: W0 = 6. 
    REAL, PARAMETER      :: PI14 =  0.7511255444649
    REAL,     INTENT(IN) :: S, FA
    COMPLEX              :: MORLET
    !--------------------------------------------------------------------------
    IF (FA == 0.) THEN 
      MORLET = (0.,0.)
    ELSE
      MORLET = CMPLX(PI14*EXP( (-(S*2.*PI*FA-W0)**2) / 2.),0.)
    ENDIF

  END FUNCTION MORLET
  !
  ! Timer implementation based on system_clock intrinsic, which is also working on IBM machines (c) R. Baader
  !  
  REAL FUNCTION dwalltime()
    INTEGER(ik) :: count
    IF (first) THEN
       first = .false.
       CALL system_clock(count, count_rate, count_max)
       conversion = 1.0d0 / dble(count_rate)
    ELSE
       CALL system_clock(count)
    END IF
    dwalltime = count * conversion
  END FUNCTION dwalltime


  FUNCTION XYinTriangle(x,y,xtri,ytri,epsilon)
     IMPLICIT NONE
    ! Function result type
     logical :: XYinTriangle
     REAL x,y,xtri(3),ytri(3),epsilon
     REAL xi, eta
     REAL volume
     INTENT(in) :: x,y,xtri,ytri,epsilon

     ! Compute volume
     volume = 0.5*( (xtri(2)-xtri(1))*(ytri(3)-ytri(1))-(xtri(3)-xtri(1))*(ytri(2)-ytri(1)) )

     xi  = 0.5/volume*( ( xtri(3)*ytri(1) - xtri(1)*ytri(3) ) + &
         x*(ytri(3)-ytri(1)) + &
         y*(xtri(1)-xtri(3)) )
     eta = 0.5/volume*( ( xtri(1)*ytri(2) - xtri(2)*ytri(1) ) + &
         x*(ytri(1)-ytri(2)) + &
         y*(xtri(2)-xtri(1)) )

     ! Because of numerical errors, it is possible that a point which lies     !
     ! exactly on the boundary is considered outside the element               !
     ! So we set a tolerance value of epsilon which allows points lo lie       !
     ! slightly out of the element. This is necessary because fluxes are       !
     ! calculated on boundary points!                                          !
     ! For epsilon we choose 1e-5                                              !

     if((xi.lt.(0.-epsilon)).or.(eta.lt.(0.-epsilon)).or. &
         (eta.gt.(1.-xi+epsilon))) then
         XYinTriangle = .false.
     else
         XYinTriangle = .true.
     endif

   END FUNCTION XYinTriangle

END MODULE COMMON_operators_mod
