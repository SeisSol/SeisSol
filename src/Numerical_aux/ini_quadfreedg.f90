!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
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

MODULE ini_QuadFreeDG_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  INTERFACE IniQuadFreeRKDG
     MODULE PROCEDURE IniQuadFreeRKDG
  END INTERFACE
  INTERFACE IniQuadFreeADERDG
     MODULE PROCEDURE IniQuadFreeADERDG
  END INTERFACE
  INTERFACE IniQuadFreeRKDG3D
     MODULE PROCEDURE IniQuadFreeRKDG3D
  END INTERFACE
  INTERFACE IniQuadFreeADERDG3D
     MODULE PROCEDURE IniQuadFreeADERDG3D
  END INTERFACE
  !---------------------------------------------------------------------------!
  PUBLIC :: IniQuadFreeRKDG
  PUBLIC :: IniQuadFreeADERDG
  PUBLIC :: IniQuadFreeRKDG3D
  PUBLIC :: IniQuadFreeADERDG3D
  !---------------------------------------------------------------------------!

CONTAINS
  
    SUBROUTINE IniQuadFreeRKDG(FMatrix,Kxi,Keta,nDegFr,nPoly,nDegFrRec,nPolyRec,GlobalElemType,IO) 
    !---------------------------------------------------------------------------!
    USE COMMON_operators_mod
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tInputOutput)   :: IO
    INTEGER              :: nDegFr, nDegFrRec
    INTEGER              :: nPoly, nPolyRec
    INTEGER              :: GlobalElemType
    REAL                 :: FMatrix(nDegFr,nDegFrRec,0:GlobalElemType,1:GlobalElemType)
    REAL                 :: Kxi(nDegFr,nDegFrRec)
    REAL                 :: Keta(nDegFr,nDegFrRec)
    ! Local variable declaration
    INTEGER              :: i, j, k, l, m, iDegFr
    INTEGER              :: STAT
    REAL                 :: dummy
    CHARACTER(LEN=600)   :: FileName, SectionString
    CHARACTER(LEN=200)   :: DGPATH
    LOGICAL              :: configexist
    !-------------------------------------------------------------------------!
    INTENT(IN)  :: nDegFr, nDegFrRec, nPoly, nPolyRec, IO
    INTENT(OUT) :: FMatrix, Kxi, Keta
    !-------------------------------------------------------------------------!

    INQUIRE(                                            & !
         FILE= 'DGPATH'                               , & !
         EXIST=configexist                              ) !
    !                                                     !
    IF (configexist) THEN                                 !
       !                                                  !       
       OPEN(                                            & !
            UNIT= IO%UNIT%FileIn                      , & !
            FILE= 'DGPATH'                          , & !
            IOSTAT = STAT                               ) !
       !                                                  !
       IF (stat.NE.0) THEN                                !
          logError(*) 'cannot open DGPATH'   !
          STOP
       END IF                                             !
       !                                                  !
       READ(IO%UNIT%FileIn,'(A)') DGPATH                  !
       !                                                  !
       CLOSE(IO%UNIT%FileIn)                              !
       !                                                  !          
       !                                                  !
       logInfo(*) 'Path to the DG directory is: ',TRIM(DGPATH)
       IF (GlobalElemType.EQ.3) THEN ! Triangles          !
          WRITE(FileName,'(a,a13,i2.2,a4)') TRIM(DGPATH), 'QFDGMatricesP',nPolyRec,'.dat'
       ELSE                          ! Rectangles         !
          WRITE(FileName,'(a,a13,i2.2,a7)') TRIM(DGPATH), 'QFDGMatricesP',nPolyRec,'qua.dat'
       ENDIF                                              !

    ELSE                                                  !
       !                                                  !
       logWarning(*) 'Configuration file DGPATH missing!'
       logWarning(*) 'Use . as default path'    !
       !
       IF (GlobalElemType.EQ.3) THEN ! Triangles          !
          WRITE(FileName,'(a13,i2.2,a4)') 'QFDGMatricesP',nPolyRec,'.dat'
       ELSE                          ! Rectangles         !
          WRITE(FileName,'(a13,i2.2,a7)') 'QFDGMatricesP',nPolyRec,'qua.dat'
       ENDIF                                              !
       !                                                  ! 
    END IF                                                !

    CALL OpenFile(IO%UNIT%FileIn,FileName,.FALSE.)

    logInfo('(a52,i2.2,a12)') 'Reading matrices for quadrature-free DG method, P', nPolyRec, ' elements.'
    IF (GlobalElemType.EQ.4) THEN ! TRectangles
          logWarning(*) 'QuadFree linear DG-Method on rectangles highly experimental!'
    ENDIF

    FMatrix(:,:,:,:)   = 0.
    Kxi(:,:)           = 0.
    Keta(:,:)          = 0.

    ! Coefficients of the basis polynomials
    !READ(IO%UNIT%FileIn,*) SectionString
    !WRITE(IO%UNIT%stdOut,*) '|    Reading section 1/7 '
    !DO iDegFr = 0, nDegFrRec-1
    ! DO j = 0, nPolyRec
    !  DO i = 0, nPolyRec
    !    READ(IO%UNIT%FileIn,*) dummy ! cPoly(i,j,iDegFr,nPoly) 
    !  ENDDO
    ! ENDDO
    !ENDDO

    ! Entries of the mass matrix
    !READ(IO%UNIT%FileIn,*) SectionString
    !WRITE(IO%UNIT%stdOut,*) '|    Reading section 2/7 '
    !DO l = 1, nDegFrRec
    ! DO m = 1, nDegFrRec
    !   READ(IO%UNIT%FileIn,*) dummy ! MassMatrix(l,m,nPoly) 
    ! ENDDO
    !ENDDO

    IF (GlobalElemType.EQ.4) THEN ! Rectangles

      ! Entries of the Kxi stiffness matrix
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 1/6 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           Kxi(l,m) = dummy
         ENDIF
       ENDDO
      ENDDO

      ! Entries of the Keta stiffness matrix
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 2/6 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           Keta(l,m) = dummy
         ENDIF 
       ENDDO
      ENDDO

      ! Entries of the fluxmatrices for side 1
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 3/6 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
        DO k = 0, 4
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           FMatrix(l,m,k,1) = dummy
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      ! Entries of the fluxmatrices for side 2
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 4/6 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
        DO k = 0, 4
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           FMatrix(l,m,k,2) = dummy
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      ! Entries of the fluxmatrices for side 3
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 5/6 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
        DO k = 0, 4
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           FMatrix(l,m,k,3) = dummy
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      ! Entries of the fluxmatrices for side 4
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 6/6 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
        DO k = 0, 4
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           FMatrix(l,m,k,4) = dummy !!! Hier ist der Wurm drin !!!
         ENDIF
        ENDDO
       ENDDO
      ENDDO

    ELSE ! Triangles

      ! Entries of the Kxi stiffness matrix
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 1/5 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           Kxi(l,m) = dummy
         ENDIF
       ENDDO
      ENDDO

      ! Entries of the Keta stiffness matrix
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 2/5 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           Keta(l,m) = dummy
         ENDIF 
       ENDDO
      ENDDO

      ! Entries of the fluxmatrices for side 1
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 3/5 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
        DO k = 0, 3
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           FMatrix(l,m,k,1) = dummy
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      ! Entries of the fluxmatrices for side 2
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 4/5 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
        DO k = 0, 3
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           FMatrix(l,m,k,2) = dummy
         ENDIF
        ENDDO
       ENDDO
      ENDDO

      ! Entries of the fluxmatrices for side 3
      READ(IO%UNIT%FileIn,*) SectionString
      logInfo(*) 'Reading section 5/5 '
      DO l = 1, nDegFrRec
       DO m = 1, nDegFrRec
        DO k = 0, 3
         READ(IO%UNIT%FileIn,*) dummy
         IF(l.LE.nDegFr) THEN
           FMatrix(l,m,k,3) = dummy
         ENDIF
        ENDDO
       ENDDO
      ENDDO
    
    ENDIF

    CLOSE(IO%UNIT%FileIn)     

    logInfo(*) '  ----- Quadrature-free DG matrices read ------'
    

 END SUBROUTINE IniQuadFreeRKDG

 SUBROUTINE IniQuadFreeADERDG(Coeff_level0,nDegFr,nPoly,GlobalElemType,IO) 
    !---------------------------------------------------------------------------!
    USE COMMON_operators_mod
    IMPLICIT NONE
    !---------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tInputOutput)   :: IO
    INTEGER              :: nDegFr
    INTEGER              :: nPoly
    INTEGER              :: GlobalElemType
    REAL                 :: Coeff_level0( 0:nPoly, 0:nPoly, 1:nDegFr, 1:nDegFr )
    ! Local variable declaration
    INTEGER              :: r, r1, l, m 
    INTEGER              :: STAT
    DOUBLE PRECISION     :: ReadVariable
    CHARACTER(LEN=600)   :: FileName, SectionString
    CHARACTER(LEN=200)   :: DGPATH
    LOGICAL              :: configexist
    !-------------------------------------------------------------------------!
    INTENT(IN)  :: nDegFr, nPoly, IO
    INTENT(OUT) :: Coeff_level0
    !-------------------------------------------------------------------------!

    INQUIRE(                                            & !
         FILE= 'DGPATH'                               , & !
         EXIST=configexist                              ) !
    !                                                     !
    IF (configexist) THEN                                 !
       !                                                  !       
       OPEN(                                            & !
            UNIT= IO%UNIT%FileIn                      , & !
            FILE= 'DGPATH'                            , & !
            IOSTAT = STAT                               ) !
       !                                                  !
       IF (stat.NE.0) THEN                                !
          logError(*) 'cannot open DGPATH'   !
          STOP
       END IF                                             !
       !                                                  !
       READ(IO%UNIT%FileIn,'(A)') DGPATH                  !
       !                                                  !
       CLOSE(IO%UNIT%FileIn)                              !
       !                                                  !          
       !                                                  !
       logInfo(*) 'Path to the DG directory is: ',TRIM(DGPATH)
       IF (GlobalElemType.EQ.7) THEN ! Triangles          !
          WRITE(FileName,'(a,a15,i2.2,a4)') TRIM(DGPATH), 'ADERDGMatricesP',nPoly,'.dat'
       ELSE                          ! Rectangles         !
          WRITE(FileName,'(a,a15,i2.2,a7)') TRIM(DGPATH), 'ADERDGMatricesP',nPoly,'qua.dat'
       ENDIF                                              !
       !                                                  !
    ELSE                                                  !
       !                                                  !
       logWarning(*) 'Configuration file DGPATH missing!'
       logWarning(*) 'Use . as default path'
       !                                                  !
       IF (GlobalElemType.EQ.3) THEN ! Triangles          !
          WRITE(FileName,'(a15,i2.2,a4)') 'ADERDGMatricesP',nPoly,'.dat'
       ELSE                          ! Rectangles         !
          WRITE(FileName,'(a15,i2.2,a7)') 'ADERDGMatricesP',nPoly,'qua.dat'
       ENDIF                                              !
       !                                                  ! 
    END IF                                                !

    CALL OpenFile(IO%UNIT%FileIn,FileName,.FALSE.)

    logInfo('(a41,i2.2,a12)') 'Reading matrices for ADER-DG method, P', nPoly, ' elements.'
  

    Coeff_level0(:,:,:,:) = 0.

    ! Entries of the ADER Clm matrices
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 1/1 '

    DO l = 1, nDegFr
     DO m = 1, nDegFr
      DO r = 0, nPoly
       DO r1 = r, 0, -1
       READ(IO%UNIT%FileIn,*) ReadVariable
       Coeff_level0(r,r1,l,m) = ReadVariable 
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    CLOSE(IO%UNIT%FileIn)     

    logInfo(*) '----- ADER-DG matrices read ------'

 END SUBROUTINE IniQuadFreeADERDG

 
 
 
 ! -------------------------------------------------------------------------- !
 !                               3D Functions                                 !
 ! -------------------------------------------------------------------------- !

    SUBROUTINE IniQuadFreeRKDG3D(FMatrix3D,Kxi,Keta,Kzeta,nDegFr,nPoly,nDegFrRec,nPolyRec,IO ) 
    !-------------------------------------------------------------------------!
    USE COMMON_operators_mod
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tInputOutput)   :: IO
    INTEGER              :: nDegFr
    INTEGER              :: nPoly
    INTEGER              :: nDegFrRec
    INTEGER              :: nPolyRec
    REAL                 :: FMatrix3D(nDegFr,nDegFrRec,0:4,1:3,1:4)
    REAL                 :: Kxi(  nDegFr,nDegFrRec)
    REAL                 :: Keta( nDegFr,nDegFrRec)
    REAL                 :: Kzeta(nDegFr,nDegFrRec)
    ! Local variable declaration
    INTEGER              :: i, j, k, l, m, iDegFr
    INTEGER              :: STAT
    REAL                 :: dummy
    CHARACTER(LEN=600)   :: FileName, SectionString
    CHARACTER(LEN=200)   :: DGPATH
    LOGICAL              :: configexist
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: nDegFr, nPoly, nDegFrRec, nPolyRec, IO
    INTENT(OUT)   :: FMatrix3D, Kxi, Keta, Kzeta
    !-------------------------------------------------------------------------!

    INQUIRE(                                            & !
         FILE= 'DGPATH'                               , & !
         EXIST=configexist                              ) !
    !                                                     !
    IF (configexist) THEN                                 !
       !                                                  !       
       OPEN(                                            & !
            UNIT= IO%UNIT%FileIn                      , & !
            FILE= 'DGPATH'                          , & !
            IOSTAT = STAT                               ) !
       !                                                  !
       IF (stat.NE.0) THEN                                !
          logError(*) 'cannot open DGPATH'   !
          STOP
       END IF                                             !
       !                                                  !
       READ(IO%UNIT%FileIn,'(A)') DGPATH                !
       !                                                  !
       CLOSE(IO%UNIT%FileIn)                              !
       !                                                  !          
       !                                                  !
       logInfo(*) 'Path to the DG directory is: ',TRIM(DGPATH)
       WRITE(FileName,'(a,a15,i2.2,a4)') TRIM(DGPATH), 'QFDG3DMatricesP',nPolyRec,'.dat'

    ELSE                                                  !
       !                                                  !
       logWarning(*) 'Configuration file DGPATH missing!'
       logWarning(*) 'Use . as default path'    !
       !                                                  !
       WRITE(FileName,'(a15,i2.2,a4)') 'QFDG3DMatricesP',nPolyRec,'.dat'
       !                                                  ! 
    END IF                                                !

    CALL OpenFile(IO%UNIT%FileIn,FileName,.FALSE.)

    logInfo('(a52,i2.2,a12)') 'Reading matrices for quadrature-free DG method, P', nPolyRec, ' elements.'

    FMatrix3D(:,:,:,:,:)    = 0.
    Kxi(:,:)                = 0.
    Keta(:,:)               = 0.
    Kzeta(:,:)              = 0.

    ! Entries of the Kxi stiffness matrix
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 1/7 '
    DO l = 1, nDegFrRec
     DO m = 1, nDegFrRec
       IF(l.GT.nDegFr) THEN
         READ(IO%UNIT%FileIn,*) dummy
       ELSE
         READ(IO%UNIT%FileIn,*) Kxi(l,m) 
       ENDIF
     ENDDO
    ENDDO

    ! Entries of the Keta stiffness matrix
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 2/7 '
    DO l = 1, nDegFrRec
     DO m = 1, nDegFrRec
       IF(l.GT.nDegFr) THEN
         READ(IO%UNIT%FileIn,*) dummy
       ELSE
         READ(IO%UNIT%FileIn,*) Keta(l,m) 
       ENDIF
     ENDDO
    ENDDO

    ! Entries of the Kzeta stiffness matrix
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 3/7 '
    DO l = 1, nDegFrRec
     DO m = 1, nDegFrRec
       IF(l.GT.nDegFr) THEN
         READ(IO%UNIT%FileIn,*) dummy
       ELSE
         READ(IO%UNIT%FileIn,*) Kzeta(l,m) 
       ENDIF
     ENDDO
    ENDDO

    ! Entries of the fluxmatrices for side 1
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 4/7 '
    ! Flux contribution of the element itself
    DO l = 1, nDegFrRec
     DO m = 1, nDegFrRec
       IF(l.GT.nDegFr) THEN
         READ(IO%UNIT%FileIn,*) dummy
       ELSE
         READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,0,1,1)
       ENDIF
     ENDDO
    ENDDO
    DO k = 1, 4
     DO j = 1, 3
       DO l = 1, nDegFrRec
        DO m = 1, nDegFrRec
           IF(l.GT.nDegFr) THEN
             READ(IO%UNIT%FileIn,*) dummy
           ELSE
             READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,k,j,1)
           ENDIF
        ENDDO
       ENDDO
     ENDDO
    ENDDO

    ! Entries of the fluxmatrices for side 2
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 5/7 '
    ! Flux contribution of the element itself
    DO l = 1, nDegFrRec
     DO m = 1, nDegFrRec
        IF(l.GT.nDegFr) THEN
          READ(IO%UNIT%FileIn,*) dummy
        ELSE
          READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,0,1,2)
        ENDIF
     ENDDO
    ENDDO
    DO k = 1, 4
     DO j = 1, 3
       DO l = 1, nDegFrRec
        DO m = 1, nDegFrRec
            IF(l.GT.nDegFr) THEN
              READ(IO%UNIT%FileIn,*) dummy
            ELSE
              READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,k,j,2)
            ENDIF
        ENDDO
       ENDDO
     ENDDO
    ENDDO

    ! Entries of the fluxmatrices for side 3
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 6/7 '
    ! Flux contribution of the element itself
    DO l = 1, nDegFrRec
     DO m = 1, nDegFrRec
        IF(l.GT.nDegFr) THEN
          READ(IO%UNIT%FileIn,*) dummy
        ELSE
          READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,0,1,3)
        ENDIF
     ENDDO
    ENDDO
    DO k = 1, 4
     DO j = 1, 3
       DO l = 1, nDegFrRec
        DO m = 1, nDegFrRec
          IF(l.GT.nDegFr) THEN
            READ(IO%UNIT%FileIn,*) dummy
          ELSE
            READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,k,j,3)
          ENDIF
        ENDDO
       ENDDO
     ENDDO
    ENDDO

    ! Entries of the fluxmatrices for side 4
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 7/7 '
    ! Flux contribution of the element itself
    DO l = 1, nDegFrRec
     DO m = 1, nDegFrRec
        IF(l.GT.nDegFr) THEN
          READ(IO%UNIT%FileIn,*) dummy
        ELSE
          READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,0,1,4)
        ENDIF
     ENDDO
    ENDDO
    DO k = 1, 4
     DO j = 1, 3
       DO l = 1, nDegFrRec
        DO m = 1, nDegFrRec
           IF(l.GT.nDegFr) THEN
             READ(IO%UNIT%FileIn,*) dummy
           ELSE
             READ(IO%UNIT%FileIn,*) FMatrix3D(l,m,k,j,4)
           ENDIF
        ENDDO
       ENDDO
     ENDDO
    ENDDO
    !
    CLOSE(IO%UNIT%FileIn)     
    !
    logInfo(*) '----- Quadrature-free DG matrices read ------'
    !
 END SUBROUTINE IniQuadFreeRKDG3D

 SUBROUTINE IniQuadFreeADERDG3D(Coeff_level03D,nDegFr,nPoly,IO) 
    !---------------------------------------------------------------------------!
    USE COMMON_operators_mod
    IMPLICIT NONE
    !---------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tInputOutput)   :: IO
    INTEGER              :: nDegFr
    INTEGER              :: nPoly
    REAL                 :: Coeff_level03D( 0:nPoly, 0:nPoly, 0:nPoly, 1:nDegFr, 1:nDegFr )
    ! Local variable declaration
    INTEGER              :: r, r1, r2, l, m 
    INTEGER              :: STAT
    DOUBLE PRECISION     :: ReadVariable
    CHARACTER(LEN=600)   :: FileName, SectionString
    CHARACTER(LEN=200)   :: DGPATH
    LOGICAL              :: configexist
    !-------------------------------------------------------------------------!
    INTENT(IN)  :: nDegFr, nPoly, IO
    INTENT(OUT) :: Coeff_level03D
    !-------------------------------------------------------------------------!

    INQUIRE(                                            & !
         FILE= 'DGPATH'                               , & !
         EXIST=configexist                              ) !
    !                                                     !
    IF (configexist) THEN                                 !
       !                                                  !       
       OPEN(                                            & !
            UNIT= IO%UNIT%FileIn                      , & !
            FILE= 'DGPATH'                            , & !
            IOSTAT = STAT                               ) !
       !                                                  !
       IF (stat.NE.0) THEN                                !
          logError(*) 'cannot open DGPATH'   !
          STOP
       END IF                                             !
       !                                                  !
       READ(IO%UNIT%FileIn,'(A)') DGPATH                  !
       !                                                  !
       CLOSE(IO%UNIT%FileIn)                              !
       !                                                  !          
       !                                                  !
       logInfo(*) 'Path to the ADER directory is: ',TRIM(DGPATH)
       WRITE(FileName,'(a,a17,i2.2,a4)') TRIM(DGPATH), 'ADERDG3DMatricesP',nPoly,'.dat'

    ELSE                                                  !
       !                                                  !
       logWarning(*) 'Configuration file DGPATH missing!'
       logWarning(*) 'Use . as default path'    !
       !                                                  !
       WRITE(FileName,'(a17,i2.2,a4)') 'ADERDG3DMatricesP',nPoly,'.dat'
       !                                                  ! 
    END IF                                                !

    CALL OpenFile(IO%UNIT%FileIn,FileName,.FALSE.)

    logInfo('(a41,i2.2,a12)') 'Reading matrices for 3D ADER-DG method, P', nPoly, ' elements.'
  

    Coeff_level03D(:,:,:,:,:) = 0.

    ! Entries of the ADER Clm matrices
    READ(IO%UNIT%FileIn,*) SectionString
    logInfo(*) 'Reading section 1/1 '

    DO l = 1, nDegFr
     DO m = 1, nDegFr
      DO r = 0, nPoly
       DO r1 = r, 0, -1
        DO r2 = r-r1, 0, -1
         READ(IO%UNIT%FileIn,*) ReadVariable
         Coeff_level03D(r,r1,r2,l,m) = ReadVariable 
        ENDDO
       ENDDO
      ENDDO
     ENDDO
    ENDDO

    CLOSE(IO%UNIT%FileIn)     

    logInfo(*) '----- ADER-DG 3D matrices read ------'

 END SUBROUTINE IniQuadFreeADERDG3D


END MODULE ini_QuadFreeDG_mod
