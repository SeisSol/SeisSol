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

PROGRAM GAMBIT2METIS_HEX
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Type declaration
    INTEGER                       :: NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL
    INTEGER                       :: uid1, uid2
    INTEGER                       :: i
    REAL,    ALLOCATABLE          :: Vertices(:,:)
    INTEGER, ALLOCATABLE          :: Elements(:,:)
    CHARACTER (LEN=350)           :: InFilename
    CHARACTER (LEN=350)           :: OutFilename
    CHARACTER (LEN=70)            :: StringLine            
    !--------------------------------------------------------------------------
    !
    WRITE(*,*) ' | -------------------------------------------------- |' 
    WRITE(*,*) ' |                  GAMBIT2METIS_HEX                  |'
    WRITE(*,*) ' | -------------------------------------------------- |' 
    !
    WRITE(*,*) '  ' 
    WRITE(*,*) '   Enter the GAMBIT NEUTRAL file name (.neu is appended!):  '
    WRITE(*,*) '  ' 

    READ (*,*) InFilename
    OutFilename = TRIM(InFilename)//'.met'
    InFilename  = TRIM(InFilename)//'.neu'
    WRITE(*,*) '   Reading from file ', TRIM(InFilename)
    WRITE(*,*) '   Writing to   file ', TRIM(OutFilename)
    
    uid1 = 10
    uid2 = 11
    
    OPEN( UNIT = uid1, FILE = TRIM(InFilename)  ) 
        DO i = 1,6
           READ(uid1,*)  StringLine
        ENDDO
        
        READ(uid1,*)  NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL

        WRITE(*,*) '  ' 
        WRITE(*,*) '   Found ',NELEM,' hexahedral elements.'
        WRITE(*,*) '  ' 

        ALLOCATE( Vertices(NUMNP,4)  )
        ALLOCATE( Elements(NELEM,11) )

        READ(uid1,*)  StringLine
        READ(uid1,*)  StringLine

        DO i = 1,NUMNP
            READ(uid1,*) Vertices(i,1:4)
        ENDDO

        READ(uid1,*)  StringLine
        READ(uid1,*)  StringLine

        DO i = 1,NELEM
            READ(uid1,*) Elements(i,1:10)
            READ(uid1,*) Elements(i,11)
        ENDDO
    CLOSE(uid1)

    OPEN( UNIT = uid2, FILE = TRIM(OutFilename), STATUS = 'UNKNOWN' )
        WRITE(uid2,*) NELEM, 3
        DO i = 1, NELEM
          WRITE(uid2,'(I10,I10,I10,I10,I10,I10,I10,I10)') Elements(i,4:11)
        ENDDO
    CLOSE(uid2)

    WRITE(*,*) ' | -------------------------------------------------- |' 
    WRITE(*,*) ' |     GAMBIT2METIS_HEX finished successfully         |'
    WRITE(*,*) ' | -------------------------------------------------- |' 
    WRITE(*,*) '  ' 
   
        
END PROGRAM GAMBIT2METIS_HEX
