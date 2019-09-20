!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2010-2016, SeisSol Group
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

MODULE receiver_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE PGM_output
     MODULE PROCEDURE PGM_output
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC  :: PGM_output
  !----------------------------------------------------------------------------

CONTAINS
  SUBROUTINE PGM_output(IO,MPI)
    !--------------------------------------------------------------------------
    ! PGM_output writes the 
    ! Peak Ground Displacement (PGD),
    ! Peak Ground Velocity (PGV),
    ! Peak Ground Acceleration (PGA),
    ! into a file together with the xyz-coordinates of the registration point.
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    ! local Variables
    INTEGER                  :: i, iErr
    INTEGER                  :: allocstat
    !--------------------------------------------------------------------------
    INTENT(INOUT)            :: IO
    !--------------------------------------------------------------------------
    !
    ! output peak ground motion data
    IF(IO%nPGMRecordPoint.GT.0)THEN
        logInfo(*) 'Writing PGM data to file PGM_map.dat ... !'

#ifdef PARALLEL
        !IF(MPI%myrank.EQ.0)THEN
           ALLOCATE(MPI%PGMarray(IO%nPGMRecordPoint,MPI%nCPU),  &
                    STAT = allocStat                              )
           IF (allocStat .NE. 0) THEN
                 logError(*) 'could not allocate',&
                 ' PGMarray for PGM output! '
                 STOP
           END IF
           MPI%PGMarray = 0.
        !ENDIF

        !Collect PGM data from all CPUs in a common array MPI%PGMarray
        CALL MPI_GATHER(IO%PGM,IO%nPGMRecordPoint,MPI%MPI_AUTO_REAL,MPI%PGMarray(:,:),IO%nPGMRecordPoint,MPI%MPI_AUTO_REAL,0,MPI%commWorld,iErr)
        
        !CPU 0 is outputting the collected PGM data
        IF(MPI%myrank.EQ.0)THEN
            MPI%PGMarray(:,1) = MAXVAL(MPI%PGMarray,DIM = 2)
            OPEN(UNIT=999,FILE='PGM_map.dat',FORM='FORMATTED',RECL=500)
            WRITE(999,'(a75)')'x                       y                       z                       PGV'
            DO i = 1, IO%nPGMRecordPoint
                WRITE(999,*) IO%UnstructRecPoint(IO%PGMstartindex+i-1)%X,    &
                             IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Y,    &
                             IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Z,    &
                             MPI%PGMarray(i,1)
            ENDDO
            CLOSE(999)
            DEALLOCATE(MPI%PGMarray)
        ENDIF
#else   
        OPEN(UNIT=999,FILE='PGM_map.dat',FORM='FORMATTED',RECL=500)
        WRITE(999,'(a75)')'x                       y                       z                       PGV'
        DO i = 1, IO%nPGMRecordPoint
           WRITE(999,*) IO%UnstructRecPoint(IO%PGMstartindex+i-1)%X,    &
                        IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Y,    &
                        IO%UnstructRecPoint(IO%PGMstartindex+i-1)%Z,    &
                        IO%PGM(i)
        ENDDO
        CLOSE(999)
#endif  
        
        DEALLOCATE(IO%PGM)
        !DEALLOCATE(IO%PGD_tmp)
        !DEALLOCATE(IO%PGA_tmp)

    ENDIF
    ! 
  END SUBROUTINE PGM_output
END MODULE receiver_mod
