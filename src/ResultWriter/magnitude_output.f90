!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2014, SeisSol Group
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
!! routine outputs the magnitude for each MPI domain that contains a subfault
!! the results need to be gathered and summarized in a postprocessing step
!! magnitude = scalar seismic moment = slip per element * element face * shear modulus
!! unit: N*m = 10^7 dyne*cm

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE magnitude_output_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
  PRIVATE
  !---------------------------------------------------------------------------!
  PUBLIC  :: magnitude_output
  PUBLIC  :: energy_rate_output
  !---------------------------------------------------------------------------!
  INTERFACE magnitude_output
     MODULE PROCEDURE magnitude_output
     MODULE PROCEDURE energy_rate_output
  END INTERFACE

CONTAINS

  SUBROUTINE magnitude_output(MaterialVal,DISC,MESH,MPI,IO,DR_comm)
    !< routine outputs the magnitude for each MPI domain that contains a subfault
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization), target   :: DISC                                              !
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    integer                         :: DR_comm !< dynamic rupture communicator
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: iElem,iSide,nSide,iFace
    INTEGER                         :: stat, UNIT_MAG, iErr, rankDR
    REAL                            :: magnitude, magnitude0
    REAL                            :: MaterialVal(:,:)
    LOGICAL                         :: exist
    CHARACTER (LEN=5)               :: cmyrank
    CHARACTER (len=200)             :: MAG_FILE
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, MESH, MPI, IO, DR_comm
    !-------------------------------------------------------------------------!
    !
    ! Compute output
    magnitude = 0.0D0
    nSide = MESH%FAULT%nSide
    DO iFace = 1,nSide
       iElem = MESH%Fault%Face(iFace,1,1)          ! Remark:
       iSide = MESH%Fault%Face(iFace,2,1)          ! iElem denotes "+" side
       IF (iElem.EQ.0) THEN
          cycle
       ENDIF
       IF (DISC%DynRup%magnitude_out(iFace)) THEN
           ! magnitude = scalar seismic moment = slip per element * element face * shear modulus
           magnitude = magnitude + DISC%DynRup%averaged_Slip(iFace)*DISC%Galerkin%geoSurfaces(iSide,iElem)*MaterialVal(iElem,2)
       ENDIF
    ENDDO
#ifdef PARALLEL
    CALL MPI_REDUCE(magnitude,magnitude0,1,MPI%MPI_AUTO_REAL,MPI_SUM,0, DR_comm,iErr)
    CALL MPI_Comm_rank(DR_comm, rankDR, iErr)
#else
    magnitude0 = magnitude
    rankDR = 0
#endif
    IF (rankDR.EQ.0) THEN

        WRITE(MAG_FILE, '(a,a4,a4)') TRIM(IO%OutputFile),'-MAG','.dat'
        UNIT_MAG = 299875
        !
        INQUIRE(FILE = MAG_FILE, EXIST = exist)
        IF(exist) THEN
            ! If file exists, then append data
            OPEN(UNIT     = UNIT_MAG                                          , & !
                 FILE     = MAG_FILE                                          , & !
                 FORM     = 'FORMATTED'                                      , & !
                 STATUS   = 'OLD'                                            , & !
                 POSITION = 'APPEND'                                         , & !
                 RECL     = 80000                                            , & !
                 IOSTAT   = stat                                               ) !
            IF(stat.NE.0) THEN                                              !
               logError(*) 'cannot open ',MAG_FILE         !
               logError(*) 'Error status: ', stat                !
               call exit(134)                                                          !
            ENDIF
        ELSE
            ! open file
            OPEN(UNIT   = UNIT_MAG                                            , & !
                 FILE     = MAG_FILE                                          , & !
                 FORM     = 'FORMATTED'                                      , & !
                 STATUS   = 'NEW'                                            , & !
                 RECL     = 80000                                            , & !
                 IOSTAT   = stat                                               ) !
            IF(stat.NE.0) THEN                                              !
               logError(*) 'cannot open ',MAG_FILE         !
               logError(*) 'Error status: ', stat                !
               call exit(134)                                                          !
            ENDIF
            !
        ENDIF
        !
        ! Write output
        WRITE(UNIT_MAG,*) magnitude0
        logInfo(*) 'seismic moment', magnitude0, 'Mw', 2./3.*log10(magnitude0)-6.07
        CLOSE( UNIT_Mag )

    ENDIF 

  END SUBROUTINE magnitude_output

  SUBROUTINE energy_rate_output(MaterialVal,time,DISC,MESH,MPI,IO)
    !< routine outputs the moment magnitude rate for each MPI domain that contains a subfault
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization), target   :: DISC                                              !
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: iElem,iSide,nSide,iFace,iBndGP
    INTEGER                         :: stat, UNIT_MAG
    REAL                            :: MomentRate,averageSR
    REAL                            :: FrictionalEnRate,averageFER
    REAL                            :: time
    REAL                            :: MaterialVal(:,:)
    LOGICAL                         :: exist
    CHARACTER (LEN=5)               :: cmyrank
    CHARACTER (len=200)             :: MAG_FILE
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, MESH, MPI, IO
    !-------------------------------------------------------------------------!

    ! generate unique name out of MPI rank
    IF (MESH%FAULT%nSide.EQ.0) THEN
       RETURN
    ENDIF

#ifdef PARALLEL
    ! pure MPI case
    WRITE(cmyrank,'(I5.5)') MPI%myrank                           ! myrank -> cmyrank
    WRITE(MAG_FILE, '(a,a7,a5,a4)') TRIM(IO%OutputFile),'-EnF_t-',TRIM(cmyrank),'.dat'
    UNIT_MAG = 399875+MPI%myrank
#else
    WRITE(MAG_FILE, '(a,a6,a4)') TRIM(IO%OutputFile),'-EnF_t','.dat'
    UNIT_MAG = 399875
#endif
    !
    INQUIRE(FILE = MAG_FILE, EXIST = exist)
    IF(exist) THEN
        ! If file exists, then append data
        OPEN(UNIT     = UNIT_MAG                                          , & !
             FILE     = MAG_FILE                                          , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'OLD'                                            , & !
             POSITION = 'APPEND'                                         , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        IF(stat.NE.0) THEN                                              !
           logError(*) 'cannot open ',MAG_FILE         !
           logError(*) 'Error status: ', stat                !
           call exit(134)                                                          !
        ENDIF
    ELSE
        ! open file
        OPEN(UNIT   = UNIT_MAG                                            , & !
             FILE     = MAG_FILE                                          , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'NEW'                                            , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        IF(stat.NE.0) THEN                                              !
           logError(*) 'cannot open ',MAG_FILE         !
           logError(*) 'Error status: ', stat                !
           call exit(134)                                                          !
        ENDIF
        WRITE(UNIT_MAG,*) '#time MomentRate FrictionalEnergyRate'
        !
    ENDIF
    !
    ! Compute output
    MomentRate = 0.0D0
    FrictionalEnRate = 0.0D0

    nSide = MESH%FAULT%nSide
    DO iFace = 1,nSide
       iElem = MESH%Fault%Face(iFace,1,1)          ! Remark:
       iSide = MESH%Fault%Face(iFace,2,1)          ! iElem denotes "+" side
       averageSR = 0d0
       averageFER = 0d0
       IF (iElem.EQ.0) THEN
          cycle
       ENDIF
       DO iBndGP=1,DISC%Galerkin%nBndGP
          averageSR = averageSR + sqrt(DISC%DynRup%SlipRate1(iBndGP,iFace)**2+DISC%DynRup%SlipRate2(iBndGP,iFace)**2)/DISC%Galerkin%nBndGP
          !frictional energy, based on the formula by Xu et al. 2012, p. 1333
          averageFER = averageFER + (DISC%DynRup%TracXY(iBndGP,iFace)*DISC%DynRup%SlipRate1(iBndGP,iFace)&
                      +DISC%DynRup%TracXZ(iBndGP,iFace)*DISC%DynRup%SlipRate2(iBndGP,iFace))/DISC%Galerkin%nBndGP
       ENDDO
       ! magnitude = scalar seismic moment = slip per element * element face * shear modulus
       MomentRate = MomentRate + averageSR*DISC%Galerkin%geoSurfaces(iSide,iElem)*MaterialVal(iElem,2)
       ! frictional energy, integrate over each element
       FrictionalEnRate = FrictionalEnRate + averageFER*DISC%Galerkin%geoSurfaces(iSide,iElem)
    ENDDO
    !
    ! Write output
    WRITE(UNIT_MAG,*) time, MomentRate, FrictionalEnRate

    CLOSE( UNIT_Mag )

  END SUBROUTINE energy_rate_output

END MODULE magnitude_output_mod
