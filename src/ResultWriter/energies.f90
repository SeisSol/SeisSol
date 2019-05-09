!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Stephanie Wollherr
!!
!! @section LICENSE
!! Copyright (c) 2016, SeisSol Group
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
!! routine outputs the total dissipated plastic energy for each MPI domain
!! the results need to be gathered and summarized in a postprocessing step
!! formulas see Xu et al. 2012
!! 

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE energies_output_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  PUBLIC  :: energies_output
  !---------------------------------------------------------------------------!
  INTERFACE energies_output
     MODULE PROCEDURE energies_output
  END INTERFACE

CONTAINS
  SUBROUTINE energies_output(i_fullUpdateTime, i_timeStepWidth, i_receiverTime, DISC, EQN, MESH, MPI,IO,time_op,dt_op)

    !< routine outputs the diisipated plastic energy for each MPI domain
    !-------------------------------------------------------------------------!
    IMPLICIT NONE
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tDiscretization), target   :: DISC                                              !
    TYPE(tUnstructMesh)             :: MESH
    TYPE(tEquations)                :: EQN
    TYPE(tMPI)                      :: MPI
    TYPE(tInputOutput)              :: IO
    REAL,OPTIONAL                   :: time_op, dt_op
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: iElem, nElem
    INTEGER                         :: stat, UNIT_ENERGY
    REAL                            :: time, dt
    REAL                            :: plast_energy, kinetic_energy, estrain_energy, fracture_energy
    !REAL                            :: MaterialVal(:,:)
    LOGICAL                         :: exist
     REAL                           :: localpicktime
    CHARACTER (LEN=5)               :: cmyrank
    CHARACTER (len=200)             :: ENERGY_FILE
#ifdef OMP
    INTEGER                         :: TID,omp_get_thread_num
    CHARACTER (LEN=2)               :: c_TID
#endif
    real*8                          :: i_fullUpdateTime
    real*8                          :: i_timeStepWidth
    real*8                          :: i_receiverTime

    !-------------------------------------------------------------------------!
    INTENT(IN)    :: DISC, EQN, MESH, MPI, time_op, dt_op
    INTENT(INOUT) :: IO
    !-------------------------------------------------------------------------!
        time          = i_fullUpdateTime
        dt            = i_timeStepWidth
        localpicktime = i_receiverTime !change that to picktime_energy!


    !only output at specific timesteps/times
    DO WHILE( (localpicktime.GE.time).AND.(localpicktime.LE.time+dt+1e-10).AND.(localpicktime.LE.DISC%EndTime+1e-10) )

    ! generate unique name out of MPI rank
#ifdef PARALLEL
#ifdef OMP
    ! hybrid
    TID = omp_get_thread_num()
    WRITE(c_TID,'(I2.2)') TID
    WRITE(cmyrank,'(I5.5)') MPI%myrank                           ! myrank -> cmyrank
    WRITE(ENERGY_FILE, '(a,a4,a5,a4)') TRIM(IO%OutputFile),'-EN-',TRIM(cmyrank),'.dat'
    UNIT_ENERGY = 99875+MPI%myrank+TID
#else
    !pure mpi case
    WRITE(cmyrank,'(I5.5)') MPI%myrank                           ! myrank -> cmyrank
    WRITE(ENERGY_FILE, '(a,a4,a5,a4)') TRIM(IO%OutputFile),'-EN-',TRIM(cmyrank),'.dat'
    UNIT_ENERGY = 99875+MPI%myrank
#endif
#else
    !no parallelization
    WRITE(ENERGY_FILE, '(a,a4,a4)') TRIM(IO%OutputFile),'-EN-','.dat'
    UNIT_ENERGY = 99875
#endif

    !
    INQUIRE(FILE = ENERGY_FILE, EXIST = exist)
    IF(exist) THEN
        ! If file exists, then append data
        OPEN(UNIT     = UNIT_ENERGY                                      , & !
             FILE     = ENERGY_FILE                                      , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'OLD'                                            , & !
             POSITION = 'APPEND'                                         , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        IF(stat.NE.0) THEN                                                   !
           logError(*) 'cannot open ',ENERGY_FILE                            !
           logError(*) 'Error status: ', stat                                !
           STOP                                                              !
        ENDIF
    ELSE
        ! open file
        OPEN(UNIT   = UNIT_ENERGY                                        , & !
             FILE     = ENERGY_FILE                                      , & !
             FORM     = 'FORMATTED'                                      , & !
             STATUS   = 'NEW'                                            , & !
             RECL     = 80000                                            , & !
             IOSTAT   = stat                                               ) !
        IF(stat.NE.0) THEN                                                   !
           logError(*) 'cannot open ',ENERGY_FILE                            !
           logError(*) 'Error status: ', stat                                !
           STOP                                                              !
        ENDIF
        IF (EQN%Plasticity .EQ. 0) THEN
            WRITE(UNIT_ENERGY,*) '#time KineticEnergy '
        ELSE
            WRITE(UNIT_ENERGY,*) '#time KineticEnergy PlasticEnergy ElasticStrainEnergy '
        ENDIF
    ENDIF
    !
    ! Compute output
    ! sum over each element in the mpi domain
    nElem = MESH%nELEM
    kinetic_energy = 0.0
    plast_energy = 0.0
    estrain_energy = 0.0

    DO iElem = 1,nElem
           kinetic_energy = kinetic_energy + EQN%Energy(1,iElem)
           plast_energy = plast_energy + EQN%Energy(2,iElem)
           estrain_energy = estrain_energy + EQN%Energy(3,iElem)
    ENDDO
    !
    ! Write output
    IF (EQN%Plasticity .EQ. 0) THEN
        WRITE(UNIT_ENERGY,*) time, kinetic_energy
    ELSE
        WRITE(UNIT_ENERGY,*) time, kinetic_energy, plast_energy, estrain_energy
    ENDIF

    CLOSE( UNIT_ENERGY )


    ! update picktime
    IF (IO%pickDtType .EQ. 2) THEN !timestep-wise output
        localpicktime = localpicktime + IO%pickdt_energy * dt
    ELSE !time-wise output
        localpicktime = localpicktime + IO%pickdt_energy
    ENDIF


   ENDDO ! END DO WHILE localpicktime

    IF( DISC%Galerkin%DGMethod.EQ.1) THEN !GTS
      IO%picktime_energy = localpicktime
    ENDIF

  END SUBROUTINE energies_output

END MODULE energies_output_mod

