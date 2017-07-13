!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Stefan Wenk (wenk AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/wenk)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!!
!! @section LICENSE
!! Copyright (c) 2013-2016, SeisSol Group
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

MODULE receiver_hdf_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------

  INTERFACE ini_receiver_hdf
    MODULE PROCEDURE ini_receiver_hdf
  END INTERFACE

  INTERFACE PGM_output_hdf
    MODULE PROCEDURE PGM_output_hdf
  END INTERFACE

  INTERFACE receiver_hdf
    MODULE PROCEDURE receiver_hdf
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC  :: ini_receiver_hdf
  PUBLIC  :: receiver_hdf
  PUBLIC  :: PGM_output_hdf
  !----------------------------------------------------------------------------

CONTAINS

  !< initialize hdf5 receiver output
  !< generate file structure with empty datasets
  !< write attributes
  !< write receiver location
  SUBROUTINE ini_receiver_hdf(EQN,MESH,DISC,SOURCE,IO,MPI,adaptive)
    !--------------------------------------------------------------------------
    ! HDF: This routine initializes the structure of the HDF5 output file. It !
    ! writes receiver information (location, nrec, dt, nsamp, t_end), but it  !
    ! doesn't enter values for the field data. This is the job of the sub     !
    ! receiver, which passes the values at the points to the output file in   !
    ! parallel buffered blocks.                                               !
    !--------------------------------------------------------------------------
    USE common_receiver_mod
    USE hdf_output_utils_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    TYPE (tEquations)        :: EQN
    TYPE (tUnstructMesh)     :: MESH
    TYPE (tDiscretization)   :: DISC
    TYPE (tSource)           :: SOURCE
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    LOGICAL, OPTIONAL        :: adaptive
    !--------------------------------------------------------------------------
    ! local Variables
    !--------------------------------------------------------------------------
    INTEGER,PARAMETER                 :: length = 5                        ! string length of receiver components
    !--------------------------------------------------------------------------
    CHARACTER(len=200)                :: file_name                         ! hdf5 file name
    CHARACTER(len=length)             :: comp(5,9), loc_comp(3), t_comp(1) ! field components, receiver location components, time components
    CHARACTER(len=length),ALLOCATABLE :: comp_sel(:,:)                     ! selected components
    CHARACTER(len=100),ALLOCATABLE    :: dset_name(:)                      ! selected hdf5 dataset names (stress, velocity, rotation, mtc, curl-divergence
    CHARACTER(len=100)                :: atrb_name(5)                      ! hdf5 root group attribute names
    REAL                              :: atrb(5)                           ! hdf5 root group attribute values
    INTEGER                           :: i, j, iRot                        ! loop counter
    INTEGER                           :: Vel_first, Vel_last               ! first/last position of velocity field components in outputmask
    INTEGER                           :: nsamp, natrb, ncomp(5)            ! number of receivers, samples, attributes, components
    INTEGER                           :: dset_count, any_rec, nranks       ! dataset counter, receiver switch, number of ranks
    INTEGER                           :: dset_sel(5)                       ! switch for dataset selection
    INTEGER                           :: dset_size(5,3)                    ! dataset size
    INTEGER,ALLOCATABLE               :: ncomp_sel(:), dset_mask(:)        ! number of selected components, dataset selection mask
    INTEGER,ALLOCATABLE               :: all_rec(:), rec_ranks(:)          ! mpi array of receiver switch, receiver ranks
    !------------------------------------------------------------------
    INTEGER                           :: world_grp, rec_grp, rec_comm
    !--------------------------------------------------------------------------
    INTENT(IN)               :: MESH, DISC, SOURCE
    INTENT(INOUT)            :: IO
    !--------------------------------------------------------------------------
    !
    ! get receiver locations
    !
    CALL common_ini_receiver(EQN,MESH,DISC,SOURCE,IO,MPI)
    !
#ifdef PARALLEL
    !
    !--------------------------------------------------------------------------
    ! setup mpi sub-communicator of processes containing receivers
    !--------------------------------------------------------------------------
    !
    ! identify ranks with receivers
    any_rec = 0
    IF(ANY(IO%UnstructRecpoint(:)%inside)) THEN
      any_rec = 1
    ENDIF
    !
    ! distribute information if a rank contains a receiver
    ALLOCATE( all_rec(0:MPI%nCPU-1))
    CALL MPI_ALLGATHER(any_rec, 1, MPI_INTEGER, all_rec, 1, MPI_INTEGER, MPI%commWorld, MPI%iErr)
    !
    ! count number of ranks with receivers and gather mpirank
    nranks = COUNT(all_rec.EQ.1)
    ALLOCATE( rec_ranks (nranks))
    j = 0
    DO i=0,MPI%nCPU-1
      IF(all_rec(i).EQ.1) THEN
        j = j+1
        rec_ranks(j) = i
      ENDIF
    ENDDO
    !
    ! generate new mpi communicator for all processes with receivers
    ! find group of mpi comm world communicator
    CALL MPI_COMM_GROUP(MPI%commWorld, world_grp, MPI%iErr)
    ! put selected ranks (with receivers) of global group in a new subgroup
    CALL MPI_GROUP_INCL(world_grp, nranks, rec_ranks, rec_grp, MPI%iErr)
    ! create communicator for new subgroup
    CALL MPI_COMM_CREATE(MPI%commWorld, rec_grp, rec_comm, MPI%iErr)
    ! deallocate
    DEALLOCATE(all_rec)
#endif
    !
    !--------------------------------------------------------------------------
    ! exclude ranks without receivers
    !--------------------------------------------------------------------------
    !
    IF ( IO%nlocalRecordPoint.EQ.0 ) THEN
      ! If no receiver lies inside the domain, tell the code the modified number.
      logInfo(*) 'No temporal signal or PGM is picked'
      IO%nRecordPoint      = 0
      IO%ntotalRecordPoint = 0
    ELSE
      logInfo(*) 'Pick temporal signal at ',IO%nlocalRecordPoint,' points.'
      ! WE GOT THE NUMBER OF RECEIVERS, THE NAME, THE LOCATION AND THE ELEMENT
      ALLOCATE( IO%igloblocRecordPoint(IO%nlocalRecordPoint) )
      IO%igloblocRecordPoint(:) = pack([(i,i=1,size(IO%UnstructRecPoint(:)%inside))],IO%UnstructRecPoint(:)%inside)
      !
      !--------------------------------------------------------------------------
      ! initialize hdf5 file
      !--------------------------------------------------------------------------
      !
      ! Define receiver output filename
      WRITE(file_name,'(a,a)')  TRIM(IO%OutputFile), '_rec.h5'
      !
      logInfo(*) 'Initialize file ', TRIM(file_name),' to save time signals'
      !
      ! allocate hdf5 receiver pointer
      ALLOCATE(IO%hd_rec)
      !
      ! assign mpi info.
      IO%hd_rec%mpi_grp = rec_grp
      IO%hd_rec%mpi_comm = rec_comm
      !
      ! number of samples starts at time=0s (nsamp+1)
      !!! ATTENTION !!! For future LTS implementation without Cauchy Kowalevski nsamp is not constant!
      nsamp = NINT(DISC%EndTime/IO%Pickdt)+1
      IF(nsamp.GT.DISC%MaxIteration) nsamp = DISC%MaxIteration
      !
      ! Define buffer size to 10% of nsamp
      !!! ATTENTION !!! cast to DOUBLE
      IO%MaxPickStore = NINT(DBLE(nsamp)/10.0)
      !
      ! log info
      logInfo(*) '          Number of time samples:', nsamp
      logInfo(*) 'Number of samples in HDF5 buffer:', IO%MaxPickStore
      !
      ALLOCATE( IO%CurrentPick(IO%nlocalRecordPoint) )
      ! in case of gts just a single integer for the start position is neccessary, since the data of all local receivers are
      ! written at once. but, for lts a vector is required.
      ALLOCATE( IO%hd_rec%start(IO%nlocalRecordPoint) )
      IO%CurrentPick(:) = 0
      IO%hd_rec%start(:) = 0
      !
      !--------------------------------------------------------------------------
      ! Generate internal buffers and output masks for datasets
      !--------------------------------------------------------------------------
      !
      ! time buffer
      ALLOCATE( IO%TmpTime( IO%nlocalRecordPoint, IO%MaxPickStore ) )
      !
      ! field datasets
      ! names
      IO%hd_rec%dset_name(1) = 'stress'
      IO%hd_rec%dset_name(2) = 'velocity'
      IO%hd_rec%dset_name(3) = 'rotation'
      IO%hd_rec%dset_name(4) = 'mtc'
      IO%hd_rec%dset_name(5) = 'curl_divergence'
      !
      dset_sel     = 0
      !
      ! Dataset: stress, velocity (leave out rho,mu,lambda in IO%Outputmask)
      EQN%nVar_sts = 0
      EQN%nVar_vel = 0
      DO i = 1, EQN%nVar
        IF (i.LE.6) THEN
          IF(IO%Outputmask(EQN%Dimension+i)) THEN
            EQN%nVar_sts = EQN%nVar_sts + 1
          ENDIF
        ELSE
          IF(IO%Outputmask(EQN%Dimension+i)) THEN
            EQN%nVar_vel = EQN%nVar_vel + 1
          ENDIF
        ENDIF
      ENDDO
      IF (EQN%nVar_sts.NE.0) ALLOCATE(IO%pickmask_sts(EQN%nVar_sts))
      IF (EQN%nVar_vel.NE.0) ALLOCATE(IO%pickmask_vel(EQN%nVar_vel))
      EQN%nVar_sts = 0
      EQN%nVar_vel = 0
      DO i = 1, EQN%nVar
        IF (i.LE.6) THEN
          IF(IO%Outputmask(EQN%Dimension+i)) THEN
            EQN%nVar_sts = EQN%nVar_sts + 1
            IO%pickmask_sts(EQN%nVar_sts) = i
          ENDIF
        ELSE
          IF(IO%Outputmask(EQN%Dimension+i)) THEN
            EQN%nVar_vel = EQN%nVar_vel + 1
            IO%pickmask_vel(EQN%nVar_vel) = i
          ENDIF
        ENDIF
      ENDDO
      ! stress buffer
      IF(EQN%nVar_sts.NE.0) THEN
        ALLOCATE( IO%TmpState_sts( IO%nlocalRecordPoint, IO%MaxPickStore, EQN%nVar_sts ) )
        dset_sel(1) = 1                           ! if selected, this dataset will be written to the hdf5 file
      ENDIF
      ! velocity buffer
      IF(EQN%nVar_vel.NE.0) THEN
        ALLOCATE( IO%TmpState_vel( IO%nlocalRecordPoint, IO%MaxPickStore, EQN%nVar_vel ) )
        dset_sel(2) = 1                           ! if selected, this dataset will be written to the hdf5 file
      ENDIF
      !
      ! Dataset: rotation (#rot_comp = #vel_comp)
      iRot = size(IO%RotationMask)
      IF(IO%Rotation.NE.0) THEN
        EQN%nVar_rot = 0
        DO i = 1,iRot
          IF(IO%RotationMask(i)) THEN
            EQN%nVar_rot = EQN%nVar_rot + 1
          ENDIF
        ENDDO
        ALLOCATE(IO%pickmask_rot(EQN%nVar_rot))
        EQN%nVar_rot = 0
        DO i = 1,iRot
          IF(IO%RotationMask(i))THEN
            EQN%nVar_rot = EQN%nVar_rot+1
            IO%pickmask_rot(EQN%nVar_rot) = i
          ENDIF
        ENDDO
        ! rotation buffer
        ALLOCATE(IO%TmpState_rot( IO%nlocalRecordPoint, IO%MaxPickStore, EQN%nVar_rot ))
        IF(IO%Rotation.EQ.1) dset_sel(3) = 1 ! if selected, rot dataset will be written to the hdf5 file
        IF(IO%Rotation.EQ.2) dset_sel(4) = 1 ! if selected, mtc dataset will be written to the hdf5 file
        IF(IO%Rotation.EQ.3) dset_sel(5) = 1 ! if selected, cd dataset will be written to the hdf5 file
      ENDIF
      !
      !--------------------------------------------------------------------------
      ! Generate output masks for components
      !--------------------------------------------------------------------------
      !
      ! stress
      ncomp(1) = 6
      comp(1,1:ncomp(1)) = (/'xx','yy','zz','xy','yz','xz'/)
      ! velocity
      ncomp(2) = 3
      comp(2,7:9) = (/'u','v','w'/)
      ! rotation
      ncomp(3) = 3
      comp(3,1:ncomp(3)) = (/'rx','ry','rz'/)
      ! mtc
      ncomp(4) = 9
      comp(4,1:ncomp(4)) = (/'Gxx','Gxy','Gxz','Gyx','Gyy','Gyz','Gzx','Gzy','Gzz'/)
      ! curl div
      ncomp(5) = 4
      comp(5,1:ncomp(5)) = (/'cx','cy','cz','dv'/)
      !
      !--------------------------------------------------------------------------
      ! Gather selected datasets and components
      !--------------------------------------------------------------------------
      !
      IO%hd_rec%ndset = COUNT(dset_sel.EQ.1)
      ALLOCATE(dset_name(IO%hd_rec%ndset))
      ALLOCATE(ncomp_sel(IO%hd_rec%ndset))
      ALLOCATE(dset_mask(IO%hd_rec%ndset))
      ALLOCATE(comp_sel(IO%hd_rec%ndset,9))
      dset_count=0
      DO i=1,5
        IF(dset_sel(i).EQ.1) THEN
          IF(i.EQ.1) THEN
            dset_count = dset_count+1
            dset_mask(dset_count) = i
            comp_sel(dset_count,1:EQN%nVar_sts) = comp(i,IO%pickmask_sts(:))
            ncomp_sel(dset_count) = EQN%nVar_sts
          ELSE IF(i.EQ.2) THEN
            dset_count = dset_count+1
            dset_mask(dset_count) = i
            comp_sel(dset_count,1:EQN%nVar_vel) = comp(i,IO%pickmask_vel(:))
            ncomp_sel(dset_count) = EQN%nVar_vel
          ELSE
            dset_count = dset_count+1
            dset_mask(dset_count) = i
            comp_sel(dset_count,1:EQN%nVar_rot) = comp(i,IO%pickmask_rot(:))
            ncomp_sel(dset_count) = EQN%nVar_rot
          ENDIF
        ENDIF
      ENDDO
      dset_name = IO%hd_rec%dset_name(dset_mask(:))
      !
      !--------------------------------------------------------------------------
      ! Determine absolute dataspace of datasets in the hdf5 file
      !--------------------------------------------------------------------------
      !
      j = 0
      DO i=1,5
        IF(dset_sel(i).EQ.1) THEN
          j=j+1
          ! global dataset size
          dset_size(j,:) = (/ IO%ntotalRecordPoint, nsamp, ncomp_sel(j)/)
        ENDIF
      ENDDO
      !
      !--------------------------------------------------------------------------
      ! Generate hdf5 file with gathered infos
      !--------------------------------------------------------------------------
      !
      IO%hd_rec%grp_name(1) = '/receiver'
      IO%hd_rec%grp_name(2) = '/source'
      !
      ! write data, which will not change in time loop, by one process, which contains receivers
#ifdef PARALLEL
      IF(MPI%myrank.EQ.rec_ranks(1)) THEN
#endif
        ! in the hdf5 forum it was written, that generating the datasets has to be done collectively by all processes which will
        ! write to the hdf5 file. but, anyway initialization of the structure and writing attributes works also with one process.
        !
        ! create hdf5 file
        CALL create_hdf5_file(file_name, 2, IO%hd_rec%grp_name, IO%hd_rec%file_id)
        ! write receiver attributes to root group
        atrb_name(1) = 'dt'
        atrb_name(2) = 'n_rec'
        atrb_name(3) = 'n_samp'
        atrb_name(4) = 't_end'
        atrb_name(5) = 'max_iter'
        !!! ATTENTION !!! cast to DOUBLE
        atrb      = (/ DBLE(IO%pickdt),DBLE(IO%ntotalRecordPoint),DBLE(nsamp),DBLE(DISC%EndTime),DBLE(DISC%MaxIteration) /)
        DO i=1,5
          CALL write_hdf5_atrb(IO%hd_rec%file_id, atrb_name(i), 1, 1, (/atrb(i)/), IO%hd_rec%grp_name(1))
        ENDDO
        ! create dataset receiver location and write components as attributes
        loc_comp = (/'x','y','z'/)
        CALL create_hdf5_dset(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), 'location', 2, (/IO%ntotalRecordPoint,3/) )
        CALL write_hdf5_atrb(IO%hd_rec%file_id, 'components', length, 3, loc_comp, IO%hd_rec%grp_name(1), 'location')
        ! create dataset time and write component as attribute
        ! 1D or 2D dataset for time (LTS without Cauchy-Kowalewski)?
        t_comp = (/'t'/)
        CALL create_hdf5_dset(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), 'time' , 2, (/IO%ntotalRecordPoint,nsamp/) )
        CALL write_hdf5_atrb(IO%hd_rec%file_id, 'components', length, 1, t_comp, IO%hd_rec%grp_name(1), 'time')
        ! create selected wavefield datasets and write components as attributes
        DO i=1,IO%hd_rec%ndset
          CALL create_hdf5_dset(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), dset_name(i), 3, dset_size(i,:) )
          CALL write_hdf5_atrb(IO%hd_rec%file_id, 'components', length, ncomp_sel(i), comp_sel(i,1:ncomp_sel(i)), &
                               IO%hd_rec%grp_name(1), dset_name(i))
        ENDDO
        ! write receiver location dataset
        CALL write_hdf5_data_serial(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), 'location', 2, (/IO%ntotalRecordPoint,3/), &
                                    (/ IO%UnstructRecPoint%x, IO%UnstructRecPoint%y, IO%UnstructRecPoint%z /) )
        ! close file
        CALL close_hdf5_file(IO%hd_rec%file_id)
        !
        !--------------------------------------------------------------------------
        ! to be continued...
        ! 1. write source info:
        !      atrb_name(1) = 'type' depending on source type, different source properties can be written as hdf5 attributes
        ! 2. read external text file and write additional info:
        !      svn revision number, user name, host name, compiler/compiler flags, background model...
        !--------------------------------------------------------------------------
        !
#ifdef PARALLEL
      ENDIF
      ! write log info
      logInfo0(*) 'HDF5 receiver initialization done by rank:', rec_ranks(1)
      !
      DEALLOCATE(rec_ranks)
#endif
      ! deallocate
      DEALLOCATE(comp_sel)
      DEALLOCATE(ncomp_sel)
      DEALLOCATE(dset_mask)
      DEALLOCATE(dset_name)
      !
      !--------------------------------------------------------------------------
      ! Open hdf5 file for parallel io and overwrite file id
      !--------------------------------------------------------------------------
      !
      CALL open_hdf5_file_parallel(file_name,IO%hd_rec%mpi_comm,IO%hd_rec%file_id)
      !
      ! write log info
      logInfo(*) 'Open HDF5 file for parallel access'
      logInfo(*) '<--------------------------------------------------------->'
      logInfo(*) ' '
      !
    ENDIF
    !
    ! Set local pick time for all processes
    IF(DISC%Galerkin%DGMethod.EQ.3) THEN
      ALLOCATE( IO%LocalPickTime(IO%ntotalRecordPoint) )
      IO%LocalPickTime(:) = 0.
    ENDIF
    !
    !------------------------------------------------
    ! Initialization of maximum values for PGM output
    !------------------------------------------------
    !
    IF(IO%nPGMRecordPoint.GT.0)THEN
        !Check, which output is required (PGD, PGV, PGA)
        ALLOCATE( IO%PGM(IO%nPGMRecordPoint)      ) !PGV
        !ALLOCATE( IO%PGD_tmp(IO%nPGMRecordPoint,3)) !PGD
        !ALLOCATE( IO%PGA_tmp(IO%nPGMRecordPoint,3)) !PGA

        IO%PGM     = 0.
        !IO%PGD_tmp = 0.
        !IO%PGA_tmp = 0.
    ENDIF

  END SUBROUTINE ini_receiver_hdf
  !
  !< hdf5 receiver output
  !<
  SUBROUTINE receiver_hdf(EQN,MESH,DISC,MPI,IO,time_op,dt_op)
    !--------------------------------------------------------------------------
    USE common_receiver_mod
    USE hdf_output_utils_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !--------------------------------------------------------------------------
    ! argument list declaration
    !--------------------------------------------------------------------------
    TYPE (tEquations)              :: EQN
    TYPE (tUnstructMesh)           :: MESH
    TYPE (tDiscretization)         :: DISC
    TYPE (tMPI)                    :: MPI
    TYPE (tInputOutput)            :: IO
    REAL, OPTIONAL                 :: time_op, dt_op
    !--------------------------------------------------------------------------
    ! local Variables
    !--------------------------------------------------------------------------
    REAL                           :: time, dt, localpicktime
    REAL                           :: state( EQN%nVarTotal), state_rot(9)
    REAL                           :: TaylorDOF(DISC%Galerkin%nDegFrRec,EQN%nVarTotal,0:DISC%Galerkin%nPolyRec)
    REAL                           :: PGD_tot, PGV_tot, PGA_tot
    INTEGER                        :: i, j
    INTEGER                        :: iElem, jGlob, chunk_size, rec_start, io_flag
    LOGICAL                        :: write_out
    !--------------------------------------------------------------------------
    INTENT(IN)                     :: EQN,MESH                                !
    INTENT(INOUT)                  :: IO,DISC                                 !
    !--------------------------------------------------------------------------
    !
    ! register epik/scorep function receiver
    EPIK_FUNC_REG("receiver")
    SCOREP_USER_FUNC_DEFINE()
    ! start epik/scorep function receiver
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN("receiver")
    !
    localpicktime = 0
    !
    ! exclude processes without receivers (nlocalRecordPoint = 0)
    ! loop over local receivers
    DO j = 1,IO%nlocalRecordPoint
      !
      jGlob = IO%igloblocRecordPoint(j)
      !
      SELECT CASE(DISC%Galerkin%DGMethod)
      CASE(3)
        ! LTS
        CALL common_receiver_ck(EQN,DISC,IO,jGlob,TaylorDOF)
        write_out     = .true.                              ! write data whenever buffer is full
        io_flag       = 0                                   ! independent io
        chunk_size    = 1
        rec_start     = j
      CASE DEFAULT
        ! GTS
        CALL common_receiver_ck(EQN,DISC,IO,jGlob,TaylorDOF)
        write_out     = .false.
        IF(j.EQ.IO%nlocalRecordPoint) THEN
          write_out   = .true.                              ! write data when last local receiver is reached
          io_flag     = 1                                   ! collective io
          chunk_size  = IO%nlocalRecordPoint
          rec_start   = 1
        ENDIF
      END SELECT
      !
      ! do while only loops if timestep sub-sampling is allowed. Currently it is prohibited, see SUBROUTINE common_receiver_ck
      ! last sample is not written, if timestep is adapted
      DO WHILE( (localpicktime.GE.time).AND.(localpicktime.LE.time+dt+1e-10).AND.(localpicktime.LE.DISC%EndTime+1e-10) )
        !
        CALL common_receiver_interp(EQN,MESH,DISC,IO,jGlob,TaylorDof,time,localpicktime,state,state_rot)
        !
        IF(j.LE.IO%nlocalRecordPoint)THEN   ! only output subsequent samples for seismogram record points but no pgm
          !
          !< buffer receiver output
          IO%CurrentPick(j) = IO%CurrentPick(j) + 1
          IO%TmpTime(j,IO%CurrentPick(j)) = localpicktime
          IF(associated(IO%pickmask_sts)) IO%TmpState_sts(j,IO%CurrentPick(j),:) = State(IO%pickmask_sts(:))
          IF(associated(IO%pickmask_vel)) IO%TmpState_vel(j,IO%CurrentPick(j),:) = State(IO%pickmask_vel(:))
          IF(associated(IO%pickmask_rot)) IO%TmpState_rot(j,IO%CurrentPick(j),:) = State_rot(IO%pickmask_rot(:))
          !
          ! if write_out is true, check if buffer is full
          ! In GTS case, all receivers are written at once.
          ! In LTS case, each receiver is written separately.
          IF(write_out) THEN
            ! if buffer is full or end time is reached, write out data
            IF(IO%CurrentPick(j).EQ.IO%MaxPickStore.OR.(localpicktime+IO%pickdt).GE.DISC%EndTime) THEN
              !
              logInfo(*)'-----------------------------------------------------------'
              logInfo(*)'    PLOT receiver data for current time interval'
              logInfo(*)'                      start time :', IO%TmpTime(1,1)
              logInfo(*)'                        end time :', IO%TmpTime(1,IO%CurrentPick(j))
              logInfo(*)'-----------------------------------------------------------'
              !
              ! write time
              CALL write_hdf5_data_parallel(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), 'time', &
                                            (/ chunk_size, IO%CurrentPick(j) /), &
                                            IO%igloblocRecordPoint(rec_start:j)-1, &
                                            IO%hd_rec%start(j), IO%TmpTime(rec_start:j,1:IO%CurrentPick(j)), io_flag )
              ! write stress
              IF(associated(IO%pickmask_sts)) THEN
                CALL write_hdf5_data_parallel(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), IO%hd_rec%dset_name(1), &
                                              (/ chunk_size, IO%CurrentPick(j), EQN%nVar_sts /), &
                                              IO%igloblocRecordPoint(rec_start:j)-1, IO%hd_rec%start(j), &
                                              IO%TmpState_sts(rec_start:j,1:IO%CurrentPick(j),:), io_flag )
              ENDIF
              ! write velocity
              IF(associated(IO%pickmask_vel)) THEN
                CALL write_hdf5_data_parallel(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), IO%hd_rec%dset_name(2), &
                                              (/ chunk_size, IO%CurrentPick(j), EQN%nVar_vel /), &
                                              IO%igloblocRecordPoint(rec_start:j)-1, IO%hd_rec%start(j), &
                                              IO%TmpState_vel(rec_start:j,1:IO%CurrentPick(j),:), io_flag )
              ENDIF
              ! write rotation
              IF(IO%Rotation.EQ.1) THEN
                CALL write_hdf5_data_parallel(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), IO%hd_rec%dset_name(3), &
                                              (/ chunk_size, IO%CurrentPick(j), EQN%nVar_rot /), &
                                              IO%igloblocRecordPoint(rec_start:j)-1, IO%hd_rec%start(j), &
                                              IO%TmpState_rot(rec_start:j,1:IO%CurrentPick(j),:), io_flag )
              ENDIF
              ! write mtc
              IF(IO%Rotation.EQ.2) THEN
                CALL write_hdf5_data_parallel(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), IO%hd_rec%dset_name(4), &
                                              (/ chunk_size, IO%CurrentPick(j), EQN%nVar_rot /), &
                                              IO%igloblocRecordPoint(rec_start:j)-1, IO%hd_rec%start(j), &
                                              IO%TmpState_rot(rec_start:j,1:IO%CurrentPick(j),:), io_flag )
              ENDIF
              ! write curl, divergence
              IF(IO%Rotation.EQ.3) THEN
                CALL write_hdf5_data_parallel(IO%hd_rec%file_id, IO%hd_rec%grp_name(1), IO%hd_rec%dset_name(5), &
                                              (/ chunk_size, IO%CurrentPick(j), EQN%nVar_rot /), &
                                              IO%igloblocRecordPoint(rec_start:j)-1, IO%hd_rec%start(j), &
                                              IO%TmpState_rot(rec_start:j,1:IO%CurrentPick(j),:), io_flag )
              ENDIF
              ! set new start position (in GTS case, just IO%hd_rec%start(end) is used) and reset buffer counter
              IO%hd_rec%start(j) = IO%hd_rec%start(j) + IO%CurrentPick(j)
              IO%CurrentPick(rec_start:j) = 0.
            ENDIF
          ENDIF
        ELSE  ! check PGM record point
          i = j-IO%nRecordPoint
          ! First compute PGD as integral of velocities
          !IO%PGD_tmp(i,1) = IO%PGD_tmp(i,1) + state(7)*IO%pickdt
          !IO%PGD_tmp(i,2) = IO%PGD_tmp(i,2) + state(8)*IO%pickdt
          !IO%PGD_tmp(i,3) = IO%PGD_tmp(i,3) + state(9)*IO%pickdt
          !PGD_tot         = SQRT(IO%PGD_tmp(i,1)**2 + IO%PGD_tmp(i,2)**2 + IO%PGD_tmp(i,3)**2)
          !IF( PGD_tot.GT.IO%PGM(i,1) )THEN
          !   IO%PGM(i,1)  = PGD_tot
          !ENDIF

          ! Second compute PGV
          PGV_tot = SQRT(state(7)**2+state(8)**2+state(9)**2)
          IF( PGV_tot.GT.IO%PGM(i) )THEN
             IO%PGM(i) = PGV_tot
          ENDIF

          ! Third compute PGA as derivative of velocities
          !IO%PGA_tmp(i,1) = (state(7)-IO%PGA_tmp(i,1)) / IO%pickdt
          !IO%PGA_tmp(i,2) = (state(8)-IO%PGA_tmp(i,2)) / IO%pickdt
          !IO%PGA_tmp(i,3) = (state(9)-IO%PGA_tmp(i,3)) / IO%pickdt
          !PGA_tot         = SQRT(IO%PGA_tmp(i,1)**2 + IO%PGA_tmp(i,2)**2 + IO%PGA_tmp(i,3)**2)
          !IF( PGA_tot.GT.IO%PGM(i,3) )THEN
          !   IO%PGM(i,3)  = PGA_tot
          !ENDIF
          ! Store current velocity values for derivation of PGA in the next time step
          !IO%PGA_tmp(i,1) = state(7)
          !IO%PGA_tmp(i,2) = state(8)
          !IO%PGA_tmp(i,3) = state(9)
        ENDIF  ! seismogram or PGM record points
        IF (IO%pickDtType .EQ. 2) THEN
          localpicktime = localpicktime + IO%pickdt * dt
        ELSE
          localpicktime = localpicktime + IO%pickdt
        ENDIF
      ENDDO ! END DO WHILE localpicktime
      !
      IF( DISC%Galerkin%DGMethod.EQ.3) THEN
        IO%LocalPickTime(jGlob) = localpicktime
      ENDIF
      !
    ENDDO ! END DO nlocalReceiverPoint
    !
    IF( DISC%Galerkin%DGMethod.EQ.1) THEN
      IO%picktime = localpicktime
    ENDIF
    ! end epik/scorep function receiver
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()

  END SUBROUTINE receiver_hdf
  !
  !< PGM_output writes the
  !< Peak Ground Displacement (PGD),
  !< Peak Ground Velocity (PGV),
  !< Peak Ground Acceleration (PGA),
  !< into a file together with the xyz-coordinates of the registration point.
  SUBROUTINE PGM_output_hdf(IO,MPI)
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
        CALL MPI_GATHER(IO%PGM,IO%nPGMRecordPoint,MPI%MPI_AUTO_REAL,MPI%PGMarray(:,:),IO%nPGMRecordPoint, &
                        MPI%MPI_AUTO_REAL,0,MPI%commWorld,iErr)

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
  END SUBROUTINE PGM_output_hdf

END MODULE receiver_hdf_mod
