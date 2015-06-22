!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Dmitry Pinaev (pinaev AT in.tum.de)
!!
!! @section LICENSE
!! Copyright (c) 2013, SeisSol Group
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
!! Routines handling HDF5 fault output

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE hdf_faultoutput_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE COMMON_operators_mod
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  ! Public procedures and functions
  !

  !
  !---------------------------------------------------------------------------!
  PUBLIC   :: ini_fault_receiver_hdf
  PUBLIC   :: write_fault_output_receiverwise_hdf

  PRIVATE  :: init_mpi_hdf_communication
  PRIVATE  :: init_hdf_fault_file
  !---------------------------------------------------------------------------!

CONTAINS

  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE init_mpi_hdf_communication(DISC, MPI, rec_grp, rec_comm, number_of_ranks, rec_ranks)

    IMPLICIT NONE

#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif

    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    TYPE(tMPI)              :: MPI
    INTEGER                 :: rec_grp, rec_comm
    INTEGER                 :: number_of_ranks

    !--------- local variables
    INTEGER                 :: any_rec, j, i
    INTEGER, ALLOCATABLE    :: all_rec(:)   ! mpi array of receiver switch, receiver ranks
    INTEGER, ALLOCATABLE    :: rec_ranks(:)
    INTEGER                 :: world_grp

    ! count if there any pickpoints are inside this domain
    any_rec = 0
    IF(DISC%DynRup%DynRup_out_atPickpoint%nDR_pick.GT.0) THEN
      any_rec = 1
    ENDIF

    ! distribute information if a rank contains a receiver
    ALLOCATE( all_rec(0 : MPI%nCPU - 1))
    CALL MPI_ALLGATHER(any_rec, 1, MPI_INTEGER, all_rec, 1, MPI_INTEGER, MPI_COMM_WORLD, MPI%iErr)

    ! count number of ranks with receivers
    number_of_ranks = COUNT(all_rec.EQ.1)

    logInfo(*) 'There are ', number_of_ranks, ' ranks with receivers available '

    ! collect number of nodes with receivers
    ALLOCATE( rec_ranks (number_of_ranks))
    j = 0
    DO i = 0, MPI%nCPU - 1
      IF(all_rec(i).EQ.1) THEN
        j = j + 1
        rec_ranks(j) = i
      ENDIF
    ENDDO

    ! generate new mpi communicator for all processes with receivers
    ! find group of mpi comm world communicator
    CALL MPI_COMM_GROUP(MPI_COMM_WORLD, world_grp, MPI%iErr)
    logInfo(*) 'After MPI_COMM_GROUP ', MPI%iErr

    ! put selected ranks (with receivers) of global group in a new subgroup
    CALL MPI_GROUP_INCL(world_grp, number_of_ranks, rec_ranks, rec_grp, MPI%iErr)
    logInfo(*) 'After MPI_GROUP_INCL ', MPI%iErr

    ! create communicator for new subgroup
    CALL MPI_COMM_CREATE(MPI_COMM_WORLD, rec_grp, rec_comm, MPI%iErr)
    logInfo(*) 'After MPI_COMM_CREATE ', MPI%iErr

    ! deallocate
    DEALLOCATE(all_rec)

  END SUBROUTINE init_mpi_hdf_communication



  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE init_hdf_fault_file(EQN, MESH, DISC, &
                                 IO, MPI, rec_grp, rec_comm, &
                                 number_of_ranks, rec_ranks)
    USE hdf_output_utils_mod

    IMPLICIT NONE

#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif

    TYPE (tEquations)        :: EQN
    TYPE (tUnstructMesh)     :: MESH
    TYPE (tDiscretization)   :: DISC
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    INTEGER                  :: rec_grp, rec_comm
    INTEGER                  :: number_of_ranks
    INTEGER                  :: rec_ranks(:)

    INTEGER,PARAMETER        :: length = 5    ! string length of receiver components
    CHARACTER(len=length)    :: location_component_names(3) ! names for data in hdf5

    CHARACTER(len=200)       :: file_name     ! hdf5 file name
    INTEGER                  :: number_of_samples, natrb, ncomp(5)            ! number of receivers, samples, attributes, components
    INTEGER                  :: i, j

    CHARACTER(len=100)       :: parameter_names(11)
    CHARACTER(len=100)       :: dataset_names(5)
    INTEGER                  :: dataset_sizes(5)
    INTEGER                  :: dataset_name_offsets(5)
    INTEGER                  :: dataset_size(5,3)                    ! dataset size
    INTEGER                  :: dataset_selected(5)


    ! all elements are zero
    dataset_selected = DISC%DynRup%DynRup_out_atPickpoint%OutputMask


    ! allocate hdf5 receiver pointer
    ALLOCATE(DISC%DynRup%hd_rec)

    ! Set MPI group ids
    DISC%DynRup%hd_rec%mpi_grp =  rec_grp
    DISC%DynRup%hd_rec%mpi_comm = rec_comm

    ! names for output variables
    parameter_names = (/'SRs', 'SRd', 'T_s', 'T_d', 'P_n', 'u_n', 'Mud', 'SRs', 'ts0', 'td0', 'p0'/)

    ! output variables are grouped
    dataset_names = (/ 'SRs_SRd', 'T_s_T_d_P_n', 'u_n', 'Mud_StV', 'ts0_td0_p0'/)
    dataset_sizes = (/ 2, 3, 1, 2, 3 /)
    dataset_name_offsets = (/ 0, 2, 5, 6, 7 /)


    ! number of time samples
    logInfo(*) '- End time ', DISC%EndTime
    logInfo(*) '- Timestep ', IO%Pickdt

    number_of_samples = (NINT(DISC%EndTime/IO%Pickdt)+1) / DISC%DynRup%DynRup_out_atPickpoint%printtimeinterval

    IF(number_of_samples.GT.DISC%MaxIteration) number_of_samples = DISC%MaxIteration
    logInfo(*) 'Number of time samples ', number_of_samples

    logInfo(*) 'Before counting dataset_size'

    j = 0
    DO i = 1, 5
      IF (DISC%DynRup%DynRup_out_atPickpoint%OutputMask(i).EQ.1) THEN
        logInfo(*) 'Pickpoint output ', i, ' is switched on '
        j = j + 1
        ! global dataset size
        dataset_size(j,:) = (/ DISC%DynRup%DynRup_out_atPickpoint%nOutPoints, number_of_samples, dataset_sizes(i) /)
      ENDIF
    ENDDO


    ! Define receiver output filename
    WRITE(file_name,'(a,a)')  TRIM(IO%OutputFile), '_fault_receiver.h5'

    DISC%DynRup%hd_rec%grp_name(1) = '/receiver'
    DISC%DynRup%hd_rec%grp_name(2) = '/source'

      ! write data, which will not change in time loop, by one process, which contains receivers
!#ifdef PARALLEL
    IF(MPI%myrank.EQ.rec_ranks(1)) THEN
      logInfo(*) 'Initialize file ', TRIM(file_name),' to save fault data'
      logInfo(*) '<--------------------------------------------------------->'
      logInfo(*) ' '

      ! create hdf5 file with two groups
      CALL create_hdf5_file(file_name, 2, DISC%DynRup%hd_rec%grp_name, DISC%DynRup%hd_rec%file_id)
!#endif
      !logInfo(*) 'Initialize HDF5 file', TRIM(file_name),' to save fault receiver signals'

      ! write receiver attributes to root group
      !atrb_name(1) = 'dt'
      !atrb_name(2) = 'n_rec'
      !atrb_name(3) = 'n_samp'
      !atrb_name(4) = 't_end'
      !atrb_name(5) = 'max_iter'

      !atrb = (/ dt, DBLE(DISC%DynRup%DynRup_out_atPickpoint%nOutPoints), &)

      ! dataset for receiver locations
      location_component_names = (/'x','y','z'/)
      CALL create_hdf5_dset(DISC%DynRup%hd_rec%file_id, &
                            DISC%DynRup%hd_rec%grp_name(1), &
                            'location',  & ! name
                            2,  &  !rank
                            (/DISC%DynRup%DynRup_out_atPickpoint%nOutPoints, 3/) ) ! dimensions

      !CALL write_hdf5_atrb(DISC%DynRup%hd_rec%file_id, &
                           !'components', &
                           !length, &
                           !3, &
                           !loc_comp, &
                           !DISC%DynRup%hd_rec%grp_name(1), 'location')

      DO i = 1, 5
        if (DISC%DynRup%DynRup_out_atPickpoint%OutputMask(1).EQ.1) THEN
          CALL create_hdf5_dset(DISC%DynRup%hd_rec%file_id, &
                                DISC%DynRup%hd_rec%grp_name(1), & ! receiver
                                dataset_names(i), &
                                3, & !number of dimensions - rec.pont, time samples, components
                                dataset_size(i,:) )
        ENDIF
      ENDDO

      ! TODO write receiver positions
      CALL write_hdf5_data_serial(DISC%DynRup%hd_rec%file_id, &
                                  DISC%DynRup%hd_rec%grp_name(1), &
                                  'location', 2, &
                                  (/DISC%DynRup%DynRup_out_atPickpoint%nOutPoints, 3/), &
                                  (/ DISC%DynRup%DynRup_out_atPickpoint%RecPoint%x, &
                                     DISC%DynRup%DynRup_out_atPickpoint%RecPoint%y, &
                                     DISC%DynRup%DynRup_out_atPickpoint%RecPoint%z /))

      CALL close_hdf5_file(DISC%DynRup%hd_rec%file_id)
    ENDIF

  END SUBROUTINE init_hdf_fault_file


  !-----------------------------------------------------------------------!
  !-----------------------------------------------------------------------!
  !< initialize hdf5 receiver output
  !< generate file structure with empty datasets
  !< write attributes
  !< write receiver location
  SUBROUTINE ini_fault_receiver_hdf(EQN,MESH,DISC,IO,MPI,adaptive)

    USE common_receiver_mod
    USE common_fault_receiver_mod
    USE hdf_output_utils_mod
    !-------------------------------------------------------------
    IMPLICIT NONE
    !-------------------------------------------------------------
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    TYPE (tEquations)        :: EQN
    TYPE (tUnstructMesh)     :: MESH
    TYPE (tDiscretization)   :: DISC
    TYPE (tInputOutput)      :: IO
    TYPE (tMPI)              :: MPI
    LOGICAL, OPTIONAL        :: adaptive
    !-------------------------------------------------------------
    ! local Variables
    !-------------------------------------------------------------
    INTEGER,PARAMETER                 :: length = 5                        ! string length of receiver components
    !-------------------------------------------------------------
    CHARACTER(len=length)             :: comp(5,9), loc_comp(3), t_comp(1) ! field components, receiver location components, time components
    CHARACTER(len=length),ALLOCATABLE :: comp_sel(:,:)                     ! selected components
    CHARACTER(len=100),ALLOCATABLE    :: dset_name(:)                      ! selected hdf5 dataset names (stress, velocity, rotation, mtc, curl-divergence
    CHARACTER(len=100)                :: atrb_name(5)                      ! hdf5 root group attribute names
    INTEGER                           :: i, j, iRot                        ! loop counter
    INTEGER                           :: Vel_first, Vel_last               ! first/last position of velocity field components in outputmask
    INTEGER                           :: dset_count, any_rec, nranks       ! dataset counter, receiver switch, number of ranks
    INTEGER,ALLOCATABLE               :: ncomp_sel(:), dset_mask(:)        ! number of selected components, dataset selection mask
    INTEGER                           :: rec_grp, rec_comm
    INTEGER, ALLOCATABLE              :: rec_ranks(:)   ! mpi array of receiver switch, receiver ranks

    INTEGER                           :: number_of_ranks

    !------------------------------------



    ! initialization common to hdf5 ascii output
    CALL ini_common_fault_receiver(DISC, MESH)

    IF (DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output .EQ. .TRUE.) THEN

      ! create mpi communicator for domains with receivers
      CALL init_mpi_hdf_communication(DISC, MPI, rec_grp, rec_comm, &
                                      number_of_ranks, rec_ranks)

      ! initialize hdf5 file structure
      CALL init_hdf_fault_file(EQN, MESH, DISC, &
                               IO, MPI, rec_grp, rec_comm, &
                               number_of_ranks, rec_ranks)
    ENDIF


  END SUBROUTINE ini_fault_receiver_hdf


  SUBROUTINE write_fault_output_receiverwise_hdf(EQN, DISC, MESH, IO, MPI, MaterialVal, BND, time, dt)

    IMPLICIT NONE

    TYPE(tEquations)              :: EQN
    TYPE(tDiscretization)         :: DISC
    TYPE(tUnstructMesh)           :: MESH
    TYPE(tInputOutput)            :: IO                                         ! IO structure
    TYPE(tMPI)                    :: MPI                                        ! MPI
    REAL                          :: MaterialVal(MESH%nElem,EQN%nBackgroundVar) ! Local Mean Values
    TYPE(tBoundary)               :: BND                                        ! BND    data structure

    REAL    :: dt, time                                                       ! Timestep and time

    INTEGER       :: number_of_receivers
    INTEGER       :: receiver_index
    LOGICAL       :: buffer_is_full
    LOGICAL       :: last_time_step


    TYPE(tDynRup_output) :: dyn_rup_output

    !---------------------------------------------------

    number_of_receivers = DISC%DynRup%DynRup_out_atPickpoint%nDR_pick

    ! stortcut for a long name
    dyn_rup_output = DISC%DynRup%DynRup_out_atPickpoint


    ! loop over number of output receivers for this domain
    DO receiver_index = 1,number_of_receivers

      buffer_is_full = dyn_rup_output%CurrentPick(receiver_index).GE.dyn_rup_output%MaxPickStore
      last_time_step = ABS(DISC%EndTime-time).LE.(dt*1.005d0)

      IF (buffer_is_full.OR.last_time_step) THEN
        logInfo(*)'-----------------------------------------------------------'
        logInfo(*)'    PLOT fault receiver data for current time interval'
        !logInfo(*)'            start time :', IO%TmpTime(1,1)
        !logInfo(*)'            end time :', IO%TmpTime(1,IO%CurrentPick(j))
        logInfo(*)'-----------------------------------------------------------'

        !CALL write_hdf5_data_parallel(DISC%DynRup%hd_rec%file_id, &
                                      !DISC%DynRup%hd_rec%grp_name(1), & ! receiver
                                      !IO%hd_rec%dset_name(1), &
                                      !(/ chunk_size, IO%CurrentPick(j), EQN%nVar_sts /), &
                                      !IO%igloblocRecordPoint(rec_start:j)-1, &
                                      !IO%hd_rec%start(j), &
                                      !IO%TmpState_sts(rec_start:j,1:IO%CurrentPick(j),:), &
                                      !io_flag )

      ENDIF
    ENDDO ! iOutPoints = 1,DISC%DynRup%nOutPoints
  END SUBROUTINE


END MODULE hdf_faultoutput_mod
