!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Luca Passone (luca.passone AT kaust.edu.sa, http://ces.kaust.edu.sa/Pages/Luca-Passone.aspx)
!!
!! @section LICENSE
!! Copyright (c) 2010, SeisSol Group
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

MODULE hdf_output_utils_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  !----------------------------------------------------------------------------
  IMPLICIT NONE

  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE write_data
    MODULE PROCEDURE write_real_data, write_integer_data
  END INTERFACE

  INTERFACE create_hdf5_file
    MODULE PROCEDURE create_hdf5_file
  END INTERFACE

  INTERFACE create_hdf5_dset
    MODULE PROCEDURE create_hdf5_dset
  END INTERFACE

  INTERFACE open_hdf5_file_parallel
    MODULE PROCEDURE open_hdf5_file_parallel
  END INTERFACE

  INTERFACE write_hdf5_atrb
    MODULE PROCEDURE write_hdf5_real_atrb, write_hdf5_char_atrb
  END INTERFACE

  INTERFACE write_hdf5_data_serial
    MODULE PROCEDURE write_hdf5_data_serial
  END INTERFACE

  INTERFACE write_hdf5_data_parallel
    MODULE PROCEDURE write_hdf5_2d_data_parallel, write_hdf5_3d_data_parallel
  END INTERFACE

  INTERFACE close_hdf5_file
    MODULE PROCEDURE close_hdf5_file
  END INTERFACE

  !----------------------------------------------------------------------------

  PUBLIC  :: write_data, create_hdf5_file, create_hdf5_dset, open_hdf5_file_parallel, write_hdf5_atrb, write_hdf5_data_serial, &
             write_hdf5_data_parallel, close_hdf5_file
  !----------------------------------------------------------------------------

CONTAINS

  !> used for outputting integer data to an HDF5 file
  !!
  !! @param file_name The desired output ﬁle name.
  !! @param file_length The total size of the output ﬁle. This is a multidimensional array of rank equal to the size of the output (usually MESH%nNode total or MESH%nElem total) and depends on the choice to write the data separate or together.
  !! @param start_position Where to start writing in the ﬁle (for example MESH%nodeOﬀset or MESH%elemOﬀset).
  !! @param data_rank The rank of the data. If writing data separately, then data rank must be equal to 1 (see write_seperate parameter).
  !! @param data_dims The dimensions of the portion of the data the process is about to write (for example MESH%nNode).
  !! @param data_names The names of the variables being written to ﬁle.
  !! @param data The data to be written to ﬁle.
  !! @param write_separate when writing a matrix we can have three options: 0 write the matrix all together.  1 write the matrix by column (e.g.. data(iOutVar, :) ).  2 write the matrix by row (e.g. data(:, iOutVar) ).
  !! @param create If create equals to 1, a new ﬁle is created and the previous one overwritten. If create equals to 0, then data will be added to the ﬁle. For example, the connectivity and geometry are written to the same ﬁle, but have diﬀerent data types, ﬁrst the connectivity is written with create equals to 1, and then geometry is written with create equals to 0.
  !<
  SUBROUTINE write_integer_data(file_name, file_length, start_position, data_rank, data_dims, data_names, data, &
                                write_separate, create)

    USE HDF5

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    CHARACTER(*)   :: file_name
    INTEGER        :: data_rank
    INTEGER        :: write_separate ! wether or not to write a matrix all together, or separate. If 0, it will write everything together. If 1 it will write column first and if 2 row first.
    INTEGER        :: create
    INTEGER(HID_T) :: file_id       ! File identifier
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id      ! Property list identifier
    INTEGER        :: ioutvar
    INTEGER        :: error, error_n  ! Error flags

    INTEGER(HSIZE_T), DIMENSION(*)                  :: file_length, start_position
    INTEGER(HSIZE_T), DIMENSION(2)                  :: data_dims ! MAX len x num_elem
    INTEGER, DIMENSION(data_dims(1),data_dims(2))   :: data ! Data to write
    CHARACTER(LEN=20), DIMENSION(*)                 :: data_names

    INTENT(IN) :: file_name, file_length, start_position, data_rank, data_dims, data_names, data, write_separate, create

    !PRINT *, 'File length ', file_length(1)
    !PRINT *, 'start position ', start_position(1)
    !PRINT *, 'data_rank ', data_rank
    !PRINT *, 'Dims = ', data_dims(1)
    !PRINT *, file_name

    !
    ! Initialize FORTRAN predefined datatypes
    !
    CALL h5open_f(error)

    !
    ! Setup file access property list with parallel I/O access.
    !
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)

    !
    ! Create the file collectively.
    !
    IF (create .EQ. 1) THEN
      CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    ELSE
      CALL h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    END IF

    CALL h5pclose_f(plist_id, error)

    IF (write_separate.EQ.0) THEN
      ! Create the data space for the  dataset.
      !
      CALL h5screate_simple_f(data_rank, file_length, filespace, error)
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, data_names(1), H5T_NATIVE_INTEGER, filespace, &
                       dset_id, error)
      CALL h5sclose_f(filespace, error)
      CALL h5screate_simple_f(data_rank, data_dims, memspace, error)
      !
      ! Select hyperslab in the file.
      !
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, start_position, data_dims, error)

      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

      !
      ! Write the dataset collectively.
      !
      CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data, data_dims, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
      !
      !Close dataspaces.
      !
      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      !
      ! Close the dataset and property list.
      !
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
    ELSE
      DO iOutVar=1, data_dims(3-write_separate)
        ! Create the data space for the  dataset.
        CALL h5screate_simple_f(data_rank, file_length, filespace, error)
        !PRINT *, "Error 5: ", error
        ! Create the dataset with default properties.
        !PRINT *, 'data_names(x): ', data_names(iOutVar)
        CALL h5dcreate_f(file_id, data_names(iOutVar), H5T_NATIVE_INTEGER, filespace, &
                         dset_id, error)
        !PRINT *, "Error 6: ", error
        CALL h5sclose_f(filespace, error)
        !
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        !
        CALL h5screate_simple_f(data_rank, data_dims(write_separate), memspace, error)
        !
        ! Select hyperslab in the file.
        !
        CALL h5dget_space_f(dset_id, filespace, error)
        CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, start_position, data_dims(write_separate), error)

        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        !
        ! Write the dataset collectively.
        !
        IF (write_separate .EQ. 1) THEN
          CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data(iOutVar, :), data_dims(write_separate:write_separate), error, &
                        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
        ELSEIF (write_separate .EQ. 2) THEN
          CALL h5dwrite_f(dset_id, H5T_NATIVE_INTEGER, data(:, iOutVar), data_dims(write_separate:write_separate), error, &
                        file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
        ELSE
          PRINT *, 'No output can be written to HDF5 output because the wrong data dimension was specified. Only 0, 1 and 2 are valid.'
        ENDIF
        !
        !Close dataspaces.
        !
        CALL h5sclose_f(filespace, error)
        CALL h5sclose_f(memspace, error)
        !
        ! Close the dataset and property list.
        !
        CALL h5dclose_f(dset_id, error)
        CALL h5pclose_f(plist_id, error)
      END DO
    END IF
    !
    ! Close the file.
    !
    CALL h5fclose_f(file_id, error)
    !
    ! Close FORTRAN predefined datatypes.
    !
    CALL h5close_f(error)

  END SUBROUTINE write_integer_data

  !> used for writing real(8) data to an HDF5 file
  !!
  !! @param file_name The desired output ﬁle name.
  !! @param file_length The total size of the output ﬁle. This is a multidimensional array of rank equal to the size of the output (usually MESH%nNode total or MESH%nElem total) and depends on the choice to write the data separate or together.
  !! @param start_position Where to start writing in the ﬁle (for example MESH%nodeOﬀset or MESH%elemOﬀset).
  !! @param data_rank The rank of the data. If writing data separately, then data rank must be equal to 1 (see write_seperate parameter).
  !! @param data_dims The dimensions of the portion of the data the process is about to write (for example MESH%nNode).
  !! @param data_names The names of the variables being written to ﬁle.
  !! @param data The data to be written to ﬁle.
  !! @param write_separate when writing a matrix we can have three options: 0 write the matrix all together.  1 write the matrix by column (e.g.. data(iOutVar, :) ).  2 write the matrix by row (e.g. data(:, iOutVar) ).
  !! @param create If create equals to 1, a new ﬁle is created and the previous one overwritten. If create equals to 0, then data will be added to the ﬁle. For example, the connectivity and geometry are written to the same ﬁle, but have diﬀerent data types, ﬁrst the connectivity is written with create equals to 1, and then geometry is written with create equals to 0.
  !<
  SUBROUTINE write_real_data(file_name, file_length, start_position, data_rank, data_dims, data_names, data, write_separate, create)

    USE HDF5

    IMPLICIT NONE

    INCLUDE 'mpif.h'

    CHARACTER(*)   :: file_name
    INTEGER        :: data_rank
    INTEGER        :: write_separate ! wether or not to write a matrix all together, or separate. If 0, it will write everything together. If 1 it will write column first and if 2 row first.
    INTEGER(HID_T) :: file_id       ! File identifier
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id      ! Property list identifier
    INTEGER        :: ioutvar, create
    INTEGER        :: error, error_n  ! Error flags

    INTEGER(HSIZE_T), DIMENSION(*)                :: file_length, start_position
    INTEGER(HSIZE_T), DIMENSION(2)                :: data_dims ! MAX len x num_elem
    CHARACTER(LEN=20), DIMENSION(*)               :: data_names
    REAL(4), DIMENSION(data_dims(2),data_dims(1)) :: data ! Data to write

    INTENT(IN) :: file_name, file_length, start_position, data_rank, data_dims, data_names, data, write_separate, create

    !
    ! Initialize FORTRAN predefined datatypes
    !
    CALL h5open_f(error)
    !
    ! Setup file access property list with parallel I/O access.
    !
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, MPI_COMM_WORLD, MPI_INFO_NULL, error)
    !
    ! Create the file collectively.
    !
    IF (create .EQ. 1) THEN
      CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, error, access_prp = plist_id)
    ELSE
      CALL h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    END IF

    CALL h5pclose_f(plist_id, error)

    IF (write_separate.EQ.0) THEN
      !
      ! Create the data space for the  dataset.
      !
      CALL h5screate_simple_f(data_rank, file_length, filespace, error)
      !
      ! Create the dataset with default properties.
      !
      CALL h5dcreate_f(file_id, data_names(1), H5T_NATIVE_REAL, filespace, &
                       dset_id, error)
      CALL h5sclose_f(filespace, error)
      CALL h5screate_simple_f(data_rank, data_dims, memspace, error)
      !
      ! Select hyperslab in the file.
      !
      CALL h5dget_space_f(dset_id, filespace, error)
      CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, start_position, data_dims, error)

      CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
      !
      ! Write the dataset collectively.
      !
      CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data, data_dims, error, &
                      file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
      !
      !Close dataspaces.
      !
      CALL h5sclose_f(filespace, error)
      CALL h5sclose_f(memspace, error)
      !
      ! Close the dataset and property list.
      !
      CALL h5dclose_f(dset_id, error)
      CALL h5pclose_f(plist_id, error)
    ELSE
      DO iOutVar=1, data_dims(3-write_separate)
        !
        ! Create the data space for the  dataset.
        !
        CALL h5screate_simple_f(data_rank, file_length, filespace, error)
        !
        ! Create the dataset with default properties.
        !
        CALL h5dcreate_f(file_id, data_names(iOutVar), H5T_NATIVE_REAL, filespace, &
                         dset_id, error)
        CALL h5sclose_f(filespace, error)
        !
        ! Each process defines dataset in memory and writes it to the hyperslab
        ! in the file.
        !
        CALL h5screate_simple_f(data_rank, data_dims(write_separate), memspace, error)
        !
        ! Select hyperslab in the file.
        !
        CALL h5dget_space_f(dset_id, filespace, error)
        CALL h5sselect_hyperslab_f (filespace, H5S_SELECT_SET_F, start_position, data_dims(write_separate), error)

        CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
        CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)

        !
        ! Write the dataset collectively.
        !
        IF (write_separate .EQ. 1) THEN
          CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data(iOutVar, :), data_dims(write_separate:write_separate), error, &
                          file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
        ELSEIF (write_separate .EQ. 2) THEN
          CALL h5dwrite_f(dset_id, H5T_NATIVE_REAL, data(:, iOutVar), data_dims(write_separate:write_separate), error, &
                          file_space_id = filespace, mem_space_id = memspace, xfer_prp = plist_id)
        ELSE
          PRINT *, 'No output can be written to HDF5 output because the wrong data dimension was specified. Only 0, 1 and 2 are valid.'
        ENDIF
        !
        !Close dataspaces.
        !
        CALL h5sclose_f(filespace, error)
        CALL h5sclose_f(memspace, error)
        !
        ! Close the dataset and property list.
        !
        CALL h5dclose_f(dset_id, error)
        CALL h5pclose_f(plist_id, error)
      END DO
    END IF
    !
    ! Close the file.
    !
    CALL h5fclose_f(file_id, error)
    !
    ! Close FORTRAN predefined datatypes.
    !
    CALL h5close_f(error)

  END SUBROUTINE write_real_data
  !-----------------------------------------------------------------------------------------
  !< receiver output
  !-----------------------------------------------------------------------------------------
  !
  !< create hdf5 file
  !< create multiple groups
  !< put out the file id
  SUBROUTINE create_hdf5_file(file_name,ngrp,grp_name,file_id)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER(HID_T)                        :: file_id       ! file identifier
    INTEGER(HID_T)                        :: grp_id        ! group identifier
    !------------------------------------------------------------------
    INTEGER                               :: ngrp          ! number of groups
    CHARACTER(*)                          :: file_name     ! file name
    CHARACTER(*),DIMENSION(ngrp)          :: grp_name      ! group name
    INTEGER                               :: error         ! error flag
    INTEGER                               :: i             ! loop counter
    !------------------------------------------------------------------
    INTENT(IN)                            :: file_name, ngrp, grp_name
    INTENT(OUT)                           :: file_id
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! create file.
    CALL h5fcreate_f(file_name, H5F_ACC_TRUNC_F, file_id, error)
    ! group
    DO i=1,ngrp
      ! create
      CALL h5gcreate_f(file_id, grp_name(i), grp_id, error)
      ! close
      CALL h5gclose_f(grp_id, error)
    ENDDO
    ! close interface.
    CALL h5close_f(error)

  END SUBROUTINE create_hdf5_file
  !
  !< create an empty dataset with predefined dataspace in a group of the hdf5 file
  !< create an attribute in this dataset
  SUBROUTINE create_hdf5_dset(file_id, grp_name, dset_name, dset_rank, dset_dims)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER                               :: dset_rank     ! dataset rank
    !------------------------------------------------------------------
    INTEGER(HID_T)                        :: file_id       ! file identifier
    INTEGER(HID_T)                        :: grp_id        ! group identifier
    INTEGER(HID_T)                        :: dset_id       ! Dataset identifier
    INTEGER(HID_T)                        :: attr_id       ! Attribute identifie
    INTEGER(HID_T)                        :: space_id      ! Dataspace identifier
    INTEGER(HID_T)                        :: type_id       ! Datatype identifier
    INTEGER(HSIZE_T),DIMENSION(dset_rank) :: dset_hdf_dims ! Dataset dimensions
    !------------------------------------------------------------------
    CHARACTER(*)                          :: grp_name      ! group name
    CHARACTER(*)                          :: dset_name     ! dataset name
    INTEGER,DIMENSION(dset_rank)          :: dset_dims     ! dataset dimensions
    INTEGER                               :: error         ! error handle
    !------------------------------------------------------------------
    INTENT(IN)                            :: file_id, grp_name, dset_name, dset_rank, dset_dims
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open existing group.
    CALL h5gopen_f(file_id, grp_name, grp_id, error)
    ! dataspace dimensions
    dset_hdf_dims = dset_dims
    ! Create data space for the dataset.
    CALL h5screate_simple_f(dset_rank, dset_hdf_dims, space_id, error)
    ! Make the dataspace empty
    CALL h5sselect_none_f(space_id, error)
    ! create dataset
    CALL h5dcreate_f(grp_id, dset_name, H5T_NATIVE_REAL, space_id, dset_id, error)
    ! close dataspace
    CALL h5sclose_f(space_id, error)
    ! Close the dataset.
    CALL h5dclose_f(dset_id, error)
    ! End access to the group.
    CALL h5gclose_f(grp_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

  END SUBROUTINE create_hdf5_dset
  !
  !< open hdf5 file for parallel file access
  !< use mpi communicator of mpi domains which contain receivers
  SUBROUTINE open_hdf5_file_parallel(file_name,comm,file_id)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER(HID_T)                        :: file_id       ! file identifier
    INTEGER(HID_T)                        :: plist_id      ! Property list identifier
    !------------------------------------------------------------------
    CHARACTER(*)                          :: file_name     ! file name
    INTEGER                               :: comm          ! mpi communicator
    INTEGER                               :: error         ! error handle
    !------------------------------------------------------------------
    INTENT(IN)                            :: file_name, comm
    INTENT(OUT)                           :: file_id
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Setup file access property list with parallel I/O access.
    CALL h5pcreate_f(H5P_FILE_ACCESS_F, plist_id, error)
    CALL h5pset_fapl_mpio_f(plist_id, comm, MPI_INFO_NULL, error)
    ! Create file.
    CALL h5fopen_f(file_name, H5F_ACC_RDWR_F, file_id, error, access_prp = plist_id)
    ! Close property list
    CALL h5pclose_f(plist_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

  END SUBROUTINE open_hdf5_file_parallel
  !
  !< write real attribute to group or dataset
  SUBROUTINE write_hdf5_real_atrb(file_id, atrb_name, atrb_len, atrb_dims, atrb_data ,grp_name, dset_name)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER(HID_T)                 :: file_id       ! file identifier
    INTEGER(HID_T)                 :: grp_id        ! group identifier
    INTEGER(HID_T)                 :: dset_id       ! dataset identifier
    INTEGER(HID_T)                 :: attr_id       ! Attribute identifier
    INTEGER(HID_T)                 :: space_id      ! Attribute Dataspace identifier
    INTEGER(HID_T)                 :: type_id       ! Attribute type identifier
    INTEGER(HSIZE_T),DIMENSION(1)  :: atrb_hdf_dims ! HDF5 data dimension
    !------------------------------------------------------------------
    CHARACTER(*),OPTIONAL          :: grp_name      ! group name
    CHARACTER(*),OPTIONAL          :: dset_name     ! dataset name
    CHARACTER(*)                   :: atrb_name     ! attribute name
    REAL,DIMENSION(:)              :: atrb_data     ! attribute value
    INTEGER                        :: atrb_dims     ! attribute size
    INTEGER                        :: atrb_len      ! attribute length !!! is not used for real data !!!
    INTEGER                        :: rank = 1      ! Attribute rank
    INTEGER                        :: error         ! error handle
    !------------------------------------------------------------------
    INTENT(IN)                     :: file_id, atrb_name, atrb_len, atrb_dims, atrb_data, grp_name, dset_name
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open group.
    IF( present( grp_name ) ) THEN
      ! Open existing group.
      CALL h5gopen_f(file_id, grp_name, grp_id, error)
    ELSE
      ! Open root group.
      CALL h5gopen_f(file_id, "/", grp_id, error)
    ENDIF
    IF( present( dset_name ) ) THEN
      ! Open existing dataset.
      CALL h5dopen_f(grp_id, dset_name, dset_id, error)
    ENDIF
    ! data dimensions
    atrb_hdf_dims(1) = atrb_dims
    ! Create scalar data space for the attribute.
    CALL h5screate_simple_f(rank, atrb_hdf_dims, space_id, error)
    ! Create datatype.
    CALL h5tcopy_f(H5T_NATIVE_DOUBLE, type_id, error)
    ! Create attribute.
    IF( present( dset_name ) ) THEN
      CALL h5acreate_f(dset_id, atrb_name, type_id, space_id, attr_id, error)
    ELSE
      CALL h5acreate_f(grp_id, atrb_name, type_id, space_id, attr_id, error)
    ENDIF
    ! Write the attribute data.
    CALL h5awrite_f(attr_id, type_id, atrb_data, atrb_hdf_dims, error)
    ! Close the attribute.
    CALL h5aclose_f(attr_id, error)
    ! Terminate access to the data space.
    CALL h5sclose_f(space_id, error)
    ! Terminate access to the datatype.
    CALL h5tclose_f(type_id, error)
    ! Close the dataset if present.
    IF( present( dset_name ) ) THEN
      CALL h5dclose_f(dset_id, error)
    ENDIF
    ! End access to the group.
    CALL h5gclose_f(grp_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

  END SUBROUTINE write_hdf5_real_atrb
  !
  !< write character attribute to group or dataset
  SUBROUTINE write_hdf5_char_atrb(file_id, atrb_name, atrb_len, atrb_dims, atrb_data ,grp_name, dset_name)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER(HID_T)                 :: file_id       ! file identifier
    INTEGER(HID_T)                 :: grp_id        ! group identifier
    INTEGER(HID_T)                 :: dset_id       ! dataset identifier
    INTEGER(HID_T)                 :: attr_id       ! Attribute identifier
    INTEGER(HID_T)                 :: space_id      ! Attribute Dataspace identifier
    INTEGER(HID_T)                 :: type_id       ! Attribute type identifier
    INTEGER(HSIZE_T),DIMENSION(1)  :: atrb_hdf_dims ! Attribute dimensions
    INTEGER(SIZE_T)                :: atrb_hdf_len  ! Attribute length
    !------------------------------------------------------------------
    CHARACTER(*),OPTIONAL          :: grp_name      ! group name
    CHARACTER(*),OPTIONAL          :: dset_name     ! dataset name
    CHARACTER(*)                   :: atrb_name     ! attribute name
    CHARACTER(*),DIMENSION(:)      :: atrb_data     ! attribute data
    INTEGER                        :: atrb_dims     ! Attribute dimension
    INTEGER                        :: atrb_len      ! Number of characters
    INTEGER                        :: rank = 1      ! Attribute rank
    INTEGER                        :: error         ! error handle
    !------------------------------------------------------------------
    INTENT(IN)                     :: file_id, atrb_name, atrb_len, atrb_dims, atrb_data, grp_name, dset_name
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open existing group.
    IF( present( grp_name ) ) THEN
      CALL h5gopen_f(file_id, grp_name, grp_id, error)
    ELSE
      ! if no group name is given, open the root group.
      CALL h5gopen_f(file_id, "/", grp_id, error)
    ENDIF
    ! Open existing dataset.
    IF( present( dset_name ) ) THEN
      CALL h5dopen_f(grp_id, dset_name, dset_id, error)
    ENDIF
    ! attribute dimensions
    atrb_hdf_dims(1) = atrb_dims
    atrb_hdf_len = atrb_len
    ! Create scalar data space for the attribute.
    CALL h5screate_simple_f(rank, atrb_hdf_dims, space_id, error)
    ! Create character datatype.
    CALL h5tcopy_f(H5T_NATIVE_CHARACTER, type_id, error)
    ! Make it the length of atrb_length.
    CALL h5tset_size_f(type_id, atrb_hdf_len, error)
    ! Create attribute.
    IF( present( dset_name ) ) THEN
      ! create a dataset attribute
      CALL h5acreate_f(dset_id, atrb_name, type_id, space_id, attr_id, error)
    ELSE
      ! if no dataset is given, create a group attribute
      CALL h5acreate_f(grp_id, atrb_name, type_id, space_id, attr_id, error)
    ENDIF
    ! Write the attribute data.
    CALL h5awrite_f(attr_id, type_id, atrb_data, atrb_hdf_dims, error)
    ! Close the attribute.
    CALL h5aclose_f(attr_id, error)
    ! Terminate access to the data space.
    CALL h5sclose_f(space_id, error)
    ! Terminate access to the datatype.
    CALL h5tclose_f(type_id, error)
    ! Close the dataset if present.
    IF( present( dset_name ) ) THEN
      CALL h5dclose_f(dset_id, error)
    ENDIF
    ! End access to the group.
    CALL h5gclose_f(grp_id, error)
    ! Close fortran interface.
    CALL h5close_f(error)

  END SUBROUTINE write_hdf5_char_atrb
  !
  !< write a dataset to one group in serial
  SUBROUTINE write_hdf5_data_serial(file_id, grp_name, dset_name, rank, dset_dims, data)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER                        :: rank          ! data rank
    !------------------------------------------------------------------
    INTEGER(HID_T)                 :: file_id       ! file identifier
    INTEGER(HID_T)                 :: grp_id        ! Group identifier
    INTEGER(HID_T)                 :: dset_id       ! Dataset identifier
    INTEGER(HID_T)                 :: space_id      ! Dataspace identifier
    INTEGER(HID_T)                 :: type_id       ! Datatype identifier
    INTEGER(HSIZE_T),DIMENSION(rank) :: dset_hdf_dims ! Dataset dimensions
    !------------------------------------------------------------------
    CHARACTER(*)                   :: grp_name      ! group name
    CHARACTER(*)                   :: dset_name     ! dataset name
    REAL,DIMENSION(*)              :: data          ! data
    INTEGER,DIMENSION(rank)        :: dset_dims     ! dataset dimensions
    INTEGER                        :: error         ! error handle
    INTEGER                        :: i             ! count
    !------------------------------------------------------------------
    INTENT(IN)                     :: file_id, grp_name, dset_name, rank, dset_dims, data
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open existing group.
    CALL h5gopen_f(file_id, grp_name, grp_id, error)
    ! Open existing dataset.
    CALL h5dopen_f(grp_id, dset_name, dset_id, error)
    ! dataspace dimensions
    DO i=1,rank
      dset_hdf_dims(i) = dset_dims(i)
    ENDDO
    ! Create scalar data space for the attribute.
    CALL h5screate_simple_f(rank, dset_hdf_dims, space_id, error)
    ! Create string datatype of length atrb_length.
    CALL h5tcopy_f(H5T_NATIVE_DOUBLE, type_id, error)
    ! Write the attribute data.
    CALL h5dwrite_f(dset_id, type_id, data, dset_hdf_dims, error)
    ! Terminate access to the data space.
    CALL h5sclose_f(space_id, error)
    ! Terminate access to the datatype.
    CALL h5tclose_f(type_id, error)
    ! Close the dataset.
    CALL h5dclose_f(dset_id, error)
    ! End access to the group.
    CALL h5gclose_f(grp_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

  END SUBROUTINE write_hdf5_data_serial
  !
  !< write a 2d dataset to one group in parallel
  SUBROUTINE write_hdf5_2d_data_parallel(file_id, grp_name, dset_name, buf, start_r1, start_r2, data, io_flag)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER                        :: rank = 2      ! data rank
    !------------------------------------------------------------------
    INTEGER(HID_T)                 :: file_id       ! file identifier
    INTEGER(HID_T)                 :: grp_id        ! group identifier
    INTEGER(HID_T)                 :: dset_id       ! Dataset identifier
    INTEGER(HID_T)                 :: plist_id      ! property list identifier
    INTEGER(HID_T)                 :: type_id       ! Datatype identifier
    INTEGER(HID_T)                 :: dset_space_id ! Dataset dataspace identifier
    INTEGER(HID_T)                 :: memspace_id  ! Dataset dataspace identifier
    INTEGER(HSIZE_T),DIMENSION(2)  :: buf_hdf       ! Buffer dimensions
    INTEGER(HSIZE_T),DIMENSION(2)  :: count_hdf     ! number of blocks
    INTEGER(HSIZE_T),DIMENSION(2)  :: start_hdf     ! start of current buffer
    !------------------------------------------------------------------
    CHARACTER(*)                   :: grp_name      ! group name
    CHARACTER(*)                   :: dset_name     ! dataset name
    REAL*8,DIMENSION(:,:)          :: data          ! data
    INTEGER,DIMENSION(2)           :: buf           ! buffer dimensions (chunk)
    INTEGER,DIMENSION(:)           :: start_r1      ! start position rank 1
    INTEGER                        :: start_r2      ! start position rank 2
    INTEGER                        :: io_flag       ! collective = 1, independent = 0
    INTEGER                        :: error         ! error handle
    INTEGER                        :: i             ! count
    !------------------------------------------------------------------
    INTENT(IN)                     :: file_id, grp_name, dset_name, buf, start_r1, start_r2, data, io_flag
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open existing group.
    CALL h5gopen_f(file_id, grp_name, grp_id, error)
    ! open existing dataset
    CALL h5dopen_f(grp_id, dset_name, dset_id, error)
    ! number of blocks
    count_hdf(1) = 1
    count_hdf(2) = buf(2)
    ! buffer dimensions
    buf_hdf(1) = buf(1)
    buf_hdf(2) = buf(2)
    ! Create memory space
    CALL h5screate_simple_f(rank, buf_hdf, memspace_id, error)
    ! Create datatype
    CALL h5tcopy_f(H5T_NATIVE_DOUBLE, type_id, error)
    ! Select whole dataspace of selected dataset
    CALL h5dget_space_f(dset_id, dset_space_id, error)
    ! select subset of that dataset by an irregular hyperslab (add multiple hyperslab selections)
    DO i=1,buf(1)
      start_hdf = (/ start_r1(i), start_r2 /)
      IF (i.EQ.1) THEN
        ! overwrite dataset selection
        CALL h5sselect_hyperslab_f(dset_space_id, H5S_SELECT_SET_F, start_hdf, count_hdf, error)
      ELSE
        ! add new hyperslab
        CALL h5sselect_hyperslab_f(dset_space_id, H5S_SELECT_OR_F, start_hdf, count_hdf, error)
      ENDIF
    ENDDO
    ! Create property list for collective dataset write. independent H5FD_MPIO_INDEPENDENT_F
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    ! set type of data output
    IF (io_flag.EQ.1) THEN
      ! collective
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    ELSE
      ! independent
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    ENDIF
    ! write data to selected hyperslab dataspace
    CALL h5dwrite_f(dset_id, type_id, data, buf_hdf, error, &
                    file_space_id = dset_space_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    ! flush data to disk
    !CALL h5fflush_f(file_id,H5F_SCOPE_GLOBAL_F,error)
    ! Close property list
    CALL h5pclose_f(plist_id, error)
    ! Terminate access to the dataspace.
    CALL h5sclose_f(dset_space_id, error)
    ! Terminate access to the memspace.
    CALL h5sclose_f(memspace_id, error)
    ! Terminate access to the datatype.
    CALL h5tclose_f(type_id, error)
    ! Close the dataset.
    CALL h5dclose_f(dset_id, error)
    ! Close access to the group.
    CALL h5gclose_f(grp_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

  END SUBROUTINE write_hdf5_2d_data_parallel
  !
  !< write a 3d dataset to one group in parallel
  SUBROUTINE write_hdf5_3d_data_parallel(file_id, grp_name, dset_name, buf, start_r1, start_r2, data, io_flag)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER                        :: rank = 3      ! data rank
    !------------------------------------------------------------------
    INTEGER(HID_T)                 :: file_id       ! file identifier
    INTEGER(HID_T)                 :: grp_id        ! group identifier
    INTEGER(HID_T)                 :: dset_id       ! Dataset identifier
    INTEGER(HID_T)                 :: plist_id      ! property list identifier
    INTEGER(HID_T)                 :: type_id       ! Datatype identifier
    INTEGER(HID_T)                 :: dset_space_id ! Dataset dataspace identifier
    INTEGER(HID_T)                 :: memspace_id  ! Dataset dataspace identifier
    INTEGER(HSIZE_T),DIMENSION(3)  :: buf_hdf       ! Buffer dimensions
    INTEGER(HSIZE_T),DIMENSION(3)  :: count_hdf     ! number of blocks
    INTEGER(HSIZE_T),DIMENSION(3)  :: start_hdf     ! start of current buffer
    !------------------------------------------------------------------
    CHARACTER(*)                   :: grp_name      ! group name
    CHARACTER(*)                   :: dset_name     ! dataset name
    REAL*8,DIMENSION(:,:,:)        :: data          ! data
    INTEGER,DIMENSION(3)           :: buf           ! buffer dimensions (chunk)
    INTEGER,DIMENSION(:)           :: start_r1      ! start position rank 1
    INTEGER                        :: start_r2      ! start position rank 2
    INTEGER                        :: io_flag       ! collective = 1, independent = 0
    INTEGER                        :: error         ! error handle
    INTEGER                        :: i             ! count
    !------------------------------------------------------------------
    INTENT(IN)                     :: file_id, grp_name, dset_name, buf, start_r1, start_r2, data, io_flag
    !------------------------------------------------------------------
    ! Initialize FORTRAN interface.
    CALL h5open_f(error)
    ! Open existing group.
    CALL h5gopen_f(file_id, grp_name, grp_id, error)
    ! open existing dataset
    CALL h5dopen_f(grp_id, dset_name, dset_id, error)
    ! number of blocks in each dimension
    count_hdf(1) = 1
    count_hdf(2) = buf(2)
    count_hdf(3) = buf(3)
    ! buffer dimensions
    buf_hdf(1) = buf(1)
    buf_hdf(2) = buf(2)
    buf_hdf(3) = buf(3)
    ! Create memory space
    CALL h5screate_simple_f(rank, buf_hdf, memspace_id, error)
    ! Create datatype
    CALL h5tcopy_f(H5T_NATIVE_DOUBLE, type_id, error)
    ! Select whole dataspace of selected dataset
    CALL h5dget_space_f(dset_id, dset_space_id, error)
    ! select subset of that dataset by an irregular hyperslab (add multiple hyperslab selections)
    DO i=1,buf(1)
      start_hdf = (/ start_r1(i), start_r2, 0 /)
      IF (i.EQ.1) THEN
        ! overwrite dataset selection
        CALL h5sselect_hyperslab_f(dset_space_id, H5S_SELECT_SET_F, start_hdf, count_hdf, error)
      ELSE
        ! add new hyperslab
        CALL h5sselect_hyperslab_f(dset_space_id, H5S_SELECT_OR_F, start_hdf, count_hdf, error)
      ENDIF
    ENDDO
    ! Create property list for collective dataset write
    CALL h5pcreate_f(H5P_DATASET_XFER_F, plist_id, error)
    ! set type of data output
    IF (io_flag.EQ.1) THEN
      ! collective
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_COLLECTIVE_F, error)
    ELSE
      ! independent
      CALL h5pset_dxpl_mpio_f(plist_id, H5FD_MPIO_INDEPENDENT_F, error)
    ENDIF
    ! write data collectively to selected hyperslab dataspace
    CALL h5dwrite_f(dset_id, type_id, data, buf_hdf, error, &
                    file_space_id = dset_space_id, mem_space_id = memspace_id, xfer_prp = plist_id)
    ! flush data to disk
    !CALL h5fflush_f(file_id,H5F_SCOPE_GLOBAL_F,error)
    ! Close property list
    CALL h5pclose_f(plist_id, error)
    ! Terminate access to the dataspace.
    CALL h5sclose_f(dset_space_id, error)
    ! Terminate access to the memspace.
    CALL h5sclose_f(memspace_id, error)
    ! Terminate access to the datatype.
    CALL h5tclose_f(type_id, error)
    ! Close the dataset.
    CALL h5dclose_f(dset_id, error)
    ! Close access to the group.
    CALL h5gclose_f(grp_id, error)
    ! Close FORTRAN interface.
    CALL h5close_f(error)

  END SUBROUTINE write_hdf5_3d_data_parallel
  !
  !< close hdf5 file
  !< put out the file id
  !< use mpi communicator of mpi domains which contain receivers
  SUBROUTINE close_hdf5_file(file_id)
    !------------------------------------------------------------------
    USE HDF5
    !------------------------------------------------------------------
    IMPLICIT NONE
    !------------------------------------------------------------------
    INCLUDE 'mpif.h'
    !------------------------------------------------------------------
    INTEGER(HID_T)                        :: file_id       ! file identifier
    !------------------------------------------------------------------
    INTEGER                               :: error         ! error flag
    !------------------------------------------------------------------
    INTENT(IN)                            :: file_id
    !------------------------------------------------------------------
    ! Initialize fortran interface.
    CALL h5open_f(error)
    ! close file.
    CALL h5fclose_f(file_id, error)
    ! close interface.
    CALL h5close_f(error)

  END SUBROUTINE close_hdf5_file

END MODULE hdf_output_utils_mod
