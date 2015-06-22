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

#include <Initializer/preProcessorMacros.fpp>

MODULE hd_output_mod
  !----------------------------------------------------------------------------
  USE TypesDef
  !----------------------------------------------------------------------------
  IMPLICIT NONE

  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE close_hdf_wavefield
    MODULE PROCEDURE close_hdf_wavefield
  END INTERFACE

  INTERFACE init
    MODULE PROCEDURE init
  END INTERFACE

  INTERFACE write_hd_data
    MODULE PROCEDURE write_hd_data
  END INTERFACE

  INTERFACE InterpolateBaryToNode
    MODULE PROCEDURE InterpolateBaryToNode
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC  :: close_hdf_wavefield
  PUBLIC  :: write_hd_data
  PUBLIC  :: init
  PRIVATE :: compute_refined_mesh
  !----------------------------------------------------------------------------

CONTAINS

!> This routine writes footer of the XDMF descriptor ﬁle.
!! @param IO variable containing Input and Output information
!<
  SUBROUTINE  close_hdf_wavefield (IO,MPI)

    USE HDF5

    IMPLICIT NONE

#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !------------------------------------------------------------------
    TYPE (tInputOutput) :: IO
    TYPE (tMPI)         :: MPI
    !------------------------------------------------------------------
    INTEGER             :: error
    CHARACTER(LEN=100)  :: filename
    !------------------------------------------------------------------
    INTENT(IN)          :: IO, MPI
    !------------------------------------------------------------------

#ifdef PARALLEL
    !Close the XDMF file.
    IF (MPI%myrank.EQ.0) THEN
#endif
      WRITE(filename,'(a,a)') TRIM(IO%OutputFile),'.xdmf'
      OPEN (UNIT = 5, FILE = filename, POSITION='APPEND')
      WRITE (5,'(A)') '      </Grid>'
      WRITE (5,'(A)') '   </Domain>'
      WRITE (5,'(A)') '</Xdmf>'
      CLOSE (5)
#ifdef PARALLEL
    ENDIF
#endif
  END SUBROUTINE close_hdf_wavefield                    !

!> This subroutine writes the node data to file and updates the XDMF descriptor with the information for loading the time step.
!! @param OptionalFields variable containing optional fields
!! @param MESH
!! @param IO
!! @param timestep the timestep parameter is used in conjunction with IO%OutputFile for constructing the output ﬁle name.
!<
  SUBROUTINE write_hd_data(OptionalFields, MESH, IO, timestep, EQN, DISC, MPI)

    USE HDF5
    USE hdf_output_utils_mod
    IMPLICIT NONE

#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !------------------------------------------------------------------
    TYPE (tUnstructOptionalFields) :: OptionalFields
    TYPE (tUnstructMesh)           :: MESH
    TYPE (tInputOutput)            :: IO
    TYPE (tEquations)              :: EQN
    TYPE (tDiscretization)         :: DISC
    TYPE (tMPI)                    :: MPI
    !------------------------------------------------------------------
    INTEGER(HSIZE_T),DIMENSION(1)  :: count = 0
    INTEGER(HSSIZE_T),DIMENSION(1) :: offset = (0)
    INTEGER(HSIZE_T),DIMENSION(1)  :: file_length, start_position
    INTEGER(HSIZE_T),DIMENSION(2)  :: data_dims ! MAX len x num_elem
    !------------------------------------------------------------------
    CHARACTER(LEN=20),DIMENSION(IO%nOutputMask) :: data_names
    CHARACTER(LEN=100)             :: filename  ! File name
    CHARACTER(LEN=20)              :: dummy
    INTEGER                        :: timestep, pos
    INTEGER                        :: ielem, iNode, iOutVar, OutVarStart
    INTEGER                        :: n_output_var = 0   ! Number of variables to output
    INTEGER                        :: rank = 1 ! Dataset rank
    REAL                           :: NodeVal, write_time_start, write_time_end
    REAL (4), ALLOCATABLE          :: data (:,:)  ! Data to write
#ifdef PARALLEL
    INTEGER                        :: info
#endif
    !------------------------------------------------------------------
    INTENT(IN)                     :: MESH, IO, timestep,OptionalFields, EQN, DISC, MPI
    !------------------------------------------------------------------

    ! register epik/scorep function
    EPIK_FUNC_REG("write_hd_data")
    SCOREP_USER_FUNC_DEFINE()

    ! start epik/scorep function
    EPIK_FUNC_START()
    SCOREP_USER_FUNC_BEGIN("write_hd_data")

    OutVarStart = 4

    !find out the number of variables
    n_output_var = 0
    DO iOutVar=OutVarStart,IO%nOutputMask
       IF(IO%OutputMask(iOutVar)) THEN
          n_output_var = n_output_var + 1
          dummy = TRIM(IO%TitleMask(iOutVar))
          data_names(n_output_var) = dummy(3:LEN_TRIM(IO%TitleMask(iOutVar))-1)
       ENDIF
    ENDDO


    IF (IO%hd_out%refinement .EQ. 0) THEN
      ALLOCATE (data(n_output_var, IO%hd_out%node_n))
      file_length(1) = IO%hd_out%node_total
      start_position(1) = IO%hd_out%node_offset
      data_dims(1) = IO%hd_out%node_n
      data_dims(2) = n_output_var
      DO inode=1,IO%hd_out%node_n
        n_output_var = 1
        DO iOutVar=OutVarStart,IO%nOutputMask
            IF(IO%OutputMask(iOutVar)) THEN
                !logInfo(*) "RANK: ", MPI%myrank, " - n_output_var: ", n_output_var, " - iNode: ", iNode
                CALL InterpolateBaryToNode(              &  !
                    OptionalFields%FieldMask(iOutVar)%PTR, & !
                    OptionalFields%weight , &
                    NodeVal               , &   !
                    iNode                 , &   !
                    MESH                    )   !
                data(n_output_var, iNode) = NodeVal
                n_output_var = n_output_var + 1
            ENDIF
        ENDDO
      ENDDO
    ELSEIF (IO%hd_out%refinement .EQ. 1) THEN !output at barycenters
      !logInfo(*) 'elem_n = ', IO%hd_out%elem_n, 'n_output_var', n_output_var, 'EQN%nVar', EQN%nVar
      ALLOCATE (data(EQN%nVar, IO%hd_out%elem_n))
      file_length(1) = IO%hd_out%elem_total
      start_position(1) = IO%hd_out%elem_offset
      data_dims(1) = IO%hd_out%elem_n
      data_dims(2) = EQN%nVar
      !logInfo(*) 'Calling evaluate solution'
      CALL evaluate_solution(data, IO%hd_out, EQN, DISC, MESH)
      !logInfo(*) 'Returned'
    ELSEIF (IO%hd_out%refinement .EQ. 2) THEN !Have 4 subdivisions, and output at their baricenters
      !logInfo(*) 'elem_n = ', IO%hd_out%elem_n, 'n_output_var', n_output_var, 'EQN%nVar', EQN%nVar
      ALLOCATE (data(EQN%nVar, IO%hd_out%elem_n))
      file_length(1) = IO%hd_out%elem_total
      start_position(1) = IO%hd_out%elem_offset
      data_dims(1) = IO%hd_out%elem_n
      data_dims(2) = EQN%nVar
      !logInfo(*) 'Calling evaluate solution'
      CALL evaluate_solution(data, IO%hd_out, EQN, DISC, MESH)
      !logInfo(*) 'Returned'
    ELSEIF (IO%hd_out%refinement .EQ. 3) THEN !Have 4 subdivisions, and output at their baricenters
      !logInfo(*) 'elem_n = ', IO%hd_out%elem_n, 'n_output_var', n_output_var, 'EQN%nVar', EQN%nVar
      ALLOCATE (data(EQN%nVar, IO%hd_out%elem_n))
      file_length(1) = IO%hd_out%elem_total
      start_position(1) = IO%hd_out%elem_offset
      data_dims(1) = IO%hd_out%elem_n
      data_dims(2) = EQN%nVar
      !logInfo(*) 'Calling evaluate solution'
      CALL evaluate_solution(data, IO%hd_out, EQN, DISC, MESH)
      !logInfo(*) 'Returned'
    ELSE
      logError(*)  'this mesh refinemente is not supported'
      STOP
    END IF
    WRITE(filename,'(a,a,i6.6,a)') TRIM(IO%OutputFile), '_', timestep, '.h5'
    ! to stay consitent with the loginfo of data_output.f90
    logInfo(*)'                               file   :  ', TRIM(filename)
    logInfo(*)'-----------------------------------------------------------'
    !logInfo(*) MPI%myrank, ' at barrier'
    !CALL MPI_Barrier(MPI_COMM_WORLD, mpierror)
    !logInfo(*) MPI%myrank, 'passed barrier'
    !CALL CPU_TIME(write_time_start)
    CALL write_data(filename, file_length, start_position, rank, data_dims, data_names, data, 1, 1)
    !CALL CPU_TIME(write_time_end)
    !logInfo(*) MPI%myrank, ' started at ', write_time_start, 'ended at', write_time_end, ' to write ', data_dims(1), ' x ', data_dims(2), ' of REAL(4)'
    !logInfo(*) 'Done'
    !
    ! Deallocate data buffer.
    !
    DEALLOCATE(data)

! Write XDMF data do go with the HDF
! First, find the name of the mesh file.
#ifdef PARALLEL
    IF (MPI%myrank.EQ.0) THEN
#endif
      !logInfo(*) 'Append snapshot to XDMF file'
      pos = INDEX(IO%OutputFile, "/", BACK = .TRUE.) + 1
      WRITE(filename,'(a,a)') TRIM(IO%OutputFile), '.xdmf'
      !logInfo(*) 'Appending to ', filename
      OPEN (UNIT = 5, FILE = filename, POSITION='APPEND')
      WRITE (5,'(a,i6.6,a)') '    <Grid Name="', timestep, '" GridType="Uniform">'
      WRITE (5,'(A)') '       <Topology Reference="/Xdmf/Domain/Topology[1]"/>'
      WRITE (5,'(A)') '       <Geometry Reference="/Xdmf/Domain/Geometry[1]"/>'
      IF (IO%hd_out%refinement .EQ. 0 ) THEN
        DO iOutVar=OutVarStart,IO%nOutputMask
         IF(IO%OutputMask(iOutVar)) THEN
           WRITE (5,'(A,A,A)') '      <Attribute Name=', TRIM(IO%TitleMask(iOutVar)), ' Center="Node">'
           WRITE (5, *)     '       <DataItem DataType="Float" Precision="4" Format="HDF" Dimensions="', IO%hd_out%node_total, '">'
           dummy = TRIM(IO%TitleMask(iOutVar))
           WRITE (5,'(A,A,i6.6,A,A)') TRIM(IO%OutputFile(pos:)), '_', timestep, '.h5:/',dummy(3:LEN_TRIM(IO%TitleMask(iOutVar))-1)
           WRITE (5,'(A)') '        </DataItem>'
           WRITE (5,'(A)') '      </Attribute>'
         END IF
       END DO
     ELSE
       DO iOutVar=OutVarStart, EQN%nVar+OutVarStart
         IF(IO%OutputMask(iOutVar)) THEN
           WRITE (5,'(A,A,A)') '      <Attribute Name=', TRIM(IO%TitleMask(iOutVar)), ' Center="Cell">'
           WRITE (5, *)     '       <DataItem DataType="Float" Precision="4" Format="HDF" Dimensions="', IO%hd_out%elem_total, '">'
           dummy = TRIM(IO%TitleMask(iOutVar))
           WRITE (5,'(A,A,i6.6,A,A)') TRIM(IO%OutputFile(pos:)), '_', timestep, '.h5:/',dummy(3:LEN_TRIM(IO%TitleMask(iOutVar))-1)
           WRITE (5,'(A)') '        </DataItem>'
           WRITE (5,'(A)') '      </Attribute>'
         END IF
       END DO
     END IF
     ! This is the partition part (to see which processor got which part of the mesh)
     WRITE (5,'(A)') '      <Attribute Name="proc_id" Center="Node">'
     WRITE (5, *)     '       <DataItem DataType="Float" Precision="4" Format="HDF" Dimensions="', IO%hd_out%node_total, '">'
     WRITE (5,'(A,A)') TRIM(IO%OutputFile(pos:)), '_partition.h5:/partition'
     WRITE (5,'(A)') '        </DataItem>'
     WRITE (5,'(A)') '      </Attribute>'
     WRITE (5,'(A)') '    </Grid>'
     CLOSE(5)
#ifdef PARALLEL
   ENDIF
#endif

    ! end epik function data_output
    EPIK_FUNC_END()
    SCOREP_USER_FUNC_END()

  END SUBROUTINE write_hd_data

!> This routine takes care of three things:
!! 1 write the geometry of the mesh to an HDF5 ﬁle.
!! 2 write the connectivity to the same HDF5 ﬁle as above.
!! 3 Initialise important HDF5 structures.
 !!
!! It is important that this subroutine is called ﬁrst before any HDF5 output is done as it is responsible for initialising the basic HDF5 structures used for output.
!! @param IO
!! @param MPI
!! @param MESH
!<
  SUBROUTINE init(IO, MPI, MESH, DISC)

    USE HDF5
    USE hdf_output_utils_mod

    IMPLICIT NONE
#ifdef PARALLEL
    INCLUDE 'mpif.h'
#endif
    !------------------------------------------------------------------
    TYPE (tInputOutput)            :: IO                                       !
    TYPE (tUnstructMesh)           :: MESH                                     !
    TYPE (tMPI)                    :: MPI                                      !
    TYPE (tDiscretization)         :: DISC
    !------------------------------------------------------------------
    CHARACTER(LEN=100)              :: file_name! File name
    CHARACTER(LEN=20), DIMENSION(1) :: dsetname   ! Dataset name
    !------------------------------------------------------------------
    INTEGER(HID_T) :: file_id       ! File identifier
    INTEGER(HID_T) :: dset_id       ! Dataset identifier
    INTEGER(HID_T) :: filespace     ! Dataspace identifier in file
    INTEGER(HID_T) :: memspace      ! Dataspace identifier in memory
    INTEGER(HID_T) :: plist_id      ! Property list identifier

    INTEGER(HSIZE_T), DIMENSION(2) :: dimsf ! Dataset dimensions.
    INTEGER(HSIZE_T), DIMENSION(2) :: count
    INTEGER(HSIZE_T), DIMENSION(2) :: offset
    !------------------------------------------------------------------
    INTEGER :: rank = 2, mpierror, pos ! Dataset rank
    INTEGER, DIMENSION(2) :: total_size
    REAL(4), ALLOCATABLE           :: real_data (:,: )  ! Data to write
    INTEGER, ALLOCATABLE           :: int_data (:,: )  ! Data to write
    !------------------------------------------------------------------
    INTENT(IN)    :: MPI, DISC
    INTENT(INOUT) :: MESH, IO
    !------------------------------------------------------------------
#ifdef PARALLEL
    ! STEP 1 Find out how many output points are ncessary
    logInfo(*) 'Initializing hdf5 variables'
    CALL MPI_SCAN(MESH%nNode, MESH%nodeOffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, MPI%iErr)
    MESH%nodeOffset = MESH%nodeOffset - MESH%nNode
    CALL MPI_SCAN(MESH%nElem, MESH%elemOffset, 1, MPI_INTEGER, MPI_SUM, MPI_COMM_WORLD, MPI%iErr)
    MESH%elemOffset = MESH%elemOffset - MESH%nElem
    !logInfo(*) MPI%myrank, " has ", MESH%nNode, " nodes, and the offset is ", MESH%nodeOffset
    !logInfo(*) MPI%myrank, " has ", MESH%nElem, " elements, and the offset is ", MESH%elemOffset
    ! Broadcast the total number of nodes, including boundry nodes.
    IF (MPI%myrank == MPI%nCPU -1) THEN
       MESH%nNode_total = MESH%nodeOffset + MESH%nNode
       MESH%nElem_total = MESH%elemOffset + MESH%nElem
    ENDIF
    total_size(1) = MESH%nNode_total
    total_size(2) = MESH%nElem_total
    CALL MPI_Bcast(total_size, 2, MPI_INTEGER, MPI%nCPU-1, MPI_COMM_WORLD, MPI%iErr)
    MESH%nNode_total = total_size(1)
    MESH%nElem_total = total_size(2)
    !CALL MPI_Bcast(MESH%nElem_total, 1, MPI_INTEGER, MPI%nCPU-1, MPI_COMM_WORLD, MPI%iErr)
#else
    MESH%nodeOffset = 0
    MESH%elemOffset = 0
#endif
    ALLOCATE(IO%hd_out)
    SELECT CASE (IO%hd_out%refinement)
      CASE(0)
        IO%hd_out%elem_n      = MESH%nElem
        IO%hd_out%elem_total  = MESH%nElem_total
        IO%hd_out%elem_offset = MESH%elemOffset
        IO%hd_out%node_n      = MESH%nNode
        IO%hd_out%node_total  = MESH%nNode_total
        IO%hd_out%node_offset = MESH%nodeOffset
      CASE(1)
        IO%hd_out%elem_n      = MESH%nElem
        IO%hd_out%elem_total  = MESH%nElem_total
        IO%hd_out%elem_offset = MESH%elemOffset
        IO%hd_out%node_n      = MESH%nNode
        IO%hd_out%node_total  = MESH%nNode_total
        IO%hd_out%node_offset = MESH%nodeOffset
        IO%hd_out%points_n    = 1
      CASE(2)
        IO%hd_out%elem_n      = MESH%nElem * 4
        IO%hd_out%elem_total  = MESH%nElem_total * 4
        IO%hd_out%elem_offset = MESH%elemOffset * 4
        IO%hd_out%node_n      = MESH%nNode + IO%hd_out%elem_n
        IO%hd_out%node_total  = MESH%nNode_total + IO%hd_out%elem_total
        IO%hd_out%node_offset = MESH%nodeOffset + MESH%elemOffset
        IO%hd_out%points_n    = 4
      CASE(3)
        IO%hd_out%elem_n      = MESH%nElem * 8
        IO%hd_out%elem_total  = MESH%nElem_total * 8
        IO%hd_out%elem_offset = MESH%elemOffset * 8
        IO%hd_out%node_n      = MESH%nNode + MESH%nElem * 6
        IO%hd_out%node_total  = MESH%nNode_total + MESH%nElem_total * 6
        IO%hd_out%node_offset = MESH%nodeOffset + MESH%elemOffset * 6
        IO%hd_out%points_n    = 8
    END SELECT
    logInfo(*) 'Done initializing hdf5 variables'

    ! STEP 2 Write connectivity and geometry
    file_name = TRIM(IO%OutputFile) // '_mesh.h5'
    logInfo(*) 'Beginning to write connectivity to ', TRIM(file_name)

    CALL compute_refined_mesh(IO, MESH, real_data, int_data)
    !dimsf = (/4, IO%hd_out%total_hd_out_points/)
    dimsf = (/4, IO%hd_out%elem_total/)

    count(1) = 4
    count(2) = IO%hd_out%elem_n
    offset(1) = 0
    offset(2) = IO%hd_out%elem_offset
    dsetname = 'connect'
    logInfo(*) 'About ot write', dimsf(2), offset(2), count(2), int_data(1,10)
    CALL write_data(file_name, dimsf, offset, 2, count, dsetname, &
                    int_data((/1,2,3,4/), :), 0, 1)

    logInfo(*) 'Finished writing connectivity and beginning to write geometry'

    ! HDF5 begin geometry
    dimsf = (/3, IO%hd_out%node_total/)
    !
    ! Create the dataset with default properties.
    !
    dsetname = 'geometry'
    !
    ! Each process defines dataset in memory and writes it to the hyperslab
    !  CALL TrafoXYZ2XiEtaZeta(IO%hd_out%recpoint(1)%xi,IO%hd_out%recpoint(1)%eta,IO%hd_out%recpoint(1)%zeta,io_x,io_y,io_z,xV,yV,zV,MESH%LocalVrtxType(iElem))
    ! in the file.
    !
    count(1) = 3
    count(2) = IO%hd_out%node_n
    offset(1) = 0
    offset(2) = IO%hd_out%node_offset

    CALL write_data(file_name, dimsf, offset, 2, count, dsetname,&
                             real_data, 0, 0)
    logInfo(*) 'Done writing geometry'
    DEALLOCATE(real_data)
    DEALLOCATE(int_data)

! STEP 3 Start writing the XDMF

    IF (MPI%myrank.EQ.0) THEN
      pos = INDEX(IO%OutputFile, '/', BACK = .TRUE.) + 1
      WRITE(file_name,'(a,a)') TRIM(IO%OutputFile),'.xdmf'
      OPEN (UNIT = 5, FILE = file_name)!, POSITION='APPEND')
      WRITE (5,'(A)') '<?xml version="1.0" ?>'
      WRITE (5,'(A)') '<Xdmf Version="2.0">'
      WRITE (5,'(A)') '  <Domain>'
      WRITE (5, *)     '     <Topology TopologyType="Tetrahedron" NumberOfElements="', IO%hd_out%elem_total, '">'
      WRITE (5, *)     '       <DataItem Format="HDF" DataType="Int" Dimensions="', IO%hd_out%elem_total, '4" >'
      WRITE (5,'(2A)') TRIM(IO%OutputFile(pos:)), '_mesh.h5:/connect'
      WRITE (5,'(A)') '        </DataItem>'
      WRITE (5,'(A)') '      </Topology>'
      WRITE (5, *)     '     <Geometry name="geo" GeometryType="XYZ" NumberOfElements="', IO%hd_out%node_total, '">'
      WRITE (5, *)     '       <DataItem Format="HDF" Dimensions="', IO%hd_out%node_total, '3">'
      WRITE (5,'(2A)') TRIM(IO%OutputFile(pos:)), '_mesh.h5:/geometry'
      WRITE (5,'(A)') '        </DataItem>'
      WRITE (5,'(A)') '      </Geometry>'
      WRITE (5,'(A)') '      <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">'
      WRITE (5,'(A)') '         <Time TimeType="Hyperslab">'
      WRITE (5,'(A)') '            <DataItem Format="XML" NumberType="Float" Dimensions="3">'
      IF (IO%outInterval%printIntervalCriterion .EQ. 1) THEN
         WRITE (5, *) '               0.0', IO%outInterval%Interval, DISC%MaxIteration
      ENDIF
      IF (IO%outInterval%printIntervalCriterion .EQ. 2) THEN
          WRITE (5, *) '               0.0', IO%outInterval%TimeInterval, DISC%EndTime
      ENDIF
      WRITE (5,'(A)') '            </DataItem>'
      WRITE (5,'(A)') '         </Time>'
      CLOSE(5)
    ENDIF


! STEP 4 write partition
#ifdef PARALLEL
    ALLOCATE (int_data(IO%hd_out%node_n,1))
    int_data=MPI%myrank
    WRITE(file_name,'(a,a)') TRIM(IO%OutputFile),  '_partition.h5'
    WRITE(dsetname, '(A)') 'partition'
    ! Write the dataset collectively.
    count = IO%hd_out%node_n
    offset(1) = IO%hd_out%node_offset
    dimsf =  (/IO%hd_out%node_total, 1/)
    CALL write_data(file_name, dimsf, offset, 1, count, dsetname, int_data, 0, 1)
    DEALLOCATE(int_data)
#endif
  END SUBROUTINE init

!> This subroutine is currently just duplicated code from plot_ﬁelds_mod due to cyclic dependencies. In the pipeline there is a version that allows for tet 10 to be used for output, and this function will have to be modiﬁed.
!! @todo add support for tet10, and find a better place for this routine to avoid cyclic dependencies.
!<
  SUBROUTINE InterpolateBaryToNode(BaryField,weight,NodeValue,iNode,MESH)
    !--------------------------------------------------------------------------

    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! argument list declaration
    TYPE (tUnstructMesh)   :: MESH
    REAL                   :: BaryField(MESH%nElem)
    REAL                   :: weight(MESH%nElem)
    REAL                   :: NodeValue
    INTEGER                :: iNode
    ! local variable declaration
    INTEGER                :: iElem
    INTEGER                :: ElemNr
    REAL                   :: vol
    REAL                   :: work
    !--------------------------------------------------------------------------
    INTENT(IN)  :: BaryField, weight, iNode, MESH
    INTENT(OUT) :: NodeValue
    !--------------------------------------------------------------------------
    ! interpolating cell average field into vertex field!
    !                                                   !
    vol  = 0.0                                          !
    work = 0.0                                          !
    !                                                   !
    DO ielem = 1,MESH%VRTX%NrOfElementsConnected(inode) !
       !
       ElemNr = MESH%VRTX%Element(inode,ielem)          !
       !
       vol  = vol                                     & !
            + MESH%ELEM%Volume(ElemNr)                & !
            * weight(ElemNr)                            !
       !                                                !
       work = work                                    & !
            + BaryField(  ElemNr)                     & !
            * MESH%ELEM%Volume(ElemNr)                & !
            * weight(     ElemNr)                       ! cartesian:
       !                                                ! weight = 1.0 = const.
       !                                                ! cylindrical:
       !                                                ! weight = radius
    ENDDO                                               !
    !                                                   !
    NodeValue = work/vol                                !
    !
    IF(ABS(NodeValue).LT.1e-20) THEN
       NodeValue = 0.
    ENDIF
    !                                                   !
  END SUBROUTINE InterpolateBaryToNode                  !

  SUBROUTINE compute_refined_mesh(IO, MESH, refined_nodes, refined_tets)
    TYPE (tUnstructMesh)           :: MESH
    TYPE (tInputOutput)            :: IO
    REAL(4), ALLOCATABLE           :: refined_nodes(:,:), t(:,:,:)
    INTEGER, ALLOCATABLE           :: refined_tets(:,:)
    INTEGER                        :: i,j,k, n_extra_nodes, n_extra_tets, counter

    INTENT(IN)  :: MESH
    INTENT(OUT) :: refined_nodes, refined_tets
    INTENT(INOUT)  :: IO

    logInfo(*) 'Computing refined mesh'
    IF (IO%hd_out%refinement .EQ. 0 ) THEN
      ALLOCATE(refined_nodes(3, MESH%nNode))
      ALLOCATE(refined_tets(4, MESH%nElem))
      refined_nodes(:,1:MESH%nNode) = MESH%VRTX%xyNode
      refined_tets(:,1:MESH%nElem) = MESH%ELEM%Vertex
    ELSEIF (IO%hd_out%refinement .EQ. 1 ) THEN
      ALLOCATE(refined_nodes(3, MESH%nNode))
      ALLOCATE(refined_tets(4, MESH%nElem))
      refined_nodes(:,1:MESH%nNode) = MESH%VRTX%xyNode
      refined_tets(:,1:MESH%nElem) = MESH%ELEM%Vertex
      ALLOCATE (IO%hd_out%recpoint(1))
      IO%hd_out%recpoint(1)%xi    = 0.25
      IO%hd_out%recpoint(1)%eta   = 0.25
      IO%hd_out%recpoint(1)%zeta  = 0.25
    ELSEIF (IO%hd_out%refinement .EQ. 2) THEN
      !If the refinement level is 2, then we need to add 1 node to the barycenter of each tetrahedron, and the number of tetrahedrons multiplies by 4
      n_extra_nodes = MESH%nElem
      n_extra_tets = MESH%nElem * 4
      ALLOCATE(refined_nodes(3, n_extra_nodes + MESH%nNode))
      ALLOCATE(refined_tets(4, n_extra_tets))
      refined_nodes(:,1:MESH%nNode) = MESH%VRTX%xyNode
      DO i=1, MESH%nElem
        refined_tets(:,(i-1)*4+1) = (/MESH%ELEM%Vertex(1,i), MESH%ELEM%Vertex(2,i), MESH%nNode + i, MESH%ELEM%Vertex(4,i)/)
        refined_tets(:,(i-1)*4+2) = (/MESH%ELEM%Vertex(1,i), MESH%ELEM%Vertex(2,i), MESH%ELEM%Vertex(3,i), MESH%nNode + i/)
        refined_tets(:,(i-1)*4+3) = (/MESH%ELEM%Vertex(4,i), MESH%ELEM%Vertex(1,i), MESH%ELEM%Vertex(3,i), MESH%nNode + i/)
        refined_tets(:,(i-1)*4+4) = (/MESH%ELEM%Vertex(2,i), MESH%ELEM%Vertex(4,i), MESH%ELEM%Vertex(3,i), MESH%nNode + i/)
        !This could be avoided, but for sake of consistency I do it
        refined_nodes(:, i+MESH%nNode) = MESH%ELEM%xybary(:,i)
      END DO
      ALLOCATE (IO%hd_out%recpoint(4))
      IO%hd_out%recpoint(1)%xi    = 0.3125
      IO%hd_out%recpoint(1)%eta   = 0.0625
      IO%hd_out%recpoint(1)%zeta  = 0.3125
      IO%hd_out%recpoint(2)%xi    = 0.3125
      IO%hd_out%recpoint(2)%eta   = 0.3125
      IO%hd_out%recpoint(2)%zeta  = 0.0625
      IO%hd_out%recpoint(3)%xi    = 0.0625
      IO%hd_out%recpoint(3)%eta   = 0.3125
      IO%hd_out%recpoint(3)%zeta  = 0.3125
      IO%hd_out%recpoint(4)%xi    = 0.3125
      IO%hd_out%recpoint(4)%eta   = 0.3125
      IO%hd_out%recpoint(4)%zeta  = 0.3125
    ELSEIF (IO%hd_out%refinement .EQ. 3) THEN
      ! Based on figure 2 from  "QUALITY LOCAL REFINEMENT OF TETRAHEDRAL MESHES BASED ON 8-SUBTETRAHEDRON SUBDIVISION"
      n_extra_nodes = MESH%nElem * 6
      n_extra_tets = MESH%nElem * 8
      ALLOCATE(refined_nodes(3, n_extra_nodes + MESH%nNode))
      ALLOCATE(refined_tets(4, n_extra_tets))
      ALLOCATE(t(3,4,4))
      refined_nodes(:,1:MESH%nNode) = MESH%VRTX%xyNode
      counter = 1
      DO i=1, MESH%nElem
        !Create the new nodes
        DO j=1,4
          DO k=1,j
            IF (k .LT. j) THEN
              refined_nodes(:,MESH%nNode+counter) = (MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(k,i)) + MESH%VRTX%xyNode(:,MESH%ELEM%Vertex(j,i))) * 0.5
              counter = counter + 1;
            ENDIF
          END DO
        END DO
        !Use the new nodes for the new tets!
        refined_tets(:,(i-1)*8+1) = (/MESH%ELEM%Vertex(1,i), &
                                      MESH%nNode + counter - 6, &
                                      MESH%nNode + counter - 5, &
                                      MESH%nNode + counter - 3/)
        refined_tets(:,(i-1)*8+2) = (/MESH%nNode + counter - 6, &
                                      MESH%ELEM%Vertex(2,i),    &
                                      MESH%nNode + counter - 4, &
                                      MESH%nNode + counter - 2/)
        refined_tets(:,(i-1)*8+3) = (/MESH%nNode + counter - 5, &
                                      MESH%nNode + counter - 4, &
                                      MESH%ELEM%Vertex(3,i),    &
                                      MESH%nNode + counter - 1/)
        refined_tets(:,(i-1)*8+4) = (/MESH%nNode + counter - 3, &
                                      MESH%nNode + counter - 2, &
                                      MESH%nNode + counter - 1, &
                                      MESH%ELEM%Vertex(4,i)/)
        refined_tets(:,(i-1)*8+5) = (/MESH%nNode + counter - 6, &
                                      MESH%nNode + counter - 2, &
                                      MESH%nNode + counter - 3, &
                                      MESH%nNode + counter - 5/)
        refined_tets(:,(i-1)*8+6) = (/MESH%nNode + counter - 6, &
                                      MESH%nNode + counter - 4, &
                                      MESH%nNode + counter - 2, &
                                      MESH%nNode + counter - 5/)
        refined_tets(:,(i-1)*8+7) = (/MESH%nNode + counter - 1, &
                                      MESH%nNode + counter - 5, &
                                      MESH%nNode + counter - 4, &
                                      MESH%nNode + counter - 2/)
        refined_tets(:,(i-1)*8+8) = (/MESH%nNode + counter - 1, &
                                      MESH%nNode + counter - 3, &
                                      MESH%nNode + counter - 5, &
                                      MESH%nNode + counter - 2/)
      END DO
      ALLOCATE (IO%hd_out%recpoint(8))
      IO%hd_out%recpoint(1)%xi    = 0.1250
      IO%hd_out%recpoint(1)%eta   = 0.1250
      IO%hd_out%recpoint(1)%zeta  = 0.1250
      IO%hd_out%recpoint(2)%xi    = 0.6250
      IO%hd_out%recpoint(2)%eta   = 0.1250
      IO%hd_out%recpoint(2)%zeta  = 0.1250
      IO%hd_out%recpoint(3)%xi    = 0.1250
      IO%hd_out%recpoint(3)%eta   = 0.6250
      IO%hd_out%recpoint(3)%zeta  = 0.1250
      IO%hd_out%recpoint(4)%xi    = 0.1250
      IO%hd_out%recpoint(4)%eta   = 0.1250
      IO%hd_out%recpoint(4)%zeta  = 0.6250
      IO%hd_out%recpoint(5)%xi    = 0.2500
      IO%hd_out%recpoint(5)%eta   = 0.1250
      IO%hd_out%recpoint(5)%zeta  = 0.2500
      IO%hd_out%recpoint(6)%xi    = 0.3750
      IO%hd_out%recpoint(6)%eta   = 0.2500
      IO%hd_out%recpoint(6)%zeta  = 0.1250
      IO%hd_out%recpoint(7)%xi    = 0.2500
      IO%hd_out%recpoint(7)%eta   = 0.3750
      IO%hd_out%recpoint(7)%zeta  = 0.2500
      IO%hd_out%recpoint(8)%xi    = 0.1250
      IO%hd_out%recpoint(8)%eta   = 0.2500
      IO%hd_out%recpoint(8)%zeta  = 0.3750
    ENDIF
    refined_tets = refined_tets((/1,2,3,4/), :) -1+ IO%hd_out%node_offset
  END SUBROUTINE compute_refined_mesh

  SUBROUTINE evaluate_solution(data, hd_out,EQN,DISC,MESH)
    !--------------------------------------------------------------------------!
    USE DGBasis_mod
    !--------------------------------------------------------------------------!
    IMPLICIT NONE
    !--------------------------------------------------------------------------!
    ! Argument list declaration
    REAL(4), ALLOCATABLE  , INTENT(INOUT)      :: data(:,:)
    TYPE (tEquations)     , INTENT(IN)         :: EQN                          !< EQN global variable
    TYPE (tDiscretization), INTENT(IN)         :: DISC                         !< DISC global variable
    TYPE(thd_output),       INTENT(IN)         :: hd_out                                         !< HD output data
    TYPE (tUnstructMesh),   INTENT(IN)         :: MESH
    !--------------------------------------------------------------------------!
    ! Local variable declaration
    INTEGER                  :: iDegFr,nDegFr,nVar,nPoly,iElem
    INTEGER                  :: ipoint, nhd_out_points
    REAL                     :: xi, eta, zeta
    REAL                     :: dgvar(DISC%Galerkin%nDegFr,EQN%nVar)
    REAL(4)                  :: state(EQN%nVar)
    REAL, POINTER            :: cPoly3D(:,:,:,:,:)         => NULL()
    INTEGER, POINTER         :: NonZeroCPoly(:,:)          => NULL()
    INTEGER, POINTER         :: NonZeroCPolyIndex(:,:,:,:) => NULL()
    REAL                     :: phi

  nVar              = EQN%nVar
  nDegFr            = DISC%Galerkin%nDegFr
  nPoly             = DISC%Galerkin%nPoly
  cPoly3D           => DISC%Galerkin%cPoly3D_Tet
  NonZeroCPoly      => DISC%Galerkin%NonZeroCPoly_Tet
  NonZeroCPolyIndex => DISC%Galerkin%NonZeroCPolyIndex_Tet

  DO iElem = 1, MESH%nElem
    DO ipoint = 1, hd_out%points_n

    ! load element ID and data
    dgvar(1:nDegFr,1:nVar) = DISC%Galerkin%dgvar(1:nDegFr,1:EQN%nVar,iElem,1)

    ! load position
    ! logInfo(*) 'loading position'
    xi   = hd_out%recpoint(ipoint)%xi
    eta  = hd_out%recpoint(ipoint)%eta
    zeta = hd_out%recpoint(ipoint)%zeta

    ! get state
    ! logInfo(*) 'Setting state'
    state(:) = 0.0D0
    DO iDegFr = 1, nDegFr
       CALL BaseFunc3D(phi,iDegFr,xi,eta,zeta,nPoly,cPoly3D,NonZeroCPoly,NonZeroCPolyIndex)
       state(1:nVar) = state(1:nVar) + phi*dgvar(iDegFr,:)
    ! logInfo(*) 'phi: ', phi, 'state: ', state(1:nVar)!, 'dgvar(',iDegFr,')', dgvar(iDegFr,:,1,1)
    ENDDO

    ! store
    ! logInfo(*) 'Storing to data to index ', (iElem-1)*4+ipoint
    data(:,(iElem-1)*hd_out%points_n+ipoint) = state(:)
   ENDDO
  ENDDO ! iElem = 1,nhd_out_tets

  END SUBROUTINE evaluate_solution
END MODULE hd_output_mod
