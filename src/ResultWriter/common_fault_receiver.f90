!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Dmitry Pinaev (pinaev AT in.tum.de)
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
!! Module contains intialization functions common to both hdf and non-hdf output

! logging macro

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE common_fault_receiver_mod

  USE TypesDef

  IMPLICIT NONE
  PRIVATE

  !---- public functions

  PUBLIC :: ini_common_fault_receiver


  !---- privatec functions

  ! checks if all the points belong to a fault
  PRIVATE :: locate_receiver_points

  ! removes points which are not in the fault
  PRIVATE :: reallocate_and_fill_recpoint_array

  ! check how many variables are to output
  PRIVATE :: count_output_variables_allocate_arrays


CONTAINS

  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  SUBROUTINE locate_receiver_points(DISC, MESH, LocalRecPoint, in)

    USE TypesDef
    USE create_fault_rotationmatrix_mod
    USE DGBasis_mod

    IMPLICIT NONE

    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!
    TYPE(tUnstructMesh)     :: MESH                   ! Mesh structure        !
    TYPE(tUnstructPoint), ALLOCATABLE       :: LocalRecPoint(:)
    INTEGER                 :: in

    !local variables
    REAL                    :: tolerance  ! tolerance if the receiver belongs to this element
    REAL                    :: xmin,xmax,ymin,ymax,zmin,zmax ! Maximum extend of computational domain
    INTEGER                 :: i, j, iFault, iElem, iSide
    REAL                    :: io_x, io_y, io_z              ! temp store of receiver location
    REAL                    :: xV(MESH%nVertexMax)
    REAL                    :: yV(MESH%nVertexMax)
    REAL                    :: zV(MESH%nVertexMax)
    REAL                    :: xi,eta,zeta

    ! tolerance if the receiver belongs to this element
    tolerance = 1.0e-5
    !
    ! Maximum extend of computational domain
    xmin = MINVAL( MESH%VRTX%xyNode(1,:) )
    xmax = MAXVAL( MESH%VRTX%xyNode(1,:) )
    ymin = MINVAL( MESH%VRTX%xyNode(2,:) )
    ymax = MAXVAL( MESH%VRTX%xyNode(2,:) )
    zmin = MINVAL( MESH%VRTX%xyNode(3,:) )
    zmax = MAXVAL( MESH%VRTX%xyNode(3,:) )
    !
    ! Assume that the point is not inside the domain
    DISC%DynRup%DynRup_out_atPickpoint%RecPoint(:)%inside =.false.
    DISC%DynRup%DynRup_out_atPickpoint%RecPoint(:)%index  = -1
    !
    ! ALLOCATE local list
    ALLOCATE(LocalRecPoint(DISC%DynRup%DynRup_out_atPickpoint%nOutPoints))
    LocalRecPoint(:)%inside = .false.
    LocalRecPoint(:)%index = -1
    LocalRecPoint(:)%globalreceiverindex = -1
    LocalRecPoint(:)%X = 0.0D0
    LocalRecPoint(:)%Y = 0.0D0
    LocalRecPoint(:)%Z = 0.0D0

    ! Loop over all fault output receivers (i)
    in = 0
    DO i=1, DISC%DynRup%DynRup_out_atPickpoint%nOutPoints
       !
       io_x = DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%X
       io_y = DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%Y
       io_z = DISC%DynRup%DynRup_out_atPickpoint%RecPoint(i)%Z
       !
       IF (     (io_x.GE.xmin).and.(io_x.LE.xmax)   &
           .and.(io_y.GE.ymin).and.(io_y.LE.ymax)   &
           .and.(io_z.GE.zmin).and.(io_z.LE.zmax)    ) THEN
          !
          DO iFault = 1, MESH%Fault%nSide
             !
             iElem               = MESH%Fault%Face(iFault,1,1)
             iSide               = MESH%Fault%Face(iFault,2,1)
             IF (iElem == 0) CYCLE   ! put receiver only in '+' elements
             !
             IF (XYZInElement(io_x, io_y, io_z, iElem,tolerance,MESH,DISC)) THEN
                LocalRecPoint(i)%inside = .true. ! Point is located in an element
                LocalRecPoint(i)%globalreceiverindex = i
                LocalRecPoint(i)%X = io_x
                LocalRecPoint(i)%Y = io_y
                LocalRecPoint(i)%Z = io_z
                LocalRecPoint(i)%index  = iFault  ! number in Fault%Face list!
                xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
                yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
                zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
                CALL TrafoXYZ2XiEtaZeta(xi,eta,zeta,io_x,io_y,io_z,xV,yV,zV,MESH%LocalVrtxType(iElem))
                LocalRecPoint(i)%xi   = xi
                LocalRecPoint(i)%eta  = eta
                LocalRecPoint(i)%zeta = zeta
                in = in + 1
                EXIT
             ENDIF
          ENDDO ! iFace loop
          !
       ENDIF ! receiver within domain
       !
    ENDDO ! Loop over all fault output receivers (i)

  END SUBROUTINE


  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  !---------------------------------------------------------------------------!
  SUBROUTINE reallocate_and_fill_recpoint_array(DISC, LocalRecPoint, number_of_inner_receivers)

    IMPLICIT NONE

    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!

    TYPE(tUnstructPoint), ALLOCATABLE       :: LocalRecPoint(:)
    INTEGER                 :: number_of_inner_receivers

    !local variables
    INTEGER                 :: j, i

    DISC%DynRup%DynRup_out_atPickpoint%nDR_pick = number_of_inner_receivers
    ! Deallocate global array
    DEALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%RecPoint )
    ALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%RecPoint(number_of_inner_receivers))

    ! store only necessary information for this domain
    j = 0
    DO i = 1, DISC%DynRup%DynRup_out_atPickpoint%nOutPoints
      IF (LocalRecPoint(i)%inside) THEN
        j = j + 1
        DISC%DynRup%DynRup_out_atPickpoint%RecPoint(j) = LocalRecPoint(i)
      ENDIF
    ENDDO
    !
    DEALLOCATE(LocalRecPoint)

  END SUBROUTINE

  !-----------------------------------------------------------!
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!
  !-----------------------------------------------------------!
  SUBROUTINE ini_common_fault_receiver(DISC, MESH)
    IMPLICIT NONE

    ! Discretization struct.
    TYPE(tDiscretization)   :: DISC
    ! Mesh structure
    TYPE(tUnstructMesh)     :: MESH

    INTEGER                 :: number_of_inner_receivers

    ! temp array for storing inner receivers
    TYPE(tUnstructPoint), &
    ALLOCATABLE             :: LocalRecPoint(:)

    !---- implementation

    ! search for element which contains fault receiver
    ! and determine reference coordinates cleaning
    ! of double pickpoints due to MPI boundaries not necessary
    ! as there is just one '+' element
    CALL locate_receiver_points(DISC, MESH, LocalRecPoint, &
                                number_of_inner_receivers)

    ! Deallocate old RecPoint array, allocate it again
    ! and fills with collected LocalRecPoint,
    CALL reallocate_and_fill_recpoint_array(DISC, LocalRecPoint, &
                                            number_of_inner_receivers)

    ! assign DR_pick_output variable
    IF ( number_of_inner_receivers .EQ. 0 ) THEN
      ! If no pickpoint lies inside the domain,
      ! switch fault pickpoint output off.
      logInfo(*) 'No fault output receivers in this MPI domain '
      DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output = .FALSE.
    ELSE
      logInfo(*) 'Pick fault output at ', number_of_inner_receivers,' points in this MPI domain.'
      DISC%DynRup%DynRup_out_atPickpoint%DR_pick_output = .TRUE.

      ! procedure checks which variables we want to output and allocates
      ! temporal buffer arrays to collect the results of calculations
      CALL count_output_variables_allocate_arrays(DISC)
    ENDIF

  END SUBROUTINE ini_common_fault_receiver


  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  !---------------------------------------------------------------------!
  SUBROUTINE count_output_variables_allocate_arrays(DISC)

    IMPLICIT NONE

    TYPE(tDiscretization)   :: DISC                   ! Discretization struct.!

    INTEGER                 :: OutVars
    INTEGER                 :: number_of_receivers

    ! options of output components
    OutVars = 0
    IF (DISC%DynRup%DynRup_out_atPickpoint%OutputMask(1).EQ.1) OutVars = OutVars + 2
    IF (DISC%DynRup%DynRup_out_atPickpoint%OutputMask(2).EQ.1) OutVars = OutVars + 3
    IF (DISC%DynRup%DynRup_out_atPickpoint%OutputMask(3).EQ.1) OutVars = OutVars + 1
    IF (DISC%DynRup%DynRup_out_atPickpoint%OutputMask(4).EQ.1) OutVars = OutVars + 2
    IF (DISC%DynRup%DynRup_out_atPickpoint%OutputMask(5).EQ.1) OutVars = OutVars + 3

    number_of_receivers = DISC%DynRup%DynRup_out_atPickpoint%nDR_pick

    ! remember how many output variables do we have
    DISC%DynRup%DynRup_out_atPickpoint%nOutVars = OutVars

    ! Intermediate storage of output to reduce IO operations
    ALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%CurrentPick(number_of_receivers))

    DISC%DynRup%DynRup_out_atPickpoint%MaxPickStore = 50 ! every 50 levels
    ALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%TmpTime(DISC%DynRup%DynRup_out_atPickpoint%MaxPickStore))

    ALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%TmpState(DISC%DynRup%DynRup_out_atPickpoint%nDR_pick, &
                                                          DISC%DynRup%DynRup_out_atPickpoint%MaxPickStore, &
                                                          OutVars))

    ALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%OutVal(number_of_receivers, 1, OutVars))

    ! store rotation matrix for fault receiver
    ALLOCATE( DISC%DynRup%DynRup_out_atPickpoint%rotmat(number_of_receivers, 1:6, 1:6))

    DISC%DynRup%DynRup_out_atPickpoint%CurrentPick(:) = 0.
    DISC%DynRup%DynRup_out_atPickpoint%OutVal         = 0.
    DISC%DynRup%DynRup_out_atPickpoint%rotmat         = 0.

  END SUBROUTINE
  !---------------------------------------------------------------------------!

END MODULE common_fault_receiver_mod

