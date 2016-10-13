!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Johannes Klicpera (johannes.klicpera AT tum.de)
!!
!! @section LICENSE
!! Copyright (c) 2015, SeisSol Group
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
!! Routine handling input of different fault parameters

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

module faultinput_mod
  !---------------------------------------------------------------------------!
  use TypesDef
  use DGBasis_mod
  use read_backgroundstress_mod
  !---------------------------------------------------------------------------!
  implicit none
  private

  ! Variables for throwing errors only once
  logical :: param_error_thrown = .false.
  logical :: dir_error_thrown = .false.
  logical :: shape_error_thrown = .false.
  !---------------------------------------------------------------------------!
  public :: faultinput
  private :: read_specfem3d
  private :: read_dist
  private :: process_dist
  private :: set_param
  private :: cnt_dist
  private :: toString

  ! Enumerators for specification of heterogeneous distribution blocks
  enum, bind(C)
    enumerator :: par_unknown=0, par_globalStress, par_localStress, par_cohesion, par_D_C, &
                    par_IniStateVar, par_RS_a, par_RS_srW, par_Mu_S, par_Mu_D, par_IniMu, par_Strength
  end enum

  enum, bind(C)
    enumerator :: dir_unknown=0, dir_normal=1, dir_along_strike=4, dir_along_dip=6
  end enum

  enum, bind(C)
    enumerator :: global_xx=1, global_yy, global_zz, global_xy, global_yz, global_xz
  end enum

  enum, bind(C)
    enumerator :: shape_unknown=0, shape_global, shape_square, shape_rectangle, shape_rectangle_taper, &
                    shape_circle, shape_ellipse, shape_xCylinder, shape_yCylinder, shape_zCylinder
  end enum

  ! Type for collecting all the information about heterogeneous distribution blocks,
  ! used in the input of fault-local stress fields
  type :: tDist
    private

    ! Basic information
    integer :: param                                                ! Parameter on which this distribution block acts
    integer :: dir													! For stress fields: Direction on which this distribution acts
    integer :: shape												! Shape of this heterogeneous distribution block

    ! Block parameters
    real :: val														! Standard stress field value
    real :: valh													! Second stress field value for the rectangle with taper
    real :: xc, yc, zc												! Center coordinates of the distribution block
    real :: l														! Side length of the square
    real :: lx, ly, lz												! Side lengths of the rectangle, distortion length of the ellipse, height of the cylinders
    real :: r														! Radius of the circle and ellipse

  end type tDist

  contains

  ! Subroutine for reading in parameters from the SPECFEM3D-style "Par_file_faults"-file
  subroutine faultinput (disc, eqn, mesh, bnd, IO)
  	!-------------------------------------------------------------------------!
	!use TypesDef
    use JacobiNormal_mod, only: RotationMatrix3D
	use COMMON_operators_mod
  	!-------------------------------------------------------------------------!
  	implicit none
	!-------------------------------------------------------------------------!
	type(tDiscretization), intent(inout), target  :: disc
	type(tEquations), intent(inout)               :: eqn
	type(tUnstructMesh), intent(in)               :: mesh
    type (tBoundary), intent(in)                  :: bnd
    type(tInputOutput), intent(in)                :: IO
  	!-------------------------------------------------------------------------!
  	! Local variable declaration

  	! index variables
  	integer :: iFault												! Iteration over fault faces
    integer :: iBndGP                                               ! Iteration over boundary Gaussian points
  	integer :: iDist												! Iteration over heterogeneous distributions
  	integer :: iElem												! Index of fault neighboring element
  	integer :: iSide												! Index of the fault neigboring side inside the neighboring element
  	integer :: i													! General purpose counter
  	integer :: iLocal					        					! Index of fault neighboring element on "-" side
  	integer :: iLocalSide   			          					! Index of the fault neigboring side inside the neighboring element on "-" side
    integer :: MPIIndex, iObject                                    ! Indices for MPI call for neighboring element in different MPI-domain

  	real	:: geoZ(3)												! base vector in z-direction (0,0,1)

  	! Temporary information about the fault face currently being processed
  	real :: geoStrike(3)											! Vector in along-strike direction
  	real :: geoDip(3)												! Vector in along-dip direction
  	real :: rotationMatrix(eqn%nVar,eqn%nVar)						! Rotation matrix from fault-local to global coordinates
  	real :: invMatrix(eqn%nVar,eqn%nVar)							! Rotation matrix from global to fault-local coordinates (unused)
  	real, dimension(mesh%GlobalVrtxType) :: xV, yV, zV              ! Coordinates of the tetragon's (or other body's) vertices
  	real :: chi, tau, xi, eta, zeta                                 ! Temporary variables for transformation of boundary Gaussian point coordinates into XYZ coordinate system
    real :: xGP, yGP, zGP                                           ! Coordinates of the currently processed boundary Gaussian point
    real :: localStress(6) 	                                        ! Fault-local stress field: 1=xx=fault-normal stress, 4=xy=along-strike shear, 6=xz=along-dip shear

  	! Variables for file input
    integer :: ios                                                  ! I/O-Status of input
    real    :: localStress_0(3)                                     ! Default values for the fault-local stress field
  	integer :: nSum													! Total number of heterogeneous distribution blocks
    integer :: nSpecfem                                             ! Number of distribution blocks needed for SPECFEM3D blocks

  	type(tDist), allocatable  :: dists(:)							! Heterogeneous distribution blocks
  	!-------------------------------------------------------------------------!

  	! set default values
  	geoZ(:) = (/0, 0, 1/)


  	open(IO%unit%FileIn_Fault, file="Par_file_faults", status='old', action="read", iostat=ios)
    if (ios /= 0) then
        logError(*) 'Error opening the local fault parameter file. File "Par_file_faults" present?'
    end if

    ! Get total number of heterogeneous distribution blocks
    nSum = cnt_dist(eqn, IO)
  	allocate(dists(nSum))

    ! Process namelists in SPECFEM3D-style
    call read_specfem3d(disc, eqn, IO, ios, dists, nSum, nSpecfem, localStress_0)

    ! Set parameters to default values
    eqn%IniStateVar(:,:) =  eqn%RS_sv0
    DISC%DynRup%cohesion(:,:) = DISC%DynRup%cohesion_0

    ! Read in other distribution blocks (parameter type and direction from namelist, SeisSol-style)
  	do i = nSpecfem + 1, nSum
  		call read_dist(eqn, IO, dists(i), ios)
        if (ios /= 0) then
            logError(*) 'Error reading in a dist2d namelist.'
        end if
  	end do

  	close(IO%unit%FileIn_Fault)


  	! Process stress input
  	FaultFacesIteration: do iFault = 1, mesh%Fault%nSide

  		! switch for rupture front output: RF
  		if (disc%DynRup%RFtime_on == 1) then
  		  ! rupture front output just for + side elements!
  		  if (mesh%fault%Face(iFault,1,1) .NE. 0) disc%DynRup%RF(iFault,:) = .TRUE.
  		end if

        ! Get index of neighbouring element and the neighbour side within that element
        iElem = mesh%Fault%Face(iFault,1,1)
        iSide = mesh%Fault%Face(iFault,2,1)
        if (iElem == 0) then ! iElem is in the neighbor domain (different MPI domain)
            iLocal = mesh%Fault%Face(iFault,1,2)   ! iLocal denotes local "-" side
            iLocalSide = mesh%Fault%Face(iFault,2,2)
            iObject  = mesh%elem%BoundaryToObject(iLocalSide,iLocal)
            MPIIndex = mesh%elem%MPINumber(iLocalSide,iLocal)

            ! Get XYZ coordinates of the element
            xV(:) = bnd%ObjMPI(iObject)%NeighborCoords(1,:,MPIIndex)
            yV(:) = bnd%ObjMPI(iObject)%NeighborCoords(2,:,MPIIndex)
            zV(:) = bnd%ObjMPI(iObject)%NeighborCoords(3,:,MPIIndex)
        else
            xV(:) = mesh%vrtx%xyNode(1,mesh%elem%Vertex(:,iElem))
            yV(:) = mesh%vrtx%xyNode(2,mesh%elem%Vertex(:,iElem))
            zV(:) = mesh%vrtx%xyNode(3,mesh%elem%Vertex(:,iElem))
        end if

        ! Get base vectors of the fault-local coordinate system
        geoStrike = mesh%Fault%geoNormals(:,iFault) .x. geoZ								! cross-product to get along-strike vector
        geoStrike = geoStrike / sqrt(geoStrike(1)**2 + geoStrike(2)**2)	                    ! normalize along-strike vector
        geoDip = mesh%Fault%geoNormals(:,iFault) .x. geoStrike								! cross-product to get along-dip vector

        ! Get rotation matrix from fault-local to global coordinate system
        call RotationMatrix3D( n1  = mesh%fault%geoNormals(:,iFault), &
                               n2  = geoStrike, &
                               n3  = geoDip, &
                               T   = rotationMatrix, &
                               iT  = invMatrix, &
                               eqn = eqn )

        ! Iterate over Gaussian points
        GaussPointIteration: do iBndGP = 1, disc%Galerkin%nBndGP

            ! constant stress tensor in fault-local coordinates
            localStress(1) = localStress_0(3) ! fault-normal
            localStress(2) = 0.0
            localStress(3) = 0.0
            localStress(4) = localStress_0(1) ! along-strike
            localStress(5) = 0.0
            localStress(6) = localStress_0(2) ! along-dip

            ! Transformation of boundary Gaussian point's coordinates into XYZ coordinate system
            chi  = mesh%elem%BndGP_Tri(1,iBndGP)
            tau  = mesh%elem%BndGP_Tri(2,iBndGP)
            call TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
            call TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)

            ! process the heterogeneous distribution blocks
            DistIteration: do iDist = 1, nSum
                call process_dist(disc, eqn, localStress, dists(iDist), iFault, iBndGP, xGP, yGP, zGP)
            end do DistIteration

            !--------- Convert fault-local stress field to global coordinates ---------!

            ! Transform local stress field to global coordinates
            localStress = matmul(rotationMatrix(:6,:6), localStress)

            ! Add local stress field to global stress field
            eqn%IniBulk_xx(iFault,iBndGP)  =  eqn%IniBulk_xx(iFault,iBndGP) + localStress(1)
            eqn%IniBulk_yy(iFault,iBndGP)  =  eqn%IniBulk_yy(iFault,iBndGP) + localStress(2)
            eqn%IniBulk_zz(iFault,iBndGP)  =  eqn%IniBulk_zz(iFault,iBndGP) + localStress(3)
            eqn%IniShearXY(iFault,iBndGP)  =  eqn%IniShearXY(iFault,iBndGP) + localStress(4)
            eqn%IniShearYZ(iFault,iBndGP)  =  eqn%IniShearYZ(iFault,iBndGP) + localStress(5)
            eqn%IniShearXZ(iFault,iBndGP)  =  eqn%IniShearXZ(iFault,iBndGP) + localStress(6)

            !--------------------------------------------------------------------------!

            end do GaussPointIteration

  	end do FaultFacesIteration

  end subroutine faultinput

  ! Subroutine for reading in namelists in SPECFEM3D format
  subroutine read_specfem3d(disc, eqn, IO, ios, dists, nSum, nSpecfem, localStress_0)
  	!-------------------------------------------------------------------------!
  	implicit none
  	!-------------------------------------------------------------------------!
	type(tDiscretization), intent(inout), target :: disc
	type(tEquations), intent(inout)                 :: eqn
    type(tInputOutput)                           :: IO
    integer, intent(out)                         :: ios             ! I/O-Status of input
  	type(tDist), allocatable, intent(inout)      :: dists(:)        ! Heterogeneous distribution blocks
    integer, intent(in)                          :: nSum            ! Total number of heterogeneous distribution blocks
    integer, intent(out)                         :: nSpecfem        ! Number of distribution blocks needed for SPECFEM3D blocks
    real, intent(out)                            :: localStress_0(3)! Constant stress field in fault-local coordinates
  	!-------------------------------------------------------------------------!
  	integer :: i													! General purpose counter
    integer :: cnt_dist2d                                           ! Number of dist2d blocks that have already been read in

  	real        :: Sigma(6)											! File input for constant stress fields in the global coordinate system
  	real        :: S1, S2, S3										! File input for constant fault-local stress fields: S1=along-strike shear, S2=along-dip shear, S3=fault-normal stress
  	integer     :: n1, n2, n3									    ! Number of heterogeneous distribution blocks for n1=along-strike shear, n2=along-dip shear, n3=fault-normal stress
  	real        :: mus, mud, dc										! File input for constant SWF: mus=static friction coefficient, mud=dynamic friction coefficient, dc=critical slip-weakening distance
  	integer     :: nmus, nmud, ndc									! Number of heterogeneous distribution blocks for SWF parameters

  	! Namelists
  	namelist /stress_tensor/ Sigma                                  ! Constant stress field in global coordinates
  	namelist /init_stress/ S1, S2, S3, n1, n2, n3                   ! Stress field in fault-local coordinates
    namelist /SWF/ mus, mud, dc, nmus, nmud, ndc                    ! Parameters for slip-weakening friction
  	!-------------------------------------------------------------------------!

    ! Set default values
  	Sigma(global_xx) = eqn%Bulk_xx_0
  	Sigma(global_yy) = eqn%Bulk_yy_0
  	Sigma(global_zz) = eqn%Bulk_zz_0
  	Sigma(global_xy) = eqn%ShearXY_0
  	Sigma(global_yz) = eqn%ShearYZ_0
  	Sigma(global_xz) = eqn%ShearXZ_0
  	S1 = 0
  	S2 = 0
  	S3 = 0
  	n1 = 0
  	n2 = 0
  	n3 = 0
    mus = disc%DynRup%Mu_S_ini
    mud = disc%DynRup%Mu_D_ini
    dc = disc%DynRup%D_C_ini
    nmus = 0
    nmud = 0
    ndc = 0


    ! Read in stress_tensor, init_stress and SWF namelists (each may not be present)
  	read(IO%unit%FileIn_Fault, nml=stress_tensor, iostat=ios)
    if (ios > 0) then
        logError(*) 'Error reading in the stress_tensor namelist.'
    end if
    rewind(IO%unit%FileIn_Fault)
    ! Save constant values
    eqn%IniBulk_xx(:,:)  =  Sigma(global_xx)
    eqn%IniBulk_yy(:,:)  =  Sigma(global_yy)
    eqn%IniBulk_zz(:,:)  =  Sigma(global_zz)
    eqn%IniShearXY(:,:)  =  Sigma(global_xy)
    eqn%IniShearYZ(:,:)  =  Sigma(global_yz)
    eqn%IniShearXZ(:,:)  =  Sigma(global_xz)

  	read(IO%unit%FileIn_Fault, nml=init_stress, iostat=ios)
    if (ios > 0) then
        logError(*) 'Error reading in the init_stress namelist.'
    end if
    rewind(IO%unit%FileIn_Fault)
    localStress_0(1) = S1 ! along-strike
    localStress_0(2) = S2 ! along-dip
    localStress_0(3) = S3 ! fault-normal

  	read(IO%unit%FileIn_Fault, nml=SWF, iostat=ios)
    if (ios > 0) then
        logError(*) 'Error reading in the SWF namelist.'
    end if
    if (ios == 0 .and. eqn%FL /= 2) then
        logError(*) 'Slip weakening parameters set in Par_file_faults, while friction type is not set to linear slip weakening.'
    end if
    rewind(IO%unit%FileIn_Fault)
    ! Save constant values
    disc%DynRup%Mu_S(:,:) = mus
    disc%DynRup%Mu_D(:,:) = mud
    disc%DynRup%D_C(:,:) = dc


    ! Check if number of heterogeneous distribution blocks is sufficient
    nSpecfem = n1 + n2 + n3 + nmus + nmud + ndc
    if (nSum < nSpecfem) then
        logError(*) 'Error reading in the dist2d namelists. Too few heterogeneous distribution blocks. ', trim(adjustl(toString(nSpecfem))), ' at least expected, only ', trim(adjustl(toString(nSum))), ' present!'
    end if



  	! Read in the heterogeneous distribution blocks specified in init_stress
  	! 1. along-strike
    cnt_dist2d = 0
  	do i = cnt_dist2d + 1, cnt_dist2d + n1
  		call read_dist(eqn, IO, dists(i), ios)
        if (ios /= 0) then
            logError(*) 'Error reading in a dist2d namelist.'
        end if
        dists(i)%param = par_localStress
  		dists(i)%dir = dir_along_strike
  	end do
    cnt_dist2d = cnt_dist2d + n1
  	! 2. along-dip
  	do i = cnt_dist2d + 1, cnt_dist2d + n2
  		call read_dist(eqn, IO, dists(i), ios)
        if (ios /= 0) then
            logError(*) 'Error reading in a dist2d namelist.'
        end if
        dists(i)%param = par_localStress
  		dists(i)%dir = dir_along_dip
  	end do
    cnt_dist2d = cnt_dist2d + n2
  	! 3. fault-normal
  	do i = cnt_dist2d + 1, cnt_dist2d + n3
  		call read_dist(eqn, IO, dists(i), ios)
        if (ios /= 0) then
            logError(*) 'Error reading in a dist2d namelist.'
        end if
        dists(i)%param = par_localStress
  		dists(i)%dir = dir_normal
  	end do
    cnt_dist2d = cnt_dist2d + n3

  	! Read in the heterogeneous distribution blocks specified in SWF
  	! 1. static friction coefficient
  	do i = cnt_dist2d + 1, cnt_dist2d + nmus
  		call read_dist(eqn, IO, dists(i), ios)
        if (ios /= 0) then
            logError(*) 'Error reading in a dist2d namelist.'
        end if
        dists(i)%param = par_Mu_S
  	end do
    cnt_dist2d = cnt_dist2d + nmus
  	! 2. dynamic friction coefficient
  	do i = cnt_dist2d + 1, cnt_dist2d + nmud
  		call read_dist(eqn, IO, dists(i), ios)
        if (ios /= 0) then
            logError(*) 'Error reading in a dist2d namelist.'
        end if
        dists(i)%param = par_Mu_D
  	end do
    cnt_dist2d = cnt_dist2d + nmud
  	! 3. critical slip-weakening distance
  	do i = cnt_dist2d + 1, cnt_dist2d + ndc
  		call read_dist(eqn, IO, dists(i), ios)
        if (ios /= 0) then
            logError(*) 'Error reading in a dist2d namelist.'
        end if
        dists(i)%param = par_D_C
  	end do
    cnt_dist2d = cnt_dist2d + ndc

  end subroutine read_specfem3d


  ! Subroutine for checking if the point is in the distribution block and setting the parameter
  subroutine process_dist(disc, eqn, localStress, dist, iFault, iBndGP, x, y, z)
  	!-------------------------------------------------------------------------!
  	implicit none
	!-------------------------------------------------------------------------!
	type(tDiscretization), intent(inout), target  :: disc
	type(tEquations), intent(inout)               :: eqn

    real, intent(out)                             :: localStress(6)  ! Local stress field at the current boundary Gaussian point
    type(tdist), intent(inout)                    :: dist            ! Current distribution block
    integer, intent(in)                           :: iFault, iBndGP  ! Iteration number over fault faces and boundary Gaussian points
    real, intent(in)                              :: x, y, z         ! Location of the current boundary Gaussian point
  	!-------------------------------------------------------------------------!

    select case (dist%shape)
        case (shape_global)
            call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
        case (shape_square)	! (actually a cube)
            if ( abs(x - dist%xc) < dist%l / 2.0 .and. &		! x in bounds?
                    abs(y - dist%yc) < dist%l / 2.0 .and. &	! y in bounds?
                    abs(z - dist%zc) < dist%l / 2.0 ) then	! z in bounds?

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
            end if
        case (shape_rectangle) ! (actually a cuboid)
            if ( abs(x - dist%xc) < dist%lx / 2.0 .and. &		! x in bounds?
                    abs(y - dist%yc) < dist%ly / 2.0 .and. &	! y in bounds?
                    abs(z - dist%zc) < dist%lz / 2.0 ) then	! z in bounds?

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
            end if
        case (shape_rectangle_taper) ! (actually a cuboid)
            if ( abs(x - dist%xc) < dist%lx / 2.0 .and. &		! x in bounds?
                    abs(y - dist%yc) < dist%ly / 2.0 .and. &	! y in bounds?
                    abs(z - dist%zc) < dist%lz / 2.0 ) then	! z in bounds?

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, &
                    ! linear interpolation in z-direction between val and valh:
                    (dist%valh + dist%val) * 0.5 &                      ! add central value
                    + (z - dist%zc) * (dist%valh - dist%val) / dist%lz)	! add interpolation value
            end if
        case (shape_circle)	! (actually a sphere)
            if ( (x - dist%xc)**2 + (y - dist%yc)**2 + (z - dist%zc)**2 & ! distance from circle center
                    < dist%r**2 ) then

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
            end if
        case (shape_ellipse)
            if ( ((x - dist%xc)/dist%lx)**2 &
                    + ((y - dist%yc)/dist%ly)**2 &
                    + ((z - dist%zc)/dist%lz)**2 &
                    < 1 ) then

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
            end if
        case (shape_xCylinder)
            if ( (y - dist%yc)**2 + (z - dist%zc)**2 < dist%r**2 & ! distance from circle center in yz-plane
                    .and. abs(x - dist%xc) < dist%lz / 2.0) then					 ! x in bounds?

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
            end if
        case (shape_yCylinder)
            if ( (x - dist%xc)**2 + (z - dist%zc)**2 < dist%r**2 & ! distance from circle center in xz-plane
                    .and. abs(y - dist%yc) < dist%lz / 2.0) then					 ! y in bounds?

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
            end if
        case (shape_zCylinder)
            if ( (x - dist%xc)**2 + (y - dist%yc)**2 < dist%r**2 & ! distance from circle center in xy-plane
                    .and. abs(z - dist%zc) < dist%lz / 2.0) then					 ! z in bounds?

                call set_param(disc, eqn, localStress, dist, iFault, iBndGP, dist%val)
            end if
        case default ! unknown shape
            if ( .not. shape_error_thrown ) then
                logError(*) 'Unknown shape in heterogeneous distribution block!'
                shape_error_thrown = .true.
            end if
    end select

  end subroutine process_dist


  ! Subroutine for choosing and setting a parameter to the specified value
  subroutine set_param(disc, eqn, localStress, dist, iFault, iBndGP, val)
  	!-------------------------------------------------------------------------!
  	implicit none
	!-------------------------------------------------------------------------!
	type(tDiscretization), intent(inout), target  :: disc
	type(tEquations), intent(inout)               :: eqn

    real, intent(out)                             :: localStress(6)  ! Local stress field at the current boundary Gaussian point
    type(tdist), intent(inout)                    :: dist            ! Current distribution block
    integer, intent(in)                           :: iFault, iBndGP  ! Iteration number over fault faces and boundary Gaussian points
    real, intent(in)                              :: val             ! Value to set the parameter to
  	!-------------------------------------------------------------------------!

    select case (dist%param)
        case (par_localStress)
            localStress(dist%dir) = val
        case (par_globalStress)
            select case (dist%dir)
                case (global_xx)
                    eqn%IniBulk_xx(iFault,iBndGP) = val
                case (global_yy)
                    eqn%IniBulk_yy(iFault,iBndGP) = val
                case (global_zz)
                    eqn%IniBulk_zz(iFault,iBndGP) = val
                case (global_xy)
                    eqn%IniShearXY(iFault,iBndGP) = val
                case (global_yz)
                    eqn%IniShearYZ(iFault,iBndGP) = val
                case (global_xz)
                    eqn%IniShearXZ(iFault,iBndGP) = val
                case default ! unknown direction
                    if ( .not. dir_error_thrown ) then
                        logError(*) 'Unknown direction in heterogeneous distribution block for the stress field in global coordinates!'
                        dir_error_thrown = .true.
                    end if
            end select
        case (par_cohesion)
            disc%DynRup%cohesion(iFault,iBndGP) = val
        case (par_D_C)
            disc%DynRup%D_C(iFault,iBndGP) = val
        case (par_IniStateVar)
            eqn%IniStateVar(iFault,iBndGP) = val
        case (par_RS_a)
            disc%DynRup%RS_a_array(iFault,iBndGP) = val
        case (par_RS_srW)
            disc%DynRup%RS_srW_array(iFault,iBndGP) = val
        case (par_Mu_S)
            disc%DynRup%Mu_S(iFault,iBndGP) = val
        case (par_Mu_D)
            disc%DynRup%Mu_D(iFault,iBndGP) = val
        case (par_IniMu)
            eqn%IniMu(iFault,iBndGP) = val
        case (par_Strength)
            disc%DynRup%Strength(iFault,iBndGP) = val
        case default ! unknown parameter
            if ( .not. param_error_thrown ) then
                logError(*) 'Unknown parameter in heterogeneous distribution block!'
                param_error_thrown = .true.
            end if
    end select

  end subroutine set_param


  ! Function for counting the number of heterogeneous distribution blocks
  integer function cnt_dist(eqn, IO)
  	!-------------------------------------------------------------------------!
  	implicit none
  	!-------------------------------------------------------------------------!
	type(tEquations), intent(in)    :: eqn
    type(tInputOutput)              :: IO
  	!-------------------------------------------------------------------------!
  	! Local variable declaration
    integer :: ios                                                  ! I/O-Status of input
    integer :: cntDist                                              ! Number of dist blocks
    type(tDist) :: tmpDist                                          ! Temporary distribution block for reading in distribution blocks
  	!-------------------------------------------------------------------------!

    ios = 0
    cntDist = 0

    rewind(IO%unit%FileIn_Fault)

    do while (ios == 0)
        call read_dist(eqn, IO, tmpDist, ios)
        if (ios == 0) then
            cntDist = cntDist + 1
        end if
    end do

    rewind(IO%unit%FileIn_Fault)

    cnt_dist = cntDist

  end function cnt_dist


  ! Subroutine for reading in heterogeneous distribution blocks
  subroutine read_dist(eqn, IO, dist, ios)
  	!-------------------------------------------------------------------------!
  	implicit none
  	!-------------------------------------------------------------------------!
	type(tEquations), intent(in)    :: eqn
    type(tInputOutput)              :: IO
  	type(tDist), intent(inout)      :: dist                         ! Heterogeneous distribution block
    integer, intent(out)            :: ios                          ! I/O-Status of input
  	!-------------------------------------------------------------------------!
  	! Local variable declaration
    character(len=20) :: param                                      ! Parameter specification string
    character(len=20) :: dir                                        ! Direction specification string
  	character(len=20) :: shapeval									! Shape specification string
  	real :: val														! Standard stress field value
  	real :: valh													! Second stress field value for the rectangle with taper
  	real :: xc, yc, zc												! Center coordinates of the distribution block
  	real :: l														! Side length of the square
  	real :: lx, ly, lz												! Side lengths of the rectangle, distortion length of the ellipse, height of the cylinders
  	real :: r														! Radius of the circle and ellipse

  	! Namelist
  	namelist /dist2d/ param, dir, shapeval, val, valh, xc, yc, zc, r, l, lx, ly, lz
  	!-------------------------------------------------------------------------!

    ! set default values
    param = ''
    dir = ''
    shapeval = ''
    val = 0.0
    valh = 0.0
    xc = 0.0
    yc = 0.0
    zc = 0.0
    l = 0.0
    lx = 0.0
    ly = 0.0
    lz = 0.0
    r = 0.0

    ! read in the heterogeneous distribution block
  	read(IO%UNIT%FileIn_Fault, nml=dist2d, iostat=ios)
    if ( ios /= 0 ) then
        return
    end if

  	! Copy read-in values into distribution struct
  	dist%val = val
  	dist%valh = valh
  	dist%xc = xc
  	dist%yc = yc
  	dist%zc = zc
  	dist%l = l
  	dist%lx = lx
  	dist%ly = ly
  	dist%lz = lz
  	dist%r = r

    ! Select parameter
    select case (param)
        case ('globalstress', 'globalStress')
            dist%param = par_globalStress
        case ('localstress', 'localStress')
            dist%param = par_localStress
        case ('cohesion')
            dist%param = par_cohesion
        case ('d_c', 'D_C')
            dist%param = par_D_C
        case ('inistatevar', 'IniStateVar')
            if ( eqn%FL == 3 ) then
                dist%param = par_IniStateVar
            else
                logError(*) 'IniStateVar specified for a friction type other than rate-and-state friction.'
            end if
        case ('rs_a', 'RS_a')
            if ( eqn%FL == 3 ) then
                dist%param = par_RS_a
            else
                logError(*) 'RS_a specified for a friction type other than rate-and-state friction.'
            end if
        case ('rs_srw', 'RS_srW')
            if ( eqn%FL == 3 ) then
                dist%param = par_RS_srW
            else
                logError(*) 'RS_srW specified for a friction type other than rate-and-state friction.'
            end if
        case ('mu_s', 'Mu_S')
            if ( eqn%FL == 2 ) then
                dist%param = par_Mu_S
            else
                logError(*) 'Mu_S specified for a friction type other than linear slip weakening.'
            end if
        case ('mu_d', 'Mu_D')
            if ( eqn%FL == 2 ) then
                dist%param = par_Mu_D
            else
                logError(*) 'Mu_D specified for a friction type other than linear slip weakening.'
            end if
        case ('inimu', 'IniMu')
            if ( eqn%FL == 3 ) then
                dist%param = par_IniMu
            else
                logError(*) 'IniMu specified for a friction type other than rate-and-state friction.'
            end if
        case ('strength', 'Strength')
            dist%param = par_Strength
        case ('')
            dist%param = par_unknown
  		case default
  			dist%param = par_unknown
            logError(*) 'Parameter type "', trim(adjustl(param)), '" could not be recognized!'
    end select

    ! Select direction (if stress field input)
    if (dist%param == par_globalStress .or. dist%param == par_localStress) then
        select case (dir)
            case ('normal')
                dist%dir = dir_normal
            case ('strike', 'along-strike')
                dist%dir = dir_along_strike
            case ('dip', 'along-dip')
                dist%dir = dir_along_dip
            case ('xx', 'XX')
                dist%dir = global_xx
            case ('yy', 'YY')
                dist%dir = global_yy
            case ('zz', 'ZZ')
                dist%dir = global_zz
            case ('xy', 'XY')
                dist%dir = global_xy
            case ('yz', 'YZ')
                dist%dir = global_yz
            case ('xz', 'XZ')
                dist%dir = global_xz
            case default
                dist%dir = dir_unknown
                logError(*) 'Direction "', trim(adjustl(dir)), '" could not be recognized!'
        end select
    end if

  	! Select shape
  	select case (shapeval)
  		case ('global')
  			dist%shape = shape_global
  		case ('square')
  			dist%shape = shape_square
  		case ('rectangle')
  			dist%shape = shape_rectangle
  		case ('rectangle-taper')
  			dist%shape = shape_rectangle_taper
  		case ('circle')
  			dist%shape = shape_circle
  		case ('ellipse')
  			dist%shape = shape_ellipse
  		case ('x-cylinder')
  			dist%shape = shape_xCylinder
  		case ('y-cylinder')
  			dist%shape = shape_yCylinder
  		case ('z-cylinder')
  			dist%shape = shape_zCylinder
  		case default
  			dist%shape = shape_unknown
            logError(*) 'Shape type "', trim(adjustl(shapeval)), '" could not be recognized!'
  	end select
  end subroutine read_dist

  ! Function for converting an integer to a string (making it possible to trim whitespace)
  character(len=20) function toString(val)
  	!-------------------------------------------------------------------------!
  	implicit none
  	!-------------------------------------------------------------------------!
	integer, intent(in)    :: val                         ! Value of integer to convert
  	!-------------------------------------------------------------------------!
  	! Local variable declaration
    write(toString, *) val
  end function toString

end module faultinput_mod
