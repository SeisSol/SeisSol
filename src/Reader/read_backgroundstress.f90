!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
!!
!! @section LICENSE
!! Copyright (c) 2012, SeisSol Group
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
!! Routine reads initial background stress for dynamic rupture simulations

#ifdef BG
#include "../Initializer/preProcessorMacros.fpp"
#else
#include "Initializer/preProcessorMacros.fpp"
#endif

MODULE read_backgroundstress_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE COMMON_operators_mod
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  PUBLIC  :: extract_backgroundstress_2dplanefault
  PUBLIC  :: extract_backgroundstress_3dfault
  PUBLIC  :: read_scec_stress
  !---------------------------------------------------------------------------!  
  INTERFACE extract_backgroundstress_2dplanefault
     MODULE PROCEDURE extract_backgroundstress_2dplanefault
  END INTERFACE
  INTERFACE extract_backgroundstress_3dfault
     MODULE PROCEDURE extract_backgroundstress_3dfault
  END INTERFACE
  INTERFACE read_scec_stress
     MODULE PROCEDURE read_scec_stress
  END INTERFACE

CONTAINS

  SUBROUTINE extract_backgroundstress_2dplanefault(EQN,MESH,IO,DISC,BND)
    !< routine extracts and projects an initial backgroundstress on a 2d planar fault
    !-------------------------------------------------------------------------!
    USE DGBasis_mod
    USE TrilinearInterpolation_mod
    !-------------------------------------------------------------------------!    
    IMPLICIT NONE    
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE (tEquations)       :: EQN                                            !< Equations             !
    TYPE (tUnstructMesh)    :: MESH                                           !< Mesh structure        !
    TYPE (tInputOutput)     :: IO                                             !< IO structure          !
    TYPE (tDiscretization)  :: DISC                                           !< Discretization struct.!
    TYPE (tBoundary)        :: BND                                            !< Boundary structure    !
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                 :: i,j,k,counter
    INTEGER                 :: nx,nz
    INTEGER                 :: ibackgroundfield, nbackgroundfield
    INTEGER                 :: iElem,iSide
    INTEGER                 :: iNeighbor,iLocalNeighborSide
    INTEGER                 :: iObject,MPIIndex
    INTEGER                 :: GPwise                                         !< switch for different background options
    INTEGER                 :: iBndGP,nsamples_total
    INTEGER                 :: VertexSide(4,3)
    
    REAL                    :: x,y,z                                          !<
    REAL                    :: dx,dy,dz
    REAL                    :: xf_max,xf_min,yf_max,yf_min,zf_max,zf_min
    REAL                    :: posx_max,posx_min,posy_max,posy_min,posz_max,posz_min
    REAL                    :: chi,tau,xi,eta,zeta
    REAL                    :: xGP,yGP,zGP
    REAL                    :: xp(MESH%GlobalElemType), yp(MESH%GlobalElemType), zp(MESH%GlobalElemType)
    REAL                    :: xV(MESH%GlobalVrtxType),yV(MESH%GlobalVrtxType),zV(MESH%GlobalVrtxType)
    REAL, ALLOCATABLE       :: posx(:),posy(:),posz(:)
    REAL, ALLOCATABLE       :: P(:,:)
    REAL, ALLOCATABLE       :: bgfield(:,:,:,:) 
    
    !-------------------------------------------------------------------------!
    INTENT(IN)              :: MESH, BND
    INTENT(INOUT)           :: EQN, IO, DISC                                  !< global variables
    !-------------------------------------------------------------------------!
    
    ! Option
    ! GPwise = 1 : ini background stress variable inside an element
    ! GPwise = 0 : ini background stress constant per entire element
    GPwise = EQN%GPwise

    logInfo(*) 'Interpolating background stress for a 2D planar fault'

    ! fault is assumed to be xz plane
    xf_max = MAXVAL(DISC%DynRup%bg_stress%strike(:))
    xf_min = MINVAL(DISC%DynRup%bg_stress%strike(:))
    yf_max = 0.0D0
    yf_min = 0.0D0
    ! dip is positive at SCEC, negative in the mesh -> sign change
    zf_max = MAXVAL(-DISC%DynRup%bg_stress%dip(:))
    zf_min = MINVAL(-DISC%DynRup%bg_stress%dip(:))
    
    ! determine distance between two nodes of background field
    dx = DISC%DynRup%bg_stress%strike(2) - DISC%DynRup%bg_stress%strike(1)
    dy = dx ! assuming uniform grid!
    dz = dx ! assuming uniform grid!

    logInfo('("Background stess field x = [",E15.5," , ",E15.5,"]")')  xf_min, xf_max
    logInfo('("Background stess field y = [",E15.5," , ",E15.5,"]")')  yf_min, yf_max
    logInfo('("Background stess field z = [",E15.5," , ",E15.5,"]")')  zf_min, zf_max
    
    ! copy backgroundfield on tmp variable bgfield
    !
    ! note:
    ! variable 'fields' contains in this order
    ! normal_stress
    ! traction_strike
    ! traction_dip
    ! mu_S
    ! mu_D
    ! D_C
    ! cohesion
    ! forced_rupture_time
    !
    nbackgroundfield = SIZE(DISC%DynRup%bg_stress%fields,DIM=1)
    nx = DISC%DynRup%bg_stress%nx
    nz = DISC%DynRup%bg_stress%nz
    
    ALLOCATE(bgfield(nbackgroundfield,nx,3,nz)) ! create dummy array in z direction
    
    DO ibackgroundfield = 1,nbackgroundfield
      DO i = 1,nx
      DO j = 1,3 ! create dummy array in y direction
      DO k = 1,nz ! reordering array according to routine 'TrilinearFromRaster' necessary! happens in target argument j
        bgfield(ibackgroundfield,i,j,nz-k+1) = DISC%DynRup%bg_stress%fields(ibackgroundfield, i+(k-1)*nx )
      ENDDO
      ENDDO
      ENDDO
    ENDDO

    ! Extract local information from global background stress type
    IF (GPwise == 0) THEN

        nsamples_total = MESH%Fault%nSide
        ALLOCATE( posx(nsamples_total), posy(nsamples_total), posz(nsamples_total) )
        ALLOCATE(P(nbackgroundfield,nsamples_total))
        
        VertexSide(1,:) =  (/ 1, 3, 2 /)   ! Local tet. vertices of tet. side I   !
        VertexSide(2,:) =  (/ 1, 2, 4 /)   ! Local tet. vertices of tet. side II  !
        VertexSide(3,:) =  (/ 1, 4, 3 /)   ! Local tet. vertices of tet. side III !
        VertexSide(4,:) =  (/ 2, 3, 4 /)   ! Local tet. vertices of tet. side IV  !
        
        DO i = 1,nsamples_total
        
          ! get coordinates of elements fault face
          IF (MESH%Fault%Face(i,1,1) .NE. 0) THEN
              iElem               = MESH%Fault%Face(i,1,1)
              iSide               = MESH%Fault%Face(i,2,1)
              !
              DO j=1,3
                xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
                yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
                zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),iElem))
              ENDDO           
          ELSEIF (MESH%Fault%Face(i,1,1) == 0) THEN ! in case "+" element is not present in the local domain
              !
              iSide = MESH%Fault%Face(i,2,2)
              DO j=1,3
                xp(j) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(VertexSide(iSide,j),MESH%Fault%Face(i,1,2)))
                yp(j) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(VertexSide(iSide,j),MESH%Fault%Face(i,1,2)))
                zp(j) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(VertexSide(iSide,j),MESH%Fault%Face(i,1,2)))
              ENDDO            
          ENDIF
          
          ! compute barycenter of elements fault face
          posx(i) = sum(xp(1:3))/3.0D0
          posy(i) = sum(yp(1:3))/3.0D0 ! should be zero in the planar fault case
          posz(i) = sum(zp(1:3))/3.0D0
          
        ENDDO ! i = 1,MESH%Fault%nSide
 

    ELSEIF (GPwise == 1) THEN
    
        nsamples_total = MESH%Fault%nSide*DISC%Galerkin%nBndGP
        ALLOCATE( posx(nsamples_total), posy(nsamples_total), posz(nsamples_total) )
        ALLOCATE(P(nbackgroundfield,nsamples_total))    
    
        counter = 0
        DO i = 1,MESH%Fault%nSide
            ! get vertices of complete tet
            iSide = MESH%Fault%Face(i,2,1)
            IF (MESH%Fault%Face(i,1,1) == 0) THEN
                ! iElem is in the neighbor domain
                ! The neighbor element belongs to a different MPI domain
                iNeighbor           = MESH%Fault%Face(i,1,2)          ! iNeighbor denotes "-" side
                iLocalNeighborSide  = MESH%Fault%Face(i,2,2)
                iObject  = MESH%ELEM%BoundaryToObject(iLocalNeighborSide,iNeighbor)
                MPIIndex = MESH%ELEM%MPINumber(iLocalNeighborSide,iNeighbor)
                !
                xV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(1,1:4,MPIIndex)
                yV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(2,1:4,MPIIndex)
                zV(1:4) = BND%ObjMPI(iObject)%NeighborCoords(3,1:4,MPIIndex)
            ELSE
                !
                iElem               = MESH%Fault%Face(i,1,1)           
                !
                ! get vertices
                xV(1:4) = MESH%VRTX%xyNode(1,MESH%ELEM%Vertex(1:4,iElem))
                yV(1:4) = MESH%VRTX%xyNode(2,MESH%ELEM%Vertex(1:4,iElem))
                zV(1:4) = MESH%VRTX%xyNode(3,MESH%ELEM%Vertex(1:4,iElem))
            ENDIF
            !
            DO iBndGP = 1,DISC%Galerkin%nBndGP
                !
                counter = counter + 1
                !
                ! Transformation of boundary GP's into XYZ coordinate system
                chi  = MESH%ELEM%BndGP_Tri(1,iBndGP)
                tau  = MESH%ELEM%BndGP_Tri(2,iBndGP)
                CALL TrafoChiTau2XiEtaZeta(xi,eta,zeta,chi,tau,iSide,0)
                CALL TetraTrafoXiEtaZeta2XYZ(xGP,yGP,zGP,xi,eta,zeta,xV,yV,zV)
                !
                ! store location of each GP
                posx(counter) = xGP
                posy(counter) = yGP ! should be zero in the planar fault case
                posz(counter) = zGP
                !
            ENDDO ! iBndGP
        ENDDO ! i = 1,MESH%Fault%nSide 
    ENDIF ! GPwise

    
    ! check if background stress field covers complete fault:
    posx_max = MAXVAL(posx(1:nsamples_total))
    posx_min = MINVAL(posx(1:nsamples_total))
    posy_max = MAXVAL(posy(1:nsamples_total))
    posy_min = MINVAL(posy(1:nsamples_total))
    posz_max = MAXVAL(posz(1:nsamples_total))
    posz_min = MINVAL(posz(1:nsamples_total))
    !
    ! \todo include tolerance here
!~     IF( posx_max.GT.xf_max .OR. posz_max.GT.zf_max .OR.  &
!~         posx_min.LT.xf_min .OR. posz_min.LT.zf_min)THEN
!~         logError(*) 'Background stress field ',TRIM(IO%FileName_BackgroundStress),   &
!~                     ' does not fully include the entire fault!'
!~         STOP
!~     ENDIF

    ! Interpolation of the background field onto the barycenter of the elements fault face
    DO ibackgroundfield = 1,nbackgroundfield
       !
       ! Linear interpolation
       CALL TrilinearFromRaster(P(ibackgroundfield,:),                      &
                                posx,posy,posz,                             &
                                bgfield(ibackgroundfield,:,:,:),            &
                                dx,dy,dz,                                   &
                                xf_min,yf_min,zf_min)
       !
    ENDDO ! ibackgroundfield = 1:nbackgroundfield


    ! OUTPUT
    IF (GPwise == 0) THEN
        DO i = 1,MESH%Fault%nSide 
            EQN%IniBulk_xx(i,:) = 0.0D0
            EQN%IniBulk_yy(i,:) = -P(1,i) ! ATTENTION: Sign change as compression is negative in SeisSol3D
            EQN%IniBulk_zz(i,:) = 0.0D0
            EQN%IniShearXY(i,:) = P(2,i) ! ATTENTION XY is strike direction for SCEC
            EQN%IniShearYZ(i,:) = P(3,i) ! ATTENTION YZ is dip direction for SCEC
            EQN%IniShearXZ(i,:) = 0.0D0
            EQN%IniStateVar(i,:) = 0.0D0

            DISC%DynRup%Mu_S(:,i) = P(4,i)
            DISC%DynRup%Mu_D(:,i) = P(5,i)
            DISC%DynRup%D_C(:,i)= P(6,i)
            DISC%DynRup%cohesion(:,i)= -P(7,i) ! ATTENTION: Sign change as compression is negative in SeisSol3D
            DISC%DynRup%forced_rupture_time(:,i) = P(8,i)
        ENDDO ! i = 1,MESH%Fault%nSide
!            ! DEVELOPMENT INFORMATION
!            open(111,FILE='field_output.dat',FORM='FORMATTED')
!            do i = 1,MESH%Fault%nSide
!            write(111,*) posx(i),posz(i),EQN%IniShearXY(i,1)
!            enddo
!            close(111)
!            stop
        
    ELSEIF (GPwise == 1) THEN
        counter = 0
        DO i = 1,MESH%Fault%nSide
            DO iBndGP = 1,DISC%Galerkin%nBndGP
                counter = counter + 1        
                EQN%IniBulk_xx(i,iBndGP) = 0.0D0
                EQN%IniBulk_yy(i,iBndGP) = -P(1,counter) ! ATTENTION: Sign change as compression is negative in SeisSol3D
                EQN%IniBulk_zz(i,iBndGP) = 0.0D0
                EQN%IniShearXY(i,iBndGP) = P(2,counter) ! ATTENTION XY is strike direction for SCEC
                EQN%IniShearYZ(i,iBndGP) = P(3,counter) ! ATTENTION YZ is dip direction for SCEC
                EQN%IniShearXZ(i,iBndGP) = 0.0D0
                EQN%IniStateVar(i,iBndGP) = 0.0D0

                DISC%DynRup%Mu_S(iBndGP,i) = P(4,counter)
                DISC%DynRup%Mu_D(iBndGP,i) = P(5,counter)
                DISC%DynRup%D_C(iBndGP,i)= P(6,counter)
                DISC%DynRup%cohesion(iBndGP,i)= -P(7,counter) ! ATTENTION: Sign change as compression is negative in SeisSol3D
                DISC%DynRup%forced_rupture_time(iBndGP,i) = P(8,counter)
            ENDDO ! iBndGP
        ENDDO ! i = 1,MESH%Fault%nSide
!            ! DEVELOPMENT INFORMATION
!            open(111,FILE='field_output.dat',FORM='FORMATTED')
!            counter = 0
!            do i = 1,MESH%Fault%nSide
!            do j = 1,DISC%Galerkin%nBndGP
!            counter = counter + 1
!            write(111,*) posx(counter),posz(counter),EQN%IniShearXY(i,j)
!            enddo
!            enddo
!            close(111)
!            stop             
    ENDIF ! GPwise


    ! DEALLOCATE global background stress information
    DEALLOCATE(DISC%DynRup%bg_stress%strike,DISC%DynRup%bg_stress%dip,DISC%DynRup%bg_stress%fields)
    DEALLOCATE(posx,posy,posz)
    DEALLOCATE(P,bgfield)
    
  END SUBROUTINE extract_backgroundstress_2dplanefault
  
  !---------------------------------------------------------------------------!  

  SUBROUTINE extract_backgroundstress_3dfault()
    !< routine extracts and projects an initial backgroundstress on an arbitrary
    !< shaped three dimension fault
    !-------------------------------------------------------------------------!
    IMPLICIT NONE    
    !-------------------------------------------------------------------------!
    ! Argument list declaration 

    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !

    ! May be assigning the data to the elements in a matlab script like material values?


    
  END SUBROUTINE extract_backgroundstress_3dfault
 
  !---------------------------------------------------------------------------!
  
  SUBROUTINE read_scec_stress(DISC,IO)
    !< subroutine reads in SCEC specific background stress
    !-------------------------------------------------------------------------!
    IMPLICIT NONE    
    !-------------------------------------------------------------------------!
    TYPE(tDiscretization)   :: DISC                                           !< Discretization struct.!
    TYPE(tInputOutput)      :: IO                                             !< IO structure          !
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: i
    INTEGER                         :: lines                                  !< lines in the file
    INTEGER, TARGET                 :: nodes_nx, nodes_nz
    INTEGER                         :: intDummy
    REAL                            :: realDummy
    !-------------------------------------------------------------------------!
    INTENT(INOUT)                   :: IO                                                 
    INTENT(INOUT)                   :: DISC
    !-------------------------------------------------------------------------!
    
    logInfo(*) 'Rupture model read from ', TRIM(IO%FileName_BackgroundStress)
    
    CALL OpenFile(                                       &                        
         UnitNr       = IO%UNIT%other01                , &                        
         Name         = IO%FileName_BackgroundStress   , &
         create       = .FALSE.                          )    
        
    ! Read header
    READ(IO%UNIT%other01,*) intDummy, intDummy
    READ(IO%UNIT%other01,*) nodes_nx, nodes_nz, realDummy, realDummy
    READ(IO%UNIT%other01,*) intDummy, intDummy, realDummy, realDummy, realDummy, realDummy
    
    DISC%DynRup%bg_stress%nx = nodes_nx+1
    DISC%DynRup%bg_stress%nz = nodes_nz+1
    lines = DISC%DynRup%bg_stress%nx * DISC%DynRup%bg_stress%nz    
    
    ALLOCATE(                                            &
    DISC%DynRup%bg_stress%strike(lines),                 &
    DISC%DynRup%bg_stress%dip(lines),                    &
    DISC%DynRup%bg_stress%fields(8,lines)                &
    )
    
    ! variable fields contains in this order
    ! normal_stress (is read in as positiv from file)
    ! traction_strike
    ! traction_dip
    ! mu_S
    ! mu_D
    ! D_C
    ! cohesion (is read in as positiv from file)
    ! forced_rupture_time
    
    ! Read data (14 columns)
    DO i = 1,lines
        READ(IO%UNIT%other01,*) intDummy, intDummy, &
        DISC%DynRup%bg_stress%strike(i), DISC%DynRup%bg_stress%dip(i), &
        DISC%DynRup%bg_stress%fields(1,i), DISC%DynRup%bg_stress%fields(2,i), &
        DISC%DynRup%bg_stress%fields(3,i), realDummy, realDummy, &
        DISC%DynRup%bg_stress%fields(4,i), DISC%DynRup%bg_stress%fields(5,i), DISC%DynRup%bg_stress%fields(6,i), &
        DISC%DynRup%bg_stress%fields(7,i), DISC%DynRup%bg_stress%fields(8,i)
    ENDDO ! i lines
  
    CLOSE(IO%UNIT%other01)
    
    logInfo(*) 'Rupture model read in successfully! '
  
  END SUBROUTINE read_scec_stress
  
END MODULE read_backgroundstress_mod    
