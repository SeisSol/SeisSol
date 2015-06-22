!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) 2010-2014, SeisSol Group
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

MODULE allocate_mesh_mod
  !--------------------------------------------------------------------------
  USE TypesDef
  !--------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE allocate_mesh_level0_1
     MODULE PROCEDURE allocate_mesh_level0_1
  END INTERFACE

  !INTERFACE allocate_ncmesh_level0_1
  !   MODULE PROCEDURE allocate_ncmesh_level0_1
  !END INTERFACE

  INTERFACE allocate_mesh_level0_2
     MODULE PROCEDURE allocate_mesh_level0_2
  END INTERFACE

  INTERFACE destruct_mesh_level0_1
     MODULE PROCEDURE destruct_mesh_level0_1
  END INTERFACE

  !INTERFACE destruct_ncmesh_level0_1
  !   MODULE PROCEDURE destruct_ncmesh_level0_1
  !END INTERFACE

  INTERFACE destruct_mesh_level0_2
     MODULE PROCEDURE destruct_mesh_level0_2
  END INTERFACE

  !----------------------------------------------------------------------------
  PUBLIC  :: allocate_mesh_level0_1
  !PUBLIC  :: allocate_ncmesh_level0_1
  PUBLIC  :: allocate_mesh_level0_2
  PUBLIC  :: destruct_mesh_level0_1
  !PUBLIC  :: destruct_ncmesh_level0_1
  PUBLIC  :: destruct_mesh_level0_2
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE allocate_mesh_level0_1(IO,MESH)
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH               
    TYPE (tInputOutput)       :: IO
    INTEGER                   :: allocStat
    !--------------------------------------------------------------------------
    INTENT(IN)                :: IO
    INTENT(INOUT)             :: MESH
    !--------------------------------------------------------------------------
    !
    logInfo0(*) 'Memory allocation for the mesh ...'
    !    
    ALLOCATE( MESH%ELEM%Vertex(           MESH%nVertexMax,   MESH%nElem), &   
              MESH%ELEM%reference(        0:MESH%nVertexMax, MESH%nElem), &  !should be nSideMax
              MESH%ELEM%MPIreference(     0:MESH%nVertexMax, MESH%nElem), &  
              MESH%ELEM%BoundaryToObject( MESH%nVertexMax,   MESH%nElem), &
              MESH%ELEM%SideNeighbor(     MESH%nVertexMax,   MESH%nElem), & 
              MESH%ELEM%LocalNeighborSide(MESH%nVertexMax,   MESH%nElem), & 
              MESH%ELEM%LocalNeighborVrtx(MESH%nVertexMax,   MESH%nElem), & 
              MESH%ELEM%xyBary(           MESH%Dimension,    MESH%nElem), & 
              MESH%ELEM%Volume(           MESH%nElem                   ), &
              MESH%ELEM%MinDistBarySide(  MESH%nElem                   ), &
              STAT=allocStat                                              )
    
    IF (allocStat .NE. 0) THEN
        logError(*) 'allocate_mesh_level0_1: could not allocate the whole mesh!'
        STOP
    END IF
    !
  END SUBROUTINE allocate_mesh_level0_1

  SUBROUTINE destruct_mesh_level0_1(MESH)
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH               
    !--------------------------------------------------------------------------
    !
    DEALLOCATE(MESH%ELEM%Vertex               , &
               MESH%ELEM%reference            , &
               MESH%ELEM%MPIreference         , &
               MESH%ELEM%BoundaryToObject     , &               
               MESH%ELEM%SideNeighbor         , &
               MESH%ELEM%LocalNeighborSide    , &
               MESH%ELEM%LocalNeighborVrtx    , & 
               MESH%ELEM%xyBary               , &
               MESH%ELEM%Volume               , &
               MESH%ELEM%MinDistBarySide         )
    !
  END SUBROUTINE destruct_mesh_level0_1


  !SUBROUTINE allocate_ncmesh_level0_1(IO,MESH,DISC)
  !  !--------------------------------------------------------------------------
  !  TYPE (tUnstructMesh)      :: MESH  
  !  TYPE (tDiscretization)    :: DISC             
  !  TYPE (tInputOutput)       :: IO
  !  INTEGER                   :: allocStat
  !  !--------------------------------------------------------------------------
  !  INTENT(IN)                :: IO, DISC
  !  INTENT(INOUT)             :: MESH
  !  !--------------------------------------------------------------------------
  !  !
  !  WRITE(IO%UNIT%stdOut,*) '|   Memory allocation for the mesh ... '
  !  !    
  !  ALLOCATE( MESH%ELEM%ncBndNeighbor(MESH%nNonConformingEdges,(DISC%Galerkin%nPoly + 2)**2) , &  
  !            MESH%ELEM%NCB_IndexList(MESH%nNonConformingEdges,2)                            , &  
  !            MESH%ELEM%ncBndGaussP(MESH%nNonConformingEdges,(DISC%Galerkin%nPoly + 2)**2,3) , & 
  !            MESH%ELEM%NC_BoundaryToObject(MESH%nNonConformingEdges,DISC%Galerkin%nBndGP)   , &
  !            STAT=allocStat                                                                   )
  !  
  !  IF (allocStat .NE. 0) THEN
  !      WRITE(IO%UNIT%errOut,*) 'ERROR allocate_mesh_level0_1: could not allocate the whole mesh!'
  !      STOP
  !  END IF
  !  !
  !END SUBROUTINE allocate_ncmesh_level0_1


  !SUBROUTINE destruct_ncmesh_level0_1(MESH)
  !  !--------------------------------------------------------------------------
  !  TYPE (tUnstructMesh)      :: MESH               
  !  !--------------------------------------------------------------------------
  !  !
  !  DEALLOCATE(MESH%ELEM%ncBndNeighbor     , &  
  !             MESH%ELEM%NCB_IndexList     , &  
  !             MESH%ELEM%ncBndGaussP       , & 
  !             MESH%ELEM%ncIndex           , &
  !             MESH%ELEM%NC_BoundaryToObject )
  !  !
  !END SUBROUTINE destruct_ncmesh_level0_1


  SUBROUTINE allocate_mesh_level0_2(IO,EQN,MESH)
    !--------------------------------------------------------------------------
    TYPE (tInputOutput)       :: IO
    TYPE (tEquations)         :: EQN
    TYPE (tUnstructMesh)      :: MESH
    INTEGER                   :: allocStat
    !--------------------------------------------------------------------------
    INTENT(IN)                :: IO, EQN
    INTENT(INOUT)             :: MESH
    !--------------------------------------------------------------------------

    ALLOCATE( MESH%VRTX%xyNode(                EQN%Dimension, MESH%nNode            ), &
              MESH%VRTX%Reference(             MESH%nNode                           ), &
              MESH%VRTX%NrOfElementsConnected( MESH%nNode                           ), &
              MESH%VRTX%Element(               MESH%nNode, MESH%MaxElementsOnVertex ), &
              MESH%VRTX%BoundaryToObject(      MESH%nNode                           ), &
              STAT=allocStat                                                           )

    IF (allocStat .NE. 0) THEN
       logError(*) 'allocate_mesh_level0_2: could not allocate the whole mesh!'
       STOP
    END IF
    
    MESH%VRTX%xyNode                = 0                                                !
    MESH%VRTX%Reference             = 0                                                !
    MESH%VRTX%NrOfElementsConnected = 0                                                !
    MESH%VRTX%Element               = 0                                                !
    MESH%VRTX%BoundaryToObject      = 0                                                !

  END SUBROUTINE allocate_mesh_level0_2

  SUBROUTINE destruct_mesh_level0_2(MESH)
    !--------------------------------------------------------------------------
    TYPE (tUnstructMesh)      :: MESH               
    !--------------------------------------------------------------------------
    DEALLOCATE(MESH%VRTX%xyNode               , &
               MESH%VRTX%Reference            , &
               MESH%VRTX%NrOfElementsConnected, &
               MESH%VRTX%Element              , &
               MESH%VRTX%BoundaryToObject        )

  END SUBROUTINE destruct_mesh_level0_2

END MODULE allocate_mesh_mod
