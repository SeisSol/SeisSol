!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Josep de la Puente Alvarez (josep.delapuente AT bsc.es, http://www.geophysik.uni-muenchen.de/Members/jdelapuente)
!!
!! @section LICENSE
!! Copyright (c) 2008, SeisSol Group
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

MODULE MPIExtractMesh_mod
  !----------------------------------------------------------------------------
  IMPLICIT NONE
  PRIVATE
  !----------------------------------------------------------------------------
  INTERFACE MPIExtractMesh
     MODULE PROCEDURE MPIExtractMesh
  END INTERFACE
  !----------------------------------------------------------------------------
  PUBLIC :: MPIExtractMesh
  !----------------------------------------------------------------------------

CONTAINS

  SUBROUTINE MPIExtractMesh(EQN,DISC,BND,MESH,IO,MPI)
    !--------------------------------------------------------------------------
    USE typesDef
    USE allocate_mesh_mod
    USE common_operators_mod
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Argument list declaration
    TYPE (tEquations)         :: EQN
    TYPE (tDiscretization)    :: DISC
    TYPE (tBoundary)          :: BND
    TYPE (tUnstructMesh)      :: MESH
    TYPE (tInputOutput)       :: IO
    TYPE (tMPI)               :: MPI
    ! Local variable declaration
    TYPE (tUnstructMesh)      :: NewMesh
    TYPE tNeighborDomain
        INTEGER               :: nElem                 ! number of conforming elements at interface of 2 particular CPUs
        INTEGER               :: nDomElem              ! number of non-conforming domain elements at interface of 2 particular CPUs
        INTEGER               :: nNeighElem            ! number of non-conforming neighbor elements at interface of 2 particular CPUs
        INTEGER               :: nSendElem             ! number of non-conforming domain elements which have to send their information
   !     INTEGER, POINTER      :: NC_NeighborOfElem(:,:)! Neighbor(s) of non-conforming boundary element at particular GP
   !     INTEGER, POINTER      :: NC_DomainElements(:)  ! list of nc elements at particular nc boundary
   !     INTEGER, POINTER      :: NC_NeighElements(:)   ! list of nc neighbor elements at particular nc boundary
   !     INTEGER, POINTER      :: NC_SendElements(:)    ! list of nc domain elements which have to be sent
        INTEGER, POINTER      :: NeighborElements(:)   ! list of elements at neighbor CPU
        INTEGER, POINTER      :: DomainElements(:)     ! list of elements at particular MPI boundary
        INTEGER               :: nElem_DR
    END TYPE tNeighborDomain
    TYPE tDomainInfo
        INTEGER               :: nCPU
        INTEGER, POINTER      :: CPU(:)
        TYPE(tNeighborDomain), POINTER :: NeighborDomain(:)
    END TYPE tDomainInfo
    TYPE(tDomainInfo), POINTER :: DomainInfo(:)
    INTEGER                   :: iElem, iNode, iSide, iNeighbor, counter, i, j
    INTEGER                   :: CPU1, CPU2, iDomain, nElem
    INTEGER                   :: openstat, allocstat
    INTEGER                   :: nNeighborElem, nNodeMax
    INTEGER, POINTER          :: NodeList(:), ElementList(:)
    INTEGER, POINTER          :: invElementList(:), invNodeList(:)
    INTEGER                   :: iCPU, jCPU, OldElementNr
    INTEGER                   :: jElem
    INTEGER                   :: MaxSizeNgbElem, MaxSizeDomElem
    INTEGER, POINTER          :: MPIElementList(:,:) 
    INTEGER, PARAMETER        :: nHist=10
    INTEGER                   :: eHist(MESH%nElem), NelHist(nHist)
    INTEGER                   :: iHist
    INTEGER, POINTER          :: CPUGlobal2Local(:), CPULocal2Global(:)
    INTEGER                   :: nNeighborCPU 
    REAL                      :: h1,h2,dh
    LOGICAL                   :: inList
    ! nonconforming boundaries
    INTEGER                   :: iBndGP, nDomainElem, nNeighElem, nSendElem
    INTEGER                   :: iNCBoundry
    INTEGER                   :: iNeighborSide, incNeighbor, neigh_index
    !INTEGER                   :: elemcount(0:MPI%nCPU-1,0:MPI%nCPU-1) ! counts nc domain elements
    !INTEGER                   :: neighcount(0:MPI%nCPU-1,0:MPI%nCPU-1)! counts nc neighbor elements
    !INTEGER                   :: ncDomElements(0:MPI%nCPU-1,0:MPI%nCPU-1,MESH%nNonConformingEdges)  ! list for nc domain elements
    !INTEGER                   :: ncNeighElements(0:MPI%nCPU-1,0:MPI%nCPU-1,MESH%nNonConformingEdges)! list for nc neighbor elements 
    !INTEGER, POINTER          :: NC_ElementList(:), inv_NC_ElementList(:)
    INTEGER                   :: naCPU 
    !--------------------------------------------------------------------------
    ! Dynamic Rupture variables
    INTEGER                   :: nNeighborElem_DR
    INTEGER                   :: DR_list
    INTEGER                   :: eType,NewnSide
    INTEGER                   :: nLocalFaultedges                              !< Nr of local fault edges in my domain 
    INTEGER                   :: iLocalNeighborSide,iFace
    REAL, ALLOCATABLE         :: NewFace(:,:,:)
    REAL, ALLOCATABLE         :: NewNormal(:,:),NewTangent1(:,:),NewTangent2(:,:)
    !--------------------------------------------------------------------------
    INTENT(IN)                :: IO
    INTENT(INOUT)             :: MESH, BND, MPI, EQN
    !--------------------------------------------------------------------------

    IF(MPI%nCPU.EQ.1) THEN
        RETURN
    ENDIF

#ifdef PARALLEL

    logDebug(*) 'entered MPIExtractMesh.'

    ALLOCATE(ElementList(MESH%nElem), invElementList(MESH%nElem+1),                       &
             NodeList(MESH%nNode), invNodeList(MESH%nNode),                               &
             !NC_ElementList(MESH%nNonConformingEdges),                                    &
             !inv_NC_ElementList(MESH%nNonConformingEdges),                                &
             !MESH%ELEM%NC_BoundaryToObject(MESH%nNonConformingEdges,DISC%Galerkin%nBndGP),&
             STAT = allocstat                                                             )
    IF(allocstat.NE.0) THEN
      logError(*) 'Allocate Error 1 in CPU ', MPI%myrank
      STOP
    ENDIF
    
    !CALL allocate_ncmesh_level0_1( & !allocating array for nc_elem
    !     IO   = IO,                & 
    !     MESH = Mesh,              &
    !     DISC = DISC               )
 
    logInfo(*) '  Starting MPIExtractMesh '
    ! ----------------------------------------- begin new part ----------------------------------    
    !
    ! Count the number of neighbor CPUs of MPI%myrank
    !
    ALLOCATE( CPULocal2Global(0:MPI%nCPU-1) ) 
    ALLOCATE( CPUGlobal2Local(0:MPI%nCPU-1) ) 
    CPULocal2Global = -1 
    CPUGlobal2Local = -1 
    CPULocal2Global(0) = MPI%myrank 
    nNeighborCPU = 1
    DO iElem = 1, MESH%nElem 
        CPU1 = MPI%CPUDistribution(MESH%ElemLocal2Global(iElem))         
        IF(CPU1.NE.MPI%myrank) THEN
            CYCLE
        ENDIF
        DO iSide = 1, MESH%LocalElemType(iElem)
            SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
            CASE(0,3,6) ! element and neighbor conform at interface
                ! iNeighbor has value MESH%nElem + 1 if there is no neighbor 
                !(or at non-conforming boundaries) 
                iNeighbor = MESH%ELEM%SideNeighbor(iSide,iElem) 
                IF(iNeighbor.GT.MESH%nElem) THEN
                    CYCLE
                ENDIF 
                CPU2 = MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor))
                ! Both CPU numbers differ, so add neighbor processor number to list 
                ! and increase the number of neighbor elements
                IF(CPU1.NE.CPU2) THEN
                    CALL AddNodeToList(CPULocal2Global,CPU2,nNeighborCPU,MPI%nCPU) 
                ENDIF
           END SELECT
       ENDDO
    ENDDO
    logInfo(*) '  CPU ', MPI%myrank, ' has ', nNeighborCPU, ' neighbor domains: '
    logInfo(*) '  CPULocal2Global: ', CPULocal2Global(0:nNeighborCPU-1)
    ! 
    DO i = 0, nNeighborCPU-1   
        CPUGlobal2Local( CPULocal2Global(i) ) = i 
    ENDDO
    ! 
    ! ----------------------------------------- end new part ----------------------------------
    !
    ! Generate communication information across domain boundaries
    !
    logInfo(*) '  Generating communication structure for CPU .... ', MPI%myrank
    !
    ALLOCATE(DomainInfo(0:nNeighborCPU-1), STAT = allocstat)
    IF(allocstat.NE.0) THEN
       logError(*) 'MPIExtractMesh ALLOCATE Error 1 in CPU ', MPI%myrank
       STOP 
    ENDIF
    !
    ! Init data structure
    !
    DO iCPU = 0, nNeighborCPU - 1
        IF(iCPU.EQ.0) THEN
            naCPU = nNeighborCPU
        ELSE
            naCPU = nNeighborCPU  
        ENDIF
        DomainInfo(iCPU)%nCPU = 0
        ALLOCATE( DomainInfo(iCPU)%CPU(naCPU) )
        ALLOCATE( DomainInfo(iCPU)%NeighborDomain(0:naCPU) )
        DO i = 0, naCPU
            DomainInfo(iCPU)%NeighborDomain(i)%nElem      = 0
            DomainInfo(iCPU)%NeighborDomain(i)%nElem_DR   = 0 ! * Dynamic Rupture
            DomainInfo(iCPU)%NeighborDomain(i)%nDomElem   = 0
            DomainInfo(iCPU)%NeighborDomain(i)%nNeighElem = 0
            DomainInfo(iCPU)%NeighborDomain(i)%nSendElem  = 0
        ENDDO
    ENDDO
    ! Init data structure for non-conforming boundaries
    !elemcount(:,:)  = 0
    !neighcount(:,:) = 0

    ! Loop to prepare allocation of variables of communicating processors 
    DO iElem = 1, MESH%nElem
        CPU1 = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iElem)) ) 
        IF(CPU1.EQ.-1) THEN
            CYCLE
        ENDIF        
        DO iSide = 1, MESH%LocalElemType(iElem)
            SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
            CASE(0,3,6) ! element and neighbor conform at interface
                ! iNeighbor has value MESH%nElem + 1 if there is no neighbor 
                !(or at non-conforming boundaries) 
                iNeighbor = MESH%ELEM%SideNeighbor(iSide,iElem) 
                IF(iNeighbor.GT.MESH%nElem) THEN
                    CYCLE
                ENDIF 
                CPU2 = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor)) ) 
                ! Both CPU numbers differ, so add neighbor processor number to list 
                ! and increase the number of neighbor elements
                IF(CPU1.NE.CPU2) THEN
                    ! However, we only add the information to the list, if one of the CPUs is myrank, i.e. it has local number 0
                    ! otherwise, the information is irrelevant and hence does not have to be stored! 
                    IF(CPU1*CPU2.EQ.0) THEN
                        CALL AddNodeToList(DomainInfo(CPU1)%CPU(:),CPU2,DomainInfo(CPU1)%nCPU,MPI%nCPU)
                        DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem = DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem + 1
                        ! * Dynamic Rupture
                        IF (MESH%ELEM%Reference(iSide,iElem).EQ.3) THEN
                            DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem_DR = DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem_DR + 1
                        ENDIF
                        ! *
                    ENDIF 
                ENDIF
            !CASE(2) ! non-conforming boundary, more than 1 neighbor possible
            !     iNCboundary = MESH%ELEM%ncIndex(iElem,iSide)
            !     DO iBndGP = 1, DISC%Galerkin%nBndGP ! search neighbor for each GP at interface
            !        incNeighbor = MESH%ELEM%ncBndNeighbor(iNCboundary,iBndGP)
            !        iNeighbor   = MESH%ELEM%NCB_IndexList(incNeighbor,1)
            !        CPU2 = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor)) ) 
            !        ! Both CPU numbers differ, so add neighbor processor number to list of neighbor CPUs.
            !        ! Count domain elements and neighbor elements (as number might differ).
            !        ! The number and list of domain elements is only needed inside this subroutine,
            !        ! as it will be MESH%nNonConformingEdges and the NCB_IndexList after partitioning.
            !        IF(CPU1.NE.CPU2) THEN
            !         IF(CPU1*CPU2.EQ.0) THEN 
            !           CALL AddNodeToList(DomainInfo(CPU1)%CPU(:),CPU2,DomainInfo(CPU1)%nCPU,MPI%nCPU)
            !            ! count elements on my side (CPU1)
            !           CALL AddNodeToList(ncDomElements(CPU1,CPU2,:),iNCboundary,elemcount(CPU1,CPU2),MESH%nNonConformingEdges)
            !           DomainInfo(CPU1)%NeighborDomain(CPU2)%nDomElem = elemcount(CPU1,CPU2)
            !            ! count elements of neighbor (CPU2) which are needed on my side (CPU1)
            !            CALL AddNodeToList(ncNeighElements(CPU1,CPU2,:),incNeighbor,neighcount(CPU1,CPU2),MESH%nNonConformingEdges)
            !            DomainInfo(CPU1)%NeighborDomain(CPU2)%nNeighElem = neighcount(CPU1,CPU2)
            !          ENDIF 
            !        ENDIF
            !     ENDDO
            END SELECT             
        ENDDO ! iSide
    ENDDO ! iElem

    ! now, allocate variables of communicating processors
    DO iCPU = 0, nNeighborCPU - 1
        DO i = 0, nNeighborCPU - 1
            IF(i.EQ.iCPU) THEN
                CYCLE
            ENDIF
            nNeighborElem = DomainInfo(iCPU)%NeighborDomain(i)%nElem
            IF(nNeighborElem.GT.0) THEN
                ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NeighborElements(nNeighborElem)) 
                ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%DomainElements(nNeighborElem)) 
            ELSE
                ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NeighborElements(1) ) 
                ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%DomainElements(1)   ) 
            ENDIF
            DomainInfo(iCPU)%NeighborDomain(i)%nElem = 0
            !
            ! * Dynamic Rupture
            !nNeighborElem_DR = DomainInfo(iCPU)%NeighborDomain(i)%nElem_DR
            !IF(nNeighborElem_DR .GT. 0) THEN
            !    ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%Domain_Fault_Elem(nNeighborElem_DR))
            !ENDIF
            !DomainInfo(iCPU)%NeighborDomain(i)%nElem_DR = 0 ! used as an index later
            ! *
            ! non-conforming lements
            nDomainElem = DomainInfo(iCPU)%NeighborDomain(i)%nDomElem
            nNeighElem  = DomainInfo(iCPU)%NeighborDOmain(i)%nNeighElem
            IF(nDomainElem.GT.0) THEN
                !ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_NeighborOfElem(nDomainElem,DISC%Galerkin%nBndGP))
                !ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_DomainElements(nDomainElem))
                !ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_NeighElements(nNeighElem))
            ELSE
                !ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_NeighborOfElem(1,1))
                !ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_DomainElements(1))
                !ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_NeighElements(1))
            ENDIF
            DomainInfo(iCPU)%NeighborDomain(i)%nDomElem              = 0
            !DomainInfo(iCPU)%NeighborDomain(i)%NC_DomainElements(:)  = 0
            !DomainInfo(iCPU)%NeighborDomain(i)%NC_NeighElements(:)   = 0
            !DomainInfo(iCPU)%NeighborDomain(i)%NC_NeighborOfElem(:,:)= -1
            !elemcount(iCPU,i)                                        = 0
            !neighcount(iCPU,i)                                       = 0
            !ncDomElements(iCPU,i,:)                                  = 0
            !ncNeighElements(iCPU,i,:)                                = 0
        ENDDO
    ENDDO
    !
    ALLOCATE(MESH%ELEM%MPINumber(MESH%nVertexMax,MESH%nElem),                      &
             MESH%ELEM%MPINumber_DR(MESH%nVertexMax,MESH%nElem),                   &
             !MESH%ELEM%MPI_NCNumber(MESH%nNonConformingEdges,DISC%Galerkin%nBndGP),&
             STAT = allocstat)
    IF(allocstat.NE.0) THEN
       logError(*) 'MPIExtractMesh ALLOCATE Error 2 in CPU ', MPI%myrank
       STOP 
    ENDIF
    !
    MESH%ELEM%MPINumber(:,:)    = -1
    !MESH%ELEM%MPI_NCNumber(:,:) = -1
    MESH%ELEM%MPINumber_DR      = -1
    !
    MESH%ELEM%MPIReference(:,:) = 0
    !
    ! Fill neighbor lists and domain lists, where index connects adjacent elements for conforming elements
    ! Note the number of the neighbored CPU in which the neighbor element lies
    DO iElem = 1, MESH%nElem
      CPU1 = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iElem)) ) 
      DO iSide = 1, MESH%LocalElemType(iElem)
         SELECT CASE(MESH%ELEM%Reference(iSide,iElem))
         CASE(0,3,6)
            iNeighbor = MESH%ELEM%SideNeighbor(iSide,iElem)
            IF(iNeighbor.GT.MESH%nElem) THEN
                CYCLE
            ENDIF 
            CPU2 = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor)) ) 
            IF((CPU1.NE.CPU2).AND.(CPU1*CPU2.EQ.0)) THEN
                ! Increase element counter and assign neighbor
                DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem = DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem + 1
                nNeighborElem = DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem
                ! Add neighbor element and domain element to list
                DomainInfo(CPU1)%NeighborDomain(CPU2)%NeighborElements(nNeighborElem) = MESH%ElemLocal2Global(iNeighbor) 
                DomainInfo(CPU1)%NeighborDomain(CPU2)%DomainElements(nNeighborElem)   = MESH%ElemLocal2Global(iElem) 
                ! * Dynamic Rupture
                IF (MESH%ELEM%Reference(iSide,iElem) .EQ. 3) THEN
                    DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem_DR = DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem_DR + 1
                    nNeighborElem_DR = DomainInfo(CPU1)%NeighborDomain(CPU2)%nElem_DR
                    !DomainInfo(CPU1)%NeighborDomain(CPU2)%Domain_Fault_Elem(nNeighborElem_DR) = MESH%ElemLocal2Global(iElem)
                ENDIF
                ! *
                ! Create interior MPI boundary (bnd. condition 7) only for non-unified version
                ! unified version uses variable MPIReference                
                IF(CPU1.EQ.0) THEN
                    MESH%Elem%MPIReference(iSide,iElem) = 1
                    ! search index of CPU2 in list of neighbored CPUs for the according element (and side) 
                    CALL FindNodeInList(iDomain,DomainInfo(CPU1)%CPU(:),CPU2,DomainInfo(CPU1)%nCPU,nNeighborCPU)
                    IF(iDomain.EQ.-1) THEN
                        logError(*) 'ERROR: CPU ', CPU2, ' not found in list for CPU ', CPU1, ' ! ', &
                           '# domains ', DomainInfo(CPU1)%nCPU, ' list: ', DomainInfo(CPU1)%CPU(:)
                        STOP
                    ELSE
                        MESH%ELEM%BoundaryToObject(iSide,iElem) = iDomain
                    ENDIF
                ENDIF
            ENDIF
         !CASE(2)
         !   iNCboundary = MESH%ELEM%ncIndex(iElem,iSide)
         !   DO iBndGP = 1, DISC%Galerkin%nBndGP
         !       incNeighbor = MESH%ELEM%ncBndNeighbor(iNCboundary,iBndGP)
         !       iNeighbor   = MESH%ELEM%NCB_IndexList(incNeighbor,1)
         !       CPU2 = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor)) ) 
         !       IF((CPU1.NE.CPU2).AND.(CPU1*CPU2.EQ.0)) THEN
         !           CALL AddNodeToList(ncDomElements(CPU1,CPU2,:),iNCboundary,elemcount(CPU1,CPU2),MESH%nNonConformingEdges)
         !           DomainInfo(CPU1)%NeighborDomain(CPU2)%nDomElem = elemcount(CPU1,CPU2)
         !           nDomainElem = elemcount(CPU1,CPU2)
         !           DomainInfo(CPU1)%NeighborDomain(CPU2)%NC_NeighborOfElem(nDomainElem,iBndGP) = incNeighbor
         !           DomainInfo(CPU1)%NeighborDomain(CPU2)%NC_DomainElements(nDomainElem) = iNCboundary
         !           ! Here the neighbors do NOT compulsory have the same index as the adjacent element
         !           CALL AddNodeToList(ncNeighElements(CPU1,CPU2,:),incNeighbor,neighcount(CPU1,CPU2),MESH%nNonConformingEdges)
         !           DomainInfo(CPU1)%NeighborDomain(CPU2)%nNeighElem = neighcount(CPU1,CPU2)
         !           DomainInfo(CPU1)%NeighborDomain(CPU2)%NC_NeighElements(neighcount(CPU1,CPU2))=incNeighbor
         !           IF(CPU1.EQ.0) THEN
         !               MESH%Elem%MPIReference(iSide,iElem) = 1
         !               CALL FindNodeInList(iDomain,DomainInfo(CPU1)%CPU(:),CPU2,DomainInfo(CPU1)%nCPU,MPI%nCPU)
         !               IF(iDomain.EQ.-1) THEN
         !                   PRINT *, 'ERROR: CPU ', CPU2, ' not found in list for CPU ', CPU1, ' ! ', & 
         !                     '# domains ', DomainInfo(CPU1)%nCPU, ' list: ', DomainInfo(CPU1)%CPU(:)
         !                   STOP
         !               ELSE
         !                   MESH%ELEM%NC_BoundaryToObject(iNCboundary,iBndGP) = iDomain
         !               ENDIF
         !           ENDIF
         !       ENDIF
         !   ENDDO
         END SELECT
      ENDDO
      IF(MESH%nElem.GT.20) THEN
          IF (MOD(iElem,FLOOR(MESH%nElem/20.)).EQ.0) THEN
              logInfo(*) ' ', iElem,' elements done...'
          END IF
      ENDIF
    ENDDO 
    
   ! It can happen that two elements of different CPUs share only a small part
   ! of one edge. That means, one element can have a GP within this small interface
   ! (here the element of CPU1), whereas the adjacent element (here of CPU 2)
   ! has no GP inside this interval!
   ! Hence, we need to send the information of the element without any GP
   ! (which might not be included in the list of domain elements) to the element 
   ! with a GP (which needs the neighbor information). Therefore, we create a "SEND"-List,
   ! which contains exactly those elements needed by the neighbor CPU.
   !
   !                                  CPU 2
   !                                   ___
   !           _______________________|___|_________________________
   !                            |      |
   !                            |______|
   !                  CPU 1            |        CPU 3
   !                                   |
   !                                   |


    DO iCPU = 0, nNeighborCPU - 1
        DO i = 0, nNeighborCPU - 1
            IF(i.EQ.iCPU) THEN
                CYCLE
            ENDIF
            nSendElem = DomainInfo(i)%NeighborDomain(iCPU)%nNeighElem
            DomainInfo(iCPU)%NeighborDomain(i)%nSendElem = nSendElem
            !ALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_SendElements(nSendElem))
            !DomainInfo(iCPU)%NeighborDomain(i)%NC_SendElements(:)=0
        ENDDO
    ENDDO

    !neighcount(:,:) = 0

    ! Now, we add the missing send-elements.
    !DO iNCboundary = 1, MESH%nNonConformingEdges
    !    iElem = MESH%ELEM%NCB_IndexList(iNCboundary,1)
    !    CPU1  = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iElem)) ) 
    !    DO iBndGP = 1, DISC%Galerkin%nBndGP
    !        incNeighbor = MESH%ELEM%ncBndNeighbor(iNCboundary,iBndGP)
    !        iNeighbor   = MESH%ELEM%NCB_IndexList(incNeighbor,1)
    !        CPU2 = CPUGlobal2Local( MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor)) ) 
    !        IF((CPU1.NE.CPU2).AND.(CPU1*CPU2.EQ.0)) THEN
    !            nSendElem   = DomainInfo(CPU2)%NeighborDomain(CPU1)%nSendElem
    !            CALL AddNodeToList(DomainInfo(CPU2)%NeighborDomain(CPU1)%NC_SendElements(:),incNeighbor,neighcount(CPU2,CPU1),nSendElem)
    !        ENDIF
    !    ENDDO ! iBndGP
    !ENDDO ! iNCboundary
   !
    logInfo(*) '    Communication structure summary for CPU ', MPI%myrank
    logInfo(*) '    Number of neighbor domains :  ', DomainInfo(0)%nCPU
    logInfo(*) '    Neighbor domains           :  ', &
                                 DomainInfo(0)%CPU(1:DomainInfo(0)%nCPU)
 
    DO i = 1, DomainInfo(0)%nCPU
        logInfo(*) '    CPU ', MPI%myrank, ' No. of regular neighbors with CPU      ', &
                                     DomainInfo(0)%CPU(i), ' : ',                  &
                                     DomainInfo(0)%NeighborDomain(DomainInfo(0)%CPU(i))%nElem
        logInfo(*) '    CPU ', MPI%myrank, ' No. of non-conforming elements at interface with CPU ', &
                                     DomainInfo(0)%CPU(i), ' : ',                  &
                                     DomainInfo(0)%NeighborDomain(DomainInfo(0)%CPU(i))%nDomElem
        logInfo(*) '    CPU ', MPI%myrank, ' No. of non-conforming neighbors at interface with CPU ', &
                                     DomainInfo(0)%CPU(i), ' : ',                  &
                                     DomainInfo(0)%NeighborDomain(DomainInfo(0)%CPU(i))%nNeighElem
        logInfo(*) '    CPU ', MPI%myrank, ' No. of needed non-conforming elements at interface with CPU ', &
                                     DomainInfo(0)%CPU(i), ' : ',                  &
                                     DomainInfo(0)%NeighborDomain(DomainInfo(0)%CPU(i))%nSendElem
        ! Dynamic Rupture
        logInfo(*) '    CPU ', MPI%myrank, ' No. of fault neighbors (each single side) with CPU      ', &
                                     DomainInfo(0)%CPU(i), ' : ',                  &
                                     DomainInfo(0)%NeighborDomain(DomainInfo(0)%CPU(i))%nElem_DR                                    
    ENDDO     
   
    !
    ! Save communication structure
    ! 

    BND%NoMPIDomains = DomainInfo(0)%nCPU
    ALLOCATE( BND%ObjMPI(BND%NoMPIDomains) ) 
    !
    DO iDomain = 1, BND%NoMPIDomains
       CPU1        = 0
       CPU2        = DomainInfo(0)%CPU(iDomain)
       nElem       = DomainInfo(0)%NeighborDomain(CPU2)%nElem
       nDomainElem = DomainInfo(0)%NeighborDomain(CPU2)%nDomElem
       nNeighElem  = DomainInfo(0)%NeighborDomain(CPU2)%nNeighElem 
       nSendElem   = DomainInfo(0)%NeighborDomain(CPU2)%nSendElem   
       !
       BND%ObjMPI(iDomain)%CPU          = CPULocal2Global(CPU2) 
       BND%ObjMPI(iDomain)%nElem        = nElem
       ! * Dynamic Rupture
       BND%ObjMPI(iDomain)%nFault_MPI   = DomainInfo(0)%NeighborDomain(CPU2)%nElem_DR
       ! *
       !BND%ObjMPI(iDomain)%nNCNeighElem = nNeighElem
       !BND%ObjMPI(iDomain)%nNCSendElem  = nSendElem
              
       ALLOCATE( BND%ObjMPI(iDomain)%DomainElements(nElem),                                           &
                 BND%ObjMPI(iDomain)%NeighborCoords(EQN%dimension,MESH%nVertexMax,nElem),             &
                 BND%ObjMPI(iDomain)%Domain_Fault_Elem(BND%ObjMPI(iDomain)%nFault_MPI),               & ! * Dynamic Rupture
                 !BND%ObjMPI(iDomain)%NC_SendElements(nSendElem),                                      &
                 !BND%ObjMPI(iDomain)%NC_NeighborCoords(EQN%Dimension,MESH%nVertexMax,nNeighElem),     & 
                 !BND%ObjMPI(iDomain)%NC_LocalElemType(nNeighElem),                                    &
                 stat = allocstat  )
       IF(allocstat.NE.0) THEN
          logError(*) 'Could not allocate all MPI boundary fields in domain ', iDomain
          STOP
       ENDIF      
       !
       ! The CPU with the lower number enforces the order of the boundary elements for conforming boundaries!
       ! Therewith, the indices of the list of adjacent CPUs refer to neighbored elements.
       ! Since the order in the list may have changed according to the read_mesh_tetra_geizig subroutine, 
       ! we must SORT all the lists to assure consistency and the same order of the lists, which are still
       ! containing GLOBAL element numbers (from the original mesh-file). 
       ! 

       call sort(DomainInfo(CPU1)%NeighborDomain(CPU2)%DomainElements,DomainInfo(CPU1)%NeighborDomain(CPU2)%NeighborElements)
       call sort(DomainInfo(CPU2)%NeighborDomain(CPU1)%DomainElements,DomainInfo(CPU2)%NeighborDomain(CPU1)%NeighborElements)
       !
       DR_list = 0
       DO i = 1, BND%ObjMPI(iDomain)%nElem
         IF(CPULocal2Global(CPU1).LT.CPULocal2Global(CPU2)) THEN
            ! Domain of CPU number MPI%myrank has the lower CPU number, so use the element ordering seen
            ! from this CPU
            iElem       = MESH%ElemGlobal2Local(DomainInfo(CPU1)%NeighborDomain(CPU2)%DomainElements(i)) 
            iNeighbor   = MESH%ElemGlobal2Local(DomainInfo(CPU1)%NeighborDomain(CPU2)%NeighborElements(i)) 
         ELSE
            ! Domain of CPU number MPI%myrank has the higher CPU number, so use the element ordering seen
            ! from the neighbor CPU, which has the lower number
            iElem     = MESH%ElemGlobal2Local(DomainInfo(CPU2)%NeighborDomain(CPU1)%NeighborElements(i)) 
            iNeighbor = MESH%ElemGlobal2Local(DomainInfo(CPU2)%NeighborDomain(CPU1)%DomainElements(i))
         ENDIF
         BND%ObjMPI(iDomain)%DomainElements(i) = iElem
         DO j = 1, MESH%LocalElemType(iElem) ! We assume conforming elements of the same type
            iNode = MESH%ELEM%Vertex(j,iNeighbor)
            BND%ObjMPI(iDomain)%NeighborCoords(:,j,i) = MESH%VRTX%xyNode(:,iNode)
         ENDDO
         DO iSide = 1, MESH%LocalElemType(iElem)
            IF(MESH%ELEM%SideNeighbor(iSide,iElem).EQ.iNeighbor) THEN
              ! counter i stands for the index of the neighbor inside the neighbored domain
              ! (because of the ordering)
              MESH%ELEM%MPINumber(iSide,iElem) = i
              ! * Dynamic Rupture
              IF (MESH%ELEM%Reference(iSide,iElem).EQ.3) THEN
                 DR_list = DR_list + 1
                 MESH%ELEM%MPINumber_DR(iSide,iElem) = DR_list
                 IF (DR_list .GT. BND%ObjMPI(iDomain)%nFault_MPI) THEN
                     PRINT*,'DR_list too large!'
                     STOP
                 ENDIF
                 BND%ObjMPI(iDomain)%Domain_Fault_Elem(DR_list) = iElem
              ENDIF
              ! *
            ENDIF
         ENDDO
         IF(MAXVAL(MESH%ELEM%MPINumber(:,iElem)).EQ.-1) THEN
            logError(*) 'ERROR in building MPI communication structure. '
            logError(*) i,iElem,iNeighbor,MESH%ELEM%SideNeighbor(:,iElem)
            STOP
         ENDIF
       ENDDO

       ! For non-conforming boundaries do not order the elements but save everything from both points of view.
       ! The elements of the send list are saved as ncelements by now. For exchanging DOFs we need the element 
       ! number inside the mesh.
       
       !DO i = 1, BND%ObjMPI(iDomain)%nNCSendElem
       !   iNCboundary = DomainInfo(CPU1)%NeighborDomain(CPU2)%NC_SendElements(i)
       !   iElem = MESH%Elem%NCB_IndexList(iNCboundary,1)
       !   BND%ObjMPI(iDomain)%NC_SendElements(i) = iElem
       !ENDDO
       
       ! From now on, we set nSendElem to the number of send elements of the neighbor
       nSendElem = DomainInfo(CPU2)%NeighborDomain(0)%nSendElem
       
       ! Save the index of the neighbor element to be able to recall all needed information (DOFs etc.)
       !DO i = 1, nDomainElem ! BND%ObjMPI(iDomain)nNCDomElem
       !   ! info from my (MPI%myrank) point of view
       !   iNCboundary = DomainInfo(CPU1)%NeighborDomain(CPU2)%NC_DomainElements(i)
       !   iElem       = MESH%ELEM%NCB_IndexList(iNCboundary,1)
       !   DO iBndGP = 1, DISC%Galerkin%nBndGP
       !      ! Element might simultaneously have GPs which need neighbors inside own and neighbor domain
       !      IF(DomainInfo(CPU1)%NeighborDomain(CPU2)%NC_NeighborOfElem(i,iBndGP) .NE. -1)THEN ! on myrank
       !         incNeighbor = DomainInfo(CPU1)%NeighborDomain(CPU2)%NC_NeighborOfElem(i,iBndGP)
       !         iNeighbor   = MESH%Elem%NCB_IndexList(incNeighbor,1)
       !         CALL FindNodeInList(neigh_index,DomainInfo(CPU2)%NeighborDomain(CPU1)%NC_SendElements(:),incNeighbor,nSendElem,nSendElem)
       !         IF(neigh_index.EQ.-1)THEN
       !             PRINT*,'ERROR in building MPI communication structure. '
       !             PRINT*,'Cannot find neighbor for element ', iElem
       !            EXIT
       !         ENDIF
       !         MESH%ELEM%MPI_NCNumber(iNCboundary,iBndGP) = neigh_index                
       !         ! Save neighbor element type
       !         BND%ObjMPI(iDomain)%NC_LocalElemType(neigh_index) = MESH%LocalElemType(iNeighbor)
       !         ! Save neighbor coordinates
       !         DO j = 1, MESH%LocalElemType(iNeighbor)
       !              iNode = MESH%ELEM%Vertex(j,iNeighbor)
       !              BND%ObjMPI(iDomain)%NC_NeighborCoords(:,j,neigh_index) = MESH%VRTX%xyNode(:,iNode)
       !         ENDDO
       !         ! mark in usual array, that neighbor is on other CPU
       !         Mesh%ELEM%ncBndNeighbor(iNCboundary,iBndGP) = MESH%nNonConformingEdges+1
       !     ENDIF
       !  ENDDO   
       !ENDDO
    ENDDO
    !
    logInfo(*) 'Chosen communication structure for CPU ', MPI%myrank, ' : '
    logInfo(*) 'Total number of neighbor domains: ', BND%NoMPIDomains
    ! 
    DO i = 1, BND%NoMPIDomains
        logInfo(*) 'No. of conforming bnd. elements with CPU      ', BND%ObjMPI(i)%CPU, ' : ', &
                   BND%ObjMPI(i)%nElem
    ENDDO
    !
    logInfo(*) 'Deallocating temp. comm fields for CPU ', MPI%myrank, ' : '
    DO iCPU = 0, nNeighborCPU - 1
        DO i = 0, nNeighborCPU - 1
            IF(i.EQ.iCPU) THEN
                CYCLE
            ENDIF
            DEALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NeighborElements) 
            DEALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%DomainElements) 
            !DEALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_DomainElements)
            !DEALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_NeighborOfElem)
            !DEALLOCATE(DomainInfo(iCPU)%NeighborDomain(i)%NC_SendElements)
        ENDDO
        DEALLOCATE(DomainInfo(iCPU)%CPU) 
        DEALLOCATE(DomainInfo(iCPU)%NeighborDomain) 
    ENDDO
    DEALLOCATE(DomainInfo)
    ! 
    logInfo(*) 'Generating new mesh for subdomain number ', MPI%myrank
    !    
    ! 
    ! Create a new mesh that only contains the domain that coincides with MPI%myrank
    !
    NewMesh%nElem                  = 0
    !NewMesh%nNonConformingEdges    = 0
    NewMesh%nNode                  = 0
    counter                        = 0
    Nodelist(:)                    = -1
    ElementList(:)                 = -1
    invElementList(:)              = -1
    invNodeList(:)                 = -1
    !NC_ElementList(:)              = -1
    !inv_NC_ElementList(:)          = -1
    ! 
    ! Find elements within myrank
    !
    DO iElem = 1, MESH%nElem
      IF(MPI%CPUDistribution(MESH%ElemLocal2Global(iElem)).EQ.MPI%myrank) THEN
        NewMesh%nElem = NewMesh%nElem + 1
        ElementList(NewMesh%nElem) = iElem  
        DO iNode = 1, MESH%LocalVrtxType(iElem) 
           CALL AddNodeToList(Nodelist,MESH%ELEM%Vertex(iNode,iElem),counter,MESH%nNode)
        ENDDO
      ENDIF
      invElementList(iElem) = NewMesh%nElem
    ENDDO   
    
    ! Find non-conforming elements obtained by myrank
    !DO iNCboundary = 1, MESH%nNonConformingEdges
    !    iElem = MESH%Elem%NCB_IndexList(iNCboundary,1)
    !    IF(MPI%CPUDistribution(MESH%ElemLocal2Global(iElem)).EQ.MPI%myrank)THEN
    !        NewMesh%nNonConformingEdges = NewMESH%nNonConformingEdges +1
    !        NC_ElementList(NewMesh%nNonConformingEdges) = iNCboundary
    !    ENDIF
    !    ! if iNCboundary belongs to another CPU, then inverse list shows counter of previous element
    !    ! (instead of -1)!
    !    inv_NC_ElementList(iNCboundary) = NewMesh%nNonConformingEdges
    !ENDDO
        
    NewMesh%nNode     = counter
    NewMesh%Dimension = MESH%Dimension
    NewMesh%ADERDG3D  = MESH%ADERDG3D

    DO iNode = 1, NewMesh%nNode
       invNodeList(NodeList(iNode)) = iNode
    ENDDO

    logInfo(*) 'Number of elements found for CPU ', MPI%myrank, ' : ', NewMesh%nElem
    logInfo(*) 'Number of nodes found for CPU    ', MPI%myrank, ' : ', NewMesh%nNode

    NewMesh%GlobalElemType      = MESH%GlobalElemType
    NewMesh%GlobalVrtxType      = MESH%GlobalVrtxType
    NewMesh%GlobalSideType      = MESH%GlobalSideType
    NewMesh%MaxElementsOnVertex = MESH%MaxElementsOnVertex
    NewMesh%nVertexMax          = MESH%nVertexMax
    NewMesh%nSideMax            = MESH%nSideMax

    CALL allocate_mesh_level0_1(  & ! allocating array for elem
         IO   = IO,               & ! allocating array for elem
         MESH = NewMesh           ) ! allocating array for elem
    
    !CALL allocate_ncmesh_level0_1( & !allocating array for nc_elem
    !     IO   = IO,                & 
    !     MESH = NewMesh,           &
    !     DISC = DISC               )    
    !ALLOCATE(newMESH%ELEM%ncIndex(newMESH%nElem,newMESH%nSideMax))
    
    CALL allocate_mesh_level0_2(                                            & ! allocating array for node
         IO   = IO                                                        , & ! allocating array for node
         EQN  = EQN                                                       , & ! allocating array for node
         MESH = NewMesh                                                     ) ! allocating array for node

    logInfo(*) 'Allocating temporary storage...'

    ALLOCATE( NewMESH%ELEM%MPINumber(NewMESH%nVertexMax,NewMESH%nElem))
    ALLOCATE( NewMESH%ELEM%MPINumber_DR(NewMESH%nVertexMax,NewMESH%nElem)) ! * Dynamic Rupture
    !ALLOCATE( NewMESH%ELEM%MPI_NCNumber(NewMESH%nNonConformingEdges,DISC%Galerkin%nBndGP))
    ALLOCATE( NewMESH%LocalElemType(NewMESH%nElem))

    logInfo(*) 'Copying node information to temp. data...'

    ! Copy node information to temporary data structure
    DO iNode = 1, NewMesh%nNode
       NewMesh%VRTX%xyNode(:,iNode) = MESH%VRTX%xyNode(:,Nodelist(iNode))
       !
       NewMESH%VRTX%NrOfElementsConnected(iNode) = 0
       DO iElem = 1, MESH%VRTX%NrOfElementsConnected(Nodelist(iNode))
         IF(MPI%CPUDistribution(MESH%ElemLocal2Global(MESH%VRTX%Element(Nodelist(iNode),iElem))).EQ.MPI%myrank) THEN
            NewMesh%VRTX%NrOfElementsConnected(iNode) = NewMesh%VRTX%NrOfElementsConnected(iNode) + 1
            NewMesh%VRTX%Element(iNode,NewMesh%VRTX%NrOfElementsConnected(iNode)) = &
               invElementList(MESH%VRTX%Element(Nodelist(iNode),iElem))
         ENDIF
       ENDDO
       ! 
    ENDDO

    logInfo(*) 'Copying element information to temp. data...'

    ! Copy element information to temporary data structure
    DO iElem = 1, NewMesh%nElem
       NewMESH%LocalElemType(iElem)=MESH%LocalElemType(ElementList(iElem))
       IF(NewMESH%LocalElemType(iElem).EQ.6)THEN !NewMesh%LocalVrtxType not declared
           nNodeMax = 8
       ELSE
           nNodeMax = NewMESH%LocalElemType(iElem)
       ENDIF
       DO iNode = 1, nNodeMax
          NewMESH%ELEM%Vertex(iNode,iElem) = invNodelist(MESH%ELEM%Vertex(iNode,Elementlist(iElem)))
          IF( (NewMESH%ELEM%Vertex(iNode,iElem).LE.0             ) .OR. &
              (NewMESH%ELEM%Vertex(iNode,iElem).GT.NewMesh%nNode )      )  THEN
             logError(*) iElem, iNode, MESH%ELEM%Vertex(iNode,Elementlist(iELem)), NewMESH%ELEM%Vertex(iNode,iElem)
             STOP
          ENDIF
       ENDDO
       !
       NewMesh%ELEM%MPIReference(:,iElem)       = MESH%ELEM%MPIReference(:,ElementList(iElem))
       NewMesh%ELEM%Reference(:,iElem)          = MESH%ELEM%Reference(:,Elementlist(iElem))
       NewMesh%ELEM%MPINumber(:,iElem)          = MESH%ELEM%MPINumber(:,Elementlist(iElem))
       NewMesh%ELEM%MPINumber_DR(:,iElem)       = MESH%ELEM%MPINumber_DR(:,Elementlist(iElem))! * Dynamic Rupture
       NewMesh%ELEM%BoundaryToObject(:,iElem)   = MESH%ELEM%BoundaryToObject(:,Elementlist(iElem))
       !NewMesh%ELEM%ncIndex(iElem,:)            = 0
       DO iSide = 1, NewMESH%LocalElemType(iElem)
          !IF(MESH%ELEM%ncIndex(ElementList(iElem),iSide) .NE. 0)THEN
          !  NewMesh%ELEM%ncIndex(iElem,iSide) = inv_NC_ElementList(MESH%Elem%ncIndex(ElementList(iElem),iSide))
          !ENDIF
          SELECT CASE(NewMesh%ELEM%Reference(iSide,iElem))
            ! If unified code is used, Reference has still old index (only MPIReference has changed)
          CASE(0,3,6)
                NewMESH%ELEM%SideNeighbor(iSide,iElem) = invElementList(MESH%ELEM%SideNeighbor(iSide,Elementlist(iElem)))
                IF (NewMesh%Elem%MPIReference(iSide,iElem).EQ.1) NewMesh%Elem%SideNeighbor(iSide,iElem) = NewMesh%nElem + 1
          !CASE(2)
          !      NewMESH%ELEM%SideNeighbor(iSide,iElem) = NewMESH%nElem + 1
          !      NewMesh%ELEM%ncIndex(iElem,iSide)      = inv_NC_ElementList( Mesh%ELEM%ncIndex(Elementlist(iElem),iSide) )
          CASE DEFAULT
                ! If there is no neighbor assign neigbor number nElement+1
                NewMESH%ELEM%SideNeighbor(iSide,iElem) = NewMESH%nElem + 1
          END SELECT
       ENDDO
       NewMesh%ELEM%LocalNeighborSide(:,iElem)  = MESH%ELEM%LocalNeighborSide(:,Elementlist(iElem))
       NewMesh%ELEM%LocalNeighborVrtx(:,iElem)  = MESH%ELEM%LocalNeighborVrtx(:,Elementlist(iElem))
       NewMesh%ELEM%Volume(iElem)               = MESH%ELEM%Volume(Elementlist(iElem))
       NewMesh%ELEM%xyBary(:,iElem)             = MESH%ELEM%xyBary(:,Elementlist(iElem))
    ENDDO

    ! for non-conforming boundaries 
    !DO iNCboundary = 1, NewMesh%nNonConformingEdges 
    !      ! ncBndNeighbor = -1 for MPI boundary 
    !      DO iBndGP = 1, DISC%Galerkin%nBndGP 
    !        IF(Mesh%ELEM%ncBndNeighbor(NC_ElementList(iNCboundary),iBndGP) .EQ. MESH%nNonConformingEdges+1)THEN
    !            NewMesh%ELEM%ncBndNeighbor(iNCboundary,iBndGP)= -1
    !        ELSE
    !            NewMesh%ELEM%ncBndNeighbor(iNCboundary,iBndGP) = inv_NC_ElementList( Mesh%ELEM%ncBndNeighbor(NC_ElementList(iNCboundary),iBndGP) ) 
    !        ENDIF
    !      ENDDO
    !      NewMesh%ELEM%NCB_IndexList(iNCboundary,1)        = invElementList(Mesh%ELEM%NCB_IndexList(NC_ElementList(iNCboundary),1) ) 
    !      NewMesh%ELEM%NCB_IndexList(iNCboundary,2)        = Mesh%ELEM%NCB_IndexList(NC_ElementList(iNCboundary),2) 
    !      NewMesh%ELEM%ncBndGaussP(iNCboundary,:,:)        = Mesh%ELEM%ncBndGaussP(NC_ElementList(iNCboundary),:,:)  
    !      NewMesh%ELEM%MPI_NCNumber(iNCboundary,:)         = MESH%ELEM%MPI_NCNumber(NC_Elementlist(iNCboundary),:)  
    !      NewMesh%ELEM%NC_BoundaryToObject(iNCboundary,:)  = MESH%ELEM%NC_BoundaryToObject(NC_ElementList(iNCboundary),:)  
    !ENDDO
   
    
    ! In case of Dynamic Rupture, extract relevant elements for myrank
    IF (EQN%DR == 1) THEN
      !
      ! Test if there is a fault in the local domain
      nLocalFaultedges = 0
      DO iElem = 1, NewMESH%nElem
        eType = NewMESH%LocalElemType(iElem)
        DO iSide = 1, eType
          IF(NewMesh%ELEM%Reference(iSide,iElem).EQ.3) THEN
            nLocalFaultedges = nLocalFaultedges+1
          !  tmp_el(count,1) = iElem
          !  tmp_el(count,2) = iSide
          ENDIF
        ENDDO
      ENDDO
      !
      !
      IF (nLocalFaultedges == 0) THEN
        !
        ! Switch off DR for this CPU
        logInfo(*) 'No Fault present in domain:',MPI%myrank
        EQN%DR = 0
        DISC%DynRup%DR_output = .FALSE.
        DEALLOCATE(MESH%Fault%Face)
        DEALLOCATE(MESH%Fault%geoNormals  , &
                   MESH%Fault%geoTangent1 , &
                   MESH%Fault%geoTangent2   )
        !
      ELSE
        !
        logInfo(*) 'Fault present in domain:',MPI%myrank
        logInfo(*) 'with',nLocalFaultedges,'fault edges.'
        !
        ALLOCATE(NewFace(nLocalFaultedges,2,2))
        ALLOCATE(NewNormal(3,nLocalFaultedges)   , &
                 NewTangent1(3,nLocalFaultedges) , &
                 NewTangent2(3,nLocalFaultedges)   )

        ! If fault element lies not in my domain, set Face to zero which means: obtain value from MPI-neighbor!
        NewFace(:,:,:) = 0
        NewnSide = 0
        MESH%Fault%nPlusSide = 0
        !
        ! Dynamic Rupture output just for domains with "+" elements
        DISC%DynRup%DR_output = .FALSE.
        !
        NewNormal = 0.0D0
        NewTangent1 = 0.0D0
        NewTangent2 = 0.0D0
        !
        DO iFace=1,MESH%Fault%nSide
          !
          iElem               = MESH%Fault%Face(iFace,1,1)
          iSide               = MESH%Fault%Face(iFace,2,1)
          ! 
          iNeighbor           = MESH%Fault%Face(iFace,1,2)
          iLocalNeighborSide  = MESH%Fault%Face(iFace,2,2)
          !
          IF (MPI%CPUDistribution(MESH%ElemLocal2Global(iElem)) == MPI%myrank) THEN
            !
            ! add "+" element to the list
            NewnSide              = NewnSide + 1
            MESH%Fault%nPlusSide  = MESH%Fault%nPlusSide + 1    ! counter for "+" elements
            DISC%DynRup%DR_output = .TRUE.                      ! output "on" for this domain
            NewFace(NewnSide,1,1) = invElementList(iElem)       ! assign to new element
            NewFace(NewnSide,2,1) = iSide                       ! iSide value remains
            !
            NewFace(NewnSide,2,2) = iLocalNeighborSide          ! iSide value remains for the "-" element in any case          
            !
            ! copy normal vector information
            NewNormal(:,NewnSide) = MESH%Fault%geoNormals(:,iFace)
            NewTangent1(:,NewnSide) = MESH%Fault%geoTangent1(:,iFace)
            NewTangent2(:,NewnSide) = MESH%Fault%geoTangent2(:,iFace)
            !
            IF (MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor)) == MPI%myrank) THEN
              !
              ! if present, add corresponding "-" element (Neighbor) to the list
              NewFace(NewnSide,1,2) = invElementList(iNeighbor)  ! assign to new element
              !
            ENDIF
            !
          ELSEIF (MPI%CPUDistribution(MESH%ElemLocal2Global(iNeighbor)) == MPI%myrank &
                 .AND. MPI%CPUDistribution(MESH%ElemLocal2Global(iElem)) .NE. MPI%myrank) THEN
            ! Second IF statement just as a double check
            !
            ! add "-" element (Neighbor) to the list in case that there is no corresponding "+" element
            NewnSide              = NewnSide + 1
            NewFace(NewnSide,1,2) = invElementList(iNeighbor)  ! assign to new element
            NewFace(NewnSide,2,2) = iLocalNeighborSide         ! iSide value remains
            !
            NewFace(NewnSide,2,1) = iSide                  ! iSide value remains for the "+" element in any case
            !
            ! copy normal vector information
            NewNormal(:,NewnSide) = MESH%Fault%geoNormals(:,iFace)
            NewTangent1(:,NewnSide) = MESH%Fault%geoTangent1(:,iFace)
            NewTangent2(:,NewnSide) = MESH%Fault%geoTangent2(:,iFace)
            !
          ENDIF          
          !
        ENDDO ! iFace
        !
        DEALLOCATE(MESH%Fault%Face)
        MESH%Fault%nSide = NewnSide

        ALLOCATE(MESH%Fault%Face(MESH%Fault%nSide,2,2))
        !
        DEALLOCATE(MESH%Fault%geoNormals  , &
                   MESH%Fault%geoTangent1 , &
                   MESH%Fault%geoTangent2   )
        ! allocate with new MESH%Fault%nSide value
        ALLOCATE(MESH%Fault%geoNormals(3,MESH%Fault%nSide)  , &
                 MESH%Fault%geoTangent1(3,MESH%Fault%nSide) , &
                 MESH%Fault%geoTangent2(3,MESH%Fault%nSide)   )
        !
        !
        DO iFace=1,MESH%Fault%nSide
          MESH%Fault%Face(iFace,:,:) = NewFace(iFace,:,:)
          MESH%Fault%geoNormals(:,iFace) = NewNormal(:,iFace)
          MESH%Fault%geoTangent1(:,iFace) = NewTangent1(:,iFace)
          MESH%Fault%geoTangent2(:,iFace) = NewTangent2(:,iFace)
        ENDDO
        !
        DEALLOCATE(NewFace)
        DEALLOCATE(NewNormal,NewTangent1,NewTangent2)
        !
      ENDIF ! nLocalFaultedges present in local domain
      !
    ENDIF ! DR
    
   
    logInfo(*) 'Destruct old (total) mesh. '

    CALL destruct_mesh_level0_1(                                            & ! deallocating array for elem
         MESH = MESH                                                        ) ! deallocating array for elem
    
    !CALL destruct_ncmesh_level0_1(                                          & ! deallocating array for elem
    !     MESH = MESH                                                        ) ! deallocating array for elem

    CALL destruct_mesh_level0_2(                                            & ! deallocating array for node
         MESH = MESH                                                        ) ! deallocating array for node

    DEALLOCATE(MESH%ELEM%MPINumber)
    DEALLOCATE(MESH%ELEM%MPINumber_DR) ! * Dynamic Rupture
    !DEALLOCATE(MESH%ELEM%MPI_NCNumber)
    DEALLOCATE(MESH%LocalElemType)
    
    IF(ASSOCIATED(MESH%ElemGlobal2Local)) THEN
        DEALLOCATE(MESH%ElemGlobal2Local) 
        DEALLOCATE(MESH%ElemLocal2Global)
        DEALLOCATE(MESH%NodeGlobal2Local)
        DEALLOCATE(MESH%NodeLocal2Global) 
    ENDIF

    MESH%nNode                  = NewMesh%nNode
    MESH%nElem                  = NewMesh%nElem
    MESH%MaxElementsOnVertex    = NewMesh%MaxElementsOnVertex
    MESH%Dimension              = NewMesh%Dimension
    MESH%ADERDG3D               = NewMesh%ADERDG3D
    MESH%nVertexMax             = NewMesh%nVertexMax
    !MESH%nNonConformingEdges    = NewMesh%nNonConformingEdges
    
    CALL allocate_mesh_level0_1(                                            & ! allocating array for elem
         IO   = IO                                                        , & ! allocating array for elem
         MESH = MESH                                                        ) ! allocating array for elem
    !CALL allocate_ncmesh_level0_1(                                          & ! allocating array for nc elem
    !     IO   = IO                                                        , & ! allocating array for nc elem
    !     MESH = MESH                                                      , & ! allocating array for nc elem
    !     DISC = DISC                                                        )
    CALL allocate_mesh_level0_2(                                            & ! allocating array for node
         IO   = IO                                                        , & ! allocating array for node
         EQN  = EQN                                                       , & ! allocating array for node
         MESH = MESH                                                        ) ! allocating array for node

    ALLOCATE(MESH%ELEM%MPINumber(MESH%nVertexMax,MESH%nElem))
    ALLOCATE(MESH%ELEM%MPINumber_DR(MESH%nVertexMax,MESH%nElem)) ! * Dynamic Rupture
    !ALLOCATE(MESH%ELEM%MPI_NCNumber(MESH%nNonConformingEdges,DISC%Galerkin%nBndGP))
    ALLOCATE(MESH%LocalElemType(MESH%nElem))

    logInfo(*) 'Copying node data to new mesh '


    ! Copy new node coordinates to original data structure
    DO iNode = 1, MESH%nNode
       MESH%VRTX%xyNode(:,iNode)                 = NewMESH%VRTX%xyNode(:,iNode)
       ! 
       MESH%VRTX%NrOfElementsConnected(iNode)    = NewMESH%VRTX%NrOfElementsConnected(iNode)
       MESH%VRTX%Element(iNode,:)                = NewMESH%VRTX%Element(iNode,:)
    ENDDO

    logInfo(*) 'Copying element data to new mesh '

    ! Copy new element information to original data structure
    DO iElem = 1, MESH%nElem
       MESH%ELEM%Vertex(:,iElem)             = NewMESH%ELEM%Vertex(:,iElem)
       MESH%ELEM%Reference(:,iElem)          = NewMESH%ELEM%Reference(:,iElem)
       MESH%ELEM%MPIReference(:,iElem)       = NewMESH%ELEM%MPIReference(:,iElem)
       MESH%ELEM%MPINumber(:,iElem)          = NewMESH%ELEM%MPINumber(:,iElem)
       MESH%ELEM%BoundaryToObject(:,iElem)   = NewMESH%ELEM%BoundaryToObject(:,iElem)
       MESH%ELEM%SideNeighbor(:,iElem)       = NewMESH%ELEM%SideNeighbor(:,iElem)
       MESH%ELEM%LocalNeighborSide(:,iElem)  = NewMESH%ELEM%LocalNeighborSide(:,iElem)
       MESH%ELEM%LocalNeighborVrtx(:,iElem)  = NewMESH%ELEM%LocalNeighborVrtx(:,iElem)
       MESH%ELEM%Volume(iElem)               = NewMESH%ELEM%Volume(iElem)
       MESH%ELEM%xyBary(:,iElem)             = NewMESH%ELEM%xyBary(:,iElem)
       MESH%LocalElemType(iElem)             = NewMESH%LocalElemType(iElem)
    ENDDO
    
    ! * Dynamic Rupture
    IF (EQN%DR == 1) THEN
       DO iElem = 1, MESH%nElem
          MESH%ELEM%MPINumber_DR(:,iElem)    = NewMESH%ELEM%MPINumber_DR(:,iElem)
       ENDDO
    ENDIF
    
    !DO iNCboundary = 1, NewMesh%nNonConformingEdges 
    !        MESH%Elem%MPI_NCNumber(iNCboundary,:)        = NewMESH%ELEM%MPI_NCNumber(iNCboundary,:) 
    !        MESH%ELEM%NC_BoundaryToObject(iNCboundary,:) = NewMESH%ELEM%NC_BoundaryToObject(iNCboundary,:) 
    !        MESH%ELEM%ncBndNeighbor(iNCboundary,:)       = NewMesh%ELEM%ncBndNeighbor(iNCboundary,:) 
    !        Mesh%ELEM%NCB_IndexList(iNCboundary,:)       = NewMesh%ELEM%NCB_IndexList(iNCboundary,:)  
    !        Mesh%ELEM%ncBndGaussP(iNCboundary,:,:)       = NewMesh%ELEM%ncBndGaussP(iNCboundary,:,:)       
    !ENDDO
    
    ! Update the MPI boundary information
    DO iDomain = 1, BND%NoMPIDomains
       DO i = 1, BND%ObjMPI(iDomain)%nElem
         BND%ObjMPI(iDomain)%DomainElements(i) = invElementList(BND%ObjMPI(iDomain)%DomainElements(i)) 
       ENDDO
       !DO i = 1, BND%ObjMPI(iDomain)%nNCSendElem
       !  BND%ObjMPI(iDomain)%NC_SendElements(i) = inv_NC_ElementList(BND%ObjMPI(iDomain)%NC_SendElements(i))
       !ENDDO
       ! WRITE(IO%UNIT%stdOut,*) '|   Updated bnd. elements with CPU     ', BND%ObjMPI(iDomain)%CPU, ' : ', &
       !                              BND%ObjMPI(iDomain)%DomainElements(:)
    ENDDO
    
    ! * Dynamic Rupture
    IF (EQN%DR == 1) THEN
        BND%ObjMPI(:)%nFault_MPI = BND%ObjMPI(:)%nFault_MPI/2 ! neighbors were also counted
        DO iDomain = 1, BND%NoMPIDomains
           DO i = 1, BND%ObjMPI(iDomain)%nFault_MPI
              BND%ObjMPI(iDomain)%Domain_Fault_Elem(i) = invElementList(BND%ObjMPI(iDomain)%Domain_Fault_Elem(i)) 
           ENDDO
        ENDDO
    ENDIF
    
    !
    CLOSE(IO%UNIT%FileIn)
    !
    CALL destruct_mesh_level0_1(                                            & ! deallocating array for elem
         MESH = NewMESH                                                     ) ! deallocating array for elem

    !CALL destruct_ncmesh_level0_1(                                          & ! deallocating array for elem
    !     MESH = NewMESH                                                     ) ! deallocating array for elem

    CALL destruct_mesh_level0_2(                                            & ! deallocating array for node
         MESH = NewMESH                                                     ) ! deallocating array for node
    !
    DEALLOCATE(NewMESH%ELEM%MPINumber)
    DEALLOCATE(NewMESH%ELEM%MPINumber_DR)
    !DEALLOCATE(NewMESH%ELEM%MPI_NCNumber)
    DEALLOCATE( CPULocal2Global ) 
    DEALLOCATE( CPUGlobal2Local ) 
    !
    logInfo(*) 'MPIExtractMesh done by CPU ', MPI%myrank
    !
    ! Recount number of tets and hexas
    ! Read element list to evaluate individual element types
    MESH%nElem_Tet = 0
    MESH%nElem_Hex = 0
    DO iElem = 1, MESH%nElem
      SELECT CASE(MESH%LocalElemType(iElem))
      CASE(4)
          MESH%nElem_Tet = MESH%nElem_Tet  + 1
          IF(MESH%nElem_Tet.EQ.1)THEN
            DEALLOCATE(MESH%LocalElemIndex_Tet)
            ALLOCATE(MESH%LocalElemIndex_Tet(MESH%nElem))
            MESH%LocalElemIndex_Tet = 0
          ENDIF
          MESH%LocalElemIndex_Tet(iElem) = MESH%nElem_Tet 
      CASE(6)
          MESH%nElem_Hex = MESH%nElem_Hex + 1
          IF(MESH%nElem_Hex.EQ.1)THEN
            DEALLOCATE(MESH%LocalElemIndex_Hex)
            ALLOCATE(MESH%LocalElemIndex_Hex(MESH%nElem))
            MESH%LocalElemIndex_Hex = 0
          ENDIF
          MESH%LocalElemIndex_Hex(iElem) = MESH%nElem_Hex
      CASE DEFAULT
          logError(*) 'Mixed mesh does not contain only tets and hexas!'
          STOP 
      END SELECT
    ENDDO
#endif

  END SUBROUTINE MPIExtractMesh
  subroutine sort(a,a2)
    ! sorts array a using heap-sort algorithm
    ! and arrays a2 correspondingly
    implicit none
    integer, parameter:: ik=kind(1)
    integer(ik) :: a(:),a2(size(a)),b,d,c,s,N,t,it
    N=size(a);if(N<2)return  
    d=N;s=N/2+1;
    do
       if(s > 1)then
          s=s-1;t=a(s);it=a2(s);
       else
          t=a(d);it=a2(d);a(d)=a(1);a2(d)=a2(1);d=d-1
          if(d == 1)then
             a(1)=t;a2(1)=it; return
          endif
       endif
       b = s; c = s + s
       do
          if(c > d)exit
          if(c < d)then
             if(a(c) < a(c+1)) c = c + 1
          endif
          if(t < a(c))then
             a(b)=a(c);a2(b)=a2(c);b=c;c=c+c
          else
             c = d + 1
          endif
       enddo
       a(b)=t;a2(b)=it 
    enddo
  end subroutine sort

END MODULE MPIExtractMesh_mod
