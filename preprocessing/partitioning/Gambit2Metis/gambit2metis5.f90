!>
!! @file
!! This file is part of SeisSol.
!!
!! @section LICENSE
!! Copyright (c) SeisSol Group
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

program gambit2metis5
  
  use mesh_mod
  implicit none

  character(:),allocatable :: infile,metisfile,filename,file
  character(len=512) :: arg
  type(meshtype)     :: mesh
  integer            :: n,np

  call get_command_argument(1,arg)
  infile=trim(adjustl(arg));n=len(infile)
  call get_command_argument(2,arg)
  read(arg,*)np
  call mesh%read_gambit(infile)
  filename=get_filename(infile);filename=filename(1:len(filename)-4)
  call do_partitioning_external(mesh,np,filename)

contains

  subroutine do_partitioning_external(mesh,np,filename)
    use ifport, only:system
    implicit none
    type(meshtype) :: mesh
    character(:),allocatable :: metisfile
    character(len=*) :: filename
    integer :: np
    integer :: i,iun
    metisfile=filename//'.met'  
    open(newunit=iun,file=metisfile,action='write')
    write(iun,*)mesh%nelem,1
    do i = 1,mesh%nelem
       write(iun,*)mesh%cells(i)%nodes
    end do
    close(iun)
    i=system('m2gmetis '//filename//'.met '//filename//'.met.dgraph ')
    i=system('gpmetis '//filename//'.met.dgraph '//int2str(np))
    !clean up of files
    i=system('rm '//filename//'.met')
    i=system('rm '//filename//'.met.dgraph')
    i=system('mv '//filename//'.met.dgraph.part.'//int2str(np)//' '//filename//'.met.epart.'//int2str(np))
  end subroutine do_partitioning_external

  function get_filename(filestring,mext)
    implicit none
    character(len=*) :: filestring
    character(:), allocatable :: get_filename
    integer, optional :: mext
    integer :: i
    do i=len(filestring),1,-1
       if(filestring(i:i)=='/')exit
    end do
    get_filename=trim(adjustl(filestring(i+1:len(filestring))))
  end function get_filename

  function int2str(i)
    integer :: i
    character(:), allocatable :: int2str
    character(len=20):: dum
    write(dum,*)i
    int2str=trim(adjustl(dum))
  end function int2str

!!$  subroutine do_partitioning_api
!!$
!!$    use metis_mod
!!$    USE, INTRINSIC :: ISO_C_BINDING
!!$
!!$    integer :: objval,icount,numflag,ncommon
!!$    integer, allocatable, dimension(:) :: cellnodes,cellindices,e_mpiid,n_mpiid,eptr,eind
!!$    integer, allocatable, dimension(:),target ::   xadj,adjncy
!!$    !integer,pointer     ::   pxadj(:)=> NULL(),padjncy(:)=> NULL()
!!$    integer, pointer :: vwgt(:)  => NULL()
!!$    integer, pointer :: vsize(:) => NULL()
!!$    integer, pointer :: mopts(:) => NULL()
!!$    real(kind(1.d0)), pointer :: tpwgts(:) => NULL()
!!$    integer::METIS_PartMeshDual
!!$  
!!$    allocate(cellnodes(1:sum(mesh%cells%nvert)))
!!$    allocate(cellindices(1:mesh%nelem+1))
!!$    icount=1
!!$    cellindices(1)=1
!!$    do j = 1,mesh%nelem
!!$       cellindices(j+1)=cellindices(j)+mesh%cells(j)%nvert
!!$       do i = 1,mesh%cells(j)%nvert
!!$          cellnodes(icount)=mesh%cells(j)%nodes(i)
!!$          icount =icount+1
!!$       enddo
!!$       if(j<10)then
!!$          write(*,*)mesh%cells(j)%nodes,cellindices(j),cellindices(j+1)-1
!!$          write(*,*)cellnodes(cellindices(j):cellindices(j+1)-1)
!!$       endif
!!$    enddo
!!$    
!!$    allocate(eind(0:sum(mesh%cells%nvert)-1))
!!$    allocate(eptr(0:mesh%nelem))
!!$    
!!$    eind=cellnodes
!!$    eptr=cellindices
!!$    
!!$    allocate(e_mpiid(mesh%nelem),n_mpiid(mesh%numnp))
!!$    write(*,*)'   partitioning mesh into ',np,' parts'
!!$    write(*,*)'   ...  calling metis_partmeshdual ...'
!!$    call METIS_SetDefaultOptions(metis_options);
!!$    call mymetis_default_options();mopts=>null()!metis_options
!!$    ncommon=1
!!$    numflag=1
!!$    i=METIS_PartMeshDual(mesh%nelem,mesh%numnp,eptr,eind, & 
!!$         & vwgt,vsize,ncommon,np,tpwgts,metis_options,objval,e_mpiid,n_mpiid)
!!$    write(*,*)'   ...          metis_partmeshdual done!'
!!$    write(*,*)'objval',objval,i
!!$    
!!$    file=filename//'.met.epart2.'//int2str(np)
!!$    open(newunit=iun,file=file)
!!$    do i = 1,size(e_mpiid,1)
!!$       write(iun,'(I1)')e_mpiid(i)
!!$    enddo;close(iun);
!!$
!!$  file=filename//'.met.epart3.'//int2str(np)
!!$  open(newunit=iun,file=file)  
!!$  do i = 1,size(n_mpiid,1)
!!$     write(iun,'(I1)')n_mpiid(i)
!!$  enddo;close(iun);
!!$
!!$  end subroutine do_partitioning_api

end program gambit2metis5

