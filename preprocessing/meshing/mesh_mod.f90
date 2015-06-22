!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Gilbert B. Brietzke (Gilbert.Brietzke AT lrz.de, http://www.geophysik.uni-muenchen.de/Members/brietzke)
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
module mesh_mod
  implicit none

  integer,parameter, private :: dp=kind(1.d0)

  type stringpair
     character(:), allocatable :: str,val
  end type stringpair

  type sectionlabel
     character(:), allocatable :: label,version,endlabel
   contains
     procedure :: read  => read_label
     procedure :: end   => read_endlabel
     procedure :: write => write_label
  end type sectionlabel

  type header
     type(sectionlabel):: seclab
     type(stringpair) :: file
     type(stringpair) :: soft
     type(stringpair) :: vers
     type(stringpair) :: date
     integer :: NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL
   contains
     procedure :: read => read_header
     procedure :: write => write_header
     procedure :: print => print_header
  end type header

  type nodetype
     integer :: numb
     real(dp) :: crds(3)
  end type nodetype

  type celltype
     integer :: numb,dims,nvert
     integer  :: ngrp
     integer, allocatable :: nodes(:)
  end type celltype

  type grouptype
     integer :: id
     integer :: ncells
     integer :: nprop
     integer :: nflags
     character(len=100) :: name
     integer, dimension(:), allocatable :: cell_id
  end type grouptype

  type bctype
     character(len=100) :: label
     integer :: n1,n,n2,n3
     integer, dimension(:), allocatable :: id
     integer, dimension(:), allocatable :: dim
     integer, dimension(:), allocatable :: side
  end type bctype

  type meshtype
     type(header) :: head
     integer,pointer :: NUMNP, NELEM, NGRPS, NBSETS, NDFCD, NDFVL
     type(nodetype) , allocatable :: nodes(:)
     type(celltype) , allocatable :: cells(:)
     type(grouptype), allocatable :: groups(:)
     type(bctype)   , allocatable :: bcs(:)
     character(len=100) :: fileinfo
   contains
     procedure :: read_gambit  => mesh_read_gambit
     procedure :: write_gambit => mesh_write_gambit
  end type meshtype

  private int2str

  type(sectionlabel):: seclab

contains

  subroutine mesh_read_gambit(mesh,meshfile)
    implicit none
    class(meshtype),target :: mesh
    character(len=*), intent(in) ::meshfile
    type(nodetype)   :: zeronode
    type(celltype)   :: zerocell
    type(grouptype)  :: zerogroup
    integer :: i,j,ic=0,ibc=0,iloop,iun,iostat
    write(*,'(a)')'opening meshfile: '//trim(adjustl(meshfile))//' ... : '
    open(newunit=iun,file=trim(adjustl(meshfile)),status='old',action='read')
    call mesh%head%read(iun)
    call mesh%head%print
    mesh%NUMNP => mesh%head%NUMNP; mesh%NELEM => mesh%head%NELEM
    mesh%NGRPS => mesh%head%NGRPS; mesh%NBSETS=> mesh%head%NBSETS
    mesh%NDFCD => mesh%head%NDFCD; mesh%NDFVL => mesh%head%NDFVL
    zeronode=nodetype(0,[0.d0,0.d0,0.d0]);
    zerocell=celltype(0,0,0,0,null())
    zerogroup=grouptype(0,0,0,0,'',null())
    allocate(mesh%cells(mesh%nelem));mesh%cells=zerocell
    allocate(mesh%groups(mesh%ngrps));mesh%groups=zerogroup
    allocate(mesh%bcs(mesh%nbsets));
    allocate(mesh%nodes(mesh%numnp));mesh%nodes=zeronode
    write(*,'(A)')'now looping sections:'
    do iloop = 1,10000
       call seclab%read(iun)
       if(seclab%label(1:17)=='NODAL COORDINATES')  call mesh_read_nodes(mesh,iun) 
       if(seclab%label(1:14)=='ELEMENTS/CELLS')     call mesh_read_elements(mesh,iun)
       if(seclab%label(1:13)=='ELEMENT GROUP')      call mesh_read_groups(mesh,iun)
       if(seclab%label(1:19)=='BOUNDARY CONDITIONS')call mesh_read_bcs(mesh,iun)
       cycle
    enddo
    close(iun)
  end subroutine mesh_read_gambit
  !-------------------------------------------------------------------------------------
  subroutine mesh_write_gambit(mesh,meshfile,filename)
    implicit none
    class(meshtype) :: mesh
    character(len=*), intent(in) ::meshfile
    character(len=*), intent(in) ::filename
    character(len=100) :: myform
    integer :: i,j,iun,ibc
    write(*,*)' opening meshfile: ',trim(adjustl(meshfile))
    open(newunit=iun,file=trim(adjustl(meshfile)),action='write')
    call mesh%head%write(iun)
    call mesh_write_nodes(mesh,iun)
    call mesh_write_elements(mesh,iun)
    call mesh_write_groups(mesh,iun)
    call mesh_write_bcs(mesh,iun)
    close(iun)
  end subroutine mesh_write_gambit

  subroutine mesh_read_nodes(mesh,iun)
    implicit none
    type(meshtype) :: mesh
    integer :: iun,i
    write(*,'("    ",a,a,i10,a)',advance='no')trim(seclab%label),',   reading: ',mesh%numnp,' mesh-nodes'
    do i = 1,mesh%numnp
       read(iun,*) mesh%nodes(i)%numb,mesh%nodes(i)%crds(1:mesh%ndfcd)
    enddo
    write(*,'(a,a)')', done! ',seclab%end(iun)       
  end subroutine mesh_read_nodes

  subroutine mesh_read_elements(mesh,iun)
    implicit none
    type(meshtype) :: mesh
    integer :: iun,i,iostat
    character(len=100) :: info
    write(*,'("    ",a,a,i10,a)',advance='no')trim(seclab%label),',      reading: ',mesh%nelem,' mesh-elements'
    do i = 1,mesh%nelem
       read(iun,'(A100)')info
       read(info,*)mesh%cells(i)%numb,mesh%cells(i)%dims,mesh%cells(i)%nvert
       allocate(mesh%cells(i)%nodes(mesh%cells(i)%nvert));mesh%cells(i)%nodes=0;             
       read(info,*,iostat=iostat) mesh%cells(i)%numb,mesh%cells(i)%dims, & 
            &                         mesh%cells(i)%nvert,mesh%cells(i)%nodes
       if(.not.(is_iostat_end(iostat).or.is_iostat_eor(iostat)))cycle
       read(info,*) mesh%cells(i)%numb,mesh%cells(i)%dims,mesh%cells(i)%nvert, & 
            &       mesh%cells(i)%nodes(1:min(mesh%cells(i)%nvert,7))
       if(mesh%cells(i)%nvert>7)read(iun,*) mesh%cells(i)%nodes(8:mesh%cells(i)%nvert)
    enddo
    write(*,'(a,a)')', done! ',seclab%end(iun)
  end subroutine mesh_read_elements
  
  subroutine mesh_read_groups(mesh,iun)
    implicit none
    type(meshtype) :: mesh
    integer :: iun,j,iostat
    integer, save :: ic=0
    character(len=100) :: char1
    ic=ic+1
    if(ic<=mesh%ngrps)then
       write(*,'("    ",a,a,i5,a,i5,a)',advance="no")trim(seclab%label),',       reading: ',ic,' of ',mesh%ngrps,' groups'
       read(iun,*)char1,mesh%groups(ic)%id   ,char1,mesh%groups(ic)%ncells, & 
            &     char1,mesh%groups(ic)%nprop,char1,mesh%groups(ic)%nflags
       read(iun,'(A100)')mesh%groups(ic)%name
       read(iun,*)
       allocate(mesh%groups(ic)%cell_id(mesh%groups(ic)%ncells))
       do j = 1,mesh%groups(ic)%ncells,10
          read(iun,*)mesh%groups(ic)%cell_id(j:min(j+9,mesh%groups(ic)%ncells))
       enddo
       mesh%cells(mesh%groups(ic)%cell_id)%ngrp = mesh%groups(ic)%id;
       write(*,'(a,a)')', done! ',seclab%end(iun)       
    endif
  end subroutine mesh_read_groups

  subroutine mesh_read_bcs(mesh,iun)
    implicit none
    type(meshtype) :: mesh
    integer :: iun,j,iostat
    integer, save :: ibc=0
    ibc=ibc+1
    if(ibc<=mesh%nbsets)then
       write(*,'("    ",a,a,i5,a,i5,a)',advance="no")trim(seclab%label),', reading: ',ibc,' of ',mesh%nbsets,' BCs'
       read(iun,*,iostat=iostat)mesh%bcs(ibc)%label,mesh%bcs(ibc)%n1,mesh%bcs(ibc)%n,mesh%bcs(ibc)%n2,mesh%bcs(ibc)%n3
       if(is_iostat_err(iostat)) stop'error during read!'
       select case (mesh%bcs(ibc)%n3)
       case(6)
          allocate(mesh%bcs(ibc)%id(mesh%bcs(ibc)%n))
          allocate(mesh%bcs(ibc)%dim(mesh%bcs(ibc)%n))
          allocate(mesh%bcs(ibc)%side(mesh%bcs(ibc)%n))
          do j = 1,mesh%bcs(ibc)%n
             read(iun,*,iostat=iostat)mesh%bcs(ibc)%id(j),mesh%bcs(ibc)%dim(j),mesh%bcs(ibc)%side(j)
             if(is_iostat_err(iostat)) stop'error during read!'
          enddo
       case(24)
          allocate(mesh%bcs(ibc)%id(mesh%bcs(ibc)%n))
          do j = 1,mesh%bcs(ibc)%n
             read(iun,*,iostat=iostat)mesh%bcs(ibc)%id(j)
             if(is_iostat_err(iostat)) stop'error during read!'
          enddo
       end select
       write(*,'(a,a)')', done! ',seclab%end(iun)       
    endif
  end subroutine mesh_read_bcs

  subroutine mesh_write_nodes(mesh,iun)
    implicit none
    type(meshtype) :: mesh
    integer :: iun,j
    write(*,*)' writing: ',mesh%numnp,' mesh-nodes'
    write(iun,'(A)')"   NODAL COORDINATES "//trim(mesh%head%vers%val)
    do j = 1,mesh%numnp
       write(iun,'(i10,3es20.11E3)')mesh%nodes(j)%numb,mesh%nodes(j)%crds(1:mesh%ndfcd)
    enddo
    write(iun,'("ENDOFSECTION")')
    write(*,*)'    ...   ',mesh%numnp,' nodes done'
  end subroutine mesh_write_nodes

  subroutine mesh_write_elements(mesh,iun)
    implicit none
    type(meshtype) :: mesh
    integer :: iun,j
    character(len=100) :: myform
    write(*,*)' writing: ',mesh%nelem,' mesh-elements'
    write(iun,'(A)')"      ELEMENTS/CELLS "//trim(mesh%head%vers%val)
    do j = 1,mesh%nelem
       myform='(i8,i3,i3," ",'//int2str(mesh%cells(j)%nvert,1,'(i1)')//'i8)'
       write(iun,trim(adjustl(myform)))mesh%cells(j)%numb,mesh%cells(j)%dims,mesh%cells(j)%nvert,mesh%cells(j)%nodes
    enddo
    write(iun,'("ENDOFSECTION")')
    write(*,*)'   ...   ',mesh%nelem,' elements done'
  end subroutine mesh_write_elements

  subroutine mesh_write_groups(mesh,iun)
    type(meshtype) :: mesh
    integer :: iun,i,j
    if(mesh%ngrps>0)then
       write(*,*)' writing: ',mesh%ngrps,' groups'
       do i =1,mesh%ngrps
          write(iun,'(A)')"       ELEMENT GROUP"//trim(mesh%head%vers%val)
          write(iun,'("GROUP:",i11," ELEMENTS:",i11," MATERIAL:",i11," NFLAGS:",i11)') & 
               & mesh%groups(i)%id,mesh%groups(i)%ncells,mesh%groups(i)%nprop,mesh%groups(i)%nflags
          write(iun,'(a40)')mesh%groups(i)%name
          write(iun,*)'       0'
          do j = 1,mesh%groups(i)%ncells,10
             write(iun,'(10i8)')mesh%groups(i)%cell_id(j:min(j+9,mesh%groups(i)%ncells))
          enddo
          write(iun,'("ENDOFSECTION")')
       enddo
       write(*,*)'   ...   ',mesh%ngrps,' groups done'
    endif
  end subroutine mesh_write_groups

  subroutine mesh_write_bcs(mesh,iun)
    implicit none
    type(meshtype) :: mesh
    integer :: iun,ibc,j
    if(mesh%nbsets>0)then
       do ibc =1,mesh%nbsets
          write(iun,'(A)')" BOUNDARY CONDITIONS "//trim(mesh%head%vers%val)
          write(*,'("    ",a,i5,a,i5,a)',advance="no")'BOUNDARY CONDITIONS, writing: ',ibc,' of ',mesh%nbsets,' BCs'
          write(iun,'("                        ",A8,4i8)')trim(adjustl(mesh%bcs(ibc)%label)),mesh%bcs(ibc)%n1, & 
               &                                       mesh%bcs(ibc)%n,mesh%bcs(ibc)%n2,mesh%bcs(ibc)%n3
          select case (mesh%bcs(ibc)%n3)
          case(6)
             do j = 1,mesh%bcs(ibc)%n
                write(iun,'(I10,I5,I5)')mesh%bcs(ibc)%id(j),mesh%bcs(ibc)%dim(j),mesh%bcs(ibc)%side(j)
             enddo
          case(24)
             do j = 1,mesh%bcs(ibc)%n
                write(iun,'(I10)')mesh%bcs(ibc)%id(j)
             enddo
          end select
          write(iun,'("ENDOFSECTION")')
          write(*,'(a)')', done! '
       end do
    endif
  end subroutine mesh_write_bcs

  subroutine read_label(lab,iun)
    class(sectionlabel) :: lab 
    integer :: iun
    character(len=20) :: dum1,dum2
    integer :: iostat
    read(iun,'(A20,A)',iostat=iostat)dum1,dum2; 
    if(is_iostat_end(iostat))return
    lab%label=trim(adjustl(dum1));lab%version=trim(adjustl(dum2))
  end subroutine read_label

  subroutine write_label(lab,iun)
    class(sectionlabel) :: lab 
    integer :: iun
    write(iun,'(A20,A,A)')adjustr(lab%label),' ',lab%version
  end subroutine write_label

  function read_endlabel(lab,iun)result(endlabel)
    class(sectionlabel) :: lab 
    character(len=12):: dum1,endlabel
    integer :: iun
    read(iun,'(A12)')dum1;
    endlabel='ENDOFSECTION'
    if(trim(adjustl(dum1))=='ENDOFSECTION')then
       lab%endlabel=endlabel
       return
    end if
    endlabel='INCONSISTENT'
  end function read_endlabel
  
  subroutine read_header(head,iun)
    class(header) :: head
    integer :: iun
    character(len=100) :: dum1,dum2,dum3,dum4
    call head%seclab%read(iun)
    read(iun,'(A100)')dum1; read(iun,'(A)')dum2; 
    head%file%str=trim(adjustl(dum1));head%file%val=trim(adjustl(dum2))
    read(iun,*)dum1,dum2,dum3,dum4
    head%soft%str=trim(adjustl(dum1));head%soft%val=trim(adjustl(dum2))
    head%vers%str=trim(adjustl(dum3));head%vers%val=trim(adjustl(dum4))
    read(iun,'(A100)')dum2;
    head%date%str='date';head%date%val=trim(adjustl(dum2))
    read(iun,*)
    read(iun,*)head%NUMNP,head%NELEM,head%NGRPS,head%NBSETS,head%NDFCD,head%NDFVL
    write(*,'(A)')'reading of section: '//head%seclab%label//': check: '//head%seclab%end(iun)
  end subroutine read_header

  subroutine print_header(head)
    class(header) :: head
    call head%write()
  end subroutine print_header

  subroutine write_header(head,unit)
    use iso_fortran_env
    class(header) :: head
    integer, optional :: unit
    integer :: iun
    iun=output_unit
    if(present(unit))iun=unit
    call head%seclab%write(iun)
    write(iun,'(A)')adjustl(head%file%str)
    write(iun,'(A)')adjustl(head%file%val); 
    write(iun,'(A8,A,A21,A,A12,A,A)') &
         & adjustr(head%soft%str),' ',adjustl(head%soft%val),' ',&
         & adjustr(head%vers%str),'  ',adjustl(head%vers%val); 
    write(iun,'(" ",A)')adjustl(head%date%val);
    write(iun,'("     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL")')
    write(iun,'(6i10)')head%NUMNP,head%NELEM,head%NGRPS,head%NBSETS,head%NDFCD,head%NDFVL
    write(iun,'(A)')'ENDOFSECTION'
  end subroutine write_header

  function int2str(n,nc,myform)
    implicit none
    integer :: n,nc
    character(len=nc)::int2str
    character(len=*)::myform
    write(int2str,myform)n
  end function int2str

  logical function is_iostat_err(iostat)
    integer :: iostat
    if(is_iostat_end(iostat))is_iostat_err=.false.
    if(is_iostat_eor(iostat))is_iostat_err=.false.
  end function is_iostat_err

end module mesh_mod


!!$program test_mesh
!!$
!!$  use mesh_mod
!!$
!!$  type(meshtype) :: mesh
!!$  character(len=512)::arg
!!$  character(:),allocatable :: infile,outfile
!!$  integer :: n
!!$
!!$  call get_command_argument(1,arg)
!!$  infile=trim(adjustl(arg));n=len(infile)
!!$  outfile=infile(1:n-4)//'_new'//infile(n-3:n)
!!$  call mesh%read_gambit(infile)
!!$  call mesh%write_gambit(outfile,'newinfo')
!!$
!!$end program test_mesh

!!$program test
!!$  
!!$  use testo
!!$  use mesh_mod
!!$
!!$  type(header) :: head
!!$
!!$  open(newunit=iun,file='scec_box40x20x30_small.neu')  
!!$  call head%read(iun)
!!$  close(iun)
!!$  call head%print()
!!$  open(newunit=iun,file='scec_box40x20x30_small_new.neu')
!!$  call head%write(iun)
!!$  close(iun)
!!$
!!$end program test
