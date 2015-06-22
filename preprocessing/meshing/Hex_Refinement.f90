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
PROGRAM Hex_refinement
    !--------------------------------------------------------------------------
    IMPLICIT NONE
    !--------------------------------------------------------------------------
    ! Type declaration
    INTEGER                       :: uid1, uid2
    INTEGER                       :: c, ix, iy, iz, nline, nrest, ng1, ng2
    INTEGER                       :: clocal(32), ind(8)
    INTEGER                       :: nVert, nElem, nl, nv, ntr
    INTEGER                       :: nx, nx2, ny, ny2, nz, nz2
    REAL                          :: xmin, xmax, ymin, ymax, zmin, zmax, zb
    REAL                          :: dx, dy, dz
    REAL,    ALLOCATABLE          :: X(:,:)
    INTEGER, ALLOCATABLE          :: HEX(:,:)
    INTEGER, ALLOCATABLE          :: member(:)
    CHARACTER (LEN=350)           :: InFilename
    CHARACTER (LEN=350)           :: OutFilename
    CHARACTER (LEN=70)            :: StringLine            
    !--------------------------------------------------------------------------
    !
    WRITE(*,*) ' | -------------------------------------------------- |' 
    WRITE(*,*) ' |                  HEX_REFINEMENT                    |'
    WRITE(*,*) ' | -------------------------------------------------- |' 
    !
    WRITE(*,*) '  ' 
    WRITE(*,*) '   Enter the block (brick) dimensions!  '
    WRITE(*,*) '  ' 
    WRITE(*,*) '   x_min   x_max  '
    READ (*,*) xmin, xmax
    WRITE(*,*) '   y_min   y_max  '
    READ (*,*) ymin, ymax
    WRITE(*,*) '   z_min   z_max  '
    READ (*,*) zmin, zmax
    WRITE(*,*) '   Enter depth of material boundary under which mesh is coarsened!  '
    WRITE(*,*) '  ' 
    READ (*,*) zb
    WRITE(*,*) '   Enter numbers of intervals in x,y,z-directions for the coarse mesh!  '
    WRITE(*,*) '  ' 
    READ (*,*) nx2, ny2, nz2
    WRITE(*,*) '   Enter numbers of intervals in z-directions for the fine mesh!  '
    WRITE(*,*) '  ' 
    READ (*,*) nz

    nx = nx2*3
    ny = ny2*3
    
    dx = (xmax-xmin)/nx
    dy = (ymax-ymin)/ny
    dz = (zmax-zb)  /nz

    !OutFilename = TRIM(InFilename)//'.met'
    !InFilename  = TRIM(InFilename)//'.neu'
    WRITE(*,*) '   Give output file name '
    READ (*,*) OutFilename
    
    !nVert is number of vertices, i.e. of fine mesh, course mesh and transition layer
    nVert  = (nx+1)*(ny+1)*(nz+1) + (nx2+1)*(ny2+1)*(nz2) + (nx2*ny2*8+nx2*2+ny2*2)
    !nVert is number of elements, i.e. of fine mesh, course mesh and transition layer
    nElem  = (nx*ny*nz)           + (nx2*ny2*(nz2-1))     + (nx2*ny2*13)
    
    
    ALLOCATE( X(nVert,3) )
    ALLOCATE( HEX(nElem,8) )

    !--------------------------------------------------------------------------
    ! Compute vertices of fine mesh 
    
    c = 0

    DO iz = 1, nz+1
      DO iy = 1, ny+1
        DO ix = 1, nx+1
             
           c = c+1

           X(c,1) = xmin + (ix-1)*dx
           X(c,2) = ymin + (iy-1)*dy
           X(c,3) = zmax - (iz-1)*dz
        
        ENDDO
      ENDDO
    ENDDO
         
    !--------------------------------------------------------------------------
    ! Compute vertices of transition mesh 
    
    dx = 3*dx
    dy = 3*dy
    dz = (zb-zmin)/nz2

    DO iy = 1, ny2
      DO ix = 1, nx2
            
           ! four inner vertices
           c = c+1

           X(c,1) = xmin + (ix-1)*dx + dx/3
           X(c,2) = ymin + (iy-1)*dy + dy/3
           X(c,3) = zb               - dz/3

           c = c+1

           X(c,1) = xmin + (ix-1)*dx + dx*2/3
           X(c,2) = ymin + (iy-1)*dy + dy  /3
           X(c,3) = zb               - dz  /3
 
           c = c+1

           X(c,1) = xmin + (ix-1)*dx + dx*2/3
           X(c,2) = ymin + (iy-1)*dy + dy*2/3
           X(c,3) = zb               - dz  /3

           c = c+1

           X(c,1) = xmin + (ix-1)*dx + dx/3
           X(c,2) = ymin + (iy-1)*dy + dy*2/3
           X(c,3) = zb               - dz  /3
  
           ! four boundary vertices
           c = c+1

           X(c,1) = xmin + (ix-1)*dx
           X(c,2) = ymin + (iy-1)*dy + dy*2/3
           X(c,3) = zb               - dz/2

           c = c+1

           X(c,1) = xmin + (ix-1)*dx 
           X(c,2) = ymin + (iy-1)*dy + dy/3
           X(c,3) = zb               - dz/2
 
           c = c+1

           X(c,1) = xmin + (ix-1)*dx + dx/3
           X(c,2) = ymin + (iy-1)*dy 
           X(c,3) = zb               - dz*2/3

           c = c+1

           X(c,1) = xmin + (ix-1)*dx + dx*2/3
           X(c,2) = ymin + (iy-1)*dy 
           X(c,3) = zb               - dz*2/3            
            
      ENDDO
   
      ! two boundary vertices at the end of a row in x-direction  
      c = c+1

      X(c,1) = xmax 
      X(c,2) = ymin + (iy-1)*dy + dy*2/3
      X(c,3) = zb               - dz/2

      c = c+1

      X(c,1) = xmax
      X(c,2) = ymin + (iy-1)*dy + dy/3
      X(c,3) = zb               - dz/2

    ENDDO

    DO ix = 1, nx2
      ! two boundary vertices at the end of a row in y-direction  
      c = c+1

      X(c,1) = xmin + (ix-1)*dx + dx/3
      X(c,2) = ymax
      X(c,3) = zb               - dz*2/3    

      c = c+1

      X(c,1) = xmin + (ix-1)*dx + dx*2/3
      X(c,2) = ymax
      X(c,3) = zb               - dz*2/3    
    
    ENDDO

    !--------------------------------------------------------------------------
    ! Compute vertices of coarse mesh 
    
    DO iz = 1, nz2
      DO iy = 1, ny2+1
        DO ix = 1, nx2+1
             
           c = c+1

           X(c,1) = xmin + (ix-1)*dx
           X(c,2) = ymin + (iy-1)*dy
           X(c,3) = zb - iz*dz
        
        ENDDO
      ENDDO
    ENDDO

    WRITE(*,*) ' | -------------------------------------------------- ' 
    WRITE(*,*) ' |',c,' vertices computed !                 '
    WRITE(*,*) ' | -------------------------------------------------- ' 
    WRITE(*,*) '  ' 
   
    !--------------------------------------------------------------------------
    ! Build connectivity of fine mesh 

    c  = 0
    nl = (nx+1)*(ny+1)

    DO iz = 1, nz
      DO iy = 1, ny
        DO ix = 1, nx
             
           c = c+1

           HEX(c,1) =     iz*nl + (iy-1)*(nx+1) + ix
           HEX(c,2) =     iz*nl + (iy-1)*(nx+1) + ix + 1
           HEX(c,3) =     iz*nl +     iy*(nx+1) + ix
           HEX(c,4) =     iz*nl +     iy*(nx+1) + ix + 1
           HEX(c,5) = (iz-1)*nl + (iy-1)*(nx+1) + ix
           HEX(c,6) = (iz-1)*nl + (iy-1)*(nx+1) + ix + 1
           HEX(c,7) = (iz-1)*nl +     iy*(nx+1) + ix
           HEX(c,8) = (iz-1)*nl +     iy*(nx+1) + ix + 1
        
        ENDDO
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Build connectivity of transition mesh 

    nv   = nl * nz
    ntr  = nv + nl + nx2*ny2*8 + (nx2+ny2)*2

    DO iy = 1, ny2
      DO ix = 1, nx2
             
           ! First determine absolute vertex numbers for the 
           ! local 32 vectices of a refined hexahedron

           clocal(1)  = nv + ((iy-1)*(nx+1) + (ix-1))*3 + 1
           clocal(2)  = clocal(1)  + 1
           clocal(3)  = clocal(2)  + 1
           clocal(4)  = clocal(3)  + 1
           clocal(5)  = nv + ((iy-1)*(nx+1) + (ix-1))*3 + nx+1 + 1
           clocal(6)  = clocal(5)  + 1
           clocal(7)  = clocal(6)  + 1
           clocal(8)  = clocal(7)  + 1
           clocal(9)  = nv + ((iy-1)*(nx+1) + (ix-1))*3 + (nx+1)*2 + 1
           clocal(10) = clocal(9)  + 1
           clocal(11) = clocal(10) + 1
           clocal(12) = clocal(11) + 1
           clocal(13) = nv + ((iy-1)*(nx+1) + (ix-1))*3 + (nx+1)*3 + 1
           clocal(14) = clocal(13) + 1
           clocal(15) = clocal(14) + 1
           clocal(16) = clocal(15) + 1
           clocal(17) = nv + nl + (iy-1)*(8*nx2+2) + (ix-1)*8 + 1
           clocal(18) = clocal(17) + 1
           clocal(19) = clocal(18) + 1
           clocal(20) = clocal(19) + 1
           clocal(21) = clocal(20) + 1
           clocal(22) = clocal(21) + 1
           clocal(23) = clocal(22) + 1
           clocal(24) = clocal(23) + 1
           ! Special treatment of boundary vertices at the end of each row in x-direction
           IF(ix.LT.nx2)THEN
              clocal(25) = nv + nl + (iy-1)*(8*nx2+2) + ix*8 + 5
           ELSE
              clocal(25) = nv + nl + (iy-1)*(8*nx2+2) + ix*8 + 1
           END IF
           clocal(26) = clocal(25) + 1 
           ! Special treatment of boundary vertices at the end of each row in y-direction
           IF(iy.LT.ny2)THEN
              clocal(27) = nv + nl + iy*(8*nx2+2) + (ix-1)*8 + 7
           ELSE
              clocal(27) = nv + nl + iy*(8*nx2+2) + (ix-1)*2 + 1
           END IF
           clocal(28) = clocal(27) + 1 
           clocal(29) = ntr + (iy-1)*(nx2+1) + ix
           clocal(30) = clocal(29) + 1
           clocal(31) = ntr +     iy*(nx2+1) + ix
           clocal(32) = clocal(31) + 1

           ! Now build the connectivity of the 13 sub-hexahedrons with the local indices "ind"
           ! Local hexahedron 1
           c = c+1
           ind = (/29,23,22,17,1,2,5,6/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 2
           c = c+1
           ind = (/23,24,17,18,2,3,6,7/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 3
           c = c+1
           ind = (/24,30,18,26,3,4,7,8/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 4
           c = c+1
           ind = (/22,17,21,20,5,6,9,10/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 5
           c = c+1
           ind = (/17,18,20,19,6,7,10,11/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 6
           c = c+1
           ind = (/18,26,19,25,7,8,11,12/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 7
           c = c+1
           ind = (/21,20,31,27,9,10,13,14/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 8
           c = c+1
           ind = (/20,19,27,28,10,11,14,15/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 9
           c = c+1
           ind = (/19,25,28,32,11,12,15,16/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 10
           c = c+1
           ind = (/23,24,27,28,17,18,20,19/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 11
           c = c+1
           ind = (/29,23,31,27,22,17,21,20/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 12
           c = c+1
           ind = (/29,30,31,32,23,24,27,28/)
           HEX(c,1:8) = clocal(ind)
           ! Local hexahedron 13
           c = c+1
           ind = (/24,30,28,32,18,26,19,25/)
           HEX(c,1:8) = clocal(ind)
             
      ENDDO
    ENDDO

    !--------------------------------------------------------------------------
    ! Build connectivity of coarse mesh 

    nl = (nx2+1)*(ny2+1)

    DO iz = 1, nz2-1
      DO iy = 1, ny2
        DO ix = 1, nx2
             
           c = c+1

           HEX(c,1) = ntr +     iz*nl + (iy-1)*(nx2+1) + ix
           HEX(c,2) = ntr +     iz*nl + (iy-1)*(nx2+1) + ix + 1
           HEX(c,3) = ntr +     iz*nl +     iy*(nx2+1) + ix
           HEX(c,4) = ntr +     iz*nl +     iy*(nx2+1) + ix + 1
           HEX(c,5) = ntr + (iz-1)*nl + (iy-1)*(nx2+1) + ix
           HEX(c,6) = ntr + (iz-1)*nl + (iy-1)*(nx2+1) + ix + 1
           HEX(c,7) = ntr + (iz-1)*nl +     iy*(nx2+1) + ix
           HEX(c,8) = ntr + (iz-1)*nl +     iy*(nx2+1) + ix + 1
        
        ENDDO
      ENDDO
    ENDDO

    WRITE(*,*) ' | -------------------------------------------------- ' 
    WRITE(*,*) ' |',c,' hexahedrons built ! '
    WRITE(*,*) ' | -------------------------------------------------- ' 
    WRITE(*,*) '  ' 


    !--------------------------------------------------------------------------
    ! Writing output
    WRITE(*,*) ' | -------------------------------------------------- ' 
    WRITE(*,*) ' | Writing output to ',TRIM(OutFilename)
    WRITE(*,*) ' | -------------------------------------------------- ' 
    WRITE(*,*) '  ' 
    uid2 = 11
    

    OPEN( UNIT = uid2, FILE = TRIM(OutFilename), STATUS = 'UNKNOWN' )
      ! Write GAMBIT header
      WRITE(uid2,'("        CONTROL INFO 2.2.30")')
      WRITE(uid2,'("** GAMBIT NEUTRAL FILE")')
      WRITE(uid2,'("default_id3556")')
      WRITE(uid2,'("PROGRAM:                Gambit     VERSION:  2.2.30")')
      WRITE(uid2,'(" Feb 2007")')
      WRITE(uid2,'("     NUMNP     NELEM     NGRPS    NBSETS     NDFCD     NDFVL")')
      WRITE(uid2,'(I10,I10,I10,I10,I10,I10)') nVert, nElem, 2,0,3,3
      WRITE(uid2,'("ENDOFSECTION")')
      WRITE(uid2,'("   NODAL COORDINATES 2.2.30")')
      ! Write vertices
      DO c = 1, nVert
         WRITE(uid2,'(I10, E20.11E3, E20.11E3, E20.11E3)') c, X(c,:)
      ENDDO
 
      WRITE(uid2,'("ENDOFSECTION")')
      WRITE(uid2,'("      ELEMENTS/CELLS 2.2.30")')
      ! Write connectivity
      DO c = 1,nElem
         WRITE(uid2,'(I8, I3, I3, I9,I8,I8,I8,I8,I8,I8)') c, 4, 8, HEX(c,1:7)
         WRITE(uid2,'(I23)') HEX(c,8)
      ENDDO
      WRITE(uid2,'("ENDOFSECTION")')

      ! Write refined layer as group 1
      ng1  = nx*ny*nz
      WRITE(uid2,'("       ELEMENT GROUP 2.2.30")')
      WRITE(uid2,'("GROUP:          1 ELEMENTS:",I11," MATERIAL:",I11," NFLAGS:",I11)') ng1,4,1
      WRITE(uid2,'("                           layer")')
      WRITE(uid2,'(I8)') 0

      ALLOCATE( member(ng1) )
      
      DO c = 1,ng1
          member(c) = c
      ENDDO
      nline = FLOOR(ng1/10.) + 1
      nrest = MOD(ng1,10)
      DO c = 1, nline-1
         WRITE(uid2,'(I8,I8,I8,I8,I8,I8,I8,I8,I8,I8)') member((c-1)*10+1:c*10)
      ENDDO
      IF(nrest.GE.1)THEN
         WRITE(uid2,'(I8,I8,I8,I8,I8,I8,I8,I8,I8,I8)') member(ng1-(nrest-1):ng1)
      END IF
      ! Write coarse halfspace as group 2
      ng2  = nx2*ny2*((nz2-1)+13)
      WRITE(uid2,'("ENDOFSECTION")')
      WRITE(uid2,'("       ELEMENT GROUP 2.2.30")')
      WRITE(uid2,'("GROUP:          2 ELEMENTS:",I11," MATERIAL:",I11," NFLAGS:",I11)') ng2,4,1
      WRITE(uid2,'("                       halfspace")')
      WRITE(uid2,'(I8)') 0
      DEALLOCATE( member )
      ALLOCATE( member(ng2) )
      DO c = 1,ng2
          member(c) = ng1 + c
      ENDDO
      nline = FLOOR(ng2/10.) + 1
      nrest = MOD(ng2,10)
      DO c = 1, nline-1
         WRITE(uid2,'(I8,I8,I8,I8,I8,I8,I8,I8,I8,I8)') member((c-1)*10+1:c*10)
      ENDDO
      IF(nrest.GE.1)THEN
         WRITE(uid2,'(I8,I8,I8,I8,I8,I8,I8,I8,I8,I8)') member(ng2-(nrest-1):ng2)
      END IF
      WRITE(uid2,'("ENDOFSECTION")')

    CLOSE(uid2)

    WRITE(*,*) ' | -------------------------------------------------- |' 
    WRITE(*,*) ' |     HEX_REFINEMENT  finished successfully          |'
    WRITE(*,*) ' | -------------------------------------------------- |' 
    WRITE(*,*) '  ' 
   
        
END PROGRAM Hex_refinement
