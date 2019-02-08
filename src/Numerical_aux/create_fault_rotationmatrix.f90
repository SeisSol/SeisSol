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
!! routine creates rotationmatrix for a particular point on the fault thata
!! rotates the tangential fault vectors in strike and dip direction

MODULE create_fault_rotationmatrix_mod
  !---------------------------------------------------------------------------!
  USE TypesDef
  USE JacobiNormal_mod  
  !---------------------------------------------------------------------------!
  IMPLICIT NONE
  PRIVATE
  !---------------------------------------------------------------------------!
  PUBLIC  :: create_fault_rotationmatrix
  !---------------------------------------------------------------------------!  
  INTERFACE create_fault_rotationmatrix
     MODULE PROCEDURE create_fault_rotationmatrix
  END INTERFACE

CONTAINS

  SUBROUTINE create_fault_rotationmatrix(rotmat,iFace,EQN,MESH,iRotmat)
    !> routine creates rotationmatrix for a particular point on the fault that
    !> rotates the tangential fault vectors in strike and dip direction 
    !-------------------------------------------------------------------------!
    IMPLICIT NONE    
    !-------------------------------------------------------------------------!
    ! Argument list declaration
    TYPE(tEquations)        :: EQN                                            !< Equation structure
    TYPE(tUnstructMesh)     :: MESH                                           !< Mesh structure
    !-------------------------------------------------------------------------!
    ! Local variable declaration                                              !
    INTEGER                         :: iFace                                  !< fault face index
    REAL                            :: rotmat(1:6,1:6)                        !< rotation matrix for individual fault receiver
    REAL,optional                   :: iRotmat(1:6,1:6)                       !
    REAL                            :: T(EQN%nVar,EQN%nVar)                   !< Transformation matrix
    REAL                            :: iT(EQN%nVar,EQN%nVar)                  !< inverse Transformation matrix
    REAL                            :: n(1:3)                                 !< normal vector to fault, pointing away from reference point
    REAL                            :: strike_vector(1:3)                     !< strike vector
    REAL                            :: dip_vector(1:3)                        !< dip vector
    REAL                            :: norm                                   !< norm
    !-------------------------------------------------------------------------!
    INTENT(IN)    :: iFace,MESH, EQN                                                !
    intent(out)   :: rotmat,iRotmat
    !-------------------------------------------------------------------------!
 
    ! Works only with following right-handed coordinate system:
    ! +x : east
    ! +y : north
    ! -z : depth
    
  
    ! Local side's normal (and tangential vectors - commented)
    ! Note, normal vector n is pointing away from reference point
    n = MESH%Fault%geoNormals(1:3,iFace)
    ! s = MESH%Fault%geoTangent1(1:3,iFace)
    ! t = MESH%Fault%geoTangent2(1:3,iFace)

    strike_vector(1) = n(2)/sqrt(n(1)**2+n(2)**2)
    strike_vector(2) = -n(1)/sqrt(n(1)**2+n(2)**2)
    strike_vector(3) = 0.0D0
   
    dip_vector(1) = -strike_vector(2)*n(3)
    dip_vector(2) = strike_vector(1)*n(3)
    dip_vector(3) = strike_vector(2)*n(1) - strike_vector(1)*n(2)
    norm = 1.0D0/sqrt(dip_vector(1)**2+dip_vector(2)**2+dip_vector(3)**2)
    dip_vector(1) = dip_vector(1)*norm
    dip_vector(2) = dip_vector(2)*norm
    dip_vector(3) = dip_vector(3)*norm
    
    CALL RotationMatrix3D(n,strike_vector,dip_vector,T(:,:),iT(:,:),EQN)
    rotmat = iT(1:6,1:6)

    if(present(iRotmat)) then
      iRotmat = T(1:6,1:6)
    endif  
  END SUBROUTINE create_fault_rotationmatrix
  
END MODULE create_fault_rotationmatrix_mod
