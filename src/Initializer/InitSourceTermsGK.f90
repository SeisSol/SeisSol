!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
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

!! @section DESCRIPTION
!! Conversion of the fortran source term data structures to the C++ source term data structures.

#ifdef GENERATEDKERNELS

#include <Initializer/preProcessorMacros.fpp>

! Supports only source type 50.
module InitSourceTermsGK_mod
  use TypesDef

  implicit none
  private
  
  interface InitSourceTermsGK
    module procedure InitSourceTermsGK
  end interface
  
  public :: InitSourceTermsGK
  
contains
  subroutine InitSourceTermsGK(DISC, MESH, SOURCE)
    !-------------------------------------------------------------------------!
    use DGBasis_mod
    use f_ftoc_bind_interoperability
    use iso_c_binding
    !-------------------------------------------------------------------------!
    implicit none
    !-------------------------------------------------------------------------!
    type (tDiscretization), target                        :: DISC
    type (tUnstructMesh)                                  :: MESH
    type (tSource)                                        :: SOURCE
    !-------------------------------------------------------------------------!
    integer                                                       :: l_source, l_elem, l_dof, l_numberOfSources, l_newSource
    integer, pointer                                              :: l_elements(:)
    integer, pointer                                              :: l_oldSourceIndex(:)
    real, dimension(NUMBER_OF_ALIGNED_BASIS_FUNCTIONS), target    :: l_mInvJInvPhisAtSources
    real, target                                                  :: l_momentTensor(3,3)
    real                                                          :: l_xi, l_eta, l_zeta
    
    select case(SOURCE%Type)
    case(0)
      ! No source terms
      ! Do nothing
    case(50)
      allocate( l_elements(SOURCE%RP%nSbfs(1)) )
      allocate( l_oldSourceIndex(SOURCE%RP%nSbfs(1)) )
      
      ! Clean invalid sources
      l_numberOfSources = 0
      do l_source = 1, SOURCE%RP%nSbfs(1)
        if (SOURCE%RP%Element(l_source) .GE. 1) then
          l_numberOfSources                   = l_numberOfSources + 1
          l_elements(l_numberOfSources)       = SOURCE%RP%Element(l_source)
          l_oldSourceIndex(l_numberOfSources) = l_source
        end if
      end do

      call c_interopability_allocatePointSources( i_meshIds = c_loc(l_elements),    &
                                                  i_numberOfPointSources = c_loc(l_numberOfSources) )
      
      do l_newSource = 1, l_numberOfSources
          l_elem = l_elements(l_newSource)
          l_source = l_oldSourceIndex(l_newSource)
          
          !-------------------------------------------------------------------------!
          ! Compute and set all phisAtSource                                        !
          !-------------------------------------------------------------------------!
          call TrafoXYZ2XiEtaZeta(xi    = l_xi,                                                                 &
                                  eta   = l_eta,                                                                &
                                  zeta  = l_zeta,                                                               &
                                  xP    = SOURCE%RP%SpacePosition(1, l_source),                                 &
                                  yP    = SOURCE%RP%SpacePosition(2, l_source),                                 &
                                  zP    = SOURCE%RP%SpacePosition(3, l_source),                                 &
                                  x     = MESH%VRTX%xyNode(1, MESH%ELEM%Vertex(1:MESH%nVertices_Tet, l_elem)),  &
                                  y     = MESH%VRTX%xyNode(2, MESH%ELEM%Vertex(1:MESH%nVertices_Tet, l_elem)),  &
                                  z     = MESH%VRTX%xyNode(3, MESH%ELEM%Vertex(1:MESH%nVertices_Tet, l_elem)),  &
                                  vType = MESH%LocalVrtxType(l_elem)                                            )
          do l_dof = 1, NUMBER_OF_BASIS_FUNCTIONS
            call BaseFunc3D(l_mInvJInvPhisAtSources(l_dof),    &
                            l_dof,                             &
                            l_xi,                              &
                            l_eta,                             &
                            l_zeta,                            &
                            DISC%Galerkin%nPoly,               &
                            DISC%Galerkin%cPoly3D_Tet,         &
                            DISC%Galerkin%NonZeroCPoly_Tet,    &
                            DISC%Galerkin%NonZeroCPolyIndex_Tet)

            l_mInvJInvPhisAtSources(l_dof) = l_mInvJInvPhisAtSources(l_dof) / ( 6.0d0 * MESH%ELEM%Volume(l_elem) * DISC%Galerkin%MassMatrix_Tet(l_dof, l_dof, DISC%Galerkin%nPoly) )
          end do
          
          l_momentTensor = SOURCE%RP%Area(l_source) * SOURCE%RP%MomentTensor
          
          call c_interopability_setupPointSource( i_source                = c_loc(l_newSource),                     &
                                                  i_mInvJInvPhisAtSources = c_loc(l_mInvJInvPhisAtSources),         &
                                                  i_localMomentTensor     = c_loc(l_momentTensor),                  &
                                                  i_strike                = c_loc(SOURCE%RP%strks(1,l_source)),     &
                                                  i_dip                   = c_loc(SOURCE%RP%dips(1,l_source)),      &
                                                  i_rake                  = c_loc(SOURCE%RP%rake(1,l_source)),      &
                                                  i_samples               = c_loc(SOURCE%RP%TimeHist(:,l_source)),  &
                                                  i_numberOfSamples       = c_loc(SOURCE%RP%nsteps),                &
                                                  i_onsetTime             = c_loc(SOURCE%RP%Tonset(l_source)),      &
                                                  i_samplingInterval      = c_loc(SOURCE%RP%t_samp)                 )
      end do

      deallocate(l_elements)
      deallocate(l_oldSourceIndex)
    case default
      logError(*)  'Generated Kernels: Unsupported source type: ', SOURCE%Type
      stop
    end select
    
  end subroutine
end module

#endif
