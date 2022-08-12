!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
!! @author Carsten Uphoff (c.uphoff AT tum.de, http://www5.in.tum.de/wiki/index.php/Carsten_Uphoff,_M.Sc.)
!!
!! @section LICENSE
!! Copyright (c) 2015-2017, SeisSol Group
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
!! C++/Fortran-interoperability: Fortran-bindings

#include "Initializer/preProcessorMacros.fpp"

module f_ctof_bind_interoperability
  implicit none

  contains
    subroutine copyDynamicRuptureState(domain, fromMeshId, toMeshId)
      use typesDef
      implicit none
      type(tUnstructDomainDescript), pointer :: domain
      integer     :: fromMeshId, toMeshId,iFace
    end subroutine
    !
    ! C to fortran bindings
    !
    subroutine f_interoperability_copyDynamicRuptureState( i_domain ) bind (c, name='f_interoperability_copyDynamicRuptureState')
      use iso_c_binding
      use typesDef
      use f_ftoc_bind_interoperability
      implicit none

      type(c_ptr), value                     :: i_domain
      type(tUnstructDomainDescript), pointer :: l_domain

      ! convert c to fortran pointers
      call c_f_pointer( i_domain, l_domain)
      
      call copyDynamicRuptureState(l_domain, 1, l_domain%mesh%Fault%nSide)
    end subroutine

    subroutine f_interoperability_fitAttenuation( domain, rho, mu, lambda, Qp, Qs, materialFitted) bind( c, name='f_interoperability_fitAttenuation')
      use iso_c_binding
      use TypesDef
      use ini_MODEL_mod

      type(c_ptr), value                     :: domain
      type(tUnstructDomainDescript), pointer :: l_domain

      real(kind=c_double), value             :: rho, mu, lambda, Qp, Qs

      type(c_ptr), value                     :: materialFitted
      real(kind=8), pointer                  :: l_materialFitted(:)

      real                                   :: material(5)

      call c_f_pointer( domain,                   l_domain                                         )
      call c_f_pointer( materialFitted,           l_materialFitted,  [l_domain%EQN%nBackgroundVar] )

      material(:) = (/ rho, mu, lambda, Qp, Qs /)
      call fitAttenuation(material, l_materialFitted, l_domain%EQN)
    end subroutine
end module
