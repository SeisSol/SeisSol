!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Alex Breuer (breuer AT mytum.de, http://www5.in.tum.de/wiki/index.php/Dipl.-Math._Alexander_Breuer)
!! @author Sebastian Rettenberger (sebastian.rettenberger @ tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger)
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
!! C++/Fortran-interoperability: Fortran-bindings

#include "Initializer/preProcessorMacros.fpp"

module f_ftoc_bind_interoperability
  implicit none
  !
  ! Fortran to C bindings
  !

  interface c_interoperability_setDomain
    subroutine c_interoperability_setDomain( i_domain ) bind( C, name='c_interoperability_setDomain' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_domain
    end subroutine
  end interface

  interface c_interoperability_setTimeStepWidth
    subroutine c_interoperability_setTimeStepWidth( i_meshId, i_timeStepWidth ) bind( C, name='c_interoperability_setTimeStepWidth' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: i_timeStepWidth
    end subroutine
  end interface

  interface c_interoperability_initializeClusteredLts
    subroutine c_interoperability_initializeClusteredLts( i_clustering ) bind( C, name='c_interoperability_initializeClusteredLts' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_clustering
    end subroutine
  end interface
  
  interface c_interopability_allocatePointSources
    subroutine c_interopability_allocatePointSources( i_meshIds, i_numberOfPointSources ) bind ( C, name='c_interopability_allocatePointSources' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshIds
      type(c_ptr), value :: i_numberOfPointSources
    end subroutine
  end interface
  
  interface c_interopability_setupPointSource
    subroutine c_interopability_setupPointSource( i_source, i_mInvJInvPhisAtSources, i_localMomentTensor, i_strike, i_dip, i_rake, i_samples, i_numberOfSamples, i_onsetTime, i_samplingInterval ) bind ( C, name='c_interopability_setupPointSource' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_source
      type(c_ptr), value :: i_mInvJInvPhisAtSources
      type(c_ptr), value :: i_localMomentTensor
      type(c_ptr), value :: i_strike
      type(c_ptr), value :: i_dip
      type(c_ptr), value :: i_rake
      type(c_ptr), value :: i_samples
      type(c_ptr), value :: i_numberOfSamples
      type(c_ptr), value :: i_onsetTime
      type(c_ptr), value :: i_samplingInterval
    end subroutine
  end interface

  interface c_interoperability_addReceiver
    subroutine c_interoperability_addReceiver( i_receiverId, i_meshId ) bind( C, name='c_interoperability_addReceiver' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_receiverId, i_meshId
    end subroutine
  end interface

  interface c_interoperability_setReceiverSampling
    subroutine c_interoperability_setReceiverSampling( i_receiverSampling ) bind( C, name='c_interoperability_setReceiverSampling' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_receiverSampling
    end subroutine
  end interface
  
  interface c_interoperability_setMaterial
    subroutine c_interoperability_setMaterial( i_elem, i_side, i_materialVal, i_numMaterialVals ) bind( C, name='c_interoperability_setMaterial' )
    use iso_c_binding, only: c_ptr
    implicit none
    type(c_ptr), value :: i_elem
    type(c_ptr), value :: i_side
    type(c_ptr), value :: i_materialVal
    type(c_ptr), value :: i_numMaterialVals
    end subroutine
  end interface

  interface c_interoperability_setInitialLoading
    subroutine c_interoperability_setInitialLoading( i_meshId, i_initialLoading ) bind( C, name='c_interoperability_setInitialLoading' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: i_initialLoading
    end subroutine
  end interface
  
  interface c_interoperability_initializeCellLocalMatrices
    subroutine c_interoperability_initializeCellLocalMatrices() bind( C, name='c_interoperability_initializeCellLocalMatrices' )
    end subroutine
  end interface

  interface c_interoperability_synchronizeMaterial
    subroutine c_interoperability_synchronizeMaterial() bind( C, name='c_interoperability_synchronizeMaterial' )
    end subroutine
  end interface

  interface c_interoperability_synchronizeCopyLayerDofs
    subroutine c_interoperability_synchronizeCopyLayerDofs() bind( C, name='c_interoperability_synchronizeCopyLayerDofs' )
    end subroutine
  end interface

  interface c_interoperability_enableDynamicRupture
    subroutine c_interoperability_enableDynamicRupture() bind( C, name='c_interoperability_enableDynamicRupture' )
    end subroutine
  end interface

  interface
    subroutine c_interoperability_enableWaveFieldOutput( i_waveFieldInterval, i_waveFieldFilename ) bind( C, name='c_interoperability_enableWaveFieldOutput' )
      use iso_c_binding, only: c_ptr, c_char
      implicit none
      type(c_ptr), value :: i_waveFieldInterval
      character(kind=c_char), dimension(*), intent(in) :: i_waveFieldFilename
    end subroutine
  end interface

  interface
    subroutine c_interoperability_enableCheckPointing( i_checkPointInterval, i_checkPointFilename, i_checkPointBackend ) bind( C, name='c_interoperability_enableCheckPointing' )
      use iso_c_binding, only: c_ptr, c_char
      implicit none
      type(c_ptr), value :: i_checkPointInterval
      character(kind=c_char), dimension(*), intent(in) :: i_checkPointFilename
      character(kind=c_char), dimension(*), intent(in) :: i_checkPointBackend
    end subroutine
  end interface

  interface c_interoperability_initializeIO
    subroutine c_interoperability_initializeIO( i_mu, i_slipRate1, i_slipRate2, i_slip, i_state, i_strength, i_numSides, i_numBndGP ) &
        bind( C, name='c_interoperability_initializeIO' )
      use iso_c_binding, only: c_ptr
      implicit none

      type(c_ptr), value                     :: i_mu
      type(c_ptr), value                     :: i_slipRate1
      type(c_ptr), value                     :: i_slipRate2
      type(c_ptr), value                     :: i_slip
      type(c_ptr), value                     :: i_state
      type(c_ptr), value                     :: i_strength
      type(c_ptr), value                     :: i_numSides
      type(c_ptr), value                     :: i_numBndGP
    end subroutine
  end interface

  interface c_interoperability_addToDofs
    subroutine c_interoperability_addToDofs( i_meshId, i_update ) bind( C, name='c_interoperability_addToDofs' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: i_update
    end subroutine
  end interface

  interface c_interoperability_getTimeDerivatives
    subroutine c_interoperability_getTimeDerivatives( i_meshId, o_timeDerivatives ) bind( C, name='c_interoperability_getTimeDerivatives' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: o_timeDerivatives
    end subroutine
  end interface

  interface c_interoperability_getFaceDerInt
    subroutine c_interoperability_getFaceDerInt( i_meshId, i_faceId, i_timeStepWidth, o_timeDerivativesCell, o_timeDerivativesNeighbor, o_timeIntegratedCell, o_timeIntegratedNeighbor ) bind( C, name='c_interoperability_getFaceDerInt' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: i_faceId
      type(c_ptr), value :: i_timeStepWidth
      type(c_ptr), value :: o_timeDerivativesCell
      type(c_ptr), value :: o_timeDerivativesNeighbor
      type(c_ptr), value :: o_timeIntegratedCell
      type(c_ptr), value :: o_timeIntegratedNeighbor
    end subroutine
  end interface

  interface c_interoperability_getDofs
    subroutine c_interoperability_getDofs( i_meshId, o_dofs ) bind( C, name='c_interoperability_getDofs' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: o_dofs
    end subroutine
  end interface

  interface c_interoperability_getDofsFromDerivatives
    subroutine c_interoperability_getDofsFromDerivatives( i_meshId, o_dofs ) bind( C, name='c_interoperability_getDofsFromDerivatives' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: o_dofs
    end subroutine
  end interface

  interface c_interoperability_getNeighborDofsFromDerivatives
    subroutine c_interoperability_getNeighborDofsFromDerivatives( i_meshId, i_faceId, o_dofs ) bind( C, name='c_interoperability_getNeighborDofsFromDerivatives' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_meshId
      type(c_ptr), value :: i_faceId
      type(c_ptr), value :: o_dofs
    end subroutine
  end interface

  interface c_interoperability_simulate
    subroutine c_interoperability_simulate( i_finalTime ) bind( C, name='c_interoperability_simulate' )
      use iso_c_binding, only: c_ptr
      implicit none
      type(c_ptr), value :: i_finalTime
    end subroutine
  end interface

  interface c_interoperability_finalizeIO
    subroutine c_interoperability_finalizeIO() bind( C, name='c_interoperability_finalizeIO' )
    end subroutine
  end interface
end module
