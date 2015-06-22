#if 0
!>
!! @file
!! This file is part of SeisSol.
!!
!! @author Sebastian Rettenberger (sebastian.rettenberger AT tum.de, http://www5.in.tum.de/wiki/index.php/Sebastian_Rettenberger,_M.Sc.)
!!
!! @section LICENSE
!! Copyright (c) 2014, SeisSol Group
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
!! Defines macros for instrumentation
#endif

#ifndef INSTRUMENTATION_FPP
#define INSTRUMENTATION_FPP

#if 0
! Manual instrumentation for Scalasca with epik.
! Override function calls if not compiled with EPIK.
#endif

#if defined(EPIK)
#if (!defined(__STDC__) || (defined(BG) && !defined(__bg__)))
#include "epik_user.inc"
#else
#include "epik_user.h"
#endif
#else
#define EPIK_FUNC_REG(str)
#define EPIK_FUNC_START()
#define EPIK_FUNC_END()
#define EPIK_USER_REG(id,str)
#define EPIK_USER_START(id)
#define EPIK_USER_END(id)
#define EPIK_TRACER(str)
#endif

#if 0
! empty score-p definitions without usage
#endif

#if defined(SCOREP_USER_ENABLE)

#if (!defined(__STDC__) || (defined(BG) && !defined(__bg__)))
#include <scorep/SCOREP_User.inc>
#else
#include <scorep/SCOREP_User.h>
#endif

#else

#if (!defined(__STDC__) || (defined(BG) && !defined(__bg__)))
#define SCOREP_USER_REGION_DEFINE( handle )
#define SCOREP_USER_REWIND_DEFINE( handle )
#define SCOREP_USER_OA_PHASE_BEGIN( handle, name, type  )
#define SCOREP_USER_OA_PHASE_END( handle )
#define SCOREP_USER_REWIND_POINT( handle, name )
#define SCOREP_USER_REGION_BEGIN( handle, name, type )
#define SCOREP_USER_REGION_INIT( handle, name, type )
#define SCOREP_USER_REGION_END( handle )
#define SCOREP_USER_REWIND_CHECK( handle, value )
#define SCOREP_USER_REGION_ENTER( handle )
#define SCOREP_USER_FUNC_DEFINE()
#define SCOREP_USER_FUNC_BEGIN( name )
#define SCOREP_USER_FUNC_END()
#define SCOREP_USER_GLOBAL_REGION_DEFINE( handle )
#define SCOREP_USER_GLOBAL_REGION_EXTERNAL( handle )
#define SCOREP_USER_PARAMETER_DEFINE( handle )
#define SCOREP_USER_PARAMETER_INT64( handle, name, value )
#define SCOREP_USER_PARAMETER_UINT64( handle, name, value )
#define SCOREP_USER_PARAMETER_STRING( handle, name, value )
#define SCOREP_USER_METRIC_GLOBAL( metricHandle )
#define SCOREP_USER_METRIC_EXTERNAL( metricHandle )
#define SCOREP_USER_METRIC_LOCAL( metricHandle )
#define SCOREP_USER_METRIC_INIT( metricHandle, name, unit, type, context )
#define SCOREP_USER_METRIC_INT64( metricHandle, value )
#define SCOREP_USER_METRIC_UINT64( metricHandle, value )
#define SCOREP_USER_METRIC_DOUBLE( metricHandle, value )
#define SCOREP_RECORDING_ON()
#define SCOREP_RECORDING_OFF()
#define SCOREP_RECORDING_IS_ON( isOn ) isOn = .false.
#else
#define SCOREP_USER_REGION( name, type )
#define SCOREP_USER_REGION_DEFINE( handle )
#define SCOREP_USER_OA_PHASE_BEGIN( handle, name, type  )
#define SCOREP_USER_OA_PHASE_END( handle )
#define SCOREP_USER_REGION_BEGIN( handle, name, type )
#define SCOREP_USER_REGION_INIT( handle, name, type )
#define SCOREP_USER_REGION_END( handle )
#define SCOREP_USER_REGION_ENTER( handle )
#define SCOREP_USER_FUNC_BEGIN()
#define SCOREP_USER_FUNC_END()
#define SCOREP_GLOBAL_REGION_DEFINE( handle )
#define SCOREP_GLOBAL_REGION_EXTERNAL( handle )
#define SCOREP_USER_PARAMETER_INT64( name, value )
#define SCOREP_USER_PARAMETER_UINT64( name, value )
#define SCOREP_USER_PARAMETER_STRING( name, value )
#define SCOREP_USER_METRIC_GLOBAL( metricHandle )
#define SCOREP_USER_METRIC_EXTERNAL( metricHandle )
#define SCOREP_USER_METRIC_LOCAL( metricHandle )
#define SCOREP_USER_METRIC_INIT( metricHandle, name, unit, type, context )
#define SCOREP_USER_METRIC_INT64( metricHandle, value )
#define SCOREP_USER_METRIC_UINT64( metricHandle, value )
#define SCOREP_USER_METRIC_DOUBLE( metricHandle, value )
#define SCOREP_RECORDING_ON()
#define SCOREP_RECORDING_OFF()
#define SCOREP_RECORDING_IS_ON() 0
#endif

#endif

#endif
