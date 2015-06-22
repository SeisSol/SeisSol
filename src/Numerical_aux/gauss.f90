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

module gauss_mod
  integer, parameter, private :: rk=kind(1.)!set to default single precision
  integer, parameter, private :: ik=kind(1) !set to default integer
contains
  subroutine gauss_jacobi(a,b,xes,weights) 
    use, intrinsic :: iso_fortran_env, only : error_unit
    real,intent(in) :: a,b 
    real,intent(out):: xes(:),weights(size(xes))
    integer,parameter :: dk=kind(1.d0)
    integer :: it,i1,i2,M
    real :: ab,a_b
    real(dk) :: o1,o2,o3,q,d,y,yd,f
    real(dk), parameter :: one=1.d0,two=2.d0,half=0.5d0,quart=0.25d0
    real(dk), parameter :: pi=3.1415926535897932384626433832795028841971694d0
    integer , parameter :: nt=20
    character(len=38), parameter ::message='gauss-jacobi: iteration limit exeeded '
    M=size(xes);ab=a+b;a_b=a-b;f=fexp();yd=two;it=1
    do i1=1,M 
       y=cos((half*a + i1 - quart)/(half*(one + a + b) + M)*pi)
       do it=1,nt
          d = two + a+b 
          o1 = (a_b + d*y)*half; o2 = one
          do i2=2,M 
             o3 = o2;o2=o1;d=a+b+2*i2 
             o1 = (-o3*(two*(a + i2 - 1)*(b + i2 - 1)*d) + o2*((d - one) & 
                  & *(a*a - b*b + d*(d - two)*y)))/(2*i2*(i2 + a+b)*(d - two))
          enddo
          q = (M*(a - b - d*y)*o1 + two*(M + a)*(M + b)*o2)/ & 
               & (d*(one - y*y)); yd = y; y = yd - o1/q 
          if(converged())exit
       enddo
       xes(i1) = y ; weights(i1) = f*d*2.**(a+b)/(q*o2)                               
       if(it>nt)write(error_unit,*)message
    enddo
    return
  contains
    logical function converged()
      converged = (abs(y - yd)<epsilon(y)*one)
    end function converged
    real(dk) function fexp()
      fexp = exp(log_gamma(a + M) + log_gamma(b + M) & 
           &   - log_gamma(M + one) - log_gamma(M + a + b + one))
    end function fexp 
  end subroutine gauss_jacobi
  subroutine gausslegendre(a,b,x,w)
    implicit none
    !dummies:
    real(rk),intent(in) :: a,b
    real(rk),intent(out):: x(:),w(size(x))
    !locals:
    real(rk), parameter:: env=epsilon(1._rk)*100._rk,pi=acos(-1.0_rk)
    integer :: j,N,NH
    real(rk):: pop,x1,x2,s
    N=size(x);NH=(N+1)/2
    x2 =(b+a)*0.5_rk; x1=(b-a)*0.5_rk
    do j=1,NH
       s = cos(pi*(j-.25_rk)/(.5_rk+N))
       call subint(s,pop,N)
       x(j)=x2-x1*s;x(N+1-j)=x2+x1*s
       w(j)=2._rk*x1/((1._rk-s*s)*pop*pop); w(N+1-j)=2._rk*x1/((1._rk-s*s)*pop*pop)
    enddo
    return
  contains
    subroutine subint(s,pop,N)
      implicit none
      real(rk),intent(inout):: s,pop
      integer(ik),intent(in) :: N
      real(rk) :: pmo,pmt,pmth,sone
      integer(ik) :: j
      real(rk), parameter:: env=epsilon(1._rk)*100._rk
      sone=2._rk
      do while(abs(s-sone) > env)
         pmo=1._rk; pmt=0._rk
         do j=1,N
            pmth=pmt; pmt =pmo;
            pmo =((2._rk*j-1._rk)*s*pmt-(j-1._rk)*pmth)/j
         enddo
         sone=s;pop=N*(s*pmo-pmt)/(s*s-1._rk); s=sone-pmo/pop
      enddo
    end subroutine subint
  end subroutine gausslegendre
  
end module gauss_mod
