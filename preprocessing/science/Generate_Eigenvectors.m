%%
% @file
% This file is part of SeisSol.
%
% @author Josep de la Puente Alvarez (josep.delapuente AT bsc.es, http://www.geophysik.uni-muenchen.de/Members/jdelapuente)
%
% @section LICENSE
% Copyright (c) 2007, SeisSol Group
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice,
%    this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder nor the names of its
%    contributors may be used to endorse or promote products derived from this
%    software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
% ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
% LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
% CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
% SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
% INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
% CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
% ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.
%
% @section DESCRIPTION
% MATLAB script for obtaning the eigenvalues and eigenvectors of elastic,
% viscoelastic, anisotropic or poroelastic materials
%
% The elastic and viscoelastic cases read the parameters directly from the
% script, whereas the anisotropic and poroelastic cases read the parameters
% from an external SeisSol material definition file (*.def). Notice that
% only the first material in the file is read.
%
% The resulting *.dat file can be read by SeisSol when using the initial
% conditions of type 'Planewave_Ricker_Puls' or 'Planewave_Gauss_Puls'.
% In this case, a plane wave will be generated using one chosen eigenvector
% from those listed in the *.dat file.
%
% The eigenvalues and eigenvectors are ordered in a descending way, to be
% as close as possible to the standard ordering in the literature. However,
% some ambiguity remains in the case of waves with identical wave speed
% of propagation. In this case, the ordering of the eigenvalues is
% arbitrary and visual inspection of the eigenvector matrix is recommended.

home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  Generates eigenvectors ans eigenvalues for a specified      %%')
disp('     %%  material and propagation direction.                         %%')
disp('     %%                                                              %%')
disp('     %%  The material can be either elastic, viscoelastic, aniso-    %%')
disp('     %%  tropic or poroelastic.                                      %%')
disp('     %%                                                              %%')
disp('     %%  The output is a file that can be read by SeisSol with the   %%')
disp('     %%  initial condition types "Planewave_Ricker_Puls" or          %%')
disp('     %%  "Planewave_Gauss_Puls".                                     %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')
clear, close all;

rheo=input('     Choose rheology type (1:ela., 2:viscoela., 3:aniso., 4:poro.):');
disp('     Write unitary vector of propagation direction:')
kx=input('       nx=');
ky=input('       ny=');
kz=input('       nz=');
fname_out=input('     Enter name of output file (*.dat): ','s');


if(rheo==1)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PARAMETER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    rho     = 1;                %density
    lambda  = 1;                %Lame constant
    mu      = 1;                %Lame constant
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    disp('-----------------------------------------------------------------')
    disp('     Chosen parameters are:');
    disp(sprintf('      Density rho             = %g',rho));
    disp(sprintf('      Lame constant (lambda)  = %g',lambda));
    disp(sprintf('      Lame constant (mu)      = %g',mu));
    
    %% Build Jacobian matrix A
    A = [0       0       0       0       0       0      -(lambda+2*mu)    0       0;
         0       0       0       0       0       0      -lambda           0       0;
         0       0       0       0       0       0      -lambda           0       0;
         0       0       0       0       0       0       0              -mu       0;
         0       0       0       0       0       0       0                0       0;
         0       0       0       0       0       0       0                0     -mu;
         -1/rho  0       0       0       0       0       0                0       0;
         0       0       0      -1/rho   0       0       0                0       0;
         0       0       0       0       0    -1/rho     0                0       0 ];


    %% Build Jacobian matrix B
    B = [0       0       0       0       0       0       0           -lambda       0;
         0       0       0       0       0       0       0    -(lambda+2*mu)       0;
         0       0       0       0       0       0       0           -lambda       0;
         0       0       0       0       0       0     -mu                0        0;
         0       0       0       0       0       0       0                0      -mu;
         0       0       0       0       0       0       0                0        0;
         0       0       0      -1/rho   0       0       0                0        0;
         0     -1/rho    0       0       0       0       0                0        0;
         0       0       0       0    -1/rho     0       0                0       0 ];

    %% Build Jacobian matrix C
    C = [0       0       0       0       0       0       0        0    -lambda;
         0       0       0       0       0       0       0        0    -lambda;
         0       0       0       0       0       0       0        0    -(lambda+2*mu);
         0       0       0       0       0       0       0        0        0;
         0       0       0       0       0       0       0       -mu       0;
         0       0       0       0       0       0      -mu       0        0;
         0       0       0       0       0    -1/rho     0        0        0; 
         0       0       0       0     -1/rho    0       0        0        0; 
         0       0    -1/rho     0       0       0       0        0        0 ];
    
elseif(rheo==2)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PARAMETER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    rho     = 1;                %density
    lambda  = 1;                %Lame constant
    mu      = 1;                %Lame constant
    QPval   = 20;               %Q-factor for P-waves
    QSval   = 10;               %Q-factor for S-waves
    n       = 3;                %Number of mechanisms, i.e. Maxwell bodies 
    freq    = 1;                %Central frequency of the absorption band (in Hertz)
    f_ratio = 100;              %The ratio between the maximum and minimum frequencies of our bandwidth
                                %Usually 10^4 is taken as good (Day 2003)                             
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%% Output chosen parameters on screen %%%%%%%%%%%%%%%%%%%%%
    disp('-----------------------------------------------------------------')
    disp('     Chosen parameters are:');
    disp(sprintf('      Density rho             = %g',rho));
    disp(sprintf('      Lame constant (lambda)  = %g',lambda));
    disp(sprintf('      Lame constant (mu)      = %g',mu));
    disp(sprintf('      Q for P-wave            = %g',QPval));
    disp(sprintf('      Q for S-wave            = %g',QSval));

    % Selection of the logarithmically equispaced frequencies
    % i.e.: log(w) = log(wmin) + (i-1)*log(f_ratio)/(2*(n-1))
    w0      = 2*pi*freq;
    wmin    = w0/sqrt(f_ratio); 
    kmax    = 2*n-1;              
    
    if n>1
      for i=1:kmax
        w(i) = exp(log(wmin)+(i-1)/(2*(n-1))*log(f_ratio));  
      end
    else
        w = w0;
    end


    %% Build Jacobian matrix A
    A = [0       0       0       0       0       0      -(lambda+2*mu)    0       0;
         0       0       0       0       0       0      -lambda           0       0;
         0       0       0       0       0       0      -lambda           0       0;
         0       0       0       0       0       0       0              -mu       0;
         0       0       0       0       0       0       0                0       0;
         0       0       0       0       0       0       0                0     -mu;
         -1/rho  0       0       0       0       0       0                0       0;
         0       0       0      -1/rho   0       0       0                0       0;
         0       0       0       0       0    -1/rho     0                0       0 ];

    %% extend by n mechanisms
    for i = 1:n
        A = [A;
             [0       0      0      0       0      0   -w(i)     0         0 ;
              0       0      0      0       0      0    0        0         0 ;
              0       0      0      0       0      0    0        0         0 ;
              0       0      0      0       0      0    0       -w(i)/2    0 ;
              0       0      0      0       0      0    0        0         0 ;          
              0       0      0      0       0      0    0        0       -w(i)/2];
            ]; 
    end
    A = [A,zeros(9+n*6,n*6)];

    %% Build Jacobian matrix B
    B = [0       0       0       0       0       0       0           -lambda       0;
         0       0       0       0       0       0       0    -(lambda+2*mu)       0;
         0       0       0       0       0       0       0           -lambda       0;
         0       0       0       0       0       0     -mu                0        0;
         0       0       0       0       0       0       0                0      -mu;
         0       0       0       0       0       0       0                0        0;
         0       0       0      -1/rho   0       0       0                0        0;
         0     -1/rho    0       0       0       0       0                0        0;
         0       0       0       0    -1/rho     0       0                0       0 ];
    %% extend by n mechanisms
    for i = 1:n
        B = [B;
             [0       0      0      0       0      0    0        0         0 ;
              0       0      0      0       0      0    0       -w(i)      0 ;
              0       0      0      0       0      0    0        0         0 ;
              0       0      0      0       0      0  -w(i)/2    0         0 ;
              0       0      0      0       0      0    0        0       -w(i)/2 ;          
              0       0      0      0       0      0    0        0         0];
            ];
    end
    B = [B,zeros(9+n*6,n*6)];

    %% Build Jacobian matrix C
    C = [0       0       0       0       0       0       0        0    -lambda;
         0       0       0       0       0       0       0        0    -lambda;
         0       0       0       0       0       0       0        0    -(lambda+2*mu);
         0       0       0       0       0       0       0        0        0;
         0       0       0       0       0       0       0       -mu       0;
         0       0       0       0       0       0      -mu       0        0;
         0       0       0       0       0    -1/rho     0        0        0; 
         0       0       0       0     -1/rho    0       0        0        0; 
         0       0    -1/rho     0       0       0       0        0        0 ];
    %% extend by n mechanisms
    for i = 1:n
        C = [C;
             [0       0      0      0       0      0    0        0         0 ;
              0       0      0      0       0      0    0        0         0 ;
              0       0      0      0       0      0    0        0        -w(i) ;
              0       0      0      0       0      0    0        0         0 ;
              0       0      0      0       0      0    0      -w(i)/2     0 ;          
              0       0      0      0       0      0   -w(i)/2   0         0];
            ];
    end
    C = [C,zeros(9+n*6,n*6)];


elseif(rheo==3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PARAMETER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname=input('     Enter the file name of the anisotropic properties (*.def):','s')
    fid = fopen([fname,'.def']);
      junk = fgetl(fid);
      data = fscanf(fid,'%g',[1,32]);
    fclose(fid);
    disp(['     Successfully read parameters from file: ',fname,'.def!'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %coefficients given in N/m^2 
    ce(1,1)     = data(3);
    ce(1,2)     = data(4);
    ce(1,3)     = data(5);
    ce(1,4)     = data(6);
    ce(1,5)     = data(7);
    ce(1,6)     = data(8);
    ce(2,2)     = data(9);  
    ce(2,3)     = data(10);
    ce(2,4)     = data(11);
    ce(2,5)     = data(12);
    ce(2,6)     = data(13);
    ce(3,3)     = data(14);
    ce(3,4)     = data(15);
    ce(3,5)     = data(16);
    ce(3,6)     = data(17);
    ce(4,4)     = data(18);
    ce(4,5)     = data(19);
    ce(4,6)     = data(20);
    ce(5,5)     = data(21);
    ce(5,6)     = data(22);
    ce(6,6)     = data(23);
    rho         = data(2);

    % direction of material symmetry axes
    nx = data(24);
    ny = data(25);
    nz = data(26);
    sx = data(27);
    sy = data(28);
    sz = data(29);
    tx = data(30);
    ty = data(31);
    tz = data(32);
    
    for i = 1:6 
    for j = 1:i-1
    ce(i,j) = ce(j,i);
    end
    end

    % Transformation matrix TO the rotated system  
    T(1,1) = nx^2;
    T(1,2) = ny^2;
    T(1,3) = nz^2;
    T(1,4) = 2*nz*ny;
    T(1,5) = 2*nz*nx;
    T(1,6) = 2*ny*nx;
    T(2,1) = sx^2;
    T(2,2) = sy^2;
    T(2,3) = sz^2;
    T(2,4) = 2*sz*sy;
    T(2,5) = 2*sz*sx;
    T(2,6) = 2*sy*sx;
    T(3,1) = tx^2;
    T(3,2) = ty^2;
    T(3,3) = tz^2;
    T(3,4) = 2*tz*ty;
    T(3,5) = 2*tz*tx;
    T(3,6) = 2*ty*tx;
    T(4,1) = sx*tx;
    T(4,2) = sy*ty;
    T(4,3) = sz*tz;
    T(4,4) = sz*ty+sy*tz;
    T(4,5) = sz*tx+sx*tz;
    T(4,6) = sy*tx+sx*ty;
    T(5,1) = nx*tx;
    T(5,2) = ny*ty;
    T(5,3) = nz*tz;
    T(5,4) = nz*ty+ny*tz;
    T(5,5) = nz*tx+nx*tz;
    T(5,6) = ny*tx+nx*ty;
    T(6,1) = nx*sx;
    T(6,2) = ny*sy;
    T(6,3) = nz*sz;
    T(6,4) = nz*sy+ny*sz;
    T(6,5) = nz*sx+nx*sz;
    T(6,6) = ny*sx+nx*sy;
    for i = 1:6 
    for j = 1:6
    TT(i,j) = T(j,i);
    end
    end

    % Initialize lower half of Voigt matrix using the symmetry property

    for i = 1:6 
    for j = 1:i-1
    ce(i,j) = ce(j,i);
    end
    end
    % 
    % % Initialize and rotate Voigt matrix to get local material properties

    Voigt_rot(:,:) =  T(:,:)*(ce(:,:)*TT(:,:));

    ce(:,:) = Voigt_rot(:,:);

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A = [0       0       0       0       0       0     0    0       0    ;
     0       0       0       0       0       0       0           0       0    ;
     0       0       0       0       0       0       0           0       0    ;
     0       0       0       0       0       0       0           0       0    ;
     0       0       0       0       0       0       0           0       0    ;
     0       0       0       0       0       0       0           0       0    ;
     -1/rho  0       0       0       0       0       0           0       0    ;
     0       0       0      -1/rho   0       0       0           0       0    ;
     0       0       0       0       0    -1/rho     0           0       0       ];

    A(1,7:9) = [ -ce(1,1), -ce(1,6), -ce(1,5) ];
    A(2,7:9) = [ -ce(1,2), -ce(2,6), -ce(2,5) ];
    A(3,7:9) = [ -ce(1,3), -ce(3,6), -ce(3,5) ];
    A(4,7:9) = [ -ce(1,6), -ce(6,6), -ce(5,6) ];
    A(5,7:9) = [ -ce(1,4), -ce(4,6), -ce(4,5) ];
    A(6,7:9) = [ -ce(1,5), -ce(5,6), -ce(5,5) ];

    %% Build Jacobian matrix B
    B = [0       0       0       0       0       0       0           0      0;
         0       0       0       0       0       0       0                0       0;
         0       0       0       0       0       0       0                0       0;
         0       0       0       0       0       0       0                0        0;
         0       0       0       0       0       0       0                0        0;
         0       0       0       0       0       0       0                0        0;
         0       0       0      -1/rho   0       0       0                0        0;
         0     -1/rho    0       0       0       0       0                0        0;
         0       0       0       0    -1/rho     0       0                0       0 ];
     
     B(1,7:9) = [ -ce(1,6), -ce(1,2), -ce(1,4) ];
     B(2,7:9) = [ -ce(2,6), -ce(2,2), -ce(2,4) ];
     B(3,7:9) = [ -ce(3,6), -ce(2,3), -ce(3,4) ];
     B(4,7:9) = [ -ce(6,6), -ce(2,6), -ce(4,6) ];
     B(5,7:9) = [ -ce(4,6), -ce(2,4), -ce(4,4) ];
     B(6,7:9) = [ -ce(5,6), -ce(2,5), -ce(4,5) ];
     
    %% Build Jacobian matrix C
    C = [0       0       0       0       0       0       0        0        0;
         0       0       0       0       0       0       0        0        0;
         0       0       0       0       0       0       0        0        0;
         0       0       0       0       0       0       0        0        0;
         0       0       0       0       0       0       0        0        0;
         0       0       0       0       0       0       0        0        0;
         0       0       0       0       0    -1/rho     0        0        0; 
         0       0       0       0     -1/rho    0       0        0        0; 
          0       0    -1/rho     0       0       0       0        0        0 ];
      
    C(1,7:9) = [ -ce(1,5), -ce(1,4), -ce(1,3) ];
    C(2,7:9) = [ -ce(2,5), -ce(2,4), -ce(2,3) ];
    C(3,7:9) = [ -ce(3,5), -ce(3,4), -ce(3,3) ];
    C(4,7:9) = [ -ce(5,6), -ce(4,6), -ce(3,6) ];
    C(5,7:9) = [ -ce(4,5), -ce(4,4), -ce(3,4) ];
    C(6,7:9) = [ -ce(5,5), -ce(4,5), -ce(3,5) ];

elseif(rheo==4)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% PARAMETER INPUT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    fname=input('     Enter the file name of the poroelastic properties (*.def):','s')
    fid = fopen([fname,'.def']);
      junk = fgetl(fid);
      data = fscanf(fid,'%g',[1,43]);
    fclose(fid);
    disp(['     Successfully read parameters from file: ',fname,'.def!'])
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %coefficients given in N/m^2 
    c(1,1)     = data(3);
    c(1,2)     = data(4);
    c(1,3)     = data(5);
    c(1,4)     = data(6);
    c(1,5)     = data(7);
    c(1,6)     = data(8);
    c(2,2)     = data(9);  
    c(2,3)     = data(10);
    c(2,4)     = data(11);
    c(2,5)     = data(12);
    c(2,6)     = data(13);
    c(3,3)     = data(14);
    c(3,4)     = data(15);
    c(3,5)     = data(16);
    c(3,6)     = data(17);
    c(4,4)     = data(18);
    c(4,5)     = data(19);
    c(4,6)     = data(20);
    c(5,5)     = data(21);
    c(5,6)     = data(22);
    c(6,6)     = data(23);
    rho_S      = data(2);
    rho_F      = data(24);
    K_F        = data(25);
    nu         = data(26);
    K_S        = data(27);
    Poro       = data(28);
    Kappa(1)   = data(29);
    Kappa(2)   = data(30);
    Kappa(3)   = data(31);
    Tor(1)     = data(32);
    Tor(2)     = data(33);
    Tor(3)     = data(34);
    rho        = rho_S * (1-Poro)+Poro*rho_F;
    K_Mean     = 1/9*(c(1,1)+c(2,2)+c(3,3)+2*(c(1,2)+c(1,3)+c(2,3)));
    MM         = K_S/((1-K_Mean/K_S)-Poro*(1-K_S/K_F));
    Alpha(1)   = 1 - (c(1,1)+c(1,2)+c(1,3))/(3*K_S);
    Alpha(2)   = 1 - (c(1,2)+c(2,2)+c(2,3))/(3*K_S);
    Alpha(3)   = 1 - (c(1,3)+c(2,3)+c(3,3))/(3*K_S);
    Alpha(4)   = -(c(1,4)+c(2,4)+c(3,4))/(3*K_S);
    Alpha(5)   = -(c(1,5)+c(2,5)+c(3,5))/(3*K_S);
    Alpha(6)   = -(c(1,6)+c(2,6)+c(3,6))/(3*K_S);
    Rho1(1:3)  = rho-(rho_F^2 ./ (rho_F * Tor(1:3)  /Poro ));
    Rho2(1:3)  = rho_F - (rho_F * Tor(1:3) / Poro) * rho / rho_F;
    Beta1(1:3) = rho_F ./ (rho_F * Tor(1:3) / Poro) ;
    Beta2(1:3) = rho / rho_F;

    %initialize lower triangular matrix symmetrically
    for i=1:5
        for j=i+1:6
            c(j,i) = c(i,j);
        end
    end

    % direction of material symmetry axes
    nx = data(35);
    ny = data(36);
    nz = data(37);
    sx = data(38);
    sy = data(39);
    sz = data(40);
    tx = data(41);
    ty = data(42);
    tz = data(43);

    % normalization of coordinate vectors
    len = sqrt(nx*nx+ny*ny+nz*nz);   nx=nx/len;  ny=ny/len;  nz=nz/len;
    len = sqrt(sx*sx+sy*sy+sz*sz);   sx=sx/len;  sy=sy/len;  sz=sz/len;
    len = sqrt(tx*tx+ty*ty+tz*tz);   tx=tx/len;  ty=ty/len;  tz=tz/len;

    % rotation matrix of anisotropic coefficients  (Bond matrix)       
    BM(1,1) = nx^2;      
    BM(1,2) = sx^2; 
    BM(1,3) = tx^2;
    BM(1,4) = 2*sx*tx;
    BM(1,5) = 2*nx*tx;
    BM(1,6) = 2*nx*sx;
    BM(2,1) = ny^2;
    BM(2,2) = sy^2;
    BM(2,3) = ty^2;
    BM(2,4) = 2*sy*ty;
    BM(2,5) = 2*ny*ty;
    BM(2,6) = 2*ny*sy;
    BM(3,1) = nz^2;
    BM(3,2) = sz^2;
    BM(3,3) = tz^2;
    BM(3,4) = 2*sz*tz;
    BM(3,5) = 2*nz*tz;
    BM(3,6) = 2*nz*sz;
    BM(4,1) = nz*ny;
    BM(4,2) = sz*sy;
    BM(4,3) = tz*ty;
    BM(4,4) = sz*ty+sy*tz;
    BM(4,5) = nz*ty+ny*tz;
    BM(4,6) = nz*sy+ny*sz;
    BM(5,1) = nz*nx;
    BM(5,2) = sz*sx;
    BM(5,3) = tz*tx;
    BM(5,4) = sz*tx+sx*tz;
    BM(5,5) = nz*tx+nx*tz;
    BM(5,6) = nz*sx+nx*sz;
    BM(6,1) = ny*nx;
    BM(6,2) = sy*sx;
    BM(6,3) = ty*tx;
    BM(6,4) = sy*tx+sx*ty;
    BM(6,5) = ny*tx+nx*ty;
    BM(6,6) = ny*sx+nx*sy;

    %rotation of anisotropic coefficient c into symmetry axes system 
    c = BM*c*BM';

    %computation of undrained c(i,j) coefficients
    for i=1:6
        for j=1:6
            c(i,j)=c(i,j)+MM*Alpha(i)*Alpha(j);
        end
    end



    %% Build Jacobian matrix A
    A = [ 0 0 0      0 0      0 -c(1,1) -c(1,6) -c(1,5) 0 -Alpha(1)*MM 0 0;
          0 0 0      0 0      0 -c(1,2) -c(2,6) -c(2,5) 0 -Alpha(2)*MM 0 0;
          0 0 0      0 0      0 -c(1,3) -c(3,6) -c(3,5) 0 -Alpha(3)*MM 0 0;
          0 0 0      0 0      0 -c(1,6) -c(6,6) -c(5,6) 0 -Alpha(6)*MM 0 0;
          0 0 0      0 0      0 -c(1,4) -c(4,6) -c(4,5) 0 -Alpha(4)*MM 0 0;
          0 0 0      0 0      0 -c(1,5) -c(5,6) -c(5,5) 0 -Alpha(5)*MM 0 0;
     -1/Rho1(1) 0 0  0 0      0      0      0      0 -Beta1(1)/Rho1(1) 0 0 0;
          0 0 0 -1/Rho1(1) 0   0      0      0      0    0       0      0 0;
          0 0 0      0 0 -1/Rho1(1)   0      0      0    0       0      0 0;
          0 0 0      0 0      0 Alpha(1)*MM Alpha(6)*MM Alpha(5)*MM 0 MM 0 0 ; 
      -1/Rho2(1) 0 0  0 0      0      0      0      0 -Beta2(1)/Rho2(1) 0 0 0;
          0 0 0 -1/Rho2(1) 0   0      0      0      0   0       0       0 0;
          0 0 0      0 0 -1/Rho2(1)   0      0      0   0       0       0 0];


    %% Build Jacobian matrix B
    B = [ 0      0 0 0      0      0 -c(1,6) -c(1,2) -c(1,4) 0 0 -Alpha(1)*MM 0;
          0      0 0 0      0      0 -c(2,6) -c(2,2) -c(2,4) 0 0 -Alpha(2)*MM 0;
          0      0 0 0      0      0 -c(3,6) -c(2,3) -c(3,4) 0 0 -Alpha(3)*MM 0;
          0      0 0 0      0      0 -c(6,6) -c(2,6) -c(4,6) 0 0 -Alpha(6)*MM 0;
          0      0 0 0      0      0 -c(4,6) -c(2,4) -c(4,4) 0 0 -Alpha(4)*MM 0;
          0      0 0 0      0      0 -c(5,6) -c(2,5) -c(4,5) 0 0 -Alpha(5)*MM 0;
          0      0 0 -1/Rho1(2) 0      0      0      0      0  0 0    0       0;
          0 -1/Rho1(2) 0 0      0      0      0      0 0 -Beta1(2)/Rho1(2)  0 0 0;
          0      0 0 0      -1/Rho1(2) 0      0      0      0  0 0    0       0;
          0      0 0 0      0      0 Alpha(6)*MM Alpha(2)*MM Alpha(4)*MM 0 0 MM 0;
          0      0 0 -1/Rho2(2) 0      0      0      0      0  0   0   0   0;
          0 -1/Rho2(2) 0 0      0      0      0      0 0 -Beta2(2)/Rho2(2) 0 0 0;
          0      0 0 0      -1/Rho2(2) 0      0      0      0  0  0  0  0];

    %% Build Jacobian matrix C
    C = [0 0      0 0      0      0 -c(1,5) -c(1,4) -c(1,3) 0 0 0 -Alpha(1)*MM;
         0 0      0 0      0      0 -c(2,5) -c(2,4) -c(2,3) 0 0 0 -Alpha(2)*MM;
         0 0      0 0      0      0 -c(3,5) -c(3,4) -c(3,3) 0 0 0 -Alpha(3)*MM;
         0 0      0 0      0      0 -c(5,6) -c(4,6) -c(3,6) 0 0 0 -Alpha(6)*MM;
         0 0      0 0      0      0 -c(4,5) -c(4,4) -c(3,4) 0 0 0 -Alpha(4)*MM;
         0 0      0 0      0      0 -c(5,5) -c(4,5) -c(3,5) 0 0 0 -Alpha(5)*MM;
         0 0      0 0      0 -1/Rho1(3)      0      0      0    0 0 0   0; 
         0 0      0 0 -1/Rho1(3)      0      0      0      0    0 0 0   0; 
         0 0 -1/Rho1(3) 0      0      0      0      0      0    -Beta1(3)/Rho1(3) 0 0 0;
         0 0      0 0      0      0 Alpha(5)*MM Alpha(4)*MM Alpha(3)*MM 0 0 0 MM;
         0 0      0 0      0 -1/Rho2(3)      0      0      0    0 0  0  0; 
         0 0      0 0 -1/Rho2(3)      0      0      0      0    0 0  0  0; 
         0 0 -1/Rho2(3) 0      0      0      0      0      0    -Beta2(3)/Rho2(3) 0 0 0];
      
end
%%%%%%%%%%%%%%%%%%%%%%% Compute analytic solution %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Formulate and solve Eigenvalue problem
clear i     %to make sure to have i as the imaginary i

M1   = -(A*kx + B*ky + C*kz);

[R1,r1] = eig(M1);     % R1     = matrix of eigenvectors
omega1  = diag(r1);                % omega1 = vector of eigenvalues

[omega,Y] = sort(real(omega1),'descend');
R=R1*0;
for j=1:length(omega)
    R(:,j)=real(R1(:,Y(j)));
end
%% Writing output data
disp(' ');
disp(sprintf('Writing data into file...'));
fid_out  = fopen([fname_out,'.dat'],'w');
    fprintf(fid_out,'Number of eigenvalues\n');
    fprintf(fid_out,'%d\n',length(omega));    
    fprintf(fid_out,'Eigenvalues\n');
    fprintf(fid_out,'%25.14f\n',[real(omega)]');
    fprintf(fid_out,'Eigenvectors\n');
    fprintf(fid_out,'%25.14f\n',[real(reshape(R,length(omega)^2,1))]');
fclose(fid_out);

disp(sprintf('Finished! '));                                  

