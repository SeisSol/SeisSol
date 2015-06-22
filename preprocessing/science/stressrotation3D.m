%%
% @file
% This file is part of SeisSol.
%
% @author Alice Gabriel (gabriel AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/gabriel)
%
% @section LICENSE
% Copyright (c) 2013, SeisSol Group
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
% 3D Rotation of the fault's stress tensor into the coordinatesystem of the main
% fault; modified from stressrotation2D

% convention:
disp('Right-handed coordinate system required');
disp('Please, make sure to follow this convention during model/mesh generation!');
disp('');
disp('z      y          free-surface     North');
disp('|    /                      |    /');
disp('|  /                        |  /');
disp('| ---> x            =       | --->  East');
disp('-z                        depth');
disp(' ');
disp('- Normal stress negative in compression');
disp('- Shear traction in strike direction right lateral: positive');
disp('- Shear traction in strike direction left lateral: negative');
disp('- Shear traction along dip direction: positive');
disp('- Shear traction anti dip direction: negative');
disp(' ');
disp('Script rotates values in fault coordinate system into global xyz coordinate system values');
disp(' ');

% Input normal stresses
disp('Input normal stresses');
sxx = input('sxx=');            %example TPV10 at z=-1m: 0
syy = input('syy=');            %example TPV10: -7378
disp(' ');

% Input shear stress
disp('Input shear traction')
sxy = input('sxy=');            %example TPV10: 0
syz = input('syz=');            %example TPV10: 0.55*7378
szx = input('szx=');            %example TPV10: 0
disp(' ');

disp('Input angles: If you transform your values counter clockwise, use a plus sign');
disp('             If you transform your values clockwise, use a minus sign.');

disp(' ');
disp('Input rotation angles about each axis')
thetax = input('theta x=');     %example TPV10: -30
thetay = input('theta y=');     %example TPV10: 0
thetaz = input('theta z=');     %example TPV10: 0

% stress tensor in fault coordinate system                                                                                                                        Koordinatensystem des Branch
A=[sxx sxy szx;
   sxy syy syz;
   szx syz 0];

% rotary tensor for rotation around x-axis
Tx=[1 0 0;
   0 cos(thetax*pi/180) -sin(thetax*pi/180) ;
   0 sin(thetax*pi/180) cos(thetax*pi/180) ];


% rotary tensor for rotation around y-axis
Ty=[cos(thetay*pi/180)    0   -sin(thetay*pi/180);
   0 1 0;
   sin(thetay*pi/180) 0 cos(thetay*pi/180) ];


% rotary tensor for rotation around z-axis
Tz=[cos(thetaz*pi/180) -sin(thetaz*pi/180) 0;
   sin(thetaz*pi/180) cos(thetaz*pi/180) 0
   0 0 1];

% stress tensor in coordinatesystem of main fault
% chain of rotations
A_tempx=Tx*A*Tx';
A_tempy=Ty*A_tempx*Ty';
Stress_in_global_xyz=Tz*A_tempy*Tz'
