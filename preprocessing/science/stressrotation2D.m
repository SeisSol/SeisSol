%%
% @file
% This file is part of SeisSol.
%
% @author Amaryllis Nerger (amaryllis.nerger AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/anerger)
% @author Christian Pelties (pelties AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/pelties)
%
% @section LICENSE
% Copyright (c) 2012-2013, SeisSol Group
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
% 2D Rotation of the fault's stress tensor into the coordinatesystem of the main
% fault

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
disp('- Shear traction right lateral positive');
disp('- Shear traction left lateral negative');
disp(' ');
disp('Script rotates values in fault coordinate system into global xyz coordinate system values');
disp(' ');

% Input normal stresses
disp('Input normal stresses');
sxx = input('sxx=');
syy = input('syy=');
disp(' ');

% Input shear stress
disp('Input shear traction')
sxy = input('sxy=');
disp(' ');

disp('Input angle: If you transform your values counter clockwise, use a plus sign');
disp('             If you transform your values clockwise, use a minus sign.');

theta = input('theta=');
disp(' ');

% stress tensor in fault coordinate system                                                                                                                        Koordinatensystem des Branch
A=[sxx sxy;
   sxy syy];

% rotary tensor
T=[cos(theta*pi/180) -sin(theta*pi/180);
   sin(theta*pi/180) cos(theta*pi/180)];
   
% stress tensor in coordinatesystem of main fault
Stress_in_global_xyz=T*A*T'
