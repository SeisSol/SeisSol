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

home;
disp(' ')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%                CONVERT_FAULTRECEIVERS               %%')
disp('    %%            from fault COS to absolute COSIONS       %%')
disp('    %%                                                     %%')
disp('    %%           addition to place_faultreceivers.m        %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')
disp(' Transform (along dip, along strike) 2D receiver coordinates ')
disp(' along arbitrary shaped faults into xyz coordinates ')
disp(' input: receiver dip and strike coordinates, fault dip and strike ')
disp(' output: receiver xyz coordinates')
disp(' '),disp(' ')
%
clear, close all;
format long
%read in receiver file
dip = input('   Dip of fault plane (angle between the fault and a horizontal plane, 0° to 90):  ','s');
dip = str2double(dip);
strike = input('   Strike of fault plane (direction of a line created by the intersection of a fault plane and a horizontal surface, 0° to 360°, relative to North:  ','s');
strike=str2double(strike);
rec_filename = input('   Filename of 2D receiver locations (along dip, along strike) on the fault (suffix ".dat" is appended):  ','s');
%
% load receiver stations
%
eval(['load ',rec_filename,'.dat']);
eval(['st = ',rec_filename,';']);
%
%convert dip
%
recs(:,3)=st(:,1).*sin(degtorad(dip));
recs(:,2)=st(:,1).*cos(degtorad(dip));
%
%convert strike
%
recs(:,1)=st(:,2).*cos(degtorad(strike));
disp('    Receiver coordinates:'); disp(' ')
disp(recs);
%
%write recs to file
%
fid_out  = fopen([rec_filename ,'_coordinates.dat'],'w');
fprintf(fid_out,'%20.12f%20.12f%20.12f\n',receivers');
fclose(fid_out);
disp('    Receiver coordinates saved!');
%
