%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2008, SeisSol Group
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

home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  Plot_SEISSOL_Snapshots reads in the TECPLOT format          %%')
disp('     %%  of the typical SEISSOL snapshot output and visualizes       %%')
disp('     %%  it on a color-coded triangulation surface plot.             %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;
filename = input('     Specify root-filename:                                ','s');
nsd      = input('     Give number of subdomain:                             ');
var      = input('     Specify variable number (1-5):                        ');
disp(' '),
disp('ATTENTION: Make sure to have 10 data columns in your .tri. file or change it in this matlab script!')
disp(' ')
disp(' '),  disp('     Creating seismogram-processor relation ...' )

for i=1:nsd
    
    %define character string of input files
    if(i<=10)
        in_file = [filename,'.000',num2str(i-1),'.tri.dat'];
    elseif(i>10 & i<=100)
        in_file = [filename,'.00',num2str(i-1),'.tri.dat'];
    else
        in_file = [filename,'.0',num2str(i-1),'.tri.dat'];
    end

    %open file and read data
    fid   = fopen(in_file);

        junk    = fgetl(fid); time = str2num(junk(24:end-1));
        junk    = fgetl(fid);
        junk    = fgetl(fid); nX = str2num(junk(9:20)); nE = str2num(junk(26:37));
        data    = fscanf(fid,'%g',[10,nX]); data    = data';    % change here if the .tri. file does not contain 10 columns!
        %extract vertices
        X       = data(:,1:2); data(:,1:2) = [];
        %extrace connectivity matrix
        con_mat = fscanf(fid,'%g',[3,nE]);  con_mat = con_mat';
        %visualize results 
        trisurf(con_mat,X(:,1),X(:,2),data(:,var)); hold on;
       
        disp(sprintf('     Visualizing subdomain %d',i));    

    fclose(fid);

    xlabel('x [m]'), ylabel('y [m]'), zlabel(['Variable ',num2str(var)]); view(0,90)
  
end
   
disp(sprintf('    \n Finished! \n'));
