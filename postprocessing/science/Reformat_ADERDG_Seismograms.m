%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2006, SeisSol Group
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
% Reads in TECPLOT format of SeisSol's seismogram output and converts it.
home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  Reformat_ADER-DG_Seismograms reads in the TECPLOT format    %%')
disp('     %%  of the typical SEISSOL seismogram output and converts       %%')
disp('     %%  it into a 3-dimensional tensor in MATLAB format.            %%')
disp('     %%                                                              %%')
disp('     %%  Dimension 1 has the length of the number of time samples    %%')
disp('     %%  Dimension 2 has the length of the number of variables + 1   %%')
disp('     %%              (in general this is 10: 1 time     component    %%')
disp('     %%                                      6 stress   components   %%')
disp('     %%                                      3 velocity components   %%')
disp('     %%  Dimension 3 has the length of the number of receivers       %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;
filename = input('     Specify root-filename (typically: *file*-receiver-...):  ','s');
ndt      = input('     Give number of time samples of the data:              ');
nseis    = input('     Give number of seismograms:                           ');
nvar     = input('     Give number of variables (including time):            ');

disp(' '),  disp('     Creating seismogram-processor relation ...' )
eval(['!ls -l ',filename,'-receiver* > tmp.dat']);
tmp_file = 'tmp.dat';
fid   = fopen(tmp_file);
  for i=1:nseis
    junk  = fgetl(fid); proc_s(i,1:25)=junk(end-24:end); 
  end
fclose(fid);
eval(['!rm tmp.dat']);

num = 1;
ind = 0;

for num = 1:nseis
     
    in_file = [filename,proc_s(num,:)];
    
    fid   = fopen(in_file);
    junk  = fgetl(fid); junk  = fgetl(fid); 
    junk  = fscanf(fid,'%c',[6,1]);
    x1    = fscanf(fid,'%g',[1,1]);
    junk  = fscanf(fid,'%c',[6,1]);
    x2    = fscanf(fid,'%g',[1,1]);
    junk  = fscanf(fid,'%c',[6,1]);
    x3    = fscanf(fid,'%g',[1,1]);
    disp(sprintf('X: %g   Y: %g  Z: %g',x1,x2,x3));
    data = fscanf(fid,'%g',[nvar,ndt]); data = data';
    
    %determine and delete double samples from a possible restart
    data(:,1) = (round(data(:,1)*1000000))/1000000;
    [tmp,k,kk] = unique(data(:,1));
    data  = data(k,:);
    
    disp(sprintf('\n Samples in seismogram nr. %i\t:   %i', num,size(data,1)));
    MATLAB_OUTPUT(:,:,num) = data;
    fclose(fid);
      
end

disp(' '),  disp('     Finished conversion!'), disp(' ')
out_filename = input('     Give filename to save MATLAB data:  ','s');
out_file = out_filename(1:end-4);
eval([out_file,'= MATLAB_OUTPUT;']);
eval(['save ',out_filename,'  ',out_file]);
disp(sprintf('    \n Saved data as :  %s\n',out_filename));
