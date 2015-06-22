%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
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
% Fixes problems in the TECPLOT seismograms caused by restarts.

home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                              %%')
disp('     %%  Restart_fix reads in the TECPLOT format of the pickpoint    %%')
disp('     %%  data of SEISSOL and fixes possible problems created         %%')
disp('     %%  through a restart.                                          %%')
disp('     %%                                                              %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')

clear, close all;
filename = input('     Specify root-filename (typically: file-pickpoints-):  ','s');
ndt      = input('     Give number of time samples of the data:              ');
nseis    = input('     Give number of seismograms:                           ');
nvar     = input('     Give number of variables (including time):            ');

disp(' '),  disp('     Creating seismogram-processor relation ...' )
eval(['!ls -l ',filename,'* > tmp.dat']);
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
    
    in_file  = [filename,proc_s(num,:)];
    out_file = ['fixed_',filename,proc_s(num,:)];

    fidin  = fopen(in_file);
    fidout = fopen(out_file,'w');

    for i =1:5
       junk = fgetl(fidin);
       fprintf(fidout,'%s\n',junk);
    end
    
    amin = -1;
    for i=1:ndt
       junk = fgetl(fidin);
       if(ischar(junk)==1)
         a = str2num(junk);
         if a(1)>amin
           if(nvar==13)
              fprintf(fidout,'%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n',a);
           elseif(nvar==10)
                 fprintf(fidout,'%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e%24.15e\n',a);            
           end
           amin = a(1);
         end
       else
         break;
       end
    end

    fclose(fidin);
    fclose(fidout);
    disp(in_file)
      
end

disp(' '), disp('Finished!'), disp(' ')
