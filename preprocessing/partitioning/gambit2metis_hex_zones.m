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

home;
disp(' ')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('     %%                                                             %%')
disp('     %% Gambit2Metis converts a GAMBIT-meshfile of                  %%')
disp('     %% hexahedrons stored as "filename.neu"                        %%')
disp('     %% into a METIS-file called "filename.met" where different     %%')
disp('     %% zones can be partitioned individually!                      %%')
disp('     %% This METIS-file can be used as input for the                %%')
disp('     %% METIS mesh partitioner:                                     %%')
disp('     %% http://www-users.cs.umn.edu/~karypis/metis/metis/index.html %%')
disp('     %% The METIS mesh partitioner uses the METIS-file through:     %%')
disp('     %%               partnmesh filename.met n                      %%')
disp('     %% where n specifies the number of desired partitions.         %%')
disp('     %%                                                             %%')
disp('     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' '),disp(' ')
clear, close all;
filename = input('     Specify GAMBIT-filename:   ','s');
tmp = 1;
zone_separate=[];
while tmp~=0
    tmp = input('     Specify zone number to partition separately or 0 to continue):   ');
    zone_separate = [zone_separate,tmp];
end
zone_separate(end) = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Read Gambit Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t = cputime;

fid   = fopen([filename,'.neu']);
    disp(' ')
    disp('-----------------------------------------------------------------------------------')
    disp(sprintf('\t Reading data from: \t\t%s',[filename,'.neu']));
    for i= 1:6
        junk  = fgetl(fid);
    end
    tmp    = fscanf(fid,'%i',[5,1]); 
    NUMNP  = tmp(1);                     %Number of Vertices
    NELEM  = tmp(2);                     %Number of Elements 
    NGRPS  = tmp(3);                     %Number of Element Groups
    NBSETS = tmp(4);                     %Number of Boundary Condition Sets
    NDFCD  = tmp(5);                     %Number of Dimensions (2 or 3)

    disp(sprintf('\t Reading \t\t\t\t\t%i vertices',NUMNP));
    disp(sprintf('\t Reading \t\t\t\t\t%i tetrahedrons',NELEM));
    junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
    X     = fscanf(fid,'%g',[4,NUMNP]);                            %Read Vertices
    junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
    hexa  = fscanf(fid,'%i',[11,NELEM]);                           %Read Hexahedrons
             %Check if the mesh is a pure hexahedral mesh
             if ( isempty(find(hexa(2,:)-6~=0))==0 & isempty(find(hexa(3,:)-8~=0))==0 )
                 disp(sprintf('\n\t ##ERROR: The GAMBIT input mesh is not purely hexahedral\n\n ##'))
                 return;
             end
    hexa(1:3,:) = []; 
    
    %read zones
    for ng = 1:NGRPS
         junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid); 
         eing(ng)   = str2num(junk(29:39));  junk  = fgetl(fid); junk  = fgetl(fid);
         nlines     = floor(eing/10); 
         nrest      = mod(eing,10);
         tmp        = fscanf(fid,'%i',[eing(ng),1]);
         zone(1:eing(ng),ng) = tmp;
    end
fclose(fid);
disp(sprintf('\t Read input successfully!'));
disp('-----------------------------------------------------------------------------------')

zone_non_separate = (1:1:NGRPS);
zone_non_separate(zone_separate) = [];
disp(sprintf('\t Separately partitioned zones are:'));
disp(zone_separate)
disp(sprintf('\t Jointly partitioned zones are:'));
disp(zone_non_separate)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Writing Metis Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Writing the separately partitioned zones
outstore = cell(length(zone_separate)+1,1);
out = [];
ns = length(zone_separate);
for i = 1:ns
    
    len    = eing(zone_separate(i));
    disp(sprintf('\t Write METIS-file: \t\t\t%s',[filename,num2str(i),'.met']));

    fid = fopen([filename,num2str(i),'.met'],'w');  
                fprintf(fid,'%12i %5i \n',len,3);
                hexa_tmp = hexa(:,zone(1:len,zone_separate(i)))';
                outstore{i} = hexa_tmp;
                tmp = reshape(hexa_tmp,len*8,1);
                [tmp2,j,k] = unique(tmp);
                newindex = (1:1:length(tmp2))';
                tmp = newindex(k);
                output = reshape(tmp,len,8);
                fprintf(fid,'%12i %12i %12i %12i \n',output');
    fclose(fid);
    out = [out;zone(1:len,zone_separate(i))];
end

%  Writing the remaining zones that are partitioned together
lentot  = sum(eing(zone_non_separate));
disp(sprintf('\t Write METIS-file: \t\t\t%s',[filename,num2str(ns+1),'.met']));

fid = fopen([filename,num2str(ns+1),'.met'],'w');
         fprintf(fid,'%12i %5i \n',lentot,3);
         hexa_tmp = hexa';
         hexa_tmp(out,:) = [];
         outstore{ns+1} = hexa_tmp;
         tmp = reshape(hexa_tmp,lentot*8,1);
         [tmp2,j,k] = unique(tmp);
         newindex = (1:1:length(tmp2))';
         tmp = newindex(k);
         output = reshape(tmp,lentot,8);
         fprintf(fid,'%12i %12i %12i %12i \n',output'); 
fclose(fid);

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t Conversion finished successfully!  (%g CPU sec)\n',cputime-t));
disp('-----------------------------------------------------------------------------------')
disp(' ')
pause(0.2);
proc = input('     Specify number of processor of the partition:   ','s');
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t METIS partition of each of the %i .met files starts!\n',i));
disp('-----------------------------------------------------------------------------------')
disp(' ')
for i = 1:length(zone_separate)+1
    eval(['!partnmesh ',filename,num2str(i),'.met ',proc])
    pause(0.5);
end
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t METIS partition done!'));
disp('-----------------------------------------------------------------------------------')
disp(' ')

met = cell(length(zone_separate)+1,1);

for i = 1:length(zone_separate)+1
    eval(['load ',filename,num2str(i),'.met.epart.',proc]);     
    eval(['met{',num2str(i),'}=',filename,num2str(i),';']);
    [tmp,i1,i2] = intersect(hexa',outstore{i},'rows');
    epart_final(i1,1) = met{i}(i2);
end

%  Writing the final .met file with separately partitioned zones
%  merged together
fid = fopen([filename,'.met.epart.',num2str(proc)],'w');
    fprintf(fid,'%3i \n',epart_final);
fclose(fid);

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Wrote finally merged METIS-file: \t\t%s',[filename,'.met.epart.',num2str(proc)]));
disp('-----------------------------------------------------------------------------------')
disp(' ')
for i = 1:length(zone_separate)+1
    eval(['!rm ',filename,num2str(i),'.met.epart.',proc])
    eval(['!rm ',filename,num2str(i),'.met'])
    eval(['!rm ',filename,num2str(i),'.met.npart.',proc])
end
