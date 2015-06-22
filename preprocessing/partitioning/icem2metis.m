%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2005, SeisSol Group
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
disp('     %% Icem2Metis converts a ICEM-meshfile (FIDAP format) of       %%')
disp('     %% tetrahedrons stored as "filename.neu"                       %%')
disp('     %% into a METIS-file called "filename.met".                    %%')           
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
filename = input('     Specify ICEM-CFD filename:   ','s');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Read Icem Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t        = cputime;
ZoneType = 4;
iZone    = 0;
tetra    = [];

fid   = fopen([filename,'.neu']);
    disp(' ')
    disp('-----------------------------------------------------------------------------------')
    disp(sprintf('\t Reading data from: \t\t%s',[filename,'.neu']));
    for i= 1:5
        junk  = fgetl(fid);
    end
    tmp    = fscanf(fid,'%i',[5,1]); 
    NUMNP  = tmp(1);                     %Number of Vertices
    NELEM  = tmp(2);                     %Number of Elements (Triangles or Tetrahedrons)
    NGRPS  = tmp(3);                     %Number of Element Groups
    %NBSETS = tmp(4);                     %Number of Boundary Condition Sets
    NDFCD  = 3;                          %Number of Dimensions (only 3)

    switch NDFCD                         %Choose between 2-D or 3-D mesh input
        %case 2
            % is not implemented
        case 3
            disp(sprintf('\t Reading \t\t\t\t\t%i vertices\n',NUMNP));
            for i= 1:8
                junk  = fgetl(fid);
            end
            X     = fscanf(fid,'%g',[4,NUMNP]);                            %Read Vertices
            for i= 1:3
                junk  = fgetl(fid);
            end
            while(ZoneType==4)
                junk  = fgetl(fid);  junk = fgetl(fid); 
                ZoneType = str2num(junk(43:56));
                if(ZoneType==4)
                    iZone = iZone + 1;
                    nZoneElem = str2num(junk(26:36));
                    junk  = fgetl(fid);
                    disp(sprintf('\t Reading \t\t\t%i \t elements in zone %i: %s',nZoneElem,iZone,junk(13:end)));
                    tmp   = fscanf(fid,'%g',[5,nZoneElem]); 
                    
                    tetra = [tetra,tmp];
                else
                    break;
                end
            end
            
            %number of tetrahedral elements
            NELEM = size(tetra,2);
            
            %reshape connectivity matrix for output
            tetra = tetra(2:5,:);
            
        otherwise
            disp(sprintf('\n\t ##ERROR: The ICEM CFD input file has to be 3-dimensional\n\n ##'))
            return;
            
    end
fclose(fid);
disp(sprintf('\n\t Read input successfully!'));
disp('-----------------------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Writing Metis Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('\t Write METIS-file: \t\t\t%s',[filename,'.met']));

fid = fopen([filename,'.met'],'w');
    switch NDFCD
        case 2
            
        case 3
            fprintf(fid,'%12i %5i \n',NELEM,2);
            fprintf(fid,'%12i %12i %12i %12i \n',tetra);
    end
fclose(fid);

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t Conversion finished successfully!  (%g CPU sec)\n',cputime-t));
disp('-----------------------------------------------------------------------------------')
disp(' ')
pause(0.2);
proc = input('     Specify number of processor of the partition:   ','s');
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t METIS partition starts!\n'));
disp('-----------------------------------------------------------------------------------')
disp(' ')
%eval(['!partnmesh ',filename,num2str(i),'.met ',proc])
eval(['!mesh2dual ',filename,'.met'])
eval(['!pmetis ',filename,'.met.dgraph ',proc])
pause(0.5);
%clean up of files
eval(['!del ',filename,'.met'])
eval(['!del ',filename,'.met.dgraph'])    
eval(['!move ',filename,'.met.dgraph.part.',proc,'  ',filename,'.met.epart.',proc])

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t METIS partition done!'));
disp('-----------------------------------------------------------------------------------')
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Wrote final METIS-file: \t\t%s',[filename,'.met.epart.',num2str(proc)]));
disp('-----------------------------------------------------------------------------------')
disp(' ')


