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
disp('     %% Gambit2Metis converts a GAMBIT-meshfile of                  %%')
disp('     %% triangles or tetrahedrons stored as "filename.neu"          %%')
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
filename = input('     Specify GAMBIT-filename:   ','s');

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
    NELEM  = tmp(2);                     %Number of Elements (Triangles or Tetrahedrons)
    %NGRPS  = tmp(3);                    %Number of Element Groups
    %NBSETS = tmp(4);                    %Number of Boundary Condition Sets
    NDFCD  = tmp(5);                     %Number of Dimensions (2 or 3)

    switch NDFCD                         %Choose between 2-D or 3-D mesh input
        case 2
            disp(sprintf('\t Reading \t\t\t\t\t%i vertices',NUMNP));
            disp(sprintf('\t Reading \t\t\t\t\t%i triangles',NELEM));
            junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
            X     = fscanf(fid,'%g',[3,NUMNP]);                            %Read Vertices
            junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
            tri   = fscanf(fid,'%i',[6,NELEM]);                            %Read Triangles
            %Check if the mesh is a pure triangular mesh
             if ( isempty(find(tri(2:3,:)-3~=0)) == 0 )
                 disp(sprintf('\n\t ##ERROR: The GAMBIT input mesh is not purely triangular\n\n ##'))
                 return;
             end
            %reshape connectivity matrix for output
            tri   = tri(4:6,:);
            
        case 3
            disp(sprintf('\t Reading \t\t\t\t\t%i vertices',NUMNP));
            disp(sprintf('\t Reading \t\t\t\t\t%i tetrahedrons',NELEM));
            junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
            X     = fscanf(fid,'%g',[4,NUMNP]);                            %Read Vertices
            junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
            for i = 1:NELEM 
               prefix(1:3) = fscanf(fid,'%i',[3]);
               if (prefix(3)==4)                    
                 tetra(1:4,i) = fscanf(fid,'%i',[4]);
               elseif(prefix(3)==8)
                 tetra(1:7,i) = fscanf(fid,'%i',[7]);
                 tetra(8,i)   = fscanf(fid,'%i',[1]);
               else
                 disp(sprintf('\n\t ##ERROR: The GAMBIT input mesh is neither tetrahedral nor hexahedral\n\n ##'))
                 return;                   
               end
            end
            %tetra = fscanf(fid,'%i',[7,NELEM]);                           %Read Tetrahedrons
            %Check if the mesh is a pure tetrahedral mesh
            % if ( isempty(find(tetra(2,:)-6~=0))==0 & isempty(find(tetra(3,:)-4~=0))==0 )
            %     disp(sprintf('\n\t ##ERROR: The GAMBIT input mesh is not purely tetrahedral\n\n ##'))
            %     return;
            % end
            %reshape connectivity matrix for output
            %tetra = tetra(4:7,:);
            
        otherwise
            disp(sprintf('\n\t ##ERROR: The GAMBIT input file has to be 2- or 3-dimensional\n\n ##'))
            return;
            
    end
fclose(fid);
disp(sprintf('\t Read input successfully!'));
disp('-----------------------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Writing Metis Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(sprintf('\t Write METIS-file: \t\t\t%s',[filename,'.met']));

fid = fopen([filename,'.met'],'w');
    switch NDFCD
        case 2
            fprintf(fid,'%12i %5i \n',NELEM,1);
            fprintf(fid,'%12i %12i %12i \n',tri);
        case 3
            if(prefix(3)==4)
                ElemCode = 2
            else
                ElemCode = 3
            end
            fprintf(fid,'%12i %5i \n',NELEM,ElemCode);
            fprintf(fid,'%12i %12i %12i %12i \n',tetra);
    end
fclose(fid);
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t Conversion finished successfully!  (%g CPU sec)\n',cputime-t));
disp('-----------------------------------------------------------------------------------')
disp(' ')

