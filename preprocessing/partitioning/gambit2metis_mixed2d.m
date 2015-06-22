%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
%
% @section LICENSE
% Copyright (c) 2009, SeisSol Group
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
disp('     %% Gambit2Metis_Mixed converts a GAMBIT-meshfile of            %%')
disp('     %% triangles and/or quadrilaterals stored as "filename.neu"    %%')
disp('     %% into a METIS-file called "filename.met" where different     %%')
disp('     %% zones are partitioned individually!                         %%')
disp('     %% This METIS-file can be used as input for the                %%')
disp('     %% METIS mesh partitioner:                                     %%')
disp('     %% http://www-users.cs.umn.edu/~karypis/metis/metis/index.html %%')
disp('     %% The METIS mesh partitioner uses the METIS-file through:     %%')
disp('     %%                                                             %%')
disp('     %%     mesh2dual   to convert the mesh to a graph              %%')
disp('     %%     pmetis      to partition the graph                      %%')
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
    NELEM  = tmp(2);                     %Number of Elements 
    NGRPS  = tmp(3);                     %Number of Element Groups
    NBSETS = tmp(4);                     %Number of Boundary Condition Sets
    NDFCD  = tmp(5);                     %Number of Dimensions (2 or 3)

    switch NDFCD                         %Choose between 2-D or 3-D mesh input
        case 2     %2D-meshes
            disp(sprintf('\t Reading \t\t\t\t\t%i vertices',NUMNP));
            disp(sprintf('\t Reading \t\t\t\t\t%i elements ...',NELEM));
            junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
            X     = fscanf(fid,'%g',[3,NUMNP]);  X(1,:)=[]; X=X';             %Read Vertices
            junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid);
            
            element_type = zeros(NELEM,1);
            connectivity = zeros(NELEM,4);
            
            for i=1:NELEM
                tmp             = fscanf(fid,'%i',[3,1]);            
                element_type(i) = tmp(3);
                connectivity(i,1:element_type(i)) = fscanf(fid,'%i',[element_type(i),1]);
            end
            tmp = find(element_type==3); n_tri  = length(tmp);
            tmp = find(element_type==4); n_quad = length(tmp);
            disp(sprintf('\t ... read \t\t\t\t\t%i triangles',n_tri));
            disp(sprintf('\t ... read \t\t\t\t\t%i quadrilaterals\n',n_quad));

            %read zones
            for ng = 1:NGRPS
                junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid); junk  = fgetl(fid); 
                %determine number of elements in this group (e_in_g)
                e_in_g(ng)   = str2num(junk(29:39));  junk  = fgetl(fid); junk  = fgetl(fid);
                tmp          = fscanf(fid,'%i',[e_in_g(ng),1]);
                zone(1:e_in_g(ng),ng) = tmp;
                grp_type(ng) = element_type(zone(1,ng));
                switch grp_type(ng)
                    case 3
                        disp(sprintf('\t Zone %i \t includes \t\t%i triangles.',ng,e_in_g(ng)));
                    case 4
                        disp(sprintf('\t Zone %i \t includes \t\t%i quadrilaterals.',ng,e_in_g(ng)));
                end
            end

        otherwise
            disp(sprintf('\t Currently only 2D meshes can be treated!'));
            
    end
fclose(fid);
disp(sprintf('\n\t Read input        successfully!'));
disp('-----------------------------------------------------------------------------------')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Writing Metis Data 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Writing the separately partitioned zones as input files for metis
for ng = 1:NGRPS
    
    disp(sprintf('\t Write METIS-file: \t\t\t%s',[filename,num2str(ng),'.met']));
    fid = fopen([filename,num2str(ng),'.met'],'w');
        len    = e_in_g(ng);
        switch grp_type(ng)
            case 3 %triangle
                output = connectivity(zone(1:e_in_g(ng),ng),1:grp_type(ng));
                fprintf(fid,'%12i %5i \n',len,1);            %Metis identifier for triangles is "1"
                fprintf(fid,'%12i %12i %12i \n',output');
            case 4 %quadrilateral
                output = connectivity(zone(1:e_in_g(ng),ng),1:grp_type(ng));
                fprintf(fid,'%12i %5i \n',len,4);             %Metis identifier for quadrilaterals is "4"
                fprintf(fid,'%12i %12i %12i %12i \n',output');
        end
    fclose(fid);

end

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t Conversion finished successfully!  (%g CPU sec)\n',cputime-t));
disp('-----------------------------------------------------------------------------------')
disp(' ')
pause(0.2);
n_proc = input('     Specify number of processor of the partition:   ','s');
disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t METIS partition of each of the %i .met files starts!\n',NGRPS));
disp('-----------------------------------------------------------------------------------')
disp(' ')
t = cputime;
for ng = 1:NGRPS
    %WARNING: partnmesh does not always create partitions with equal element number 
    %eval(['!partnmesh ',filename,num2str(i),'.met ',proc]) 
    %converting the mesh-file to a graph-file and partitioning the graph works better
    eval(['!mesh2dual ',filename,num2str(ng),'.met'])
    eval(['!pmetis ',filename,num2str(ng),'.met.dgraph ',n_proc])
    pause(0.5);
end
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Separate METIS partition done!  (%g CPU sec)',cputime-t'));
disp('-----------------------------------------------------------------------------------')
for ng = 1:NGRPS
   eval(['!mv ',filename,num2str(ng),'.met.dgraph.part.',n_proc,' ',filename,num2str(ng),'.met.epart.',n_proc]); 
end

%merging all metis files of different zones to one long vector in the second column
%attaching the order of the elements in the different zones in the first column 
met   = [];
proc  = str2num(n_proc);
n_e   = zeros(proc,NGRPS);
for ng = 1:NGRPS
    eval(['load ',filename,num2str(ng),'.met.epart.',n_proc]);     
    eval(['vec = ',filename,num2str(ng),';']);
    met   = [met;[zone(1:e_in_g(ng),ng),vec]];
    %check how many elements of each zone go to each processor
    for p = 1:proc
       tmp = find(vec==p-1);       %processor number in METIS starts at 0 and goes to n_proc-1  
       n_e(p,ng) = length(tmp);  
    end
end
%sort processor flags in the order of elements of the entire mesh
met = sortrows(met);
met = met(:,2);

disp(sprintf('\t Element distribution of each zones to processors!\n'));
disp(['       Zones:       ',num2str(1:NGRPS)]); disp(' ');
for p = 1:proc
    disp(['       Proc',num2str(p-1),':      ',num2str(n_e(p,:))]);
end
disp(' ');

%  Writing the final .met file with separately partitioned zones
%  merged together
fid = fopen([filename,'.met.epart.',n_proc],'w');
    fprintf(fid,'%5i \n',met');
fclose(fid);
%clean up of files
for ng = 1:NGRPS
    eval(['!rm ',filename,num2str(ng),'.met.epart.',n_proc])
    eval(['!rm ',filename,num2str(ng),'.met'])
    eval(['!rm ',filename,num2str(ng),'.met.dgraph'])    
    %eval(['!del ',filename,num2str(ng),'.met.npart.',proc])
end
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Wrote finally merged METIS-file: \t\t%s',[filename,'.met.epart.',n_proc]));
disp('-----------------------------------------------------------------------------------')
disp(' ')


%Visualization
C = colormap; c = floor(length(C)/proc);
color = C(1:c:end,:);

ct = 0;
cq = 0;
for i=1:NELEM
    C = color(met(i)+1,:);
    fill(X(connectivity(i,1:element_type(i)),1),X(connectivity(i,1:element_type(i)),2),C);  hold on
end
axis equal
