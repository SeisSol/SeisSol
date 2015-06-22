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
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp('    %%                                                     %%')
disp('    %%                   GAMBIT_HEX2TET                    %%')
disp('    %%           Creates regular tetrahedral mesh          %%')
disp('    %%           on the base of an hexahedral mesh         %%')
disp('    %%                                                     %%')
disp('    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
disp(' ')

clear, close all;

filename = input('    Filename(suffix ".neu" assumed)        :  ','s');
disp(' ')
nx       = input('    Give number of elements in x-dimension :  ');
ny       = input('    Give number of elements in y-dimension :  ');
nz       = input('    Give number of elements in z-dimension :  ');

disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Reading data from: \t\t%s',[filename,'.neu']));
    
fid_in   = fopen([filename ,'.neu']);
fid_out  = fopen([filename ,'_tetra.neu'],'w');
junk  = fgetl(fid_in);      fprintf(fid_out,'%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
param = fscanf(fid_in,'%i',[6,1]);  NX = param(1);  NT = param(2); 
param(2) = param(2)*5;      fprintf(fid_out,'\n%10.0f%10.0f%10.0f%10.0f%10.0f%10.0f',param);
junk  = fgetl(fid_in);      
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
X     = fscanf(fid_in,'%g',[4,NX]); 
                            fprintf(fid_out,'\n%10.0f%20.10e%20.10e%20.10e',X);
X = X'; X(:,1) = []; 
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
junk  = fgetl(fid_in);      fprintf(fid_out,'\n%s',junk);
hex   = fscanf(fid_in,'%g',[11,NT]); hex = hex'; hex(:,1:3) = [];
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
junk  = fgetl(fid_in);
fclose(fid_in);

disp(sprintf('\t Read \t\t\t\t\t%i vertices',NX));
disp(sprintf('\t Read \t\t\t\t\t%i hexahedral elements',NT));

xmin = min(X(:,1)); xmax = max(X(:,1));
ymin = min(X(:,2)); ymax = max(X(:,2));
zmin = min(X(:,3)); zmax = max(X(:,3));

dx = (xmax-xmin)/nx;
dy = (ymax-ymin)/ny;
dz = (zmax-zmin)/nz;

disp(' ')
disp('-----------------------------------------------------------------------------------')
disp(sprintf('\t Creating tetraheral elemenst ... \n'));
t = cputime;

tet = zeros(NT*5,7);
count = 0;
for i = 1:NT
    
    center = sum([X(hex(i,:),1), X(hex(i,:),2), X(hex(i,:),3)])/8;
    
    xpos = ceil( (center(1)-xmin)/dx );
    ypos = ceil( (center(2)-ymin)/dy );
    zpos = ceil( (center(3)-zmin)/dz );
    
    if( (mod(xpos,2)==1 & mod(ypos,2)==1 & mod(zpos,2)==1) | ...
        (mod(xpos,2)==0 & mod(ypos,2)==0 & mod(zpos,2)==1) | ...
        (mod(xpos,2)==1 & mod(ypos,2)==0 & mod(zpos,2)==0) | ...
        (mod(xpos,2)==0 & mod(ypos,2)==1 & mod(zpos,2)==0) )
    
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,1),hex(i,2),hex(i,3),hex(i,5)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,6),hex(i,2),hex(i,5),hex(i,8)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,4),hex(i,3),hex(i,2),hex(i,8)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,7),hex(i,5),hex(i,3),hex(i,8)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,3),hex(i,5),hex(i,2),hex(i,8)];
    
    else    
        
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,8),hex(i,7),hex(i,4),hex(i,6)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,3),hex(i,1),hex(i,4),hex(i,7)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,2),hex(i,1),hex(i,6),hex(i,4)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,5),hex(i,6),hex(i,1),hex(i,7)];
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,4),hex(i,6),hex(i,7),hex(i,1)];
        
    end
   
end

disp(sprintf('\t Writing output-file: \t\t\t%s',[filename,'_tetra.neu']));

fprintf(fid_out,'\n%8.0f%3.0f%3.0f%9.0f%8.0f%8.0f%8.0f',tet');
fprintf(fid_out,'\n%s','ENDOFSECTION');
fprintf(fid_out,'\n%s','       ELEMENT GROUP 2.0.0');
fprintf(fid_out,'\n%s','GROUP:          1 ELEMENTS:');
fprintf(fid_out,'%11.0f%',NT*5);
fprintf(fid_out,'%s',' MATERIAL:          2 NFLAGS:          1');
fprintf(fid_out,'\n%s','                           fluid');
fprintf(fid_out,'\n%s','       0');
fprintf(fid_out,'\n%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f%8.0f',(1:1:NT*5));
fprintf(fid_out,'\n%s\n','ENDOFSECTION');

fclose(fid_out);

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t Conversion finished successfully!  (%g CPU sec)\n',cputime-t));
disp('-----------------------------------------------------------------------------------')
disp(' ')
