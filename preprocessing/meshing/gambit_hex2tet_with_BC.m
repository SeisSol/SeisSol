%%
% @file
% This file is part of SeisSol.
%
% @author Martin Kaeser (martin.kaeser AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/kaeser)
% @author Thomas Ulrich (thomas.ulrich AT geophysik.uni-muenchen.de, http://www.geophysik.uni-muenchen.de/Members/ulrich)
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
disp('    %%                   GAMBIT_HEX2TET_with_BC            %%')
disp('    %%           Creates regular tetrahedral mesh          %%')
disp('    %%           on the base of an hexahedral mesh         %%')
disp('    %%              Then add boundary conditions           %%')
disp('    %%              useful for complex geometries          %%')
disp('    %%  i.e. when gambit cannot recreate the geometry      %%') 
disp('    %%                     from the mesh                   %%')
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
param(4)=3;
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
disp(sprintf('\t Creating tetraheral elements ... \n'));
t = cputime;

tet = zeros(NT*5,7);
count = 0;

b101=[];
b103=[];
b105=[];

for i = 1:NT
    
    center = sum([X(hex(i,:),1), X(hex(i,:),2), X(hex(i,:),3)])/8;
    
    xpos = ceil( (center(1)-xmin)/dx );
    ypos = ceil( (center(2)-ymin)/dy );
    zpos = ceil( (center(3)-zmin)/dz );

    %mx (="minus x" taken as an example) will contains 
    %(1) the id of the tetra elements having a face in x=xmin.  
    %(2) the face id of these faces (between 1 and 4). 
    %This information is used for setting up the boundary conditions
    
    mx=[];
    my=[];
    mz=[];
    px=[];
    py=[];
    pz=[];

    %hex2 is a permutation of the index of the hexahedron in order to put
    %it in a(n) (arbitrary) situation of reference
    %   4---8
    % 3---7 !
    % !   !
    %   2---6
    % 1---5
    %
    %then we have the face id: top=pz=3, down=mz=2, my=6, py=5 mx=4,px=2
    
    
    hex2=zeros(8,2);
    for j = 1:8
        hex2(j,2)=j;
        xj=X(hex(i,j),1);
        yj=X(hex(i,j),2);
        zj=X(hex(i,j),3);
        
        if (xj-center(1)<0)
            if (yj-center(2)<0)
                if (zj-center(3)<0)
                    hex2(j,1)=5;
                else
                    hex2(j,1)=7;
                end
            else
                if (zj-center(3)<0)
                    hex2(j,1)=1;
                else
                    hex2(j,1)=3;
                end
            end
        else
            if (yj-center(2)<0)
                if (zj-center(3)<0)
                    hex2(j,1)=6;
                else
                    hex2(j,1)=8;
                end
            else
                if (zj-center(3)<0)
                    hex2(j,1)=2;
                else
                    hex2(j,1)=4;
                end
            end
        end
    end
    hex2 = sortrows(hex2,1);
    hex2=hex2(:,2);
    
    if( (mod(xpos,2)==1 & mod(ypos,2)==1 & mod(zpos,2)==1) | ...
        (mod(xpos,2)==0 & mod(ypos,2)==0 & mod(zpos,2)==1) | ...
        (mod(xpos,2)==1 & mod(ypos,2)==0 & mod(zpos,2)==0) | ...
        (mod(xpos,2)==0 & mod(ypos,2)==1 & mod(zpos,2)==0) )

        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(1)),hex(i,hex2(2)),hex(i,hex2(3)),hex(i,hex2(5))];
        mx=[mx;[count,surfid('134')]];
        py=[py;[count,surfid('123')]];
        mz=[mz;[count,surfid('124')]];
        
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(6)),hex(i,hex2(2)),hex(i,hex2(5)),hex(i,hex2(8))];
        px=[px;[count,surfid('124')]];
        my=[my;[count,surfid('134')]];
        mz=[mz;[count,surfid('123')]];

        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(4)),hex(i,hex2(3)),hex(i,hex2(2)),hex(i,hex2(8))];
        px=[px;[count,surfid('134')]];
        py=[py;[count,surfid('123')]];
        pz=[pz;[count,surfid('124')]];

        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(7)),hex(i,hex2(5)),hex(i,hex2(3)),hex(i,hex2(8))];
        mx=[mx;[count,surfid('123')]];
        my=[my;[count,surfid('124')]];
        pz=[pz;[count,surfid('134')]];
        
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(3)),hex(i,hex2(5)),hex(i,hex2(2)),hex(i,hex2(8))];
        
    else    

        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(8)),hex(i,hex2(7)),hex(i,hex2(4)),hex(i,hex2(6))];
        px=[px;[count,surfid('134')]];
        my=[my;[count,surfid('124')]];
        pz=[pz;[count,surfid('123')]];
        
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(3)),hex(i,hex2(1)),hex(i,hex2(4)),hex(i,hex2(7))];
        mx=[mx;[count,surfid('124')]];
        py=[py;[count,surfid('123')]];
        pz=[pz;[count,surfid('134')]];
        
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(2)),hex(i,hex2(1)),hex(i,hex2(6)),hex(i,hex2(4))];
        px=[px;[count,surfid('134')]];
        py=[py;[count,surfid('124')]];
        mz=[mz;[count,surfid('123')]];
        
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(5)),hex(i,hex2(6)),hex(i,hex2(1)),hex(i,hex2(7))];
        mx=[mx;[count,surfid('134')]];
        my=[my;[count,surfid('124')]];
        mz=[mz;[count,surfid('123')]];
        
        count = count+1;
        tet(count,:) = [count,6,4,hex(i,hex2(4)),hex(i,hex2(6)),hex(i,hex2(7)),hex(i,hex2(1))];
        
    end
    %free surface
    if (zpos==nz)
       b101=vertcat(b101,pz);
    end
    
    %rupture
    %if (zpos>5)&&(xpos>5)&&(xpos<=45)
        if (ypos==ny/2)
           b103=vertcat(b103,py);
        end
        if (ypos==(ny/2+1))
           b103=vertcat(b103,my);
        end
    %end
    %absorbing surfaces
    if (ypos==1)
       b105=vertcat(b105,my);
    end
    if (ypos==ny)
       b105=vertcat(b105,py);
    end
    if (xpos==1)
       b105=vertcat(b105,mx);
    end
    if (xpos==nx)
       b105=vertcat(b105,px);
    end
    if (zpos==1)
       b105=vertcat(b105,mz);
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

fprintf(fid_out,' BOUNDARY CONDITIONS 2.3.16\n');
aSize=size(b101);
nbound=aSize(1);
fprintf(fid_out,'%s%8.0f%8.0f%8.0f%8.0f\n','                             101',1,nbound,0,6);
for i = 1:nbound
    fprintf(fid_out,'%10.0f%5.0f%5.0f\n', b101(i,1),6,b101(i,2));
end
fprintf(fid_out,'%s\n','ENDOFSECTION');

fprintf(fid_out,' BOUNDARY CONDITIONS 2.3.16\n');
aSize=size(b103);
nbound=aSize(1);
fprintf(fid_out,'%s%8.0f%8.0f%8.0f%8.0f\n','                             103',1,nbound,0,6);
for i = 1:nbound
    fprintf(fid_out,'%10.0f%5.0f%5.0f\n', b103(i,1),6,b103(i,2));
end
fprintf(fid_out,'%s\n','ENDOFSECTION');

fprintf(fid_out,' BOUNDARY CONDITIONS 2.3.16\n');
aSize=size(b105);
nbound=aSize(1);
fprintf(fid_out,'%s%8.0f%8.0f%8.0f%8.0f\n','                             105',1,nbound,0,6);
for i = 1:nbound
    fprintf(fid_out,'%10.0f%5.0f%5.0f\n', b105(i,1),6,b105(i,2));
end
fprintf(fid_out,'%s\n','ENDOFSECTION');

fclose(fid_out);

disp('-----------------------------------------------------------------------------------')
disp(sprintf('\n\t Conversion finished successfully!  (%g CPU sec)\n',cputime-t));
disp('-----------------------------------------------------------------------------------')
disp(' ')


%here we can plot the boundary condition to be sure that everything is
%alright

plotmesh='n';
tet=tet(:,4:7);
s_vert(1,:) = [1,3,2];   s_vert(2,:) = [1,2,4];   s_vert(3,:) = [2,3,4];   s_vert(4,:) = [1,4,3];

aSize=size(b105);
NS=aSize(1);

tri=zeros(NS,3);
for i = 1:NS
    tri(i,:)=tet(b105(i,1),s_vert(b105(i,2),:))';
end

if plotmesh == 'y'
    figure; hold on
    trimesh(tri,X(:,1),X(:,2),X(:,3),X(:,3));
    axis equal,xlabel('x'),ylabel('y'),zlabel('z')
end

