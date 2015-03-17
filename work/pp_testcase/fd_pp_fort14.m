% fd_pp_fort14.m   
% Jessica Meixner 
% March 11, 2015 

% Making a regular, rectangular fort.14 mesh for using with PP 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% Input information (also need to change depth!!!!)%%%%%%%%%%%%%%%
%%%%%%%%%%                    (and maybe boundary?)         %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NumRows = 6; %num of nodes in y direction                                  
NumCols = 101; %number of nodes in x direction

LengthY = 5; 
LengthX = 10; 

AGRID = 'PPGrid'; %string identifying grid <=24 characters

filename = 'fort.14'; 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp(['NumRows = ' num2str(NumRows)])
disp(['NumCols = ' num2str(NumCols)])



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out beginning info: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

fileID = fopen(filename,'w');
fprintf(fileID,'%s\n',AGRID);


NP = NumRows*NumCols; 
NE = 2*(NumRows-1)*(NumCols-1);
fprintf(fileID,'%i %i \n',NE,NP);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out node info: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = linspace(0,LengthX,NumCols); 
y = linspace(0,LengthY,NumRows); 
inode = 1;
for i = 1:NumRows
    for j = 1:NumCols
        
        %%%%%%%%%%%%%%%%%%%%% DEPTH, MIGHT NEED TO BE CHANGED %%%%%%%%%%%%%
        depth = 2; 
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        fprintf(fileID,'%4i  %8.8f  %8.8f  %8.8f\n',inode,x(j),y(i),depth);
        inode = inode + 1;
    end %i
end% !j

disp(['DX = ' num2str(x(2)-x(1))])
disp(['DY = ' num2str(y(2)-y(1))])
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Write out element connectivity info: 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


ielem = 0;
for j = 1:(NumRows-1)
    for i = 1:(NumCols-1)
        
        %for each box we have 2 elements 
        %  node3 = (j+1-1)*NumCols+i     node4 = (j+1-1)*NumCols+(i+1)
        %  (x(i),y(j+1))         (x(i+1),y(j+1))
        %        .------------------.
        %        |                 /|  
        %        |                / | 
        %        |               /  |
        %        |              /   |
        %        |     1       /    |
        %        |            /     |
        %        |           /      |
        %        |          /       |
        %        |         /        |
        %        |        /         |
        %        |       /          |
        %        |      /           |
        %        |     /            |
        %        |    /             |
        %        |   /        2     |
        %        |  /               |
        %        | /                |
        %        .------------------.
        %  (x(i),y(j))          (x(i+1),y(j))
        %  node1=(j-1)*NumCols+i     node2 = (j-1)*NumCols+(i+1) 
        %
        %   element 1:   node1 node4 node3 
        %   element 2:   node4 node1 node2 
        
        node1 = (j-1)*NumCols+i;
        node2 = (j-1)*NumCols+(i+1);
        node3 = (j+1-1)*NumCols+i;
        node4 = (j+1-1)*NumCols+(i+1);
        
        %elemen 1: 
        ielem = ielem+1; 
        fprintf(fileID,'%i \t 3 \t %i \t %i  \t %i\n', ielem, node4,node3,node1); 
        
        %elemen 2: 
        ielem = ielem+1; 
        fprintf(fileID,'%i \t 3 \t %i \t  %i \t  %i\n', ielem, node4,node1,node2); 
        
    end %i
end% !j

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%  Write out boundary info 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%                 boundary 3
%      .------------.-------------.--------.
%      |                                   |
%      |                                   |
%      |                                   |
%      .                                   .
%      |                                   |
%      | boundary 4                        | boundary 2
%      |                                   |
%      .                                   .
%      |                                   |
%      |                                   |
%      |                                   |
%      .------------.-------------.--------.
%               boundary 1 
%
%

NOPE = 0; %number of elevation specified boundary forcing segments.
NETA = 0; %total number of elevation specified boundary nodes

fprintf(fileID,'%i  = Number of open boundaries \n',NOPE);
fprintf(fileID,'%i  = Total number of open boundary nodes  \n',NETA);



NBOU = 4; % number of normal flow (discharge) specified boundary segments. 
          %These include zero normal flow (land) boundaries
NVEL = NumRows*2+NumCols*2-4;  % total number of normal flow specified 
                               %boundary nodes including both the front 
                               %and back nodes on internal barrier boundaries.

fprintf(fileID,'%i  = Number of land boundaries \n',NBOU);
fprintf(fileID,'%i  = Total number of land boundary nodes  \n',NVEL);


NVELL(1:4) = [NumCols; NumRows; NumCols; NumRows];
%number of nodes in normal flow specified boundary segment k
IBTYPE(1:4) = [10; 10; 10; 12]; %could want this to be 0?

                               
for k = 1:NBOU
    fprintf(fileID,'%i %i = Number of nodes for land boundary %i  \n', NVELL(k), IBTYPE(k), k);
    for j= 1:NVELL(k)
        if k==1
            fprintf(fileID,'%i\n',j);
        elseif k==2
            fprintf(fileID,'%i\n',(j-1)*NumCols+NumCols);    
        elseif k==3
            fprintf(fileID,'%i\n',(NumRows-1)*NumCols+(NumCols-j+1));
        elseif k==4
            fprintf(fileID,'%i\n',((NumRows-j+1)-1)*NumCols+1); 
        end
    end
end





%from adcirc website fort.14 explanation: 

% AGRID
% NE, NP
% 
% for k=1 to NP
%     JN, X(JN), Y(JN), DP(JN)
% end k loop
% 
% for k=1 to NE
%     JE, NHY, NM(JE,1), NM(JE,2), NM(JE,3)
% end k loop
% 
% NOPE
% NETA
% 
% for k=1 to NOPE
%     NVDLL(k), IBTYPEE(k)
%     for j=1 to NVDLL(k)
%         NBDV(k,j)
%     end j loop
% end k loop
% 
% NBOU
% NVEL
% 
% for k=1 to NBOU
%     NVELL(k), IBTYPE(k)
%     for j=1,NVELL(k)
%         NBVV(k,j) – include this line only if IBTYPE(k) = 0, 1, 2, 10, 11, 12, 20, 21, 22, 30
%         NBVV(k,j), BARLANHT(k,j), BARLANCFSP(k,j) – include this line only if IBTYPE(k) = 3, 13, 23
%         NBVV(k,j), IBCONN(k,j), BARINHT(k,j), BARINCFSB(k,j), BARINCFSP(k,j) – include this line only if IBTYPE(k) = 4, 24
%     end j loop
% end k loop


