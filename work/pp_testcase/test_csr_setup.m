%% Testing CSR setup in DGSWEM

col = dlmread('col.crs');
ptr = dlmread('ptr.crs');
val = dlmread('val.crs');

mcg = dlmread('MCGlocal.txt');

tab = dlmread('elem_table.txt');

% Build using CSR format
Acsr = zeros(length(ptr)-1);
row = 1;
colnew = 0*col;
valnew = 0*val;
for i = 1:length(ptr)-1
    nnz = ptr(i):(ptr(i+1)-1);
    Acsr(row,col(nnz)) = val(nnz);
    row = row+1;
    
    coltmp = col(nnz);
    valtmp = val(nnz);
    
    [sorti,sortj] = sort(coltmp);
    colnew(nnz) = col(nnz(sortj));
    valnew(nnz) = val(nnz(sortj));
end

fid1 = fopen('colnew.crs','w');
fid2 = fopen('valnew.crs','w');
for i = 1:length(colnew)
    fprintf(fid1,'%d\n',colnew(i));
    fprintf(fid2,'%f\n',valnew(i));
end

    

figure(1)
spy(Acsr)

% Build using standard approach
Amcg = zeros(max(max(tab)));
for l = 1:length(tab(:,1))
    for i = 1:3
        for j = 1:3
            Amcg(tab(l,i),tab(l,j)) = Amcg(tab(l,i),tab(l,j))+mcg(i,j);
        end         
    end    
end

figure(2)
spy(Amcg)

figure(3)
spy(Acsr-Amcg)

diff = Acsr-Amcg;
fprintf('Max |diff| = %g\n',max(abs(diff(:))));