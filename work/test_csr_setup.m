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
            Amcg(tab(l,i),tab(l,j)) = Amcg(tab(l,i),tab(l,j))+mcg(i,j)*0.05/2;
        end         
    end    
end

figure(2)
spy(Amcg)

figure(3)
spy(Acsr-Amcg)

diff = Acsr-Amcg;
fprintf('Max |diff| = %g\n',max(abs(diff(:))));


%%
% Build using CSR format

colp1 = dlmread('colp1.crs');
ptrp1 = dlmread('ptrp1.crs');
valp1 = dlmread('valp1.crs');
rhsp1 = dlmread('rhsp1.crs');

Ap1 = zeros(length(ptrp1)-1);
row = 1;
for i = 1:length(ptrp1(:,1))-1
    nnz = ptrp1(i,1):(ptrp1(i+1,1)-1);
    Ap1(row,colp1(nnz,1)) = valp1(nnz,1);
    row = row+1;        
end

p1coeffmat = dlmread('p1coeff.crs');
for k = 1:length(p1coeffmat(:,1))
    el = p1coeffmat(k,1);
    i = p1coeffmat(k,2);
    j = p1coeffmat(k,3);
    p1xcfmat{el,2}(i,j) = p1coeffmat(k,4);
    p1ycfmat{el,2}(i,j) = p1coeffmat(k,5);
    p1cfmat{el,2}(i,j) = p1coeffmat(k,6);
    p1rhsmat{el,2}(j) = p1coeffmat(k,7);
end

phixymat = dlmread('phixphiy.crs');
for k = 1:length(phixymat(:,1))
    el = phixymat(k,1);
    i  = phixymat(k,2);
    phixmat{el,2}(i) = phixymat(k,3);
    phiymat{el,2}(i) = phixymat(k,4);
end