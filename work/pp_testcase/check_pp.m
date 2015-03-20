etaver1 = dlmread('dgfd_ver1.txt');
etaver2 = dlmread('dgfd_ver2.txt');
for k = 1:length(etaver1(:,1))
    eta1(etaver1(k,1),etaver1(k,2)) = etaver1(k,3);
    eta2(etaver2(k,1),etaver2(k,2)) = etaver2(k,3);
end

figure(1)
subplot(1,2,1)
surf(eta1)
subplot(1,2,2)
surf(eta2)

return
%%
tmp = dlmread('prestest.txt');

for k = 1:length(tmp(:,1))
    eta_tst(tmp(k,1),tmp(k,2)) = tmp(k,3);
    u_tst(tmp(k,1),tmp(k,2)) = tmp(k,4);
    v_tst(tmp(k,1),tmp(k,2)) = tmp(k,5);
    pd_tst(tmp(k,1),tmp(k,2)) = tmp(k,6);
    pb_tst(tmp(k,1),tmp(k,2)) = tmp(k,7);
end

figure(1)
surf(eta_tst);shading interp

figure(2)
surf(-pd_tst/9.81);shading interp

figure(3)
surf(u_tst); shading interp

figure(4)
surf(v_tst); shading interp


etalng = [eta_tst(:,1) eta_tst eta_tst(:,end)];
eta_tst_x = (etalng(:,3:end)-etalng(:,1:end-2))/2/.1;
eta_tst_xx = (etalng(:,3:end)-2*etalng(:,2:end-1)+etalng(:,1:end-2))/.1/.1;

figure(5)
surf(eta_tst_x); shading interp
figure(6)
surf(eta_tst_xx); shading interp