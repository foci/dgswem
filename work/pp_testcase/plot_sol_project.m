restoredefaultpath
clear all
close all
clc

makemovie = 0;
nsnap = 100;
%direc = '../curvilinear/';
%direc = '../ptr_test2/';
% direc = '../ptr_test3/';
% direc = '../quad_loop/';
%direc = '../LDG/';
%direc = '../openmp/';
direc = '~/dgswe/';
grd_direc = '~/data-drive/galveston/dgswe/';
sol_direc = '~/data-drive/galveston/dgswe/';
grd_name = 'galveston.grd';

% direc = '..';

% FramesFolder = strcat(direc,'output/velsnap') ;
% if ( exist(FramesFolder,'dir') == 0 ) 
%     mkdir(FramesFolder) ; 
% end

% f1 = figure('Position',[1360 667 560 420]);
% f2 = figure('Position',[798 667 560 420]);

h = figure('Position',[100 100 1500 1000]) ;

% Open water column height, x momentum, and y momentum modal solution files
fid_H = fopen([sol_direc,'solution_H.d']);
fid_Qx = fopen([sol_direc,'solution_Qx.d']);
fid_Qy = fopen([sol_direc,'solution_Qy.d']);

% Read in grid name
% name = fscanf(fid_H,' %s ', 1);
% name = fscanf(fid_Qx,' %s ', 1);
% name = fscanf(fid_Qy,' %s ', 1);
% name = name(3:end); %get rid of ..

line = fgetl(fid_H);
line = fgetl(fid_Qx);
line = fgetl(fid_Qy);

% Read in grid file
[EToV,VX,HB,~,~,~] = readfort14([grd_direc,grd_name]);

% Open L2 projection information file
fid = fopen([sol_direc,'projection.d']);

N = fscanf(fid, '%g', 1); % number of degrees of freedom
mL2 = fscanf(fid,' %g ',[N 3])'; % L2 projection matrix mL2(i,j) = integral of (phi_DG_basis(j)*phi_linear_nodal(i))
mm = fscanf(fid,' %g ',[3 3])'; % linear nodal mass matrix
fclose(fid);

p = -3/2 + sqrt(9+8*N-8)/2;

Nm = (p+1)^2;

[ne ~] = size(EToV); % number of elements
[nn ~] = size(VX); % number of nodes

area = zeros(ne,1);
for el = 1:ne
    x1 = VX(EToV(el,1),1);
    x2 = VX(EToV(el,2),1);
    x3 = VX(EToV(el,3),1);
    
    y1 = VX(EToV(el,1),2);
    y2 = VX(EToV(el,2),2);
    y3 = VX(EToV(el,3),2);
    
    area(el) = .5*((x2*y3-x3*y2)+(x3*y1-x1*y3)+(x1*y2-x2*y1));
end

% Assemble global linear mass matrix

% M = zeros(nn);
% for el = 1:ne
%     for i = 1:3
%         for j = 1:3
%             M(EToV(el,i),EToV(el,j)) = M(EToV(el,i),EToV(el,j)) + .5*area(el)*mm(i,j);
%         end
%     end
% end
% 
% M = sparse(M);

iM = zeros(9*nn,1);
jM = zeros(9*nn,1);
sM = zeros(9*nn,1);
i = 0;
for el = 1:ne
    
    iM(i+1) = EToV(el,1);  jM(i+1) = EToV(el,1);
    iM(i+2) = EToV(el,1);  jM(i+2) = EToV(el,2);
    iM(i+3) = EToV(el,1);  jM(i+3) = EToV(el,3);
    
    iM(i+4) = EToV(el,2);  jM(i+4) = EToV(el,1);
    iM(i+5) = EToV(el,2);  jM(i+5) = EToV(el,2);
    iM(i+6) = EToV(el,2);  jM(i+6) = EToV(el,3);
    
    iM(i+7) = EToV(el,3);  jM(i+7) = EToV(el,1);
    iM(i+8) = EToV(el,3);  jM(i+8) = EToV(el,2);
    iM(i+9) = EToV(el,3);  jM(i+9) = EToV(el,3); 
    
    sM(i+1) = .5*area(el)*mm(1,1);
    sM(i+2) = .5*area(el)*mm(1,2);
    sM(i+3) = .5*area(el)*mm(1,3);
    
    sM(i+4) = .5*area(el)*mm(2,1);
    sM(i+5) = .5*area(el)*mm(2,2);
    sM(i+6) = .5*area(el)*mm(2,3);
   
    sM(i+7) = .5*area(el)*mm(3,1);
    sM(i+8) = .5*area(el)*mm(3,2);
    sM(i+9) = .5*area(el)*mm(3,3);
    
    i = i + 9;
end

M = sparse(iM,jM,sM,nn,nn);

zeta = zeros(nn,nsnap);
u = zeros(nn,nsnap);
v = zeros(nn,nsnap);
vel = zeros(nn,nsnap);
t = zeros(nsnap);

snap = 0;
while ~feof(fid_H) && snap < nsnap

snap = snap + 1;
    
th = fscanf(fid_H,' %g ', 1); % read in time
H = fscanf(fid_H,' %g ', [ne Nm])'; % read in H solution at time t

tqx = fscanf(fid_Qx,' %g ', 1); % read in time
Qx = fscanf(fid_Qx,' %g ', [ne Nm])'; % read in Qx solution at time t

t(snap) = fscanf(fid_Qy,' %g ', 1); % read in time
Qy = fscanf(fid_Qy,' %g ', [ne Nm])'; % read in Qy solution at time t

prod_H = mL2*H(1:N,:); % compute product of elemental degrees of freedom and projection matrix
prod_Qx = mL2*Qx(1:N,:);
prod_Qy = mL2*Qy(1:N,:);

% Assemble RHS vector
B_H = zeros(nn,1);
B_Qx = zeros(nn,1);
B_Qy = zeros(nn,1);
for el = 1:ne
    for i = 1:3
       B_H(EToV(el,i)) = B_H(EToV(el,i)) + .5*area(el)*prod_H(i,el);
       B_Qx(EToV(el,i)) = B_Qx(EToV(el,i)) + .5*area(el)*prod_Qx(i,el); 
       B_Qy(EToV(el,i)) = B_Qy(EToV(el,i)) + .5*area(el)*prod_Qy(i,el); 
    end
end

HL2 = M\B_H; % Compute L2 projected solution
QxL2 = M\B_Qx;
QyL2 = M\B_Qy;

zeta(:,snap) = HL2-HB;
u(:,snap) = QxL2./HL2;
v(:,snap) = QyL2./HL2;
vel(:,snap) = sqrt(u(:,snap).^2 + v(:,snap).^2);
end

projsol = zeros(nn,4,nsnap);
projsol(:,1,:) = zeta;
projsol(:,2,:) = u;
projsol(:,3,:) = v;
projsol(:,4,:) = vel;

velmax = max(max(vel));
zmax = max(max(zeta));

for tsnap = 1:snap
%     figure 
    
    disp(['Time snap: ',num2str(tsnap),'/',num2str(snap)])
    [day hr min sec] = s2dhms(t(tsnap));
    
%     subplot(2,1,1)
% %     figure(f1)
% %    pdeplot( VX', [], EToV','xydata',zeta(:,tsnap),'zdata',zeta(:,tsnap),'xystyle','interp','colormap','jet','mesh','off') ;
%     pdeplot( VX', [], EToV','xydata',zeta(:,tsnap),'xystyle','interp','colormap','jet','mesh','off') ;
%     grid on
%     ttext = ['Elevation solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(min),', Second:  ',num2str(sec),')'] ;
%     title(ttext)
%     xlabel('x')
%     ylabel('y')
%     axis image
%     %zlim([-.5 1])
%     %caxis([3.5 5])
    
%     subplot(2,1,2)
%     figure(f2)
    pdeplot( VX', [], EToV','xydata',vel(:,tsnap),'xystyle','interp','colormap','jet','mesh','off') ;
    hold on
    quiver( VX(:,1), VX(:,2),u(:,tsnap), v(:,tsnap), '-k' ) ;
    grid on
    hold off
    ttext = ['Velocity solution: t = ',num2str(t(tsnap)),' (Day:  ',num2str(day),', Hour:  ',num2str(hr),', Minute:  ',num2str(min),', Second:  ',num2str(sec),')'] ;
    title(ttext)
    xlabel('x')
    ylabel('y')
    axis image
    %zlim([-.5 1])
    %caxis([0 1])

    if makemovie == 1
       if tsnap == 1
           set(gca,'nextplot','replacechildren','visible','off')
            f = getframe(2);
            [im,map] = rgb2ind(f.cdata,256,'nodither');
            im(1,1,1,20) = 0;
       else
            f = getframe(2);
            im(:,:,1,tsnap) = rgb2ind(f.cdata,map,'nodither');
       end
    end

    pause(.01)

%     print('-r300','-dpng',sprintf('%s/surf%04d',FramesFolder,tsnap)) ;
    
end

if makemovie == 1
    imwrite(im,map,'Inlet.gif','DelayTime',0,'LoopCount',inf)
end
