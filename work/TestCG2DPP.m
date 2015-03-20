%% Test CG 2d PP
clear all

g = 9.81;


if 0
x = linspace(-5,5,81);
y = linspace(-5,5,81);
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

xvec = X(:);
yvec = Y(:);

elem = delaunay(xvec,yvec);
ne = length(elem(:,1));

eta = e0*exp(-(xvec.^2+yvec.^2)/sig);
etax = -2/sig*e0*exp(-(xvec.^2+yvec.^2)/sig).*xvec;
etay = -2/sig*e0*exp(-(xvec.^2+yvec.^2)/sig).*yvec;
etaxx = e0*(4*exp(-(xvec.^2+yvec.^2)/sig).*xvec.^2/sig^2-2*exp(-(xvec.^2+yvec.^2)/sig)/sig);
etayy = e0*(4*exp(-(xvec.^2+yvec.^2)/sig).*yvec.^2/sig^2-2*exp(-(xvec.^2+yvec.^2)/sig)/sig);

e0 = 0.02;
h0 = 1;
h1 = 0.2;
sig = 2;
p0 = g;
h = 0*xvec+1;

else
    tmp = dlmread('dgcgsoln.txt');
    nn = max(max(tmp(:,2:4)));
    ne = max(tmp(:,1));
    xvec = zeros(nn,1);
    yvec = xvec;
    eta = xvec;
    h = eta;
    u = eta;
    v = eta;
    pcg = eta;
    p2 = eta;
    pd = eta;
    pb = eta;
    elem = zeros(ne,3);
    for i = 1:length(tmp(:,1))
        elem(tmp(i,1),1:3)  = tmp(i,2:4);
        xvec(elem(i,:)) = tmp(i,5:7);
        yvec(elem(i,:)) = tmp(i,8:10);
        eta(elem(i,:))  = tmp(i,11:13);
        u(elem(i,:))  = tmp(i,14:16);
        v(elem(i,:))  = tmp(i,17:19);
        h(elem(i,:))  = tmp(i,20:22);
        pcg(elem(i,:))  = tmp(i,23:25);
        p2(elem(i,:))  = tmp(i,26:28);
        pd(elem(i,:))  = tmp(i,29:31);
        pb(elem(i,:))  = tmp(i,32:34);
    end
    u = u./(eta+h);
    v = v./(eta+h);
    
%     rhscg = dlmread('rhsp1.crs');
    rhscg = 0*eta;
    
    tmp = dlmread('dgfd_ver1.txt');
    for k = 1:length(tmp(:,1))
        eta2(tmp(k,1),tmp(k,2)) = tmp(k,3);
    end
    
    
end


triplot(elem,xvec,yvec);


%%
manyes = 0;

% mansoln = -1/12*p0*(8*(etay.^2+etax.^2)+10*eta.^2.*(etayy+etaxx)+6*(etayy+etaxx)+...
%     eta.*(-15+16*(etaxx+etayy)));

mansoln = 0*eta;

if manyes == 1
% Linear Flatbed
mansoln = mansoln + 1/4*p0*(5*eta-2*h0*(etaxx+etayy));
% mansoln = mansoln + (-1/12*g*(15*eta+4*h0^2*(etayy+etaxx)));

% Linear Varying
mansoln = mansoln + (-1/4*h0*h1*p0*(etay.^2+etax.^2+2*eta.*(etayy+etaxx)));
% mansoln = mansoln + (-1/12*g*h1*(3*(4+h0+5*h1*eta).*etay.^2+12*h0*etax.^2+...
%                       10*h1*eta.^2.*(etayy+etaxx)+eta.*(14*h0*(etayy+etaxx)+15*h1*etax.^2)));

% Nonlinear Flat
mansoln = mansoln + (-1/6*p0*(4*h0*etay.^2+4*h0*etax.^2+8*h0*eta.*...
    (etayy+etaxx)+5*eta.^2.*(etaxx+etayy)));
% mansoln = mansoln + (1/3*((h0*(2*g+5*h0*v0^2)+10*h0*v0^2*eta+5*v0^2*eta.^2).*etax.^2+...
%     10*u0*v0*(h0+eta).^2.*etax.*etay+h0*(2*g+5*h0*u0^2)*etax.^2+...
%     5*u0^2*eta.^2*etax.^2-h0*eta.*(g*etayy-10*u0^2*etax.^2+g*etaxx)));

% Nonlinear Varying
mansoln = mansoln + (-5/12*h1*p0*eta.*(3*etay.^2+3*etax.^2+2*eta.*(etaxx+etayy)));
% mansoln = mansoln + (5/12*h1*eta.*(4*v0^2*(2*h0+(h1+2)*eta).*etay.^2+...
%     8*u0*v0*(2*h0+(h1+2)*eta).*etay.*etax+8*h0*u0^2*etax.^2+eta.*((-2*g+3*h0*v0^2)*etayy+...
%        4*u0^2*(h1+2)*etax.^2+6*h0*u0*v0*etaxy-2*g*etaxx+3*h0*u0^2*etaxx)+...
%     3*(h1+1)*eta.^2.*(v0^2*etayy+u0*(2*v0*etaxy+u0*etaxx))));
end


%%
% find quad points
cub2d( 1, :) = [ -5.014265096581790e-01,  -5.014265096581790e-01, 2.335725514527590e-01 ] ; 
cub2d( 2, :) = [ 2.853019316358000e-03,  -5.014265096581790e-01, 2.335725514527590e-01 ] ; 
cub2d( 3, :) = [ -5.014265096581790e-01,  2.853019316358000e-03, 2.335725514527590e-01 ] ; 
cub2d( 4, :) = [ -8.738219710169960e-01,  -8.738219710169960e-01, 1.016898127404140e-01 ] ; 
cub2d( 5, :) = [ 7.476439420339910e-01,  -8.738219710169960e-01, 1.016898127404140e-01 ] ; 
cub2d( 6, :) = [ -8.738219710169960e-01,  7.476439420339910e-01, 1.016898127404140e-01 ] ; 
cub2d( 7, :) = [ -3.792950979324310e-01,  -8.937099003103660e-01, 1.657021512367470e-01 ] ; 
cub2d( 8, :) = [ -8.937099003103660e-01,  -3.792950979324310e-01, 1.657021512367470e-01 ] ; 
cub2d( 9, :) = [ 2.730049982427970e-01,  -8.937099003103660e-01, 1.657021512367470e-01 ] ; 
cub2d( 10, :) = [ -8.937099003103660e-01,  2.730049982427970e-01, 1.657021512367470e-01 ] ; 
cub2d( 11, :) = [ 2.730049982427970e-01,  -3.792950979324310e-01, 1.657021512367470e-01 ] ; 
cub2d( 12, :) = [ -3.792950979324310e-01,  2.730049982427970e-01, 1.657021512367470e-01 ] ; 
phi = zeros(3,length(cub2d(:,1)));
for pt = 1:length(cub2d(:,1))
    phi(1,pt) = -1/2*(cub2d(pt,1)+cub2d(pt,2));
    phi(2,pt) =  1/2*(1+cub2d(pt,1));
    phi(3,pt) =  1/2*(1+cub2d(pt,2));
end
phix = 0*phi;
phiy = 0*phi;
LHS = zeros(length(xvec));
RHS = LHS(:,1);

p1xcfmat = cell(ne,2);
p1ycfmat = cell(ne,2);
p1cfmat  = cell(ne,2);
p1rhsmat = cell(ne,2);
phixmat = cell(ne,2);
phiymat = cell(ne,2);

for l = 1:ne
    xloc = xvec(elem(l,:));
    yloc = yvec(elem(l,:));
    A = 1/2*det([1 1 1;xloc';yloc']');
    Al(l) = A;
    
    
    phix(1,:) = 1/2/A*(yloc(2)-yloc(3));
    phix(2,:) = 1/2/A*(yloc(3)-yloc(1));
    phix(3,:) = 1/2/A*(yloc(1)-yloc(2));
    
    phiy(1,:) = 1/2/A*(xloc(3)-xloc(2));
    phiy(2,:) = 1/2/A*(xloc(1)-xloc(3));
    phiy(3,:) = 1/2/A*(xloc(2)-xloc(1));    
    
    phixmat{l,1} = phix(:,1)';
    phiymat{l,1} = phiy(:,1)';
    
    x1 = xvec(elem(l,1));
    x2 = xvec(elem(l,2));
    x3 = xvec(elem(l,3));
    
    y1 = yvec(elem(l,1));
    y2 = yvec(elem(l,2));
    y3 = yvec(elem(l,3));
    
    e1 = eta(elem(l,1));
    e2 = eta(elem(l,2));
    e3 = eta(elem(l,3));
    
    h1 = eta(elem(l,1));
    h2 = eta(elem(l,2));
    h3 = eta(elem(l,3));
    
    % Compute Integrals
    P1xcfLin = zeros(3);
    P1ycfLin = zeros(3);
    P1cfLin  = zeros(3);
    P1xcfNon = zeros(3);
    P1ycfNon = zeros(3);
    P1cfNon  = zeros(3);
    P1rhsLin = zeros(3,1);
    P1rhsNon = zeros(3,1);
    P1rhsMan = zeros(3,1);
    
    
    for pt = 1:length(cub2d(:,1))
        zi = 0;
        hi = 0;
        zx = 0;
        hx = 0;
        zy = 0;
        hy = 0;
        
        hxy = 0;
        
        ui = 0;
        vi = 0;
        ux = 0;
        vx = 0;
        uy = 0;
        vy = 0;
        
        mi = 0;
        for dof = 1:3
            zi = zi + eta(elem(l,dof))*phi(dof,pt);
            hi = hi + h(elem(l,dof))*phi(dof,pt);
            zx = zx + eta(elem(l,dof))*phix(dof,pt);
            hx = hx + h(elem(l,dof))*phix(dof,pt);
            zy = zy + eta(elem(l,dof))*phiy(dof,pt);
            hy = hy + h(elem(l,dof))*phiy(dof,pt);
            
            ui = ui + u(elem(l,dof))*phi(dof,pt);
            ux = ux + u(elem(l,dof))*phix(dof,pt);
            uy = uy + u(elem(l,dof))*phiy(dof,pt);
            
            vi = vi + v(elem(l,dof))*phi(dof,pt);
            vx = vx + v(elem(l,dof))*phix(dof,pt);
            vy = vy + v(elem(l,dof))*phiy(dof,pt);
            
            mi = mi + mansoln(elem(l,dof))*phi(dof,pt);
        end
        for i = 1:3
            for j = 1:3                
                % Linear Flatbed
                P1xcfLin(i,j) = P1xcfLin(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/2*hi^2*phix(i,pt) ) *phix(j,pt);
                P1ycfLin(i,j) = P1ycfLin(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/2*hi^2*phiy(i,pt) ) *phiy(j,pt);
                P1cfLin(i,j)  = P1cfLin(i,j)  + A/2*cub2d(pt,3)* ...
                    ( 5/4*phi(i,pt) ) *phi(j,pt);
                
%                 % Linear Varying Bathy
%                 P1xcfLin(i,j) = P1xcfLin(i,j) + A/2*cub2d(pt,3)* ...
%                     ( 1/4*hi*hx*phi(i,pt) )  *phix(j,pt);
%                 P1ycfLin(i,j) = P1ycfLin(i,j) + A/2*cub2d(pt,3)* ...
%                     ( 1/4*hi*hy*phi(i,pt) )  *phiy(j,pt);
%                 P1cfLin(i,j)  = P1cfLin(i,j)  + A/2*cub2d(pt,3)* ...
%                     ( -1/4*((hx^2+hy^2)*phi(i,pt)+2*hi*...
%                              (hy*phiy(i,pt)+hx*phix(i,pt))) )  *phi(j,pt);
%                 
%                 % Nonlinear Flat
%                 P1xcfNon(i,j) = P1xcfNon(i,j) + A/2*cub2d(pt,3)* ...
%                     ( 1/6*(3*zi^2*phix(i,pt)+4*hi*zx*phi(i,pt)+...
%                            zi*(6*hi*phix(i,pt)+4*zx*phi(i,pt))) ) *phix(j,pt);
%                 P1ycfNon(i,j) = P1ycfNon(i,j) + A/2*cub2d(pt,3)* ...
%                     ( 1/6*(3*zi^2*phiy(i,pt)+4*hi*zy*phi(i,pt)+...
%                            zi*(6*hi*phiy(i,pt)+4*zy*phi(i,pt))) ) *phiy(j,pt);
%                 P1cfNon(i,j)  = P1cfNon(i,j)  + A/2*cub2d(pt,3)* ...
%                     ( 1/3*((hi+zi)*(zx*phix(i,pt)+zy*phiy(i,pt))+...
%                     3*phi(i,pt)*(zx^2+zy^2)) )   *phi(j,pt);
%                 
%                 % Nonlinear Varying
%                 P1xcfNon(i,j) = P1xcfNon(i,j) + A/2*cub2d(pt,3)* ...
%                     ( 1/4*zi*hx*phi(i,pt) ) *phix(j,pt);
%                 P1ycfNon(i,j) = P1ycfNon(i,j) + A/2*cub2d(pt,3)* ...
%                     ( 1/4*zi*hy*phi(i,pt) ) *phiy(j,pt);
%                 P1cfNon(i,j)  = P1cfNon(i,j)  + A/2*cub2d(pt,3)* ...
%                     ( -1/2*(zi*(hx*phix(i,pt)+hy*phiy(i,pt))+...
%                             phi(i,pt)*(hy*zy+hx*zx)) )   *phi(j,pt);
            end
%             P1rhsLin(i,1) = P1rhsLin(i,1) - A/2*cub2d(pt,3)* ...
%                 ( zi ); 

            % Linear Flat
            P1rhsLin(i,1) = P1rhsLin(i,1) - A/2*cub2d(pt,3)* ...
                ( 1/12*g*(-15*zi*phi(i,pt)+4*hi^2*(zx*phix(i,pt)+zy*phiy(i,pt))) );
%             % Linear Varying
%             P1rhsLin(i,1) = P1rhsLin(i,1) - A/2*cub2d(pt,3)* ...
%                 ( 1/12*g*(3*zi*(phi(i,pt)*(hy^2+hx^2)+2*hi*(hy*phiy(i,pt)+hx*phix(i,pt)))+...
%                           2*hi*phi(i,pt)*(hy*zy+hx*zx)) );
%             % Nonlinear Flat
%             P1rhsNon(i,1) = P1rhsNon(i,1) - A/2*cub2d(pt,3)* ...
%                 ( 1/6*(10*zi^2*phi(i,pt)*(vy^2+ux*vy+ux^2+uy*vx)+...
%                        hi*phi(i,pt)*(15*ui*vi*hxy+2*(3*g*(zy^2+zx^2)+...
%                           5*hi*(vy^2+vy*ux+ux^2+uy*vx)))+...
%                        zi*(15*ui*vi*hxy*phi(i,pt)+2*hi*g*(zx*phix(i,pt)+zy*phiy(i,pt))+...
%                           10*phi(i,pt)*(vy^2+vy*ux+ux^2+uy*vx))) );
%             % Nonlinear Varying
%             P1rhsNon(i,1) = P1rhsNon(i,1) - A/2*cub2d(pt,3)* ...
%                 ( -1/12*(-6*g*zi^2*(hy*phiy(i,pt)+hx*phix(i,pt))+...
%                 15*vi^2*hy*((zi+hi)*phiy(i,pt)+(zy+hy)*phi(i,pt))+...
%                 30*(zi+hi)*vi*vy*hy+...
%                 15*ui*hx*(ui*(hi*phix(i,pt)+phi(i,pt)*(hx+zx))+2*hi*phi(i,pt)*ux)+...
%                 zi*(15*ui^2*hx*phix(i,pt)-8*g*phi(i,pt)*(hy*zy+hx*zx)+...
%                   30*ui*phi(i,pt)*hx*ux)) );
%             

            P1rhsMan(i,1) = P1rhsMan(i,1) + A/2*cub2d(pt,3)*mi*phi(i,pt);
            
        end
        P1xcf = P1xcfLin+P1xcfNon;
        P1ycf = P1ycfLin+P1ycfNon;
        P1cf  = P1cfLin+P1cfNon;
        P1rhs = (P1rhsLin+P1rhsNon)*(1-manyes)+P1rhsMan*manyes;
    end
    LHS(elem(l,:),elem(l,:)) = LHS(elem(l,:),elem(l,:)) + (P1xcf+P1ycf+P1cf);
    RHS(elem(l,:)) = RHS(elem(l,:)) + P1rhs;
    
    p1xcfmat{l,1} = P1xcf;
    p1ycfmat{l,1} = P1ycf;
    p1cfmat{l,1}  = P1cf;
    p1rhsmat{l,1} = P1rhs;
end

LHS = sparse(LHS);
P1tst = LHS\RHS;


figure(1)
trisurf(elem,xvec,yvec,eta)
view([6.5 15])
title('eta')

figure(2)
trisurf(elem,xvec,yvec,P1tst)
view([6.5 15])
title('P1tst')

figure(3)
trisurf(elem,xvec,yvec,-3/4*(eta*g-P1tst))
view([6.5 15])
title('g*eta-P1tst')

figure(4)
trisurf(elem,xvec,yvec,pcg)
view([6.5 15])
title('pcg')

figure(5)
trisurf(elem,xvec,yvec,pcg-P1tst)
view([6.5 15])
title('pcg-P1tst')

figure(6)
trisurf(elem,xvec,yvec,p2)
view([6.5 15])
title('p2')

figure(7)
surf(eta2)
view([6.5 15])
title('eta2')

pdtmp= -1/2*(eta+h).*(g*eta-P1tst);
pbtmp = -1/3*(3*g*eta-3*P1tst+p2);

figure(8)
trisurf(elem,xvec,yvec,pd)
view([6.5 15])
title('pd')

figure(9)
trisurf(elem,xvec,yvec,pb)
view([6.5 15])
title('pb')

figure(10)
trisurf(elem,xvec,yvec,pdtmp)
view([6.5 15])
title('pdtmp')

figure(11)
trisurf(elem,xvec,yvec,pbtmp)
view([6.5 15])
title('pbtmp')

% figure(6)
% subplot(2,1,1)
% plot(rhscg)
% subplot(2,1,2)
% plot(RHS/2)