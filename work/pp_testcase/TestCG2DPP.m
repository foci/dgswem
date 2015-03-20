%% Test CG 2d PP
clear all

g = 9.81;

x = linspace(-5,5,81);
y = linspace(-5,5,81);
[X,Y] = meshgrid(x,y);
X = X';
Y = Y';

xvec = X(:);
yvec = Y(:);

elem = delaunay(xvec,yvec);
ne = length(elem(:,1));


triplot(elem,xvec,yvec);

%%
e0 = 0.02;
h0 = 1;
h1 = 0.2;
sig = 2;
p0 = g;
eta = e0*exp(-(xvec.^2+yvec.^2)/sig);
etax = -2/sig*e0*exp(-(xvec.^2+yvec.^2)/sig).*xvec;
etay = -2/sig*e0*exp(-(xvec.^2+yvec.^2)/sig).*yvec;
etaxx = e0*(4*exp(-(xvec.^2+yvec.^2)/sig).*xvec.^2/sig^2-2*exp(-(xvec.^2+yvec.^2)/sig)/sig);
etayy = e0*(4*exp(-(xvec.^2+yvec.^2)/sig).*yvec.^2/sig^2-2*exp(-(xvec.^2+yvec.^2)/sig)/sig);

h = 0*xvec+1;

manyes = 1;

% mansoln = -1/12*p0*(8*(etay.^2+etax.^2)+10*eta.^2.*(etayy+etaxx)+6*(etayy+etaxx)+...
%     eta.*(-15+16*(etaxx+etayy)));

mansoln = 0;
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

figure(1)
trisurf(elem,xvec,yvec,eta)

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

for l = 1:ne
    xloc = xvec(elem(l,:));
    yloc = yvec(elem(l,:));
    A = 1/2*det([1 1 1;xloc';yloc']');
    
    
    phix(1,:) = 1/2/A*(yloc(2)-yloc(3));
    phix(2,:) = 1/2/A*(yloc(3)-yloc(1));
    phix(3,:) = 1/2/A*(yloc(1)-yloc(2));
    
    phiy(1,:) = 1/2/A*(xloc(3)-xloc(2));
    phiy(2,:) = 1/2/A*(xloc(1)-xloc(3));
    phiy(3,:) = 1/2/A*(xloc(2)-xloc(1));    
    
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
        
        mi = 0;
        for dof = 1:3
            zi = zi + eta(elem(l,dof))*phi(dof,pt);
            hi = hi + h(elem(l,dof))*phi(dof,pt);
            zx = zx + eta(elem(l,dof))*phix(dof,pt);
            hx = hx + h(elem(l,dof))*phix(dof,pt);
            zy = zy + eta(elem(l,dof))*phiy(dof,pt);
            hy = hy + h(elem(l,dof))*phiy(dof,pt);
            
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
                
                % Linear Varying Bathy
                P1xcfLin(i,j) = P1xcfLin(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/4*hi*hx*phi(i,pt) )  *phix(j,pt);
                P1ycfLin(i,j) = P1ycfLin(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/4*hi*hy*phi(i,pt) )  *phiy(j,pt);
                P1cfLin(i,j)  = P1cfLin(i,j)  + A/2*cub2d(pt,3)* ...
                    ( -1/4*((hx^2+hy^2)*phi(i,pt)+2*hi*...
                             (hy*phiy(i,pt)+hx*phix(i,pt))) )  *phi(j,pt);
                
                % Nonlinear Flat
                P1xcfNon(i,j) = P1xcfNon(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/6*(3*zi^2*phix(i,pt)+4*hi*zx*phi(i,pt)+...
                           zi*(6*hi*phix(i,pt)+4*zx*phi(i,pt))) ) *phix(j,pt);
                P1ycfNon(i,j) = P1ycfNon(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/6*(3*zi^2*phiy(i,pt)+4*hi*zy*phi(i,pt)+...
                           zi*(6*hi*phiy(i,pt)+4*zy*phi(i,pt))) ) *phiy(j,pt);
                P1cfNon(i,j)  = P1cfNon(i,j)  + A/2*cub2d(pt,3)* ...
                    ( 1/3*((hi+zi)*(zx*phix(i,pt)+zy*phiy(i,pt))+...
                    3*phi(i,pt)*(zx^2+zy^2)) )   *phi(j,pt);
                
                % Nonlinear Varying
                P1xcfNon(i,j) = P1xcfNon(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/4*zi*hx*phi(i,pt) ) *phix(j,pt);
                P1ycfNon(i,j) = P1ycfNon(i,j) + A/2*cub2d(pt,3)* ...
                    ( 1/4*zi*hy*phi(i,pt) ) *phiy(j,pt);
                P1cfNon(i,j)  = P1cfNon(i,j)  + A/2*cub2d(pt,3)* ...
                    ( -1/2*(zi*(hx*phix(i,pt)+hy*phiy(i,pt))+...
                            phi(i,pt)*(hy*zy+hx*zx)) )   *phi(j,pt);
            end
            % Linear Flat
            P1rhsLin(i,1) = P1rhsLin(i,1) - A/2*cub2d(pt,3)* ...
                ( 1/12*g*(-15*zi*phi(i,pt)+4*hi^2*(zx*phix(i,pt)+zy*phiy(i,pt))) );
            % Linear Varying
%             P1rhsLin(i,1) = P1rhsLin(i,1) - A/2*cub2d(pt,3)* ...
%                 ( 1/12*g*(3*zi*(phi(i,pt)*(hy^2+hx^2)+2*hi*(hy*phiy(i,pt)+hx*phix(i,pt)))+...
%                           2*hi*phi(i,pt)*(hy*zy+hx*zx)) );
            % Nonlinear Flat
%             P1rhsNon(i,1) = P1rhsNon(i,1) - A/2*cub2d(pt,3)* ...
%                 ( 1/6*(10*zi^2*phi(i,pt)*(vy^2+ux*vy+ux^2+uy*vx)+...
%                        hi*phi(i,pt)*(15*ui*vi*hxy+2*(3*g*(zy^2+zx^2)+...
%                           5*hi*(vy^2+vy*ux+ux^2+uy*vx)))+...
%                        zi*(15*ui*vi*hxy*phi(i,pt)+2*hi*g*(zx*phix(i,pt)+zy*phiy(i,pt))+...
%                           10*phi(i,pt)*(vy^2+vy*ux+ux^2+uy*vx))) );
            % Nonlinear Varying
%             P1rhsNon(i,1) = P1rhsNon(i,1) - A/2*cub2d(pt,3)* ...
%                 ( -1/12*(-6*g*zi^2*(hy*phiy(i,pt)+hx*phix(i,pt))+...
%                 15*vi^2*hy*((zi+hi)*phiy(i,pt)+(zy+hy)*phi(i,pt))+...
%                 30*(zi+hi)*vi*vy*hy+...
%                 15*ui*hx*(ui*(hi*phix(i,pt)+phi(i,pt)*(hx+zx))+2*hi*phi(i,pt)*ux)+...
%                 zi*(15*ui^2*hx*phix(i,pt)-8*g*phi(i,pt)*(hy*zy+hx*zx)+...
%                   30*ui*phi(i,pt)*hx*ux)) );
            

            P1rhsMan(i,1) = P1rhsMan(i,1) + A/2*cub2d(pt,3)*mi*phi(i,pt);
            
        end
        P1xcf = P1xcfLin+P1xcfNon;
        P1ycf = P1ycfLin+P1ycfNon;
        P1cf  = P1cfLin+P1cfNon;
        P1rhs = (P1rhsLin+P1rhsNon)*(1-manyes)+P1rhsMan*manyes;
    end
    LHS(elem(l,:),elem(l,:)) = LHS(elem(l,:),elem(l,:)) + (P1xcf+P1ycf+P1cf);
    RHS(elem(l,:)) = RHS(elem(l,:)) + P1rhs;
end

LHS = sparse(LHS);
P1tst = LHS\RHS;
figure(2)
trisurf(elem,xvec,yvec,P1tst)

figure(3)
trisurf(elem,xvec,yvec,p0*eta)

figure(4)
trisurf(elem,xvec,yvec,p0*eta-P1tst)