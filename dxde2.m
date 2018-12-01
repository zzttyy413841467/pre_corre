function xdot=dxde2(e,x,sigma)
constant_sim;

% H0=121900;
% Hf=30480;
% V0=7600;
% Vf=910;
% r0=(Re+H0)/Re;
% rf=(Re+Hf)/Re;
% v0=V0/Vc;
% vf=Vf/Vc;
% 
% e0=1/r0-v0^2/2;
% ef=1/rf-vf^2/2;
alpha0=45/180*pi;
% ef=evalin('base','ef');
% e0=evalin('base','e0');
s=x(1);
r=x(2);
gamma=x(3);
V=sqrt(2*(1./r-e));

ma=V*Vc/340;
alpha=alpha0.*(ma>=10)+((45-0.612*(ma-10).^2)/180*pi).*(ma<10);

Cl = cl0+cl1*alpha+cl2.*alpha.^2;
Cd = cd0+cd1*Cl+cd2.*Cl.^2;

h=r*Re-Re;
rho=rho0.*exp(-h./hs);
q=1/2*rho.*V.^2;
L=q.*Cl.*S/m*Re;
D=q.*Cd.*S/m*Re;
% g=mu./(Re+h).^2;

rdot=sin(gamma)./D;
gammadot=1./(D.*V.^2).*(L.*cos(sigma)+(V.^2./r-1./r.^2).*cos(gamma));
sdot=-cos(gamma)./(r.*D);
taudot=1./(D.*V);
xdot=[sdot;rdot;gammadot;taudot];
end